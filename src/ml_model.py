import pandas as pd
import numpy as np
import umap
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras import backend as K
from feature_extract import smiles_to_fp
from sklearn.neighbors import KNeighborsClassifier
from collections import Counter
from scipy.spatial.distance import cdist


# ------------------- Data Loading -------------------

def load_data(filepath="processed_data_for_ml.csv"):
    df = pd.read_csv(filepath)
    print("Columns in the loaded dataframe:", df.columns)
    df["saved_smiles"] = df["canonical_smiles"]
    return df



# ------------------- UMAP Plotting -------------------

def plot_umap(df_umap):
    if not {'UMAP1', 'UMAP2', 'pIC50', 'kinase_name'}.issubset(df_umap.columns):
        print("Missing necessary columns for UMAP plotting.")
        return

    fig, ax = plt.subplots(figsize=(10, 8))
    df_umap["alpha"] = df_umap["pIC50"].rank(pct=True)
    df_umap["size"] = 20 + (df_umap["pIC50"] - df_umap["pIC50"].min()) / (
            df_umap["pIC50"].max() - df_umap["pIC50"].min()) * 180
    kinase_types = df_umap["kinase_name"].unique()
    color_map = {kinase: color for kinase, color in zip(kinase_types, sns.color_palette("Set1", len(kinase_types)))}

    for kinase in kinase_types:
        subset = df_umap[df_umap["kinase_name"] == kinase]
        ax.scatter(subset["UMAP1"], subset["UMAP2"], color=color_map[kinase], alpha=subset["alpha"], s=subset["size"],
                   label=kinase)

    ax.set_xlabel("UMAP Dimension 1")
    ax.set_ylabel("UMAP Dimension 2")
    ax.set_title("Molecular Similarity Map")
    ax.legend(title="Kinase Type", loc="best")
    plt.show()

# ------------------- SMILES to Fingerprints -------------------

def convert_fingerprints_to_numeric(df):
    fingerprints = []
    for smiles in df["canonical_smiles"]:
        fp = smiles_to_fp(smiles)
        if fp is not None:
            fingerprints.append(fp)
        else:
            fingerprints.append(np.zeros(1024))  # Placeholder for invalid SMILES
    df["fingerprint"] = fingerprints
    return df


# ------------------- Triplet Loss Function -------------------

def triplet_loss(margin=1.0):
    def loss(y_true, y_pred):
        anchor, positive, negative = y_pred[0], y_pred[1], y_pred[2]
        pos_dist = K.sum(K.square(anchor - positive), axis=1)
        neg_dist = K.sum(K.square(anchor - negative), axis=1)
        return K.mean(K.maximum(pos_dist - neg_dist + margin, 0.0))  # Triplet margin loss

    return loss


# ------------------- Building Network -------------------

def build_base_network(input_dim):
    input_layer = layers.Input(shape=(input_dim,))
    x = layers.Dense(128, activation='relu')(input_layer)
    x = layers.Dense(64, activation='relu')(x)
    x = layers.Dense(32, activation='relu')(x)
    return keras.Model(input_layer, x)


def build_triplet_network(input_dim):
    base_network = build_base_network(input_dim)

    input_anchor = layers.Input(shape=(input_dim,))
    input_positive = layers.Input(shape=(input_dim,))
    input_negative = layers.Input(shape=(input_dim,))

    processed_anchor = base_network(input_anchor)
    processed_positive = base_network(input_positive)
    processed_negative = base_network(input_negative)

    model = keras.Model(inputs=[input_anchor, input_positive, input_negative],
                        outputs=[processed_anchor, processed_positive, processed_negative])
    model.compile(optimizer='adam', loss=triplet_loss(), metrics=['accuracy'])
    return model


# ------------------- Creating Triplets -------------------

def create_triplets(fingerprints, kinase_names, num_triplets):
    triplets = []
    kinase_dict = {kinase: [] for kinase in set(kinase_names)}

    # Group samples by kinase name
    for idx, kinase in enumerate(kinase_names):
        kinase_dict[kinase].append(fingerprints[idx])

    # Compute distances between all pairs of fingerprints
    dist_matrix = cdist(fingerprints, fingerprints, metric='euclidean')
    np.fill_diagonal(dist_matrix, np.inf)  # Ignore the diagonal (same sample)

    # Mask the lower triangle to remove duplicate distances
    mask = np.triu(np.ones_like(dist_matrix, dtype=bool), k=1)
    dist_matrix_half = dist_matrix.copy()
    dist_matrix_half[~mask] = np.inf  # Set the lower triangle and diagonal to infinity

    # For sampling the closest points, prioritize smaller distances
    dist_flat = dist_matrix_half.flatten()
    dist_flat_idx = np.argsort(dist_flat) # Get the indices of the closest pairs



    positive = None
    negative = None

    # Convert the list to a numpy array for efficient comparison
    kinase_names_array = np.array(kinase_names)

    # Use broadcasting to compare the array to itself
    kinase_matrix_mask = kinase_names_array[:, None] == kinase_names_array

    pos_matrix_dist = dist_matrix_half.copy()
    pos_matrix_dist[~kinase_matrix_mask] = np.inf

    pos_matrix_dist_flat = pos_matrix_dist.flatten()
    pos_dist_idx = np.argsort(pos_matrix_dist_flat)


    neg_matrix_dist = dist_matrix_half.copy()
    neg_matrix_dist[kinase_matrix_mask] = np.inf

    neg_matrix_dist_flat = neg_matrix_dist.flatten()
    neg_dist_idx = np.argsort(neg_matrix_dist_flat)


    for idx in dist_flat_idx:
        i, j = divmod(idx, len(fingerprints))  # Get the pair indices (i, j)

        anchor = fingerprints[i]

        anchor_distances = dist_matrix[i]  # distances from anchor to all others
        # Step 1: Create a boolean mask for entries that match the kinase at index i
        mask = (kinase_names == kinase_names[i])

        # Step 2: Extract the corresponding distances for the same kinase
        same_kinase_distances = anchor_distances[mask]

        # Step 3: Get the indices of the same-kinase distances, sorted by distance
        sorted_indices = np.argsort(same_kinase_distances)

        # To get the actual indices in the original `anchor_distances` array (i.e., indices that match the mask)
        sorted_indices_in_original_pos = np.where(mask)[0][sorted_indices]



        negative_distances[np.array(kinase_names) == kinase_names[i]] = np.inf  # Ignore the same kinase

        # Check for the positive sample
        if kinase_names[i] == kinase_names[j]:
            positive = fingerprints[j]
        else:
            # Find the closest negative from a different kinase
            negative_distances = dist_matrix[i]  # distances from anchor to all others
            negative_distances[np.array(kinase_names) == kinase_names[i]] = np.inf  # Ignore the same kinase
            negative_idx = np.argmin(negative_distances)  # Get the index of the closest different kinase
            negative = fingerprints[negative_idx]

        # If we have both a valid positive and a valid negative, form the triplet
        if positive is not None and negative is not None:
            triplets.append([np.array(anchor), np.array(positive), np.array(negative)])
            last_positive = None  # Reset last positive once triplet is formed
            last_negative = None
            positive = None
            negative = None

        if len(triplets) >= num_triplets:
            break

    return np.array(triplets)



# ------------------- Main Execution -------------------

def main():
    df = load_data()
    df = convert_fingerprints_to_numeric(df)

    fingerprints = df['fingerprint'].tolist()
    kinase_names = df['kinase_name'].tolist()

    # Split into train/test
    fingerprints_train, fingerprints_test, kinase_names_train, kinase_names_test = train_test_split(fingerprints,
                                                                                                    kinase_names,
                                                                                                    test_size=0.2,
                                                                                                    random_state=42)

    # Create triplets for training and testing
    X_train = create_triplets(fingerprints_train, kinase_names_train, num_triplets=10)
    X_test = create_triplets(fingerprints_test, kinase_names_test, num_triplets=10)

    # Build and train the triplet network
    triplet_model = build_triplet_network(1024)
    triplet_model.fit([X_train[:, 0], X_train[:, 1], X_train[:, 2]], np.zeros(len(X_train)), epochs=1, batch_size=32)

    loss, accuracy = triplet_model.evaluate([X_test[:, 0], X_test[:, 1], X_test[:, 2]], np.zeros(len(X_test)))
    print(f"Triplet Network Accuracy: {accuracy * 100:.2f}%")

    if accuracy > 0.9:
        print("Accuracy is high, plotting UMAP...")

        # Extract the base network from the Triplet model
        base_network = triplet_model.get_layer(index=2)  # Assuming index=1 is the base network

        # Create a new model that takes the same input and outputs the intermediate representation
        intermediate_model = keras.Model(inputs=base_network.input, outputs=base_network.output)

        # Get transformed representations
        X_test_original = np.vstack(fingerprints_test)

        X_test_transformed = intermediate_model.predict(X_test_original)

        # Plot UMAP before and after transformation
        plot_umap(X_test_original, X_test_transformed, kinase_names_test)


if __name__ == "__main__":
    main()
