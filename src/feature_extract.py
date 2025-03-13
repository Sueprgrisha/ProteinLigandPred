import pandas as pd
import numpy as np
import mysql.connector
import umap
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors
import matplotlib.pyplot as plt
import seaborn as sns
from mysql.connector import Error

# MySQL Database Config
MYSQL_CONFIG = {
    "host": "localhost",  # Change to your MySQL server
    "user": "root",
    "password": "Grisha@123!",
    "database": "chembl_db"
}


# Fetch data from MySQL
def fetch_data_from_mysql():
    """Fetches data from the MySQL database and returns it as a DataFrame."""
    try:
        # Establish the connection
        conn = mysql.connector.connect(**MYSQL_CONFIG)
        cursor = conn.cursor()

        # Execute the query to retrieve data from the bioactivity table
        cursor.execute(
            "SELECT molecule_chembl_id, canonical_smiles, standard_type, standard_value, standard_units, kinase_name FROM bioactivity")

        # Fetch all rows from the executed query
        rows = cursor.fetchall()

        # Convert the rows into a pandas DataFrame
        df = pd.DataFrame(rows, columns=["molecule_chembl_id", "canonical_smiles", "standard_type",
                                         "standard_value", "standard_units", "kinase_name"])

        # Close the cursor and connection
        cursor.close()
        conn.close()

        # Return the DataFrame
        return df

    except Error as e:
        print(f"❌ Error fetching data from MySQL: {e}")
        return pd.DataFrame()  # Return an empty DataFrame if there is an error


# Function to get the molecular weights for a list of SMILES in a vectorized way
def get_molecular_weights(smiles_list):
    # Convert SMILES list to a pandas series
    smiles_series = pd.Series(smiles_list)

    # Calculate the molecular weight for each SMILES
    molecular_weights = smiles_series.apply(
        lambda smiles: Descriptors.MolWt(Chem.MolFromSmiles(smiles)) if Chem.MolFromSmiles(smiles) else np.nan)

    return molecular_weights


# Vectorized conversion of µg/mL to nM using molecular weights
def convert_to_nM_vectorized(values, molecular_weights):
    # Convert µg/mL to nM using the formula
    concentration_in_mg_per_ml = values / 1000  # Convert from µg to mg
    concentration_in_mol_per_l = concentration_in_mg_per_ml / molecular_weights * 1000  # Convert to mol/L
    concentration_in_nM = concentration_in_mol_per_l * 1e9  # Convert to nM (1 mol/L = 1e9 nM)

    return concentration_in_nM


# Convert to pIC50 (vectorized)
def convert_to_pic50_vectorized(df):

    # Step 1: Get the molecular weights for all SMILES in the dataset
    df['mol_weight'] = get_molecular_weights(df['canonical_smiles'])

    # Step 2:
    # Handle µg/mL to nM conversion only for valid cases
    # The function takes only rows where the units are ug.ml-1 and convert the values for these lines to nM
    df['standard_value'] = np.where(df['standard_units'] == "ug.mL-1",
                                    convert_to_nM_vectorized(df['standard_value'], df['mol_weight']),
                                    df['standard_value'])








    # Next task: Maybe use only ic50 and pic50?











    # Convert all relevant standard types to pIC50
    mask = df['standard_type'].isin(["IC50", "pIC50", "Ki", "Ki low", "Ki high", "MIC50", "Kd"])
    df.loc[mask, 'pIC50'] = -np.log10(df.loc[mask, 'standard_value'] * 1e-9)

    # For rows where we can't convert, set pIC50 to NaN
    df.loc[~mask, 'pIC50'] = np.nan

    return df

# Convert SMILES to Morgan fingerprints
def smiles_to_fp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)
    arr = np.zeros((1, 1024), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr.flatten()  # Ensure the fingerprint is flattened for consistency


# Process data
def process_data(df):
    # Filter rows where 'standard_type' is one of the convertible types
    df = df[df['standard_type'].isin(['IC50', 'pIC50', 'Ki', 'Ki low', 'Ki high', 'MIC50', 'Kd'])]

    df['standard_value'] = pd.to_numeric(df['standard_value'], errors='coerce')
    df = df.dropna()

    # Apply the vectorized conversion to pIC50
    df = convert_to_pic50_vectorized(df)

    # Drop rows where the pIC50 could not be calculated (invalid values)
    df = df.dropna(subset=["pIC50"]).reset_index(drop=True)

    # Drop rows where fingerprint is None (invalid SMILES)
    df["fingerprint"] = df["canonical_smiles"].apply(smiles_to_fp)
    df = df.dropna(subset=["fingerprint"]).reset_index(drop=True)

    # Perform UMAP dimensionality reduction on the fingerprints
    fingerprints = np.vstack(df["fingerprint"])  # Stack fingerprints into a 2D array
    umap_embedding = umap.UMAP(n_components=2).fit_transform(fingerprints)

    # Add UMAP components to the dataframe
    df_umap = pd.DataFrame(umap_embedding, columns=["UMAP1", "UMAP2"])
    df_umap["pIC50"] = df["pIC50"].values
    df_umap["kinase_name"] = df["kinase_name"].values
    df_umap["canonical_smiles"] = df["canonical_smiles"].values
    return df_umap


# Visualize UMAP results (for visualization purposes)
def plot_umap(df_umap):
    fig, ax = plt.subplots(figsize=(10, 8))

    # Normalize pIC50 values between 0 and 1 for transparency
    df_umap["alpha"] = df_umap["pIC50"].rank(pct=True)

    # Normalize pIC50 values for point size scaling (e.g., size between 20 and 200)
    df_umap["size"] = 20 + (df_umap["pIC50"] - df_umap["pIC50"].min()) / (
            df_umap["pIC50"].max() - df_umap["pIC50"].min()) * 180

    # Assign colors by kinase type
    kinase_types = df_umap["kinase_name"].unique()
    color_map = {kinase: color for kinase, color in zip(kinase_types, sns.color_palette("Set1", len(kinase_types)))}

    # Scatter plot using Matplotlib
    for kinase in kinase_types:
        subset = df_umap[df_umap["kinase_name"] == kinase]
        ax.scatter(
            subset["UMAP1"],
            subset["UMAP2"],
            color=color_map[kinase],
            alpha=subset["alpha"],  # Individual transparency per point
            s=subset["size"],  # Point size varies with pIC50
            label=kinase
        )

    ax.set_xlabel("UMAP Dimension 1")
    ax.set_ylabel("UMAP Dimension 2")
    ax.set_title("Molecular Similarity Map (Color by Kinase, Transparency & Size by pIC50)")
    ax.legend(title="Kinase Type", loc="best")
    plt.show()


# Main function
def main():
    # Fetch data from MySQL
    df = fetch_data_from_mysql()

    if df.empty:
        print("❌ No data fetched from MySQL.")
        return

    # Process the data (filtering, feature selection, UMAP)
    df_umap = process_data(df)

    # Save the final dataset for ML training
    df_umap.to_csv("processed_data_for_ml.csv", index=False)
    print("✅ Processed dataset saved as 'processed_data_for_ml.csv'!")

    # Generate UMAP visualization of the processed data
    plot_umap(df_umap)


if __name__ == "__main__":
    main()
