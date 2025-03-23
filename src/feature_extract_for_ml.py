import pandas as pd
import numpy as np
import mysql.connector
import umap
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors
import matplotlib.pyplot as plt
import seaborn as sns
from mysql.connector import Error
import sys


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
            "SELECT molecule_chembl_id, canonical_smiles, standard_type, standard_value, standard_units, kinase_name, protein_sequence FROM bioactivity")
        # Fetch all rows from the executed query
        rows = cursor.fetchall()

        # Convert the rows into a pandas DataFrame
        df = pd.DataFrame(rows, columns=["molecule_chembl_id", "canonical_smiles", "standard_type",
                                         "standard_value", "standard_units", "kinase_name", "protein_sequence"])

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



# Convert to pIC50 (vectorized)
def convert_to_pic50_vectorized(df):

    # Step 1: Get the molecular weights for all SMILES in the dataset
    df['mol_weight'] = get_molecular_weights(df['canonical_smiles'])


    # Step 2:
    # Handle µg/mL ang mg/kg to nM conversion only for valid cases
    # Apply the conversion using np.where for both unit types
    df['standard_value'] = np.where(
        df['standard_units'] == "ug.mL-1",  # If the unit is µg/mL
        df['standard_value'] / df['mol_weight'] * 10 ** 9,  # Directly calculate nM for µg/mL
        np.where(
            df['standard_units'] == "mg kg-1",  # If the unit is mg/kg
            df['standard_value'] * 1000 / df['mol_weight'] * 10 ** 9,  # Convert mg to µg first, then calculate nM
            df['standard_value']  # If the unit is neither µg/mL nor mg/kg, leave it unchanged
        )
    )

    # Convert all relevant standard types to pIC50


    df['pIC50'] = -np.log10(df['standard_value'] * 1e-9)


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
    df = df[df['standard_type'].isin(['IC50', 'pIC50'])]

    df['standard_value'] = pd.to_numeric(df['standard_value'], errors='coerce')
    df = df.dropna()


    # List of allowed values for standard_units
    allowed_values = ['mg kg-1', 'ug', 'ug.mL-1', 'nM']

    # Check if any value in the set is not in the allowed values
    invalid_values = {unit for unit in set(df['standard_units']) if unit not in allowed_values}

    # Check if there are any 'pIC50' values in standard_type column
    invalid_pIC50 = df['standard_type'].isin(['pIC50'])

    # If there are invalid values in 'standard_units' or 'pIC50' in 'standard_type', stop the program
    if invalid_values:
        print("Invalid values found in 'standard_units':", invalid_values)
        sys.exit("Stopping program due to invalid 'standard_units' values.")
    elif invalid_pIC50.any():  # Check if there's at least one 'pIC50' value
        print("Found 'pIC50' in 'standard_type' column, which is not allowed.")
        sys.exit("Stopping program due to 'pIC50' values in 'standard_type'.")
    else:
        print("All values are valid.")


    df = df[df['standard_units'].isin(['mg kg-1', 'ug.mL-1', 'nM'])]
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
    # Increase the width of the figure to give more space for the legend
    fig, ax = plt.subplots(figsize=(14, 8))  # Increase width to 14

    # Normalize pIC50 values between 0 and 1 for transparency
    df_umap["alpha"] = df_umap["pIC50"].rank(pct=True)

    # Define sizes for each pIC50 range
    conditions = [
        (df_umap["pIC50"] > 6),  # Strong binding
        (df_umap["pIC50"] <= 6) & (df_umap["pIC50"] >= 5),  # Moderate binding
        (df_umap["pIC50"] < 5)  # Weak binding
    ]
    sizes = [150, 75, 25]  # Assign size: 300 for strong, 150 for moderate, 50 for weak

    # Assign sizes based on pIC50 ranges
    df_umap["size"] = np.select(conditions, sizes, default=50)

    # Assign colors by kinase type
    kinase_types = df_umap["kinase_name"].unique()
    color_map = {kinase: color for kinase, color in zip(kinase_types, sns.color_palette("Set1", len(kinase_types)))}

    # Scatter plot using Matplotlib
    for kinase in kinase_types:
        subset = df_umap[df_umap["kinase_name"] == kinase]
        ax.scatter(
            subset["UMAP1"],
            subset["UMAP2"],
            color=color_map[kinase],  # Color by kinase type (not pIC50)
            alpha=subset["alpha"],  # Individual transparency per point
            s=subset["size"],  # Point size varies with pIC50 ranges (not color)
            label=kinase
        )

    # Add a legend for kinase types, placed outside the plot
    ax.set_xlabel("UMAP Dimension 1")
    ax.set_ylabel("UMAP Dimension 2")
    ax.set_title("Molecular Similarity Map (Color by Kinase, Transparency & Size by pIC50)")

    # Move the legend outside the plot area to avoid overlap with points
    ax.legend(title="Kinase Type", loc="upper left", bbox_to_anchor=(1.05, 1), borderaxespad=0.)

    # Adjust the plot to give more space on the right for the legend
    plt.subplots_adjust(right=0.6)  # Adjust this value to your preference

    # Apply tight layout for better spacing inside the plot
    plt.tight_layout()

    # Show the plot
    plt.show()


# Main function
def main():


    # to be returned once we use the sql on the ubuntu comptuer
    """
    # Fetch data from MySQL
    df = fetch_data_from_mysql()

    if df.empty:
        print("❌ No data fetched from MySQL.")
        return
    """

    # to be commented out , once we use ubuntu and not windows
    df = pd.read_csv('../data/processed_data_for_feature_extract.csv')


    # Process the data (filtering, feature selection, UMAP)
    df_umap = process_data(df)

    # Save the final dataset for ML training
    df_umap.to_csv("../data/processed_data_for_ml.csv", index=False)
    print("✅ Processed dataset saved as 'processed_data_for_ml.csv'!")

    # Generate UMAP visualization of the processed data
    plot_umap(df_umap)


if __name__ == "__main__":
    main()
