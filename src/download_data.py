import pandas as pd
import requests
import mysql.connector
import asyncio
import aiohttp
from mysql.connector import Error

# ChEMBL API URLs
ASSAY_URL = "https://www.ebi.ac.uk/chembl/api/data/assay/"
TARGET_URL = "https://www.ebi.ac.uk/chembl/api/data/target/"

# MySQL Database Config
MYSQL_CONFIG = {
    "host": "localhost",  # Change to your MySQL server
    "user": "root",
    "password": "Grisha@123!",
    "database": "chembl_db"
}

# Constants
LIMIT = 1000  # Number of records per request
MAX_RECORDS = 1000  # Total records to retrieve
MIN_LIGAND_PER_KINASE = 2

# Fetching data from ChEMBL API asynchronously
async def get_kinase_info_bulk(assay_ids):
    """Fetches kinase names for multiple assay IDs in bulk using async calls."""
    kinase_map = {}
    async with aiohttp.ClientSession() as session:
        tasks = []

        for assay_id in assay_ids:
            tasks.append(fetch_assay_data(session, assay_id, kinase_map))

        # Run all tasks asynchronously
        await asyncio.gather(*tasks)

    return kinase_map


async def fetch_assay_data(session, assay_id, kinase_map):
    """Fetches the assay data and target data for a given assay ID."""
    assay_url = f"{ASSAY_URL}{assay_id}.json"
    async with session.get(assay_url) as response:
        if response.status == 200:
            assay_data = await response.json()
            target_id = assay_data.get("target_chembl_id", None)

            if target_id:
                target_url = f"{TARGET_URL}{target_id}.json"
                async with session.get(target_url) as target_response:
                    if target_response.status == 200:
                        target_data = await target_response.json()
                        kinase_map[assay_id] = target_data.get("pref_name", "Unknown")
                    else:
                        kinase_map[assay_id] = "Unknown"
            else:
                kinase_map[assay_id] = "Unknown"
        else:
            kinase_map[assay_id] = "Unknown"


# Main function to fetch and process data
def fetch_data():
    """Fetches bioactivity data from ChEMBL and filters for IC50/pIC50 values."""
    all_data = pd.DataFrame(
        columns=["molecule_chembl_id", "canonical_smiles", "assay_chembl_id", "standard_type", "standard_value",
                 "standard_units"])

    a, b = 0, 0
    offset = 0
    # Get at least MAX_RECORDS amount of points of data, and make sure at least 2 proteins have "MIN_LIGAND_PER_KINASE" amount of ligand pairing examples
    while not (len(all_data) >= MAX_RECORDS and a >= MIN_LIGAND_PER_KINASE and b >= MIN_LIGAND_PER_KINASE):
        url = f"https://www.ebi.ac.uk/chembl/api/data/activity.json?limit={LIMIT}&offset={offset}"
        response = requests.get(url)

        if response.status_code != 200:
            print(f"❌ Failed to fetch data (status code: {response.status_code})")
            break

        data = response.json()
        activities = data.get("activities", [])

        if not activities:
            break  # Stop when no more records are available

        df = pd.DataFrame(activities)[
            ["molecule_chembl_id", "canonical_smiles", "assay_chembl_id", "standard_type", "standard_value",
             "standard_units"]]

        filtered_df = df[df['standard_type'].isin(['IC50', 'pIC50', 'Ki', 'Ki low', 'Ki high', 'MIC50', 'Kd'])]

        all_data = pd.concat([all_data, filtered_df], ignore_index=True)

        most_common_assays = all_data['assay_chembl_id'].value_counts().head(2)
        a, b = most_common_assays[0], most_common_assays[1]

        print(most_common_assays)
        offset += LIMIT
        print(f"✅ Downloaded {len(all_data)} records...", end="\r")

    # Get kinase names in bulk asynchronously
    assay_ids = list(set(all_data["assay_chembl_id"]))
    kinase_map = asyncio.run(get_kinase_info_bulk(assay_ids))

    # Add kinase names to DataFrame
    all_data["kinase_name"] = all_data["assay_chembl_id"].map(kinase_map)

    return all_data


# MySQL functions and main script as before
def create_mysql_table():
    """Creates the MySQL table if it doesn't exist."""
    try:
        conn = mysql.connector.connect(**MYSQL_CONFIG)
        cursor = conn.cursor()

        cursor.execute("""
        CREATE TABLE IF NOT EXISTS bioactivity (
            id INT AUTO_INCREMENT PRIMARY KEY,
            molecule_chembl_id VARCHAR(50),
            canonical_smiles TEXT,
            standard_type VARCHAR(20),
            standard_value FLOAT,
            standard_units VARCHAR(10),
            kinase_name VARCHAR(255)
        );
        """)

        conn.commit()
        cursor.close()
        conn.close()
        print("✅ MySQL Table Ready!")
    except Error as e:
        print(f"❌ Error creating MySQL table: {e}")


def store_data(df):
    """Stores the fetched data in a MySQL database using batch insert."""
    try:
        conn = mysql.connector.connect(**MYSQL_CONFIG)
        cursor = conn.cursor()

        query = """
        INSERT INTO bioactivity (molecule_chembl_id, canonical_smiles, standard_type, standard_value, standard_units, kinase_name)
        VALUES (%s, %s, %s, %s, %s, %s);
        """
        data = df.values.tolist()

        cursor.executemany(query, data)  # Bulk Insert
        conn.commit()

        cursor.close()
        conn.close()
        print("\n✅ Data stored in MySQL database!")
    except Error as e:
        print(f"❌ Error inserting data into MySQL: {e}")


def main():
    df = fetch_data()
    if not df.empty:
        df = df[["molecule_chembl_id", "canonical_smiles", "standard_type", "standard_value", "standard_units",
                 "kinase_name"]]
        store_data(df)


if __name__ == "__main__":
    main()
