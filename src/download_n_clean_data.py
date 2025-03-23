import pandas as pd
import requests
import mysql.connector
import asyncio
import aiohttp
import time
from mysql.connector import Error
from data_cleaning.data_cleaning import *
import json

def load_settings(file_path="settings.json"):
    with open(file_path, "r") as f:
        return json.load(f)



# MySQL Database Config
MYSQL_CONFIG = {
    "host": "localhost",  # Change to your MySQL server
    "user": "root",
    "password": "Grisha@123!",
    "database": "chembl_db"
}



# Fetching data from ChEMBL API asynchronously
async def get_kinase_info_bulk(assay_ids):
    """Fetches kinase names and protein sequences for multiple assay IDs in bulk using async calls."""
    kinase_map = {}
    async with aiohttp.ClientSession() as session:
        tasks = []

        for assay_id in assay_ids:
            tasks.append(fetch_assay_data(session, assay_id, kinase_map))

        # Run all tasks asynchronously
        await asyncio.gather(*tasks)

    return kinase_map


async def fetch_assay_data(session, assay_id, kinase_map):
    """Fetches the assay data, target data, and target component data for a given assay ID, including kinase name and protein sequence."""
    assay_url = f"{ASSAY_URL}{assay_id}.json"
    retries = 5  # Number of retries
    delay = 2  # Initial delay in seconds

    for attempt in range(retries):
        try:
            async with session.get(assay_url, timeout=30) as response:  # Set timeout to 30 seconds
                if response.status == 200:
                    assay_data = await response.json()
                    target_id = assay_data.get("target_chembl_id", None)

                    if target_id:
                        # Fetch the target data
                        target_url = f"{TARGET_URL}{target_id}.json"
                        async with session.get(target_url, timeout=30) as target_response:
                            if target_response.status == 200:
                                target_data = await target_response.json()

                                # Fetch the kinase name from target data (pref_name)
                                kinase_name = target_data.get("pref_name", "Unknown")

                                # Fetch the target component numeric ID
                                try:
                                    target_numeric_id = target_data['target_components'][0]['component_id']
                                except (KeyError, IndexError) as e:
                                    target_numeric_id = None
                                    print(f"❌ Error extracting target_numeric_id for assay {assay_id}: {e}")

                                if target_numeric_id:
                                    # Now query the target component endpoint for protein sequence
                                    target_component_url = f"{TARGET_COMPONENT_URL}{target_numeric_id}.json"
                                    async with session.get(target_component_url, timeout=30) as component_response:
                                        if component_response.status == 200:
                                            component_data = await component_response.json()

                                            # Extract the protein sequence from the component data
                                            protein_sequence = component_data.get('sequence', "No protein sequence found")

                                            # Store both kinase name and protein sequence in the map
                                            kinase_map[assay_id] = {
                                                "kinase_name": kinase_name,
                                                "protein_sequence": protein_sequence
                                            }
                                        else:
                                            kinase_map[assay_id] = {"kinase_name": kinase_name, "protein_sequence": "Target component data not found"}
                                else:
                                    print(assay_url)
                                    print(target_data)
                                    kinase_map[assay_id] = {"kinase_name": kinase_name, "protein_sequence": "No numeric target ID found"}
                            else:
                                kinase_map[assay_id] = {"kinase_name": kinase_name, "protein_sequence": "Target data not found"}
                    else:
                        kinase_map[assay_id] = {"kinase_name": "Unknown", "protein_sequence": "No target ID found"}
                else:
                    kinase_map[assay_id] = {"kinase_name": "Unknown", "protein_sequence": "Assay data not found"}

                return  # Exit on successful request
        except asyncio.TimeoutError:
            print(f"❌ Timeout occurred for assay ID {assay_id}. Retrying...")
        except aiohttp.ClientError as e:
            print(f"❌ Request failed for assay ID {assay_id}: {str(e)}. Retrying...")

        # Exponential backoff
        time.sleep(delay)
        delay *= 2  # Double the delay for each retry attempt


# Main function to fetch and process data
def fetch_data(std_type_list):
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

        filtered_df = df[df['standard_type'].isin(std_type_list)] # delete?

        all_data = pd.concat([all_data, filtered_df], ignore_index=True)

        most_common_assays = all_data['assay_chembl_id'].value_counts().head(2)
        a, b = most_common_assays[0], most_common_assays[1]

        print(most_common_assays)
        offset += LIMIT
        print(f"✅ Downloaded {len(all_data)} records...", end="\r")

    # Get kinase names and protein sequences in bulk asynchronously
    assay_ids = list(set(all_data["assay_chembl_id"]))
    kinase_map = asyncio.run(get_kinase_info_bulk(assay_ids))

    # Add kinase names and protein sequences to DataFrame
    all_data["kinase_name"] = all_data["assay_chembl_id"].map(lambda x: kinase_map.get(x, {}).get("kinase_name", "Unknown"))
    all_data["protein_sequence"] = all_data["assay_chembl_id"].map(lambda x: kinase_map.get(x, {}).get("protein_sequence", "Unknown"))

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
            kinase_name VARCHAR(255),
            protein_sequence TEXT
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
        INSERT INTO bioactivity (molecule_chembl_id, canonical_smiles, standard_type, standard_value, standard_units,
         kinase_name, protein_sequence)
        VALUES (%s, %s, %s, %s, %s, %s, %s);
        """
        data = df.values.tolist()
        cursor.executemany(query, data)  # Bulk Insert
        conn.commit()
        cursor.close()
        conn.close()
        print("\n✅ Data stored in MySQL database!")
    except Error as e:
        print(f"❌ Error inserting data into MySQL: {e}")




def protein_proc(df):
    possible_unaccepted_values = ["Assay data not found", "No target ID found", "Target data not found",
                                  "No numeric target ID found", "Target component data not found",
                                  "No protein sequence found"]
    filtered_array_prot = check_prot_seq(possible_unaccepted_values, np.array(list(set(df['protein_sequence']))))
    try:
        if filtered_array_prot:
            raise ValueError(f"Found invalid values in column 'protein_sequence':\n{filtered_array_prot}")
    except ValueError as e:
        print(f"Error: {e}")
        input("Press Enter to continue after debugging...")
        stop_point_for_debug=0
    df_clean_prot_seq = df[~df['protein_sequence'].isin(filtered_array_prot+possible_unaccepted_values)]

    try:
        if df_clean_prot_seq.empty:
            raise ValueError("Empty data table!")
    except ValueError as e:
        print(f"Error: {e}")
        input("Press Enter to continue after debugging...")
        stop_point_for_debug = 0
    table_to_store= df_clean_prot_seq[["molecule_chembl_id", "canonical_smiles", "standard_type", "standard_value",
                                       "standard_units", "kinase_name", "protein_sequence"]]

def data_full_proc(df):

    protein_proc(df)
    actval_proc()


def main():
    settings = load_settings()
    std_type_list = settings.get("std_type_list", [])

    add_protein_sequence_column()
    df = fetch_data(std_type_list)
    data_full_proc(df)
    store_data(table_to_store)


if __name__ == "__main__":
    main()
