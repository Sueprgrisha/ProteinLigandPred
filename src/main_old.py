import json
import os
from src.db_utils import add_protein_sequence_column
import pandas as pd
import requests
import mysql.connector
import asyncio
import aiohttp
import time
from mysql.connector import Error
from data_cleaning.data_cleaning import *
import json


def process_activity_url(settings_to_process):
    # Get the limit and offset from the data
    activity_url_raw = settings_to_process.get("activity_url", [])
    limit = settings_to_process.get("limit", [])
    offset = settings_to_process.get("offset", [])
    # Replace the placeholders with the actual values
    activity_url = activity_url_raw.format(limit=limit, offset=offset)

    # Now you can use the activity_url with the actual values
    print(activity_url)
    return activity_url

def load_settings():
    settings_path = os.path.join(os.path.dirname(__file__), "..", "config", "settings.json")
    activity_url = process_activity_url(settings_path)

    settings_path = os.path.join(os.path.dirname(__file__), "..", "config", "processed_settings.json")
    with open(settings_path, "r") as f:
        return json.load(f)

def setup_database(settings):
    """Setup database if needed based on settings."""
    if settings.get("enable_protein_sequence_column", False):
        add_protein_sequence_column()


def fetch_data(settings):
    """Fetches bioactivity data from ChEMBL and cleans the data before saving it"""

    # Get the relevant settings
    fetch_columns_names = settings.get("fetch_columns", [])
    max_records = settings.get("max_records", 1000)
    min_ligand_per_kinase = settings.get("min_ligand_per_kinase", 2)
    first_url = settings.get("activity_url", "URL not found") #qw



    # Initialise dataframe
    all_data = pd.DataFrame(columns=fetch_columns_names)

    a, b = 0, 0 # check


    # Get at least MAX_RECORDS amount of points of data, and make sure at least 2 proteins have "MIN_LIGAND_PER_KINASE" amount of ligand pairing examples
    # while not (len(all_data) >= max_records and a >= min_ligand_per_kinase and b >= min_ligand_per_kinase): # check

    while not (len(all_data) >= max_records):  # check

        first_url =


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


def main():
    main_pipeline
    process_setting
    db_setup_settings = load_settings()

    setup_database(db_setup_settings)  # Load database changes before main execution

    # std_type_list = settings.get("std_type_list", []) - to put inside fetch data
    df = fetch_data(db_setup_settings)

    data_full_proc(df)

    store_data(table_to_store)

if __name__ == "__main__":
    main()
