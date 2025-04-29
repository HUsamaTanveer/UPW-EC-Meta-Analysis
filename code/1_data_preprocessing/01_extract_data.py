01_extract_data.py
"""
Script to extract data from the Meta-Analysis_Complete_Dataset.xlsx file.
This script reads the Excel file and converts sheets to pandas DataFrames for further analysis.
"""

import pandas as pd
import numpy as np
import os

# Set up directories
DATA_DIR = "data/excel"
OUTPUT_DIR = "data/processed"

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

def load_data():
    """Load data from Excel file into pandas DataFrames"""
    file_path = os.path.join(DATA_DIR, "Meta-Analysis_Complete_Dataset.xlsx")
    
    # Load each sheet into a separate DataFrame
    study_chars = pd.read_excel(file_path, sheet_name="Study Characteristics")
    rejection_data = pd.read_excel(file_path, sheet_name="Rejection Data")
    flux_data = pd.read_excel(file_path, sheet_name="Flux Decline Data")
    op_conditions = pd.read_excel(file_path, sheet_name="Operating Conditions")
    
    print(f"Loaded {len(study_chars)} studies")
    print(f"Loaded {len(rejection_data)} rejection data points")
    print(f"Loaded {len(flux_data)} flux decline data points")
    print(f"Loaded {len(op_conditions)} operating condition records")
    
    return study_chars, rejection_data, flux_data, op_conditions

def save_processed_data(study_chars, rejection_data, flux_data, op_conditions):
    """Save processed DataFrames to CSV files"""
    study_chars.to_csv(os.path.join(OUTPUT_DIR, "study_metadata.csv"), index=False)
    rejection_data.to_csv(os.path.join(OUTPUT_DIR, "rejection_data.csv"), index=False)
    flux_data.to_csv(os.path.join(OUTPUT_DIR, "flux_decline_data.csv"), index=False)
    op_conditions.to_csv(os.path.join(OUTPUT_DIR, "operating_conditions.csv"), index=False)
    
    print("Data saved to CSV files in", OUTPUT_DIR)

if __name__ == "__main__":
    print("Extracting data from Excel...")
    study_chars, rejection_data, flux_data, op_conditions = load_data()
    
    # Basic data validation
    assert not study_chars["study_id"].duplicated().any(), "Duplicate study IDs found"
    assert rejection_data["study_id"].isin(study_chars["study_id"]).all(), "Invalid study IDs in rejection data"
    assert flux_data["study_id"].isin(study_chars["study_id"]).all(), "Invalid study IDs in flux data"
    assert op_conditions["study_id"].isin(study_chars["study_id"]).all(), "Invalid study IDs in operating conditions"
    
    save_processed_data(study_chars, rejection_data, flux_data, op_conditions)
    print("Data extraction complete.")

