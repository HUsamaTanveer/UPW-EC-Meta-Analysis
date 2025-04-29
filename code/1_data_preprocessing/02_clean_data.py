02_clean_data.py
"""
Script to clean and preprocess data for meta-analysis.
This script handles missing values, outliers, and data transformations.
"""

import pandas as pd
import numpy as np
import os
from scipy import stats

# Set up directories
PROCESSED_DIR = "data/processed"
CLEAN_DIR = "data/processed/cleaned"

if not os.path.exists(CLEAN_DIR):
    os.makedirs(CLEAN_DIR)

def load_processed_data():
    """Load processed data from CSV files"""
    study_metadata = pd.read_csv(os.path.join(PROCESSED_DIR, "study_metadata.csv"))
    rejection_data = pd.read_csv(os.path.join(PROCESSED_DIR, "rejection_data.csv"))
    flux_data = pd.read_csv(os.path.join(PROCESSED_DIR, "flux_decline_data.csv"))
    op_conditions = pd.read_csv(os.path.join(PROCESSED_DIR, "operating_conditions.csv"))
    
    return study_metadata, rejection_data, flux_data, op_conditions

def clean_study_metadata(df):
    """Clean study metadata DataFrame"""
    # Convert publication year to integer
    df["year"] = df["year"].astype(int)
    
    # Standardize contaminant class names
    contaminant_map = {
        "PFAS": "PFAS",
        "PFASs": "PFAS",
        "Per- and polyfluoroalkyl substances": "PFAS",
        "Pharmaceuticals": "Pharmaceuticals",
        "Pharmaceutical": "Pharmaceuticals",
        "PPCP": "Pharmaceuticals",
        "PPCPs": "Pharmaceuticals",
        "CMP Organics": "CMP Organics",
        "Chemical mechanical planarization": "CMP Organics",
        "Lithography organics": "CMP Organics",
        "Micro/nanoplastics": "Micro/nanoplastics",
        "Microplastics": "Micro/nanoplastics",
        "Nanoplastics": "Micro/nanoplastics"
    }
    df["contaminant_class"] = df["contaminant_class"].replace(contaminant_map)
    
    # Standardize study scale
    scale_map = {
        "Laboratory": "Laboratory",
        "Lab": "Laboratory",
        "lab-scale": "Laboratory",
        "Pilot": "Pilot",
        "pilot-scale": "Pilot",
        "Full-scale": "Full-scale",
        "full scale": "Full-scale",
        "Full scale": "Full-scale",
        "Review": "Review",
        "Meta-review": "Review"
    }
    df["study_scale"] = df["study_scale"].replace(scale_map)
    
    return df

def clean_rejection_data(df):
    """Clean rejection data DataFrame"""
    # Calculate log rejection if missing
    mask = df["log_rejection"].isna()
    if mask.any():
        df.loc[mask, "log_rejection"] = np.log10(1 / (1 - df.loc[mask, "rejection_percent"] / 100))
    
    # Handle outliers (values beyond 3 standard deviations)
    z_scores = np.abs(stats.zscore(df["log_rejection"].dropna()))
    outliers = z_scores > 3
    
    if outliers.any():
        print(f"Found {outliers.sum()} outliers in log_rejection values")
        # Flag outliers but keep them in the dataset
        df["is_outlier"] = False
        df.loc[df["log_rejection"].notna(), "is_outlier"] = outliers
    
    # Ensure numeric types
    numeric_cols = ["molecular_weight", "log_kow", "pressure", "ph", "rejection_percent", 
                   "log_rejection", "std_dev", "n_replicates"]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    
    return df

def clean_flux_data(df):
    """Clean flux decline data DataFrame"""
    # Calculate flux decline percent if missing
    mask = df["flux_decline_percent"].isna()
    if mask.any():
        df.loc[mask, "flux_decline_percent"] = (
            (df.loc[mask, "initial_flux"] - df.loc[mask, "final_flux"]) / 
            df.loc[mask, "initial_flux"] * 100
        )
    
    # Convert operating_time to numeric hours
    if "operating_time" in df.columns:
        df["operating_time"] = pd.to_numeric(df["operating_time"], errors='coerce')
    
    # Ensure percentage columns are within valid range (0-100%)
    for col in ["flux_decline_percent", "recovery_physical", "recovery_chemical"]:
        if col in df.columns:
            df[col] = df[col].clip(0, 100)
    
    return df

def save_cleaned_data(study_metadata, rejection_data, flux_data, op_conditions):
    """Save cleaned DataFrames to CSV files"""
    study_metadata.to_csv(os.path.join(CLEAN_DIR, "study_metadata_clean.csv"), index=False)
    rejection_data.to_csv(os.path.join(CLEAN_DIR, "rejection_data_clean.csv"), index=False)
    flux_data.to_csv(os.path.join(CLEAN_DIR, "flux_decline_data_clean.csv"), index=False)
    op_conditions.to_csv(os.path.join(CLEAN_DIR, "operating_conditions_clean.csv"), index=False)
    
    print("Cleaned data saved to CSV files in", CLEAN_DIR)

if __name__ == "__main__":
    print("Cleaning data for meta-analysis...")
    study_metadata, rejection_data, flux_data, op_conditions = load_processed_data()
    
    # Apply cleaning functions
    study_metadata = clean_study_metadata(study_metadata)
    rejection_data = clean_rejection_data(rejection_data)
    flux_data = clean_flux_data(flux_data)
    
    # Basic data validation after cleaning
    assert not study_metadata["study_id"].duplicated().any(), "Duplicate study IDs found"
    assert not rejection_data["rejection_id"].duplicated().any(), "Duplicate rejection IDs found"
    
    save_cleaned_data(study_metadata, rejection_data, flux_data, op_conditions)
    print("Data cleaning complete.")
