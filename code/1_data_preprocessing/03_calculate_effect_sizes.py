"""
Script to calculate effect sizes and prepare data for meta-analysis.
This script calculates log rejection values, standard errors, and confidence intervals.
"""

import pandas as pd
import numpy as np
import os
from scipy import stats

# Set up directories
CLEAN_DIR = "data/processed/cleaned"
OUTPUT_DIR = "data/processed/meta_analysis"

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

def load_cleaned_data():
    """Load cleaned data from CSV files"""
    study_metadata = pd.read_csv(os.path.join(CLEAN_DIR, "study_metadata_clean.csv"))
    rejection_data = pd.read_csv(os.path.join(CLEAN_DIR, "rejection_data_clean.csv"))
    flux_data = pd.read_csv(os.path.join(CLEAN_DIR, "flux_decline_data_clean.csv"))
    
    return study_metadata, rejection_data, flux_data

def calculate_effect_sizes(rejection_data):
    """Calculate effect sizes (log rejection values) and confidence intervals"""
    # Create a new DataFrame for effect sizes
    effect_sizes = rejection_data.copy()
    
    # Calculate standard error based on std_dev and n_replicates
    effect_sizes["standard_error"] = effect_sizes["std_dev"] / np.sqrt(effect_sizes["n_replicates"])
    
    # Calculate 95% confidence intervals
    # For log rejection values, we use t-distribution for small sample sizes
    effect_sizes["ci_lower"] = effect_sizes.apply(
        lambda row: row["log_rejection"] - stats.t.ppf(0.975, row["n_replicates"]-1) * row["standard_error"] 
        if row["n_replicates"] > 1 else np.nan, axis=1
    )
    
    effect_sizes["ci_upper"] = effect_sizes.apply(
        lambda row: row["log_rejection"] + stats.t.ppf(0.975, row["n_replicates"]-1) * row["standard_error"]
        if row["n_replicates"] > 1 else np.nan, axis=1
    )
    
    # Calculate study weight (inverse variance)
    effect_sizes["weight"] = 1 / (effect_sizes["standard_error"]**2)
    
    # Normalize weights to sum to 100
    effect_sizes["weight_normalized"] = effect_sizes["weight"] / effect_sizes["weight"].sum() * 100
    
    return effect_sizes

def prepare_subgroup_data(effect_sizes, study_metadata):
    """Prepare data for subgroup analysis"""
    # Merge effect sizes with study metadata to get subgroup variables
    analysis_data = effect_sizes.merge(
        study_metadata[["study_id", "contaminant_class", "study_scale"]], 
        on="study_id", 
        how="left"
    )
    
    # Create additional subgroup variables
    # 1. Categorize PFAS by chain length
    analysis_data["pfas_chain_length"] = np.nan
    
    # Mapping for PFAS compounds to carbon chain length
    pfas_chain_map = {
        "TFA": 2,
        "PFBA": 4,
        "PFBS": 4,
        "PFHxA": 6,
        "PFHxS": 6,
        "GenX": 6,
        "PFOA": 8,
        "PFOS": 8,
        "6:2 FTS": 8,
        "PFNA": 9
    }
    
    # Apply mapping for PFAS compounds
    pfas_mask = analysis_data["compound_class"] == "PFAS"
    for compound, chain_length in pfas_chain_map.items():
        compound_mask = analysis_data["compound"].str.contains(compound, case=False, na=False)
        analysis_data.loc[pfas_mask & compound_mask, "pfas_chain_length"] = chain_length
    
    # Create chain length category
    analysis_data["chain_length_category"] = np.nan
    analysis_data.loc[analysis_data["pfas_chain_length"] <= 6, "chain_length_category"] = "Short-chain PFAS"
    analysis_data.loc[analysis_data["pfas_chain_length"] >= 7, "chain_length_category"] = "Long-chain PFAS"
    
    # 2. Create molecular weight categories
    mw_bins = [0, 200, 300, 400, 500, np.inf]
    mw_labels = ["<200 Da", "200-300 Da", "300-400 Da", "400-500 Da", ">500 Da"]
    analysis_data["mw_category"] = pd.cut(
        analysis_data["molecular_weight"], 
        bins=mw_bins, 
        labels=mw_labels,
        right=False
    )
    
    return analysis_data

def save_analysis_data(effect_sizes, subgroup_data):
    """Save analysis data to CSV files"""
    effect_sizes.to_csv(os.path.join(OUTPUT_DIR, "effect_sizes.csv"), index=False)
    subgroup_data.to_csv(os.path.join(OUTPUT_DIR, "subgroup_data.csv"), index=False)
    
    print("Analysis data saved to CSV files in", OUTPUT_DIR)

if __name__ == "__main__":
    print("Calculating effect sizes for meta-analysis...")
    study_metadata, rejection_data, flux_data = load_cleaned_data()
    
    effect_sizes = calculate_effect_sizes(rejection_data)
    subgroup_data = prepare_subgroup_data(effect_sizes, study_metadata)
    
    save_analysis_data(effect_sizes, subgroup_data)
    print("Effect size calculation complete.")
