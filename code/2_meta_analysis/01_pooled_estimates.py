01_pooled_estimates.py
"""
Script to calculate pooled estimates for meta-analysis.
This script implements random-effects meta-analysis using the DerSimonian and Laird method.
"""

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy import stats

# Set up directories
DATA_DIR = "data/processed/meta_analysis"
RESULTS_DIR = "results/meta_analysis"
FIGURE_DIR = "figures/meta_analysis"

for dir_path in [RESULTS_DIR, FIGURE_DIR]:
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

def load_analysis_data():
    """Load data prepared for meta-analysis"""
    effect_sizes = pd.read_csv(os.path.join(DATA_DIR, "effect_sizes.csv"))
    subgroup_data = pd.read_csv(os.path.join(DATA_DIR, "subgroup_data.csv"))
    
    return effect_sizes, subgroup_data

def random_effects_meta_analysis(data, effect_column, se_column, group_column=None):
    """
    Perform random-effects meta-analysis using DerSimonian and Laird method
    
    Parameters:
    -----------
    data : pandas DataFrame
        DataFrame containing effect sizes and standard errors
    effect_column : str
        Column name for effect sizes
    se_column : str
        Column name for standard errors
    group_column : str, optional
        Column name for grouping variable (for subgroup analysis)
        
    Returns:
    --------
    results : pandas DataFrame
        DataFrame with pooled estimates, confidence intervals, and heterogeneity statistics
    """
    if group_column is None:
        # Overall meta-analysis
        groups = [("Overall", data)]
    else:
        # Subgroup analysis
        groups = [(group, data[data[group_column] == group]) for group in data[group_column].unique()]
    
    results = []
    
    for group_name, group_data in groups:
        # Filter out missing values
        valid_data = group_data[[effect_column, se_column]].dropna()
        
        if len(valid_data) < 2:
            # Skip groups with insufficient data
            print(f"Skipping group {group_name} with insufficient data")
            continue
        
        # Calculate fixed-effects weights (inverse variance)
        weights = 1 / (valid_data[se_column]**2)
        
        # Calculate fixed-effects pooled estimate
        pooled_estimate_fe = np.sum(valid_data[effect_column] * weights) / np.sum(weights)
        
        # Calculate Q statistic (measure of heterogeneity)
        Q = np.sum(weights * (valid_data[effect_column] - pooled_estimate_fe)**2)
        
        # Degrees of freedom
        df = len(valid_data) - 1
        
        # Calculate between-study variance (tau²) using DerSimonian and Laird method
        tau_squared = max(0, (Q - df) / (np.sum(weights) - np.sum(weights**2) / np.sum(weights)))
        
        # Calculate I² statistic (percentage of variation due to heterogeneity)
        I_squared = max(0, 100 * (Q - df) / Q) if Q > 0 else 0
        
        # Calculate random-effects weights
        weights_re = 1 / (valid_data[se_column]**2 + tau_squared)
        
        # Calculate random-effects pooled estimate
        pooled_estimate_re = np.sum(valid_data[effect_column] * weights_re) / np.sum(weights_re)
        
        # Calculate standard error of pooled estimate
        se_pooled = np.sqrt(1 / np.sum(weights_re))
        
        # Calculate 95% confidence interval
        ci_lower = pooled_estimate_re - 1.96 * se_pooled
        ci_upper = pooled_estimate_re + 1.96 * se_pooled
        
        # Calculate p-value for pooled estimate
        z_score = pooled_estimate_re / se_pooled
        p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))
        
        # Calculate p-value for heterogeneity
        p_heterogeneity = 1 - stats.chi2.cdf(Q, df)
        
        # Convert log rejection to rejection percentage
        rejection_percent = (1 - 10**(-pooled_estimate_re)) * 100
        rejection_percent_lower = (1 - 10**(-ci_lower)) * 100
        rejection_percent_upper = (1 - 10**(-ci_upper)) * 100
        
        # Store results
        results.append({
            "Group": group_name,
            "N_studies": len(valid_data),
            "Pooled_estimate": pooled_estimate_re,
            "CI_lower": ci_lower,
            "CI_upper": ci_upper,
            "SE": se_pooled,
            "p_value": p_value,
            "Rejection_percent": rejection_percent,
            "Rejection_percent_lower": rejection_percent_lower,
            "Rejection_percent_upper": rejection_percent_upper,
            "Q": Q,
            "df": df,
            "p_heterogeneity": p_heterogeneity,
            "I_squared": I_squared,
            "tau_squared": tau_squared
        })
    
    return pd.DataFrame(results)

def run_analyses(subgroup_data):
    """Run overall and subgroup meta-analyses"""
    # Overall meta-analysis
    overall_results = random_effects_meta_analysis(
        subgroup_data, 
        effect_column="log_rejection", 
        se_column="standard_error"
    )
    
    # Subgroup analyses
    # 1. By membrane type
    membrane_results = random_effects_meta_analysis(
        subgroup_data, 
        effect_column="log_rejection", 
        se_column="standard_error",
        group_column="membrane_type"
    )
    
    # 2. By contaminant class
    contaminant_results = random_effects_meta_analysis(
        subgroup_data, 
        effect_column="log_rejection", 
        se_column="standard_error",
        group_column="contaminant_class"
    )
    
    # 3. By PFAS chain length category (for PFAS compounds only)
    pfas_data = subgroup_data[subgroup_data["contaminant_class"] == "PFAS"]
    chain_length_results = random_effects_meta_analysis(
        pfas_data, 
        effect_column="log_rejection", 
        se_column="standard_error",
        group_column="chain_length_category"
    )
    
    # 4. By study scale
    scale_results = random_effects_meta_analysis(
        subgroup_data, 
        effect_column="log_rejection", 
        se_column="standard_error",
        group_column="study_scale"
    )
    
    # Save results
    overall_results.to_csv(os.path.join(RESULTS_DIR, "overall_meta_analysis.csv"), index=False)
    membrane_results.to_csv(os.path.join(RESULTS_DIR, "membrane_subgroup_analysis.csv"), index=False)
    contaminant_results.to_csv(os.path.join(RESULTS_DIR, "contaminant_subgroup_analysis.csv"), index=False)
    chain_length_results.to_csv(os.path.join(RESULTS_DIR, "chain_length_subgroup_analysis.csv"), index=False)
    scale_results.to_csv(os.path.join(RESULTS_DIR, "scale_subgroup_analysis.csv"), index=False)
    
    return {
        "overall": overall_results,
        "membrane": membrane_results,
        "contaminant": contaminant_results,
        "chain_length": chain_length_results,
        "scale": scale_results
    }

if __name__ == "__main__":
    print("Running meta-analysis...")
    effect_sizes, subgroup_data = load_analysis_data()
    
    results = run_analyses(subgroup_data)
    
    # Print summary of results
    for analysis_name, result_df in results.items():
        print(f"\n{analysis_name.capitalize()} meta-analysis results:")
        print(result_df[["Group", "N_studies", "Pooled_estimate", "CI_lower", "CI_upper", 
                          "Rejection_percent", "I_squared"]].to_string(index=False))
    
    print("\nMeta-analysis complete. Results saved to", RESULTS_DIR)

