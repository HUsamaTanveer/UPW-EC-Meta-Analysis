02_heterogeneity.py
"""
#Script to assess heterogeneity in meta-analysis results.
#This script calculates and visualizes heterogeneity statistics (Q, I², tau²).
"""

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Set up directories
DATA_DIR = "data/processed/meta_analysis"
RESULTS_DIR = "results/meta_analysis"
FIGURE_DIR = "figures/meta_analysis"

# Set plot styles
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("colorblind")

def load_data():
    """Load meta-analysis data and results"""
    subgroup_data = pd.read_csv(os.path.join(DATA_DIR, "subgroup_data.csv"))
    
    # Load meta-analysis results
    membrane_results = pd.read_csv(os.path.join(RESULTS_DIR, "membrane_subgroup_analysis.csv"))
    contaminant_results = pd.read_csv(os.path.join(RESULTS_DIR, "contaminant_subgroup_analysis.csv"))
    
    return subgroup_data, membrane_results, contaminant_results

def plot_heterogeneity(membrane_results, contaminant_results):
    """Create heterogeneity plots"""
    # Plot I-squared values by subgroup
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    
    # Plot for membrane types
    membrane_results = membrane_results.sort_values("I_squared", ascending=False)
    ax1.barh(membrane_results["Group"], membrane_results["I_squared"], color="skyblue")
    ax1.set_xlim(0, 100)
    ax1.set_xlabel("I² (%)")
    ax1.set_ylabel("Membrane Type")
    ax1.set_title("Heterogeneity by Membrane Type")
    
    # Add reference lines for I-squared interpretation
    ax1.axvline(25, color='green', linestyle='--', alpha=0.7, label='Low (25%)')
    ax1.axvline(50, color='orange', linestyle='--', alpha=0.7, label='Moderate (50%)')
    ax1.axvline(75, color='red', linestyle='--', alpha=0.7, label='High (75%)')
    ax1.legend()
    
    # Plot for contaminant classes
    contaminant_results = contaminant_results.sort_values("I_squared", ascending=False)
    ax2.barh(contaminant_results["Group"], contaminant_results["I_squared"], color="lightcoral")
    ax2.set_xlim(0, 100)
    ax2.set_xlabel("I² (%)")
    ax2.set_ylabel("Contaminant Class")
    ax2.set_title("Heterogeneity by Contaminant Class")
    
    # Add reference lines
    ax2.axvline(25, color='green', linestyle='--', alpha=0.7, label='Low (25%)')
    ax2.axvline(50, color='orange', linestyle='--', alpha=0.7, label='Moderate (50%)')
    ax2.axvline(75, color='red', linestyle='--', alpha=0.7, label='High (75%)')
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURE_DIR, "heterogeneity_by_subgroup.png"), dpi=300)
    plt.close()

def analyze_moderators(subgroup_data):
    """Analyze potential moderators of heterogeneity"""
    # Identify potential moderators
    moderators = [
        "molecular_weight", 
        "log_kow", 
        "pressure", 
        "ph", 
        "pfas_chain_length"
    ]
    
    results = []
    
    # Analyze relationship between each moderator and effect size
    for moderator in moderators:
        # Skip moderators with too many missing values
        if subgroup_data[moderator].isna().sum() > len(subgroup_data) * 0.5:
            continue
        
        # Filter valid data
        valid_data = subgroup_data[[moderator, "log_rejection"]].dropna()
        
        if len(valid_data) < 10:
            continue
        
        # Calculate correlation
        corr, p_value = stats.pearsonr(valid_data[moderator], valid_data["log_rejection"])
        
        # Store results
        results.append({
            "Moderator": moderator,
            "Correlation": corr,
            "p_value": p_value,
            "N": len(valid_data)
        })
    
    # Create DataFrame of results
    moderator_df = pd.DataFrame(results)
    
    # Save results
    moderator_df.to_csv(os.path.join(RESULTS_DIR, "moderator_analysis.csv"), index=False)
    
    return moderator_df

def plot_moderator_effects(subgroup_data, moderator_df):
    """Plot effects of significant moderators"""
    # Plot significant moderators (p < 0.05)
    significant_moderators = moderator_df[moderator_df["p_value"] < 0.05]["Moderator"].tolist()
    
    for moderator in significant_moderators:
        # Filter valid data
        valid_data = subgroup_data[[moderator, "log_rejection", "membrane_type"]].dropna()
        
        # Create scatter plot
        plt.figure(figsize=(8, 6))
        
        # Color points by membrane type
        membrane_types = valid_data["membrane_type"].unique()
        colors = sns.color_palette("colorblind", n_colors=len(membrane_types))
        
        for i, membrane in enumerate(membrane_types):
            mask = valid_data["membrane_type"] == membrane
            plt.scatter(
                valid_data.loc[mask, moderator], 
                valid_data.loc[mask, "log_rejection"],
                color=colors[i],
                alpha=0.7,
                label=membrane
            )
        
        # Add best fit line
        x = valid_data[moderator]
        y = valid_data["log_rejection"]
        
        # Calculate regression line
        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
        x_line = np.linspace(x.min(), x.max(), 100)
        y_line = slope * x_line + intercept
        
        plt.plot(x_line, y_line, 'k--', alpha=0.7)
        
        # Add text with correlation information
        plt.text(
            0.05, 0.95, 
            f"r = {r_value:.2f}, p = {p_value:.4f}\ny = {slope:.3f}x + {intercept:.3f}",
            transform=plt.gca().transAxes,
            bbox=dict(facecolor='white', alpha=0.8)
        )
        
        # Format plot
        plt.xlabel(moderator.replace("_", " ").title())
        plt.ylabel("Log Rejection Value")
        plt.title(f"Effect of {moderator.replace('_', ' ').title()} on Rejection Performance")
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(os.path.join(FIGURE_DIR, f"moderator_{moderator}.png"), dpi=300)
        plt.close()

if __name__ == "__main__":
    print("Analyzing heterogeneity in meta-analysis results...")
    subgroup_data, membrane_results, contaminant_results = load_data()
    
    plot_heter
if __name__ == "__main__":
    print("Analyzing heterogeneity in meta-analysis results...")
    subgroup_data, membrane_results, contaminant_results = load_data()
    
    # Plot heterogeneity statistics
    plot_heterogeneity(membrane_results, contaminant_results)
    
    # Analyze moderators
    moderator_df = analyze_moderators(subgroup_data)
    print("Significant moderators of heterogeneity:")
    print(moderator_df[moderator_df["p_value"] < 0.05].to_string(index=False))
    
    # Plot significant moderators
    plot_moderator_effects(subgroup_data, moderator_df)
    
    print("Heterogeneity analysis complete. Results saved to", RESULTS_DIR)

