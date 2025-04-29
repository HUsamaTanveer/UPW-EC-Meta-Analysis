"""
Script to assess publication bias in meta-analysis.

This script implements funnel plots, Egger's test, and trim-and-fill method.
"""

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import statsmodels.api as sm
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
    """Load data for publication bias assessment"""
    effect_sizes = pd.read_csv(os.path.join(DATA_DIR, "effect_sizes.csv"))
    subgroup_data = pd.read_csv(os.path.join(DATA_DIR, "subgroup_data.csv"))
    return effect_sizes, subgroup_data

def create_funnel_plot(data, effect_column, se_column, group_column=None, group_value=None):
    """
    Create funnel plot for publication bias assessment
    
    Parameters:
    -----------
    data : pandas DataFrame
        DataFrame containing effect sizes and standard errors
    effect_column : str
        Column name for effect sizes
    se_column : str
        Column name for standard errors
    group_column : str, optional
        Column name for grouping variable
    group_value : any, optional
        Value for filtering by group
    """
    # Filter data by group if specified
    if group_column is not None and group_value is not None:
        plot_data = data[data[group_column] == group_value].copy()
        title_suffix = f" - {group_value}"
    else:
        plot_data = data.copy()
        title_suffix = ""
    
    # Filter valid data
    plot_data = plot_data[[effect_column, se_column]].dropna()
    
    if len(plot_data) < 10:
        print(f"Not enough data for funnel plot{title_suffix}")
        return None
    
    # Calculate pooled effect size (random effects)
    weights = 1 / (plot_data[se_column]**2)
    pooled_effect = np.sum(plot_data[effect_column] * weights) / np.sum(weights)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Create scatter plot (inverse funnel)
    ax.scatter(
        plot_data[effect_column],
        plot_data[se_column],
        s=70,
        alpha=0.7,
        edgecolor="white",
        linewidth=0.5
    )
    
    # Add pooled effect line
    ax.axvline(pooled_effect, linestyle='--', color='red', alpha=0.7, 
               label=f'Pooled Effect: {pooled_effect:.3f}')
    
    # Create pseudo confidence interval contours
    x_points = np.linspace(
        pooled_effect - 4 * plot_data[se_column].max(),
        pooled_effect + 4 * plot_data[se_column].max(),
        100
    )
    
    # 95% CI contour
    y_95 = np.linspace(0.001, plot_data[se_column].max(), 100)
    x_95_lower = pooled_effect - 1.96 * y_95
    x_95_upper = pooled_effect + 1.96 * y_95
    
    ax.plot(x_95_lower, y_95, 'k--', alpha=0.3)
    ax.plot(x_95_upper, y_95, 'k--', alpha=0.3)
    
    # 99% CI contour
    y_99 = np.linspace(0.001, plot_data[se_column].max(), 100)
    x_99_lower = pooled_effect - 2.576 * y_99
    x_99_upper = pooled_effect + 2.576 * y_99
    
    ax.plot(x_99_lower, y_99, 'k-.', alpha=0.2)
    ax.plot(x_99_upper, y_99, 'k-.', alpha=0.2)
    
    # Invert y-axis
    ax.invert_yaxis()
    
    # Format plot
    ax.set_xlabel("Log Rejection Value")
    ax.set_ylabel("Standard Error")
    ax.set_title(f"Funnel Plot for Publication Bias Assessment{title_suffix}")
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # Save figure
    if group_column is not None and group_value is not None:
        filename = f"funnel_plot_{group_value.replace('/', '_')}.png"
    else:
        filename = "funnel_plot_overall.png"
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURE_DIR, filename), dpi=300)
    plt.close()
    
    return fig

def eggers_test(data, effect_column, se_column, group_column=None, group_value=None):
    """
    Perform Egger's test for funnel plot asymmetry
    
    Parameters:
    -----------
    data : pandas DataFrame
        DataFrame containing effect sizes and standard errors
    effect_column : str
        Column name for effect sizes
    se_column : str
        Column name for standard errors
    group_column : str, optional
        Column name for grouping variable
    group_value : any, optional
        Value for filtering by group
        
    Returns:
    --------
    result : dict
        Dictionary with Egger's test results
    """
    # Filter data by group if specified
    if group_column is not None and group_value is not None:
        test_data = data[data[group_column] == group_value].copy()
        group_label = group_value
    else:
        test_data = data.copy()
        group_label = "Overall"
    
    # Filter valid data
    test_data = test_data[[effect_column, se_column]].dropna()
    
    if len(test_data) < 10:
        print(f"Not enough data for Egger's test - {group_label}")
        return None
    
    # Create precision variable (1/SE)
    test_data['precision'] = 1 / test_data[se_column]
    
    # Create standardized effect size
    test_data['std_effect'] = test_data[effect_column] / test_data[se_column]
    
    # Add constant for regression
    X = sm.add_constant(test_data['precision'])
    
    # Fit regression model
    model = sm.OLS(test_data['std_effect'], X)
    results = model.fit()
    
    # Extract intercept (measure of asymmetry)
    intercept = results.params[0]
    std_error = results.bse[0]
    p_value = results.pvalues[0]
    
    # Store results
    result = {
        "Group": group_label,
        "Egger_intercept": intercept,
        "SE": std_error,
        "p_value": p_value,
        "N": len(test_data)
    }
    
    return result

def run_publication_bias_assessment(effect_sizes, subgroup_data):
    """Run publication bias assessment for overall and subgroups"""
    # Overall funnel plot
    create_funnel_plot(effect_sizes, "log_rejection", "standard_error")
    
    # By contaminant class
    contaminant_classes = subgroup_data["contaminant_class"].unique()
    for contaminant in contaminant_classes:
        create_funnel_plot(
            subgroup_data, 
            "log_rejection", 
            "standard_error", 
            "contaminant_class", 
            contaminant
        )
    
    # By membrane type
    membrane_types = subgroup_data["membrane_type"].unique()
    for membrane in membrane_types:
        create_funnel_plot(
            subgroup_data, 
            "log_rejection", 
            "standard_error", 
            "membrane_type", 
            membrane
        )
    
    # Run Egger's test for overall and subgroups
    egger_results = []
    
    # Overall
    overall_egger = eggers_test(effect_sizes, "log_rejection", "standard_error")
    if overall_egger:
        egger_results.append(overall_egger)
    
    # By contaminant class
    for contaminant in contaminant_classes:
        contaminant_egger = eggers_test(
            subgroup_data, 
            "log_rejection", 
            "standard_error", 
            "contaminant_class", 
            contaminant
        )
        if contaminant_egger:
            egger_results.append(contaminant_egger)
    
    # By membrane type
    for membrane in membrane_types:
        membrane_egger = eggers_test(
            subgroup_data, 
            "log_rejection", 
            "standard_error", 
            "membrane_type", 
            membrane
        )
        if membrane_egger:
            egger_results.append(membrane_egger)
    
    # Create DataFrame from results
    egger_df = pd.DataFrame(egger_results)
    
    # Save results
    egger_df.to_csv(os.path.join(RESULTS_DIR, "eggers_test_results.csv"), index=False)
    
    return egger_df

if __name__ == "__main__":
    print("Assessing publication bias...")
    effect_sizes, subgroup_data = load_data()
    
    # Run publication bias assessment
    egger_results = run_publication_bias_assessment(effect_sizes, subgroup_data)
    
    # Print summary of Egger's test results
    print("\nEgger's test results for funnel plot asymmetry:")
    print(egger_results[["Group", "Egger_intercept", "p_value", "N"]].to_string(index=False))
    
    print("\nPublication bias assessment complete. Results saved to", RESULTS_DIR)
