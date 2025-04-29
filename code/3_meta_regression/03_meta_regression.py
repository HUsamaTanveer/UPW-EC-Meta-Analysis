"""
Script to perform meta-regression analysis.

This script implements meta-regression to examine the influence of 
study-level covariates on effect sizes.
"""

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import statsmodels.api as sm
import statsmodels.formula.api as smf
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
    """Load data for meta-regression"""
    subgroup_data = pd.read_csv(os.path.join(DATA_DIR, "subgroup_data.csv"))
    return subgroup_data

def simple_meta_regression(data, outcome, predictor, weights="weight"):
    """
    Perform simple meta-regression with one predictor
    
    Parameters:
    -----------
    data : pandas DataFrame
        DataFrame containing outcome, predictor, and weights
    outcome : str
        Column name for outcome variable
    predictor : str
        Column name for predictor variable
    weights : str
        Column name for weights (default: "weight")
        
    Returns:
    --------
    results : statsmodels RegressionResults
        Regression results object
    """
    # Filter valid data
    valid_data = data[[outcome, predictor, weights]].dropna()
    
    if len(valid_data) < 10:
        print(f"Not enough data for meta-regression with {predictor}")
        return None
    
    # Create regression formula
    formula = f"{outcome} ~ {predictor}"
    
    # Fit weighted least squares model
    model = smf.wls(formula=formula, data=valid_data, weights=valid_data[weights])
    results = model.fit()
    
    return results

def multiple_meta_regression(data, outcome, predictors, weights="weight"):
    """
    Perform multiple meta-regression with multiple predictors
    
    Parameters:
    -----------
    data : pandas DataFrame
        DataFrame containing outcome, predictors, and weights
    outcome : str
        Column name for outcome variable
    predictors : list
        List of column names for predictor variables
    weights : str
        Column name for weights (default: "weight")
        
    Returns:
    --------
    results : statsmodels RegressionResults
        Regression results object
    """
    # Filter valid data
    valid_data = data[[outcome] + predictors + [weights]].dropna()
    
    if len(valid_data) < len(predictors) + 10:
        print(f"Not enough data for meta-regression with {len(predictors)} predictors")
        return None
    
    # Create regression formula
    formula = f"{outcome} ~ " + " + ".join(predictors)
    
    # Fit weighted least squares model
    model = smf.wls(formula=formula, data=valid_data, weights=valid_data[weights])
    results = model.fit()
    
    return results

def run_pfas_regression(data):
    """Run meta-regression for PFAS compounds"""
    # Filter PFAS data
    pfas_data = data[data["contaminant_class"] == "PFAS"].copy()
    
    # Skip if not enough data
    if len(pfas_data) < 10:
        print("Not enough PFAS data for meta-regression")
        return None
    
    # Define predictors
    predictors = [
        "pfas_chain_length",
        "membrane_salt_rejection",
        "membrane_mwco",
        "pressure",
        "feed_concentration",
        "ph",
        "ionic_strength",
        "temperature",
        "nom_presence"
    ]
    
    # Initialize results dictionary
    results_dict = {}
    
    # Run simple meta-regression for each predictor
    for predictor in predictors:
        if pfas_data[predictor].notna().sum() < 10:
            continue
            
        result = simple_meta_regression(pfas_data, "log_rejection", predictor)
        
        if result is not None:
            results_dict[predictor] = {
                "coefficient": result.params[predictor],
                "std_error": result.bse[predictor],
                "p_value": result.pvalues[predictor],
                "r_squared": result.rsquared,
                "n": len(result.fittedvalues)
            }
    
    # Create DataFrame from results
    results_df = pd.DataFrame.from_dict(results_dict, orient="index")
    results_df.index.name = "Predictor"
    results_df.reset_index(inplace=True)
    
    # Save results
    results_df.to_csv(os.path.join(RESULTS_DIR, "pfas_meta_regression.csv"), index=False)
    
    # Try multiple meta-regression with significant predictors (p < 0.05)
    significant_predictors = results_df[results_df["p_value"] < 0.05]["Predictor"].tolist()
    
    if len(significant_predictors) >= 2:
        multi_result = multiple_meta_regression(pfas_data, "log_rejection", significant_predictors)
        
        if multi_result is not None:
            # Save multiple regression results
            with open(os.path.join(RESULTS_DIR, "pfas_multiple_regression.txt"), "w") as f:
                f.write(multi_result.summary().as_text())
    
    return results_df

def run_pharmaceutical_regression(data):
    """Run meta-regression for pharmaceutical compounds"""
    # Filter pharmaceutical data
    pharma_data = data[data["contaminant_class"] == "Pharmaceuticals"].copy()
    
    # Skip if not enough data
    if len(pharma_data) < 10:
        print("Not enough pharmaceutical data for meta-regression")
        return None
    
    # Define predictors
    predictors = [
        "molecular_weight",
        "log_kow",
        "membrane_salt_rejection",
        "membrane_mwco",
        "pressure",
        "feed_concentration",
        "ph",
        "ionic_strength",
        "temperature",
        "nom_presence"
    ]
    
    # Initialize results dictionary
    results_dict = {}
    
    # Run simple meta-regression for each predictor
    for predictor in predictors:
        if pharma_data[predictor].notna().sum() < 10:
            continue
            
        result = simple_meta_regression(pharma_data, "log_rejection", predictor)
        
        if result is not None:
            results_dict[predictor] = {
                "coefficient": result.params[predictor],
                "std_error": result.bse[predictor],
                "p_value": result.pvalues[predictor],
                "r_squared": result.rsquared,
                "n": len(result.fittedvalues)
            }
    
    # Create DataFrame from results
    results_df = pd.DataFrame.from_dict(results_dict, orient="index")
    results_df.index.name = "Predictor"
    results_df.reset_index(inplace=True)
    
    # Save results
    results_df.to_csv(os.path.join(RESULTS_DIR, "pharma_meta_regression.csv"), index=False)
    
    return results_df

def run_cmp_regression(data):
    """Run meta-regression for CMP and lithography organics"""
    # Filter CMP data
    cmp_data = data[data["contaminant_class"] == "CMP Organics"].copy()
    
    # Skip if not enough data
    if len(cmp_data) < 10:
        print("Not enough CMP data for meta-regression")
        return None
    
    # Define predictors
    predictors = [
        "molecular_weight",
        "log_kow",
        "membrane_salt_rejection",
        "membrane_mwco",
        "pressure",
        "feed_concentration",
        "ph",
        "ionic_strength",
        "temperature",
        "nom_presence"
    ]
    
    # Initialize results dictionary
    results_dict = {}
    
    # Run simple meta-regression for each predictor
    for predictor in predictors:
        if cmp_data[predictor].notna().sum() < 10:
            continue
            
        result = simple_meta_regression(cmp_data, "log_rejection", predictor)
        
        if result is not None:
            results_dict[predictor] = {
                "coefficient": result.params[predictor],
                "std_error": result.bse[predictor],
                "p_value": result.pvalues[predictor],
                "r_squared": result.rsquared,
                "n": len(result.fittedvalues)
            }
    
    # Create DataFrame from results
    results_df = pd.DataFrame.from_dict(results_dict, orient="index")
    results_df.index.name = "Predictor"
    results_df.reset_index(inplace=True)
    
    # Save results
    results_df.to_csv(os.path.join(RESULTS_DIR, "cmp_meta_regression.csv"), index=False)
    
    return results_df

def run_nanoplastics_regression(data):
    """Run meta-regression for micro/nanoplastics"""
    # Filter nanoplastics data
    np_data = data[data["contaminant_class"] == "Micro/nanoplastics"].copy()
    
    # Skip if not enough data
    if len(np_data) < 10:
        print("Not enough nanoplastics data for meta-regression")
        return None
    
    # Define predictors
    predictors = [
        "particle_size",
        "particle_charge",
        "membrane_pore_size",
        "pressure",
        "feed_concentration",
        "ph",
        "ionic_strength",
        "temperature",
        "nom_presence"
    ]
    
    # Initialize results dictionary
    results_dict = {}
    
    # Run simple meta-regression for each predictor
    for predictor in predictors:
        if np_data[predictor].notna().sum() < 10:
            continue
            
        result = simple_meta_regression(np_data, "log_rejection", predictor)
        
        if result is not None:
            results_dict[predictor] = {
                "coefficient": result.params[predictor],
                "std_error": result.bse[predictor],
                "p_value": result.pvalues[predictor],
                "r_squared": result.rsquared,
                "n": len(result.fittedvalues)
            }
    
    # Create DataFrame from results
    results_df = pd.DataFrame.from_dict(results_dict, orient="index")
    results_df.index.name = "Predictor"
    results_df.reset_index(inplace=True)
    
    # Save results
    results_df.to_csv(os.path.join(RESULTS_DIR, "nanoplastics_meta_regression.csv"), index=False)
    
    return results_df

def plot_regression_bubble(data, x_var, y_var="log_rejection", color_var="membrane_type", size_var="weight_normalized"):
    """
    Create bubble plot for regression visualization
    
    Parameters:
    -----------
    data : pandas DataFrame
        DataFrame containing variables for plot
    x_var : str
        Column name for x-axis variable
    y_var : str
        Column name for y-axis variable (default: "log_rejection")
    color_var : str
        Column name for color variable (default: "membrane_type")
    size_var : str
        Column name for bubble size variable (default: "weight_normalized")
    """
    # Filter valid data
    valid_data = data[[x_var, y_var, color_var, size_var]].dropna()
    
    if len(valid_data) < 10:
        print(f"Not enough data for bubble plot with {x_var}")
        return
    
    # Create figure
    plt.figure(figsize=(10, 7))
    
    # Get unique categories for color
    categories = valid_data[color_var].unique()
    colors = sns.color_palette("colorblind", n_colors=len(categories))
    
    # Create scatter plot
    for i, category in enumerate(categories):
        mask = valid_data[color_var] == category
        plt.scatter(
            valid_data.loc[mask, x_var],
            valid_data.loc[mask, y_var],
            s=valid_data.loc[mask, size_var] * 5,  # Scale bubble size
            color=colors[i],
            alpha=0.7,
            label=category,
            edgecolor="white",
            linewidth=0.5
        )
    
    # Calculate regression line
    slope, intercept, r_value, p_value, std_err = stats.linregress(
        valid_data[x_var], valid_data[y_var]
    )
    
    # Plot regression line
    x_range = np.linspace(valid_data[x_var].min(), valid_data[x_var].max(), 100)
    plt.plot(x_range, intercept + slope * x_range, 'k--', alpha=0.7)
    
    # Add regression equation
    plt.text(
        0.05, 0.95,
        f"y = {slope:.3f}x + {intercept:.3f}\nr = {r_value:.3f}, p = {p_value:.4f}",
        transform=plt.gca().transAxes,
        bbox=dict(facecolor='white', alpha=0.8)
    )
    
    # Format plot
    plt.xlabel(x_var.replace("_", " ").title())
    plt.ylabel("Log Rejection Value")
    plt.title(f"Relationship Between {x_var.replace('_', ' ').title()} and Rejection Performance")
    plt.legend(title=color_var.replace("_", " ").title())
    plt.grid(True, alpha=0.3)
    
    # Save figure
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURE_DIR, f"bubble_plot_{x_var}.png"), dpi=300)
    plt.close()

if __name__ == "__main__":
    print("Running meta-regression analysis...")
    data = load_data()
    
    # Run meta-regression for each contaminant class
    pfas_results = run_pfas_regression(data)
    pharma_results = run_pharmaceutical_regression(data)
    cmp_results = run_cmp_regression(data)
    np_results = run_nanoplastics_regression(data)
    
    # Create bubble plots for key relationships
    
    # PFAS chain length plot
    pfas_data = data[data["contaminant_class"] == "PFAS"].copy()
    if len(pfas_data) >= 10:
        plot_regression_bubble(pfas_data, "pfas_chain_length")
    
    # Molecular weight plot for pharmaceuticals
    pharma_data = data[data["contaminant_class"] == "Pharmaceuticals"].copy()
    if len(pharma_data) >= 10:
        plot_regression_bubble(pharma_data, "molecular_weight")
    
    # Log Kow plot for CMP organics
    cmp_data = data[data["contaminant_class"] == "CMP Organics"].copy()
    if len(cmp_data) >= 10:
        plot_regression_bubble(cmp_data, "log_kow")
    
    # Particle size plot for nanoplastics
    np_data = data[data["contaminant_class"] == "Micro/nanoplastics"].copy()
    if len(np_data) >= 10:
        plot_regression_bubble(np_data, "particle_size")
    
    print("Meta-regression analysis complete. Results saved to", RESULTS_DIR)
