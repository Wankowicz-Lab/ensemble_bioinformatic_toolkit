import argparse
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
from scipy.stats import ttest_ind
from scipy import stats
import matplotlib.pyplot as plt

def get_args():
    parser = argparse.ArgumentParser(description="Bootstrap test to compare residue flexibility (s2calc_diff) between two groups of PDBs.")
    parser.add_argument('--order_all_A', type=str, required=True, help="Path to the CSV file containing flexibility metrics.")
    parser.add_argument('--group_1_pdbs', type=str, required=True, help="Comma-separated list of PDBs for Group 1.")
    parser.add_argument("--n_iterations", type=int, default=1000, help="Number of bootstrap iterations (default 1000)")
    parser.add_argument("--alpha", type=float, default=0.05, help="Significance level for the confidence interval (default 0.05)")
    
    args = parser.parse_args()
    group_1_pdbs = args.group_1_pdbs.split(",")  # Convert comma-separated string into list
    
    return args.order_all_A, group_1_pdbs, args.n_iterations, args.alpha


def load_data(file_path):
    """Load the CSV containing the flexibility metrics for all PDBs."""
    return pd.read_csv(file_path)

def split_groups(data, group_1_pdbs):
    """Split data into Group 1 (specified PDBs) and Group B (all other PDBs)."""
    group_1_pdbs = set(group_1_pdbs)
    group_1_data = data[data['PDB'].isin(group_1_pdbs)]
    group_b_data = data[~data['PDB'].isin(group_1_pdbs)]
    return group_1_data, group_b_data

def bootstrap_resampling(group_1_residue_data, group_b_residue_data, n_iterations):
    """Perform bootstrap resampling to calculate confidence intervals."""
    np.random.seed(42)
    
    group1 = np.array(group_1_residue_data)
    group2 = np.array(group_b_residue_data)

    observed_t_stat, real_p_value = stats.ttest_ind(group1, group2, equal_var=False)
    print(f"Actual t-test p-value: {real_p_value}")

    num_bootstraps = n_iterations
    bootstrap_tstats = []

    combined_data = np.concatenate([group1, group2])
    n1 = len(group1)

    for _ in range(num_bootstraps):
        np.random.shuffle(combined_data)  # Shuffle data for permutation test
        random_group1 = combined_data[:n1]
        random_group2 = combined_data[n1:]
        
        t_stat, _ = stats.ttest_ind(random_group1, random_group2, equal_var=False)
        bootstrap_tstats.append(t_stat)

    lower_ci = np.percentile(bootstrap_tstats, 2.5)
    upper_ci = np.percentile(bootstrap_tstats, 97.5)

    # Compute bootstrap p-value
    bootstrapped_p_value = np.mean(np.abs(bootstrap_tstats) >= np.abs(observed_t_stat))

    print(f"Bootstrapped p-value: {bootstrapped_p_value}")
    print(f"95% Confidence Interval for t-statistics: ({lower_ci}, {upper_ci})")

    return lower_ci, upper_ci, bootstrapped_p_value

def create_file(residues, p_values, corrected_p_values, reject, output_file, group_1_stats, group_b_stats, ttest_p_values):
    """Save the results to a text file."""
    # Create a pandas DataFrame
    results_df = pd.DataFrame({
        'Residue': residues,
        'Lower CI': [p[0] for p in p_values],
        'Upper CI': [p[1] for p in p_values],
        'Non-corrected P-value': [p[2] for p in p_values],
        'Corrected P-value': corrected_p_values,
        'T-test P-value': ttest_p_values,
        'Significant': ['Yes' if sig else 'No' for sig in reject],
        'Group1_Median': [stat['median'] for stat in group_1_stats],
        'Group1_SD': [stat['sd'] for stat in group_1_stats],
        'GroupB_Median': [stat['median'] for stat in group_b_stats],
        'GroupB_SD': [stat['sd'] for stat in group_b_stats]
    })

    # Sort the DataFrame by significance
    results_df = results_df.sort_values(by='Significant', ascending=False)

    # Output the DataFrame to a CSV file
    results_df.to_csv(output_file, index=False)

def calculate_p_value(bootstrap_differences, observed_difference):
    """Calculate the p-value based on bootstrap differences."""
    more_extreme = np.sum(np.abs(observed_difference) >= np.abs(bootstrap_differences))
    p_value = more_extreme / len(bootstrap_differences)
    return p_value

def main():
    # Step 1: Parse arguments
    order_all_A, group_1_pdbs, n_iterations, alpha = get_args()
    
    # Step 2: Load the data
    data = load_data(order_all_A)

    # Step 3: Split data into Group 1 and Group B
    group_1_data, group_b_data = split_groups(data, group_1_pdbs)

    residues = []
    p_values = []
    
    group_1_stats = []
    group_b_stats = []
    
    ttest_p_values = []  # List to store t-test p-values
    
    # Step 4: Perform bootstrap resampling and calculate the confidence intervals for each residue
    for residue in group_1_data['resi'].unique():
        group_1_residue_data = group_1_data[group_1_data['resi'] == residue]['s2calc_diff']
        group_b_residue_data = group_b_data[group_b_data['resi'] == residue]['s2calc_diff']
        
        if len(group_1_residue_data) > 0 and len(group_b_residue_data) > 0:
            lower, upper, p_value = bootstrap_resampling(group_1_residue_data, group_b_residue_data, n_iterations, residue, p_values)
            observed_tstat = stats.ttest_ind(group_1_residue_data, group_b_residue_data, equal_var=False).statistic
            #p_value = calculate_p_value(bootstrap_tstats, observed_tstat)
            residues.append(residue)
            p_values.append([lower, upper, p_value])
            
            # Calculate t-test p-value
            t_stat, ttest_p_value = ttest_ind(group_1_residue_data, group_b_residue_data)
            ttest_p_values.append(ttest_p_value)
            
            # Calculate median and standard deviation for each group
            group_1_stats.append({
                'median': np.median(group_1_residue_data),
                'sd': np.std(group_1_residue_data)
            })
            group_b_stats.append({
                'median': np.median(group_b_residue_data),
                'sd': np.std(group_b_residue_data)
            })
    
    # Step 5: Apply FDR correction
    all_p_values = [p[2] for p in p_values]  # Use the calculated p-values
    reject, corrected_p_values, _, _ = multipletests(all_p_values, method='fdr_bh')
    
    # Step 6: Output results to a file
    output_file = "series5_comparison_results_bootstrap.csv"
    create_file(residues, p_values, corrected_p_values, reject, output_file, group_1_stats, group_b_stats, ttest_p_values)
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    main()
