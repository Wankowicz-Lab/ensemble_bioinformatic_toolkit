import argparse
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
from scipy.stats import ttest_ind
from scipy import stats
import matplotlib.pyplot as plt

def get_args():
    parser = argparse.ArgumentParser(description="Bootstrap test to compare residue flexibility (s2calc_diff) between two groups of PDBs.")
    series2 = [ 'x3677', 'x3675', 'x3707', 'x3694']
    openers = [ 'x3631',  'x3460', 'x4393',
    'x3458', 'x3439', 'x3466', 'x3476', 'x3623', 'x3422', 'x3924', 'x3406',
    'x3410', 'x3459', 'x3481', 'x3233', 'x3844', 'x4037', 'x3402', 'x3465'
]
    series3 = ['x3638', 'x3769']
    series4 = ['x3682', 'x2498', 'x4022', 'x4239', 'x3401', 'x4393', 'x3476', 'x3952']
    series5 = ['x3871', 'x3928', 'x4393', 'x3870', 'x4034', 'x3753', 'x3928', 'x4371', 'x3885', 'x3867', 'x3522', 'x3888', 'x4160']
    group_1_pdbs = series4
    parser.add_argument('--order_all_A')
    parser.add_argument("--n_iterations", type=int, default=1000, help="Number of bootstrap iterations (default 10000)")
    parser.add_argument("--alpha", type=float, default=0.05, help="Significance level for the confidence interval (default 0.05)")
    args = parser.parse_args()
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

def bootstrap_resampling(group_1_residue_data, group_b_residue_data, n_iterations, residue, p_values):
    """Perform bootstrap resampling to calculate confidence intervals and plot graphics."""
    # === Step 1: Define Data ===
    np.random.seed(42)
    group1 = np.array(group_1_residue_data)
    group2 = np.array(group_b_residue_data)

    # === Step 2: Compute Actual T-Test ===
    t_stat, real_p_value = stats.ttest_ind(group1, group2)
    print(f"Actual t-test p-value: {real_p_value}")

    # === Step 3: Bootstrap Random T-Tests ===
    num_bootstraps = n_iterations  # Number of resampling iterations
    significant_count = 0  # Count of t-tests with p < 0.05
    bootstrap_tstats = []  # Store t-statistics for each bootstrap

    combined_data = np.concatenate([group1, group2])
    n1 = len(group1)

    for _ in range(num_bootstraps):
        np.random.shuffle(combined_data)  # Shuffle labels
        random_group1 = combined_data[:n1]
        random_group2 = combined_data[n1:]
        
        t_stat, p_value = stats.ttest_ind(random_group1, random_group2)
        bootstrap_tstats.append(t_stat)
        if p_value < 0.05:
            significant_count += 1

    # === Step 4: Compute Bootstrapped P-Value ===
    bootstrapped_p_value = significant_count / num_bootstraps
    print(f"Bootstrapped p-value: {bootstrapped_p_value}")

    # === Step 5: Calculate Confidence Intervals ===
    lower_ci = np.percentile(bootstrap_tstats, 2.5)
    upper_ci = np.percentile(bootstrap_tstats, 97.5)
    print(f"95% Confidence Interval for t-statistics: ({lower_ci}, {upper_ci})")

    # Return the bootstrapped p-value and confidence intervals
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
