
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import numpy as np

# Read in the CSV files
mac1_affinity_op_subsets = pd.read_csv('mac1_affinity_op_subsets.csv')
merged_ligand_smiles_clusters = pd.read_csv('merged_ligand_smiles_clusters.csv')
close_OP = pd.read_csv('avg_close_OP.csv')


# Merge the dataframes on a common column, assuming 'PDB' is the common column
merged_data = pd.merge(mac1_affinity_op_subsets, merged_ligand_smiles_clusters, on='PDB')

print("Columns in merged_data:", merged_data.columns)

# Define the columns to plot
columns_to_plot = ['MW', 'logP', 'Hbond_donor', 'Hbond_acceptor',
                    'LogD', 'Polar_SA', 'rotatable', 'Hbond_MW']

# Set the style for the plots
sns.set(style="whitegrid")
# Loop through each column and create a separate scatter plot
for column in columns_to_plot:
    plt.figure(figsize=(8, 6))  # Create a new figure for each column
    sns.scatterplot(data=merged_data, x='s2calc_diff', y=column)
    regplot = sns.regplot(data=merged_data, x='s2calc_diff', y=column, scatter=False, color='red')  # Add best fit line
    
    # Calculate and print regression line information
    slope, intercept = np.polyfit(merged_data['s2calc_diff'], merged_data[column], 1)
    correlation_matrix = np.corrcoef(merged_data['s2calc_diff'], merged_data[column])
    correlation_xy = correlation_matrix[0, 1]
    r_squared = correlation_xy**2
    p_value = stats.pearsonr(merged_data['s2calc_diff'], merged_data[column])[1]
    print(f'Regression line for {column}: slope = {slope}, intercept = {intercept}, r^2 = {r_squared}, p-value = {p_value}')
    
    plt.xlabel('OP Difference', fontsize=24)
    plt.ylabel(column, fontsize=24)
    plt.xticks(rotation=45, fontsize=24)
    plt.yticks(fontsize=24)
    plt.tight_layout()
    plt.savefig(f'figure_out/scatter_s2calc_{column}.png', dpi=300)
    plt.close()


# Loop through each column and create a separate boxplot
for column in columns_to_plot:
    plt.figure(figsize=(8, 6))  # Create a new figure for each column
    sns.violinplot(data=merged_data, x='Cluster', y=column, inner=None)
    plt.xticks(rotation=45, fontsize=24)
    plt.yticks(fontsize=24)
    plt.tight_layout()
    plt.savefig(f'figure_out/violin_{column}.png', dpi=300)
    plt.close()

    for column in columns_to_plot:
        print(column)
        # Calculate the 1st and 4th quartiles for the current column
        first_quartile = merged_data[column].quantile(0.25)
        fourth_quartile = merged_data[column].quantile(0.75)

        # Filter the PDBs in the first and fourth quartiles
        first_quartile_pdbs = merged_data[merged_data[column] <= first_quartile]['PDB']
        fourth_quartile_pdbs = merged_data[merged_data[column] >= fourth_quartile]['PDB']

        # Get the s2calc_diff values for PDBs in the first and fourth quartiles from close_OP
        first_quartile_s2calc_diff = close_OP[close_OP['PDB'].isin(first_quartile_pdbs)]['s2calc_diff']
        fourth_quartile_s2calc_diff = close_OP[close_OP['PDB'].isin(fourth_quartile_pdbs)]['s2calc_diff']

        # Combine the data for plotting
        quartile_data = pd.DataFrame({
            's2calc_diff': pd.concat([first_quartile_s2calc_diff, fourth_quartile_s2calc_diff]),
            'Quartile': ['1st Quartile'] * len(first_quartile_s2calc_diff) + ['4th Quartile'] * len(fourth_quartile_s2calc_diff)
        })

        # Plot the violin plot
        plt.figure(figsize=(8, 6))
        sns.violinplot(data=quartile_data, x='Quartile', y='s2calc_diff', inner=None)
        plt.ylabel('OP Differences', fontsize=24)
        plt.xlabel(column, fontsize=24)
        plt.xticks(rotation=45, fontsize=24)
        plt.yticks(fontsize=24)
        plt.tight_layout()
        plt.savefig(f'figure_out/violin_s2calc_diff_{column}.png', dpi=300)
        plt.close()

        # Calculate and output statistics
        stats_data = quartile_data.groupby('Quartile')['s2calc_diff'].describe()
        print(stats_data)

        # Perform Mann-Whitney U test
        u_statistic, p_value = stats.mannwhitneyu(first_quartile_s2calc_diff, fourth_quartile_s2calc_diff, alternative='two-sided')
        print(f'Mann-Whitney U test statistic: {u_statistic}, p-value: {p_value}')
        # stats_data.to_csv(f'figure_out/stats_s2calc_diff_{column}.csv')
