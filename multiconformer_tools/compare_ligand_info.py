
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import numpy as np
from statsmodels.stats.multitest import multipletests

# Set font and figure params
plt.rcParams.update({
    'font.size': 24,
    'axes.titlesize': 24,
    'axes.labelsize': 24,
    'xtick.labelsize': 24,
    'ytick.labelsize': 24,
    'legend.fontsize': 24
})

boxprops = {'edgecolor': 'k', 'linewidth': 2}
lineprops = {'color': 'k', 'linewidth': 2}

boxplot_kwargs = dict({'boxprops': boxprops, 'medianprops': lineprops,
                       'whiskerprops': lineprops, 'capprops': lineprops,
                       'width': 0.75})

#read in affinity csv with variables, op or rmsf csv, resi close, with argpase

# Read in the CSV files
mac1_affinity_op_subsets = pd.read_csv('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/mac1_affinity_op_subsets.csv')
merged_ligand_smiles_clusters = pd.read_csv('merged_ligand_smiles_clusters.csv')
close_OP = pd.read_csv('avg_close_delta_rmsf.csv')
# Read in the avg_OP.csv file
all_OP = pd.read_csv('avg_delta_rmsf.csv')

all_OP.rename(columns={'delta_RMSF': 's2calc_diff_all'}, inplace=True)
close_OP.rename(columns={'delta_RMSF': 's2calc_diff_close'}, inplace=True)
all_OP.rename(columns={'PDB_ensemble': 'PDB'}, inplace=True)
close_OP.rename(columns={'PDB_ensemble': 'PDB'}, inplace=True)

# Merge the avg_OP data with the existing merged_data on a common column, assuming 'PDB' is the common column
merged_data = pd.merge(merged_ligand_smiles_clusters, all_OP, on='PDB')
# Merge the close_OP data with the existing merged_data_avg on a common column, assuming 'PDB' is the common column
merged_data = pd.read_csv('all_ligand_properties.csv') #pd.merge(merged_data, close_OP, on='PDB')

column_name_map = {
    'MW': 'Molecular Weight',
    'Hbond_donor': 'Hydrogen Bond Donors',
    'Hbond_acceptor': 'Hydrogen Bond Acceptors',
    'Polar_SA': 'Polar Surface Area',
    'rotatable': 'Rotatable Bonds',
    'Hbond_MW': 'Hydrogen Bonds/Molecular Weight'
}

merged_data.rename(columns=column_name_map, inplace=True) 

print(merged_data.head())
# # Define the columns to plot
columns_to_plot = ['Molecular Weight', 'logP', 'Hydrogen Bond Donors', 'Hydrogen Bond Acceptors',
                    'LogD', 'Polar Surface Area', 'Rotatable Bonds', 'Hydrogen Bonds/Molecular Weight']

#Loop through each column and create a separate scatter plot for close_OP
# for column in columns_to_plot:
#     plt.figure(figsize=(8, 6))  # Create a new figure for each column
#     sns.scatterplot(data=merged_data, x='s2calc_diff_close', y=column)
#     regplot = sns.regplot(data=merged_data, x='s2calc_diff_close', y=column, scatter=False, color='red')  # Add best fit line
    
#     # Calculate and print regression line information
#     slope, intercept = np.polyfit(merged_data['s2calc_diff_close'], merged_data[column], 1)
#     correlation_matrix = np.corrcoef(merged_data['s2calc_diff_close'], merged_data[column])
#     correlation_xy = correlation_matrix[0, 1]
#     r_squared = correlation_xy**2
#     p_value = stats.pearsonr(merged_data['s2calc_diff_close'], merged_data[column])[1]
#     print(f'Regression line for {column} (close_OP): slope = {slope}, intercept = {intercept}, r^2 = {r_squared}, p-value = {p_value}')
    
#     plt.xlabel('OP Difference', fontsize=20)
#     plt.ylabel(column, fontsize=20)
#     plt.xticks(rotation=45, fontsize=20)
#     plt.yticks(fontsize=20)
#     plt.tight_layout()
#     safe_col = column.replace(" ", "")
#     plt.savefig(f'figure_out/scatter_s2calc_{safe_col}_close.png', dpi=300)
#     plt.close()


#Loop through each column and create a separate scatter plot for avg_OP
# for column in columns_to_plot:
#     plt.figure(figsize=(8, 6))  # Create a new figure for each column
#     sns.scatterplot(data=merged_data, x='s2calc_diff_all', y=column)
#     regplot = sns.regplot(data=merged_data, x='s2calc_diff_all', y=column, scatter=False, color='red')  # Add best fit line
    
#     # Calculate and print regression line information
#     slope, intercept = np.polyfit(merged_data['s2calc_diff_all'], merged_data[column], 1)
#     correlation_matrix = np.corrcoef(merged_data['s2calc_diff_all'], merged_data[column])
#     correlation_xy = correlation_matrix[0, 1]
#     r_squared = correlation_xy**2
#     p_value = stats.pearsonr(merged_data['s2calc_diff_all'], merged_data[column])[1]
#     print(f'Regression line for {column} (avg_OP): slope = {slope}, intercept = {intercept}, r^2 = {r_squared}, p-value = {p_value}')
    
#     plt.xlabel('OP Difference', fontsize=24)
#     plt.ylabel(column, fontsize=24)
#     plt.xticks(rotation=45, fontsize=20)
#     plt.yticks(fontsize=20)
#     plt.tight_layout()
#     safe_col = column.replace(" ", "")
#     plt.savefig(f'figure_out/scatter_s2calc_{safe_col}_all.png', dpi=300)
#     plt.close()


sns.set(style="whitegrid")
# Loop through each column and create a separate scatter plot
# for column in columns_to_plot:
#     plt.figure(figsize=(8, 6))  # Create a new figure for each column
#     sns.scatterplot(data=merged_data, x='s2calc_diff_close', y=column)
#     regplot = sns.regplot(data=merged_data, x='s2calc_diff_close', y=column, scatter=False, color='red')  # Add best fit line
    
#     # Calculate and print regression line information
#     slope, intercept = np.polyfit(merged_data['s2calc_diff_close'], merged_data[column], 1)
#     correlation_matrix = np.corrcoef(merged_data['s2calc_diff_close'], merged_data[column])
#     correlation_xy = correlation_matrix[0, 1]
#     r_squared = correlation_xy**2
#     p_value = stats.pearsonr(merged_data['s2calc_diff_close'], merged_data[column])[1]
#     print(f'Regression line for {column}: slope = {slope}, intercept = {intercept}, r^2 = {r_squared}, p-value = {p_value}')
    
#     plt.xlabel('OP Difference', fontsize=24)
#     plt.ylabel(column, fontsize=24)
#     plt.xticks(rotation=45, fontsize=24)
#     plt.yticks(fontsize=24)
#     plt.tight_layout()
#     plt.savefig(f'figure_out/scatter_s2calc_{column}_close.png', dpi=300)
#     plt.close()


for column in columns_to_plot:
    print(column)
    plt.figure(figsize=(8, 6))  # Create a new figure for each column
    sns.violinplot(data=merged_data, x='Cluster', y=column, inner=None)
    plt.xticks(rotation=45, fontsize=16)
    plt.yticks(fontsize=20)
    plt.tight_layout()
    safe_col = column.replace(" ", "")
    print(safe_col)
    safe_col = safe_col.replace("/", "")
    print('safe col:')
    print(safe_col)
    plt.savefig(f'figure_out/violin_{safe_col}_close.png', dpi=300)
    plt.close()

    print(column)
    print('CLOSE')
    # Calculate the 1st and 4th quartiles for the current column
    first_quartile = merged_data[column].quantile(0.25)
    print('first quartile:')
    print(first_quartile)
    fourth_quartile = merged_data[column].quantile(0.75)
    print('fourth quartile:')
    print(fourth_quartile)

    # Filter the PDBs in the first and fourth quartiles
    first_quartile_pdbs = merged_data[merged_data[column] <= first_quartile]['PDB']
    fourth_quartile_pdbs = merged_data[merged_data[column] >= fourth_quartile]['PDB']

    # Get the s2calc_diff values for PDBs in the first and fourth quartiles from close_OP
    first_quartile_s2calc_diff = close_OP[close_OP['PDB'].isin(first_quartile_pdbs)]['s2calc_diff_close']
    fourth_quartile_s2calc_diff = close_OP[close_OP['PDB'].isin(fourth_quartile_pdbs)]['s2calc_diff_close']

    # Combine the data for plotting
    quartile_data = pd.DataFrame({
        's2calc_diff_close': pd.concat([first_quartile_s2calc_diff, fourth_quartile_s2calc_diff]),
        'Quartile': ['1st Quartile'] * len(first_quartile_s2calc_diff) + ['4th Quartile'] * len(fourth_quartile_s2calc_diff)
    }) 

    # Plot the violin plot
    plt.figure(figsize=(8, 6))
    sns.violinplot(data=quartile_data, x='Quartile', y='s2calc_diff_close', inner=None)
    plt.ylabel('OP Differences', fontsize=24)
    plt.xlabel(column, fontsize=24)
    plt.xticks(rotation=45, fontsize=24)
    plt.yticks(fontsize=24)
    plt.tight_layout()
    plt.savefig(f'figure_out/violin_s2calc_diff_close_{safe_col}.png', dpi=300)
    plt.close()

    # Calculate and output statistics
    stats_data = quartile_data.groupby('Quartile')['s2calc_diff_close'].describe()
    print(stats_data)

    # Perform Mann-Whitney U test
    u_statistic, p_value = stats.mannwhitneyu(first_quartile_s2calc_diff, fourth_quartile_s2calc_diff, alternative='two-sided')
    print(f'Mann-Whitney U test statistic: {u_statistic}, p-value: {p_value}')

    corrected_p_values = multipletests([p_value], method='fdr_bh')[1]
    print(f'Corrected p-value: {corrected_p_values[0]}')


# Loop through each column and create a separate boxplot
for column in columns_to_plot:
    plt.figure(figsize=(8, 6))  # Create a new figure for each column
    sns.violinplot(data=merged_data, x='Cluster', y=column, inner=None)
    plt.xticks(rotation=45, fontsize=20)
    plt.yticks(fontsize=20)
    plt.tight_layout()
    safe_col = column.replace(" ", "")
    print(safe_col)
    safe_col = safe_col.replace("/", "")
    print('safe col:')
    print(safe_col)
    plt.savefig(f'figure_out/violin_{safe_col}_all.png', dpi=300)
    plt.close()

    print(column)
    print('ALL')
    # Calculate the 1st and 4th quartiles for the current column
    first_quartile = merged_data[column].quantile(0.25)
    fourth_quartile = merged_data[column].quantile(0.75)

    # Filter the PDBs in the first and fourth quartiles
    first_quartile_pdbs = merged_data[merged_data[column] <= first_quartile]['PDB']
    fourth_quartile_pdbs = merged_data[merged_data[column] >= fourth_quartile]['PDB']

    # Get the s2calc_diff values for PDBs in the first and fourth quartiles from close_OP
    first_quartile_s2calc_diff = all_OP[all_OP['PDB'].isin(first_quartile_pdbs)]['s2calc_diff_all']
    fourth_quartile_s2calc_diff = all_OP[all_OP['PDB'].isin(fourth_quartile_pdbs)]['s2calc_diff_all']

    # Combine the data for plotting
    quartile_data = pd.DataFrame({
        's2calc_diff_all': pd.concat([first_quartile_s2calc_diff, fourth_quartile_s2calc_diff]),
        'Quartile': ['1st Quartile'] * len(first_quartile_s2calc_diff) + ['4th Quartile'] * len(fourth_quartile_s2calc_diff)
    })

    # Plot the violin plot
    plt.figure(figsize=(8, 6))
    sns.violinplot(data=quartile_data, x='Quartile', y='s2calc_diff_all', inner=None)
    plt.ylabel('OP Differences', fontsize=24)
    plt.xlabel(column, fontsize=20)
    plt.xticks(rotation=45, fontsize=20)
    plt.yticks(fontsize=20)
    plt.tight_layout()
    plt.savefig(f'figure_out/violin_s2calc_diff_all_{safe_col}.png', dpi=300)
    plt.close()

    # Calculate and output statistics
    stats_data = quartile_data.groupby('Quartile')['s2calc_diff_all'].describe()
    print(stats_data)

    # Perform Mann-Whitney U test
    u_statistic, p_value = stats.mannwhitneyu(first_quartile_s2calc_diff, fourth_quartile_s2calc_diff, alternative='two-sided')
    print(f'Mann-Whitney U test statistic: {u_statistic}, p-value: {p_value}')
    # stats_data.to_csv(f'figure_out/stats_s2calc_diff_{column}.csv')


# Plot log(IC50) distribution per cluster
merged_data['log_IC50'] = np.log(merged_data['Affinity'])  # Convert IC50 to log(IC50)
plt.figure(figsize=(12, 6))
sns.boxplot(data=merged_data, x='Cluster', y='log_IC50')
plt.ylabel('log(IC50)', fontsize=24)
plt.xlabel('Cluster', fontsize=24)
plt.xticks(rotation=45, fontsize=24)
plt.yticks(fontsize=24)
plt.tight_layout()
plt.savefig('figure_out/log_ic50_per_cluster.png', dpi=300)
plt.close()


plt.figure(figsize=(12, 6))
sns.boxplot(data=merged_data, x='Cluster', y='Affinity')
plt.ylabel('IC50', fontsize=24)
plt.xlabel('Cluster', fontsize=24)
plt.xticks(rotation=45, fontsize=24)
plt.yticks(fontsize=24)
plt.tight_layout()
plt.savefig('figure_out/ic50_per_cluster.png', dpi=300)
plt.close()

merged_data_lig = merged_data

# Define clusters for altloc and noaltloc
altloc_clusters = [0, 3, 5, 6, 9, 10, 13]
merged_data_lig['altloc_label'] = merged_data_lig['Cluster'].apply(lambda x: 'Alt Loc' if x in altloc_clusters else 'No Alt Loc')
merged_data['altloc_label'] = merged_data['Cluster'].apply(lambda x: 'Alt Loc' if x in altloc_clusters else 'No Alt Loc')
alt_loc_count = merged_data_lig[merged_data_lig['altloc_label'] == 'Alt Loc'].shape[0]
print(f"Number of entries with altloc_label == 'Alt Loc': {alt_loc_count}")


# # Define columns to plot
# columns_to_plot = ['s2calc_diff_all','s2calc_diff_close', 'MW', 'logP', 'Hbond_donor', 'Hbond_acceptor',
#                     'LogD', 'Polar_SA', 'rotatable', 'Hbond_MW']  # Example columns, adjust as needed

# # Plot violin plots for each column
# for column in columns_to_plot:
#     plt.figure(figsize=(12, 6))
#     sns.violinplot(data=merged_data, x='altloc_label', y=column, inner=None)
#     plt.ylabel(column, fontsize=24)
#     plt.xlabel('')
#     plt.xticks(rotation=45, fontsize=24)
#     plt.yticks(fontsize=24)
#     plt.tight_layout()
#     plt.savefig(f'figure_out/violin_{column}_altloc_vs_noaltloc.png', dpi=300)
#     plt.close()
# # Perform statistical tests to compare 'altloc' vs 'noaltloc' for each column
# for column in columns_to_plot:
#     altloc_data = merged_data[merged_data['altloc_label'] == 'Alt Loc'][column]
#     noaltloc_data = merged_data[merged_data['altloc_label'] == 'No Alt Loc'][column]
#     # Calculate and print the mean and median for 'altloc' and 'noaltloc' data for each column
#     altloc_mean = altloc_data.mean()
#     altloc_median = altloc_data.median()
#     noaltloc_mean = noaltloc_data.mean()
#     noaltloc_median = noaltloc_data.median()
    
#     print(f'{column} - Altloc: Mean = {altloc_mean}, Median = {altloc_median}')
#     print(f'{column} - No Altloc: Mean = {noaltloc_mean}, Median = {noaltloc_median}')
    
#     # Perform Mann-Whitney U test
#     u_statistic, p_value = stats.mannwhitneyu(altloc_data, noaltloc_data, alternative='two-sided')
#     print(f'Mann-Whitney U test for {column}: U statistic = {u_statistic}, p-value = {p_value}')


# columns_to_plot = ['log_IC50', 'IC50']
# for column in columns_to_plot:
#     plt.figure(figsize=(12, 6))
#     sns.violinplot(data=merged_data_lig, x='altloc_label', y=column, inner=None)
#     plt.ylabel(column, fontsize=24)
#     plt.xlabel('Altloc vs No Altloc', fontsize=24)
#     plt.xticks(rotation=45, fontsize=24)
#     plt.yticks(fontsize=24)
#     plt.tight_layout()
#     plt.savefig(f'figure_out/violin_{column}_altloc_vs_noaltloc.png', dpi=300)
#     plt.close()

# for column in columns_to_plot:
#     altloc_data = merged_data_lig[merged_data_lig['altloc_label'] == 'Alt Loc'][column]
#     noaltloc_data = merged_data_lig[merged_data_lig['altloc_label'] == 'No Alt Loc'][column]

#     altloc_mean = altloc_data.mean()
#     altloc_median = altloc_data.median()
#     noaltloc_mean = noaltloc_data.mean()
#     noaltloc_median = noaltloc_data.median()
    
#     print(f'{column} - Altloc: Mean = {altloc_mean}, Median = {altloc_median}')
#     print(f'{column} - No Altloc: Mean = {noaltloc_mean}, Median = {noaltloc_median}')
    
#     # Perform Mann-Whitney U test
#     u_statistic, p_value = stats.mannwhitneyu(altloc_data, noaltloc_data, alternative='two-sided')
#     print(f'Mann-Whitney U test for {column}: U statistic = {u_statistic}, p-value = {p_value}')



# Multiply every value in all_OP by 165
all_OP_scaled = all_OP.copy()
all_OP_scaled.iloc[:, 1:] = all_OP_scaled.iloc[:, 1:] * 165

# Merge the scaled all_OP with merged_data_lig
merged_data_scaled = pd.merge(merged_data_lig, all_OP_scaled, on='PDB', suffixes=('', '_scaled'))

# Determine the top 25% and bottom 25% of IC50
ic50_75th_percentile = merged_data_scaled['IC50'].quantile(0.75)
ic50_25th_percentile = merged_data_scaled['IC50'].quantile(0.25)

# Filter the data for top 25% and bottom 25% of IC50
top_25_ic50 = merged_data_scaled[merged_data_scaled['IC50'] >= ic50_75th_percentile]
bottom_25_ic50 = merged_data_scaled[merged_data_scaled['IC50'] <= ic50_25th_percentile]

# Plot the distribution of all_OP_scaled for top 25% and bottom 25% of IC50
plt.figure(figsize=(14, 8))
sns.kdeplot(data=top_25_ic50, x='s2calc_diff_all_scaled', label='Top 25% IC50', fill=True, color='royalblue', alpha=0.7)
sns.kdeplot(data=bottom_25_ic50, x='s2calc_diff_all_scaled', label='Bottom 25% IC50', fill=True, color='salmon', alpha=0.7)
plt.xlabel('Total Protein Conformational Heterogeneity Difference', fontsize=26, fontweight='bold')
plt.ylabel('Density', fontsize=26, fontweight='bold')
plt.xticks(fontsize=22, fontweight='bold')
plt.yticks(fontsize=22, fontweight='bold')
plt.legend(fontsize=22, title='IC50 Percentile', title_fontsize='22', loc='upper left', frameon=True, shadow=True)
plt.grid(True, linestyle='--', linewidth=0.5)
plt.tight_layout()
plt.savefig('figure_out/scaled_OP_distribution_top_vs_bottom_IC50.png', dpi=300, bbox_inches='tight')
plt.close()
# Calculate and print statistics for the scaled OP differences in top and bottom 25% IC50 groups

# Calculate means
top_25_mean = top_25_ic50['s2calc_diff_all_scaled'].mean()
bottom_25_mean = bottom_25_ic50['s2calc_diff_all_scaled'].mean()

# Calculate medians
top_25_median = top_25_ic50['s2calc_diff_all_scaled'].median()
bottom_25_median = bottom_25_ic50['s2calc_diff_all_scaled'].median()

# Calculate standard deviations
top_25_std = top_25_ic50['s2calc_diff_all_scaled'].std()
bottom_25_std = bottom_25_ic50['s2calc_diff_all_scaled'].std()

# Print the statistics
print(f"Top 25% IC50 - Mean: {top_25_mean}, Median: {top_25_median}, Std Dev: {top_25_std}")
print(f"Bottom 25% IC50 - Mean: {bottom_25_mean}, Median: {bottom_25_median}, Std Dev: {bottom_25_std}")

# Perform a t-test to compare the means of the two groups
t_statistic, p_value = stats.ttest_ind(top_25_ic50['s2calc_diff_all_scaled'], bottom_25_ic50['s2calc_diff_all_scaled'], equal_var=False)
print(f"T-test between top 25% and bottom 25% IC50: t-statistic = {t_statistic}, p-value = {p_value}")

# Perform a Mann-Whitney U test as a non-parametric alternative
u_statistic, p_value_mannwhitney = stats.mannwhitneyu(top_25_ic50['s2calc_diff_all_scaled'], bottom_25_ic50['s2calc_diff_all_scaled'], alternative='two-sided')
print(f"Mann-Whitney U test between top 25% and bottom 25% IC50: U statistic = {u_statistic}, p-value = {p_value_mannwhitney}")

