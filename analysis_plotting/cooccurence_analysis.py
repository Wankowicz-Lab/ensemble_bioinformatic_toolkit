import pandas as pd 

# =========================== INPUT PARAMETERS =========================== #
# Define file paths for input CSVs
RESIDUE_NETWORK_WEIGHTS_FILE = '/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/network/residue_network_weights.csv'
OPENER_NETWORK_WEIGHTS_FILE = '/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/network/series3_residue_network_weights.csv'
ORDER_ALL_A_FILE = '/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/order_all_A.csv'

# Define output file paths
MERGED_WEIGHTS_OUTPUT_FILE = '/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/series2_merged_weights.csv'
S2CALC_COMPARISON_OUTPUT_FILE = '/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/series2_s2calc_comparison.csv'

# Define PDB series for comparison
SERIES3_PDBS = ['x3638', 'x3769']
SERIES2_PDBS = ['x3677', 'x3675', 'x3707', 'x3694']

# Select which PDB series to use for analysis
PDB_IDS = SERIES2_PDBS  # Change this to SERIES3_PDBS if needed

# =========================== PROCESS NETWORK WEIGHTS =========================== #
# Load residue network weight data for two different conditions
residue_network_weights = pd.read_csv(RESIDUE_NETWORK_WEIGHTS_FILE)
opener_residue_network_weights = pd.read_csv(OPENER_NETWORK_WEIGHTS_FILE)

# Merge the two datasets on 'Residue1' and 'Residue2' to ensure direct comparison
merged_weights = pd.merge(
    residue_network_weights, 
    opener_residue_network_weights, 
    on=['Residue1', 'Residue2'], 
    suffixes=('_original', '_opener')
)

# Compute the difference in network weights between the two conditions
merged_weights['Weight_Difference'] = merged_weights['Weight_original'] - merged_weights['Weight_opener']

# Save the merged data with differences to a CSV file for further analysis
merged_weights.to_csv(MERGED_WEIGHTS_OUTPUT_FILE, index=False)

# Identify the 10 residue pairs with the largest weight differences
largest_differences = merged_weights.nlargest(10, 'Weight_Difference')

# Print the largest weight differences
print("Largest differences in weights between the two files:")
print(largest_differences[['Residue1', 'Residue2', 'Weight_Difference']])

# Identify the 10 residue pairs with the smallest weight differences
smallest_differences = merged_weights.nsmallest(10, 'Weight_Difference')

# Print the smallest weight differences
print("Smallest differences in weights between the two files:")
print(smallest_differences[['Residue1', 'Residue2', 'Weight_Difference']])


