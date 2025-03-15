import pandas as pd 
residue_network_weights = pd.read_csv('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/network/residue_network_weights.csv')
opener_residue_network_weights = pd.read_csv('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/network/series3_residue_network_weights.csv')


# Ensure that 'Residue1' and 'Residue2' in each file are the same before calculating differences
merged_weights = pd.merge(
    residue_network_weights, 
    opener_residue_network_weights, 
    on=['Residue1', 'Residue2'], 
    suffixes=('_original', '_opener')
)

# Calculate the absolute differences in weights between the two files
merged_weights['Weight_Difference'] = merged_weights['Weight_original'] - merged_weights['Weight_opener']

df = merged_weights
df.to_csv('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/series2_merged_weights.csv', index=False)


# Identify the largest differences in weights
largest_differences = merged_weights.nlargest(10, 'Weight_Difference')

# Print the largest differences
print("Largest differences in weights between the two files:")
print(largest_differences[['Residue1', 'Residue2', 'Weight_Difference']])

# Print the smallest differences
smallest_differences = merged_weights.nsmallest(10, 'Weight_Difference')

print("Smallest differences in weights between the two files:")
print(smallest_differences[['Residue1', 'Residue2', 'Weight_Difference']])



order_all_A = pd.read_csv('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/order_all_A.csv')
# Function to compare average s2calc values of every residue to the average s2calc value of specified PDBs
def compare_s2calc_values(order_all_A, pdb_ids):
    # Calculate the average s2calc for each residue
    print(order_all_A.head())
    avg_s2calc_residues = order_all_A.groupby('resi')['s2calc_diff'].median().reset_index()
    avg_s2calc_residues = avg_s2calc_residues.rename(columns={'s2calc_diff': 's2calc_all_pdbs'})

    # Filter the data for the specified PDBs
    filtered_data = order_all_A[order_all_A['PDB'].isin(pdb_ids)]

    # Calculate the average s2calc for the specified PDBs
    avg_s2calc_pdbs = filtered_data.groupby('resi')['s2calc_diff'].median().reset_index()
    avg_s2calc_pdbs = avg_s2calc_pdbs.rename(columns={'s2calc_diff': 's2calc_pdbs'})

    # Merge the average s2calc values of residues with the average s2calc values of the specified PDBs
    comparison_df = pd.merge(avg_s2calc_residues, avg_s2calc_pdbs, on='resi', how='left')
    print(comparison_df.head())

    # Calculate the absolute differences in s2calc values
    comparison_df['s2calc_difference'] = abs(comparison_df['s2calc_all_pdbs'] - comparison_df['s2calc_pdbs'])


    return comparison_df


series3 = ['x3638', 'x3769']
series2 = ['x3677', 'x3675', 'x3707', 'x3694']
# openers = [ 'x3631',  'x3460', 'x4393',
#     'x3458', 'x3439', 'x3466', 'x3476', 'x3623', 'x3422', 'x3924', 'x3406',
#     'x3410', 'x3459', 'x3481', 'x3233', 'x3844', 'x4037', 'x3402', 'x3465'
# ]
pdb_ids = series2

# comparison_df = compare_s2calc_values(order_all_A, pdb_ids)

# # Print the comparison dataframe and the largest differences
# print("Comparison of average s2calc values:")
# print(comparison_df.sort_values(by='s2calc_difference', ascending=False))
# print(comparison_df.sort_values(by='s2calc_difference', ascending=True))


# # Save the comparison dataframe and the largest differences to CSV files
# comparison_df.to_csv('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/series2_s2calc_comparison.csv', index=False)
