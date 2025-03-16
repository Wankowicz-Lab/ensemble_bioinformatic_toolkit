import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import zscore
import seaborn as sns

def plot_s2calc_diff_distribution(data, group_1_pdbs):
    # Filter data to only include rows with PDB in group_1
    group_1_data = data[data['PDB'].isin(group_1_pdbs)]
    print(group_1_data.head())
    
    # Initialize a dictionary to store z-scores
    z_scores = {pdb: {} for pdb in group_1_pdbs}
    
    # Iterate over each PDB in group_1
    for pdb in group_1_pdbs:
        print(pdb)
        # Filter data for the current PDB
        pdb_data = group_1_data[group_1_data['PDB'] == pdb]
        
        # Iterate over each residue in the current PDB
        for residue in pdb_data['resi'].unique():
            #print(residue)
            # Filter data for the current residue
            residue_data = pdb_data[pdb_data['resi'] == residue]['s2calc_diff']
            residue_data_all = group_1_data[group_1_data['resi'] == residue]['s2calc_diff']
            
            # Plot the s2calc_diff distribution for the entire group_1_data
            plt.figure()
            plt.hist(residue_data_all, bins=30, alpha=0.7, label='Group 1 Data')
            plt.axvline(x=residue_data.mean(), color='r', linestyle='--', label=f'Residue {residue} Mean')
            plt.xlabel('s2calc_diff')
            plt.ylabel('Frequency')
            plt.title('s2calc_diff Distribution for Group 1 Data')
            plt.legend()
            plt.savefig(f's2calc_diff_distribution_group_1_data.png')
            plt.close()
        
            # Calculate the z-score for the residue and PDB removed
            all_residue_data = group_1_data[group_1_data['resi'] == residue]['s2calc_diff']
            if not residue_data.empty:
                z_scores_array = zscore(all_residue_data)
                #print(z_scores_array)
                # Find the index of the current residue in the full dataset
                # Use the first index of residue_data to find the corresponding z-score
                current_index = residue_data.index[0]
                # Ensure the index is within bounds
                if current_index in all_residue_data.index:
                    z_scores[pdb][residue] = z_scores_array[current_index]
            else:
                z_scores[pdb][residue] = np.nan  # or handle the empty case as needed
    
    # Convert z-scores dictionary to a DataFrame
    z_scores_df = pd.DataFrame(z_scores).T
    z_scores_df.index.name = 'PDB'
    z_scores_df.columns.name = 'Residue'
    
    # Save the z-scores table to a CSV file
    z_scores_df.to_csv('z_scores_table.csv')
    
    return z_scores_df

def plot_zscore_heatmap(z_scores_df):
    """Plot a clustermap of the z-scores."""
    plt.figure(figsize=(12, 8))
    sns.clustermap(z_scores_df, cmap="coolwarm", fmt=".2f")
    plt.xlabel("Residue")
    plt.ylabel("PDB")
    plt.savefig("z_scores_clustermap.png")
    plt.close()

# Load the data
data = pd.read_csv('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/order_all_A.csv')

# Define group_1_pdbs
series5 = ['x3871', 'x3928', 'x4393', 'x3870', 'x4034', 'x3753', 'x3928', 'x4371', 'x3885', 'x3867', 'x3522', 'x3888', 'x4160']
group_1_pdbs = series5
# Plot s2calc_diff distribution and calculate z-scores
z_scores_df = plot_s2calc_diff_distribution(data, group_1_pdbs)

# Plot the z-scores heatmap
plot_zscore_heatmap(z_scores_df)

# Print the z-scores table
print(z_scores_df)
