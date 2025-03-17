import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd


def generate_op_heatmap(OP_df, output_path):
    """
    Generate a heatmap of OP and save it to the specified output path.

    Parameters:
    order_all_A (DataFrame): The input dataframe containing OP data.
    output_path (str): The path where the heatmap image will be saved.
    """
    # Calculate standard deviations
    s2_sd = OP_df.groupby(['resi', 'chain'])['s2calc'].std().reset_index()
    s2_sd_diff = OP_df.groupby(['resi', 'chain'])['s2calc_diff'].std().reset_index()
    s2_sd_clean = s2_sd.dropna(subset=['s2calc'])

    # Create a pivot table for the clustermap
    pivot_op = OP_df.pivot_table(index=['PDB'], columns='resi', values='s2calc_diff').fillna(0)
    pivot_op = pivot_op.drop(columns=[150], errors='ignore')

    # Generate the clustermap
    plt.figure()
    sns.color_palette("coolwarm", as_cmap=True)
    g = sns.clustermap(pivot_op, cmap='coolwarm_r', figsize=(10, 8), col_cluster=True)

    # Extract the row and column order from the clustermap
    row_order = g.dendrogram_row.reordered_ind
    col_order = g.dendrogram_col.reordered_ind

    # Reorder the pivot_op DataFrame based on the clustermap order
    pivot_op_reordered = pivot_op.iloc[row_order, col_order]

    # Output the reordered pivot_op DataFrame
    pivot_op_reordered.to_csv(output_path.replace('.png', '_pivot_op_reordered.csv'))

    # Save the heatmap
    plt.savefig(output_path)


def plot_op_distribution(OP_df, output_path):
    """
    Plot the distribution of OP differences and save it to the specified output path.

    Parameters:
    OP_df (DataFrame): The input dataframe containing OP data.
    output_path (str): The path where the plot image will be saved.
    """
    plt.figure(figsize=(20, 6))
    sns.boxplot(x='resi', y='s2calc', data=OP_df, flierprops=dict(marker='o', color='gray', markersize=1))
    plt.xlabel('resi')
    plt.ylabel('s2calc')
    plt.xticks(rotation=90)
    plt.savefig(output_path)

def plot_stddev_s2calc(OP_df, output_path):
    """
    Plot the standard deviation of 's2calc' for each residue and save it to the specified output path.

    Parameters:
    OP_df (DataFrame): The input dataframe containing OP data.
    output_path (str): The path where the plot image will be saved.
    """
    s2_sd = OP_df.groupby(['resi', 'chain'])['s2calc'].std().reset_index()
    
    plt.figure(figsize=(20, 6))
    plt.plot(s2_sd['resi'], s2_sd['s2calc'], marker='o', color='b')
    plt.xlabel('Residue (resi)')
    plt.ylabel('Standard Deviation of s2calc')
    plt.xticks(s2_sd['resi'], rotation=45)
    plt.savefig(output_path)



def estimate_order_parameter(resn, s2):
    """Estimate the entropy contribution for a given residue and order parameter S^2. Taken from Li & Brüschweiler 2009
    Example usage: op_df[['Estimated_OP', 'Max_Entropy', 'Min_Entropy']] = op_df.apply(lambda row: pd.Series(estimate_order_parameter(row['resn'], row['s2calc'])), axis=1)
    """
   
    # Constants from Table 1 of the paper
    AA_OP_PARAMS = {
    'V': {'A': 2.19, 'B': 1.32, 'M': 1, 'f': 'linear', 'max_entropy': 0.00697437, 'min_entropy': 0.00435153},
    'S': {'A': 2.19, 'B': 1.32, 'M': 1, 'f': 'linear', 'max_entropy': 0.00697437, 'min_entropy': 0.00435153},
    'T': {'A': 2.19, 'B': 1.32, 'M': 1, 'f': 'linear', 'max_entropy': 0.00697437, 'min_entropy': 0.00435153},
    'I': {'A': 1.95, 'B': 1.55, 'M': 2, 'f': 'linear', 'max_entropy': 0.01390900, 'min_entropy': 0.00774930},
    'L': {'A': 1.95, 'B': 1.55, 'M': 2, 'f': 'linear', 'max_entropy': 0.01390900, 'min_entropy': 0.00774930},
    'M': {'A': 2.73, 'B': 0.77, 'M': 3, 'f': 'linear', 'max_entropy': 0.02086350, 'min_entropy': 0.01627353},
    'N': {'A': 2.06, 'B': 2.08, 'M': 2, 'f': 'log', 'max_entropy': 0.00818644, 'min_entropy': 0.00818644},
    'Q': {'A': 2.16, 'B': 1.60, 'M': 3, 'f': 'log', 'max_entropy': 0.01287576, 'min_entropy': 0.01287576},
    'F': {'A': 2.07, 'B': 1.51, 'M': 2, 'f': 'linear', 'max_entropy': 0.01422692, 'min_entropy': 0.00822618},
    'H': {'A': 2.07, 'B': 1.51, 'M': 2, 'f': 'linear', 'max_entropy': 0.01422692, 'min_entropy': 0.00822618},
    'Y': {'A': 2.07, 'B': 1.51, 'M': 2, 'f': 'linear', 'max_entropy': 0.01422692, 'min_entropy': 0.00822618},
    'P': {'A': 1.90, 'B': 1.15, 'M': 1, 'f': 'linear', 'max_entropy': 0.00606035, 'min_entropy': 0.00377530},
    'K': {'A': 2.20, 'B': 1.22, 'M': 4, 'f': 'linear', 'max_entropy': 0.02718216, 'min_entropy': 0.01748560},
    'R': {'A': 2.22, 'B': 1.23, 'M': 5, 'f': 'linear', 'max_entropy': 0.03427575, 'min_entropy': 0.02205570},
    'D': {'A': 3.69, 'B': 0.44, 'M': 2, 'f': 'log', 'max_entropy': 0.01466406, 'min_entropy': 0.01466406},
    'E': {'A': 3.66, 'B': 0.64, 'M': 3, 'f': 'log', 'max_entropy': 0.02181726, 'min_entropy': 0.02181726},
    'backbone': {'A': 3.42, 'B': 0.50, 'M': 2, 'f': 'log', 'max_entropy': 0.01359108, 'min_entropy': 0.01359108}
    }

    # Boltzmann constant in kcal/(mol·K)
    k_B = 1.987e-3  
    
    # Dictionary to convert three-letter amino acid codes to one-letter codes
    three_to_one = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }

    # Convert three-letter code to one-letter code if necessary
    if len(resn) == 3:
        resn = three_to_one.get(resn.upper(), resn)
    resn = resn.upper()
    
    if resn not in AA_OP_PARAMS:
        return np.nan, np.nan, np.nan  # Skip unknown residues

    params = AA_OP_PARAMS[resn]
    M, A, B, f_type = params['M'], params['A'], params['B'], params['f']
    
    op_term = (1 - s2)
    
    if f_type == 'log':
        op_term = np.log(op_term) if op_term > 0 else 0  # Avoid log(0)

    entropy = k_B * M * (A + B * op_term)
    max_entropy = params['max_entropy']
    min_entropy = params['min_entropy']
    
    return entropy, max_entropy, min_entropy
