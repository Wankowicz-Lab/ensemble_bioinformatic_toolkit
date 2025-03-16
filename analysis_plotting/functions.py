import matplotlib.pyplot as plt
import seaborn as sns



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
