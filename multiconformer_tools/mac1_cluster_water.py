import os
import glob
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt

# Define the folder path for the PDB files
pdb_folder_path = '/Users/stephaniewanko/Downloads/vanderbilt/mac1/water/*.pdb'

# Initialize a list to store results for each PDB file
results = []
water_coords = []
# Loop through each PDB file
for file_path in glob.glob(pdb_folder_path):
    
    # Extract water coordinates from the current PDB file
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('HETATM') and 'HOH' in line:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                water_coords.append([x, y, z])
    
# Convert the list of water coordinates to a NumPy array
water_coords = np.array(water_coords)

# Perform DBSCAN clustering on the water coordinates
dbscan = DBSCAN(eps=2.0, min_samples=10)
clusters = dbscan.fit_predict(water_coords)

# Create a DataFrame to store the water coordinates and their cluster labels
water_df = pd.DataFrame(water_coords, columns=['x', 'y', 'z'])
water_df['cluster'] = clusters

# Filter clusters to ensure each cluster's radius is within 2.0 angstroms
valid_clusters = []
cluster_centers = []  # List to store cluster centers
for cluster_id in water_df['cluster'].unique():
    if cluster_id == -1:
        continue  # Skip noise points
    cluster_points = water_df[water_df['cluster'] == cluster_id][['x', 'y', 'z']].values
    cluster_center = np.mean(cluster_points, axis=0)
    distances = np.linalg.norm(cluster_points - cluster_center, axis=1)
    if np.all(distances <= 2.0):
        valid_clusters.append(cluster_id)
        cluster_centers.append((cluster_center, len(cluster_points)))  # Save cluster center and size
# Output each cluster center to PDB format with the number of waters as the b-factor
with open(f'cluster_centers_{os.path.basename(file_path)}.pdb', 'w') as pdb_file:
    for i, (center, size) in enumerate(cluster_centers):
        pdb_file.write(f"HETATM{i+1:5d}  O   HOH A{i+1:4d}    {center[0]:8.3f}{center[1]:8.3f}{center[2]:8.3f}  1.00{size:6.2f}           O\n")

# Count waters in valid clusters and outliers
in_cluster_count = water_df[water_df['cluster'].isin(valid_clusters)].shape[0]
out_of_cluster_count = water_df[~water_df['cluster'].isin(valid_clusters)].shape[0]

# Store the results for the current PDB file
results.append({
    'file': os.path.basename(file_path),
    'in_cluster': in_cluster_count,
    'out_of_cluster': out_of_cluster_count
})

# Define a function to count waters within 1.0 angstrom from cluster centers
def count_waters_within_distance(pdb_file, cluster_centers, distance_threshold=1.0):
    water_coords = []
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith('HETATM') and 'HOH' in line:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                water_coords.append([x, y, z])
    
    water_coords = np.array(water_coords)
    counts = []
    for center, _ in cluster_centers:
        distances = np.linalg.norm(water_coords - center, axis=1)
        count = np.sum(distances <= distance_threshold)
        counts.append(count)
    
    return counts, water_coords

# Example usage
all_percent_in_cluster_centers = []
pdb_folder_path = '/Users/stephaniewanko/Downloads/vanderbilt/mac1/water/*.pdb'
for pdb_file_path in glob.glob(pdb_folder_path):
    print(pdb_file_path)
    counts, water_coords = count_waters_within_distance(pdb_file_path, cluster_centers)
    total_waters = len(water_coords)
    print(total_waters)
    percent_in_cluster_center = (sum(counts) / total_waters) * 100
    all_percent_in_cluster_centers.append(percent_in_cluster_center)


plt.figure(figsize=(10, 6))
plt.hist(all_percent_in_cluster_centers, bins=20, edgecolor='black')
plt.title('Distribution of Percent Waters in Cluster Centers')
plt.xlabel('Percent of Waters in Cluster Centers')
plt.ylabel('Frequency')
plt.grid(True)
plt.savefig('cluster_centers_distribution.png')

def count_waters_between_structures_within_distance(pdb_folder_path, distance_threshold=1.0):
    pdb_files = glob.glob(pdb_folder_path)
    results = []

    for i in range(len(pdb_files)):
        for j in range(i + 1, len(pdb_files)):
            pdb_file_1 = pdb_files[i]
            pdb_file_2 = pdb_files[j]

            water_coords_1 = []
            with open(pdb_file_1, 'r') as file:
                for line in file:
                    if line.startswith('HETATM') and 'HOH' in line:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        water_coords_1.append([x, y, z])

            water_coords_2 = []
            with open(pdb_file_2, 'r') as file:
                for line in file:
                    if line.startswith('HETATM') and 'HOH' in line:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        water_coords_2.append([x, y, z])

            water_coords_1 = np.array(water_coords_1)
            water_coords_2 = np.array(water_coords_2)

            count = 0
            for coord_1 in water_coords_1:
                distances = np.linalg.norm(water_coords_2 - coord_1, axis=1)
                count += np.sum(distances <= distance_threshold)

            # Calculate the percentage
            total_waters = (len(water_coords_1) + len(water_coords_2)) / 2
            percent_within_distance = (count / total_waters) * 100

            results.append({
                'file_1': os.path.basename(pdb_file_1),
                'file_2': os.path.basename(pdb_file_2),
                'count_within_distance': count,
                'percent_within_distance': percent_within_distance
            })

    return results

# Example usage
pdb_folder_path = '/Users/stephaniewanko/Downloads/vanderbilt/mac1/water/*.pdb'
results = count_waters_between_structures_within_distance(pdb_folder_path)
for result in results:
    print(f"File 1: {result['file_1']}, File 2: {result['file_2']}, Count within 1.0 Å: {result['count_within_distance']}, Percent within 1.0 Å: {result['percent_within_distance']}")


