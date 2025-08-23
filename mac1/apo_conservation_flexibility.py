import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

file_path = '/Users/stephaniewanko/Downloads/vanderbilt/mac1/oP/7TX0_B_consurf_grades.txt'

# Read the file, skipping potential header lines
with open(file_path, 'r') as file:
    # Find the header line by checking for specific keywords
    for i, line in enumerate(file):
        if line.startswith("POS"):
            header_line = i
            break


consurf_grades = pd.read_csv(
    file_path,
    sep='\t',
    skiprows=header_line,
    engine='python',
    skipinitialspace=True
)

# Trim whitespace from column names
consurf_grades.columns = consurf_grades.columns.str.strip()

# Filter out non-numeric POS values
consurf_grades = consurf_grades[pd.to_numeric(consurf_grades['POS'], errors='coerce').notnull()]

# Convert POS to integer
consurf_grades['POS'] = consurf_grades['POS'].astype(int)

# Ensure COLOR is treated as a string before checking for asterisks
#consurf_grades['CONFIDENT'] = ~consurf_grades['COLOR'].astype(str).str.contains('\*')
consurf_grades['COLOR'] = consurf_grades['COLOR'].str.replace('*', '', regex=False).fillna(0).astype(int)
consurf_grades['SCORE'] = consurf_grades['SCORE'].astype(float)

apo_mac1 = pd.read_csv('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/apo_mac1.csv')
apo_rmsf_df = pd.read_csv('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/other_macro_output/ensemble/rmsf_output_7TX0.updated_refine_001_ensemble.csv')


# Calculate the average s2calc per residue
avg_s2calc_per_residue = apo_mac1.groupby('resi')['s2calc'].mean().reset_index()
avg_s2calc_per_residue.columns = ['Residue', 'Avg_s2calc']


high_flex_resi = [101, 103, 104, 105, 157, 159, 165, 30, 20, 161, 102, 158, 162, 29, 100, 160, 163, 106, 57, 49, 63, 78, 99, 100, 120, 121, 122, 246, 18]


# Add a new column to indicate if a residue is in high_flex_resi
avg_s2calc_per_residue['High Flexibility'] = avg_s2calc_per_residue['Residue'].apply(lambda x: 'Increased Flexibility' if x in high_flex_resi else 'No Increased Flexibility')

# Create a boxplot of Avg_s2calc separated by high flexibility status
plt.figure(figsize=(10, 6))
sns.violinplot(data=avg_s2calc_per_residue, x='High Flexibility', y='Avg_s2calc', palette=['#1F77B4', '#FF4500'], inner=None)
plt.ylabel('Average Apo Order Parameter', fontsize=17, fontweight='bold')
plt.xlabel('')
plt.xticks(rotation=30, fontsize=17, fontweight='bold')
plt.yticks(fontsize=17, fontweight='bold')
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/avg_s2calc_flexibility_violinplot.png', dpi=300)
plt.close()

# Calculate the average s2calc for residues with increased flexibility
avg_s2calc_high_flex = avg_s2calc_per_residue[avg_s2calc_per_residue['High Flexibility'] == 'Increased Flexibility upon Ligand Binding']['Avg_s2calc'].mean()

# Calculate the average s2calc for residues without increased flexibility
avg_s2calc_no_high_flex = avg_s2calc_per_residue[avg_s2calc_per_residue['High Flexibility'] == 'No Increased Flexibility upon Ligand Binding']['Avg_s2calc'].mean()

print(f"Average s2calc for residues with increased flexibility: {avg_s2calc_high_flex}")
print(f"Average s2calc for residues without increased flexibility: {avg_s2calc_no_high_flex}")

# Add High Flexibility column to apo_rmsf_df
apo_rmsf_df['High Flexibility'] = apo_rmsf_df['resi'].apply(lambda x: 'Increased Flexibility upon Ligand Binding' if x in high_flex_resi else 'No Increased Flexibility upon Ligand Binding')

# Create a boxplot of RMSF separated by high flexibility status
plt.figure(figsize=(10, 6))
sns.boxplot(data=apo_rmsf_df, x='High Flexibility', y='RMSF', palette='Set2')
plt.xlabel('Change in Flexibility', fontsize=14, fontweight='bold')
plt.ylabel('RMSF', fontsize=14, fontweight='bold')
plt.xticks(rotation=45)
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig('/Users/stephaniewanko/Downloads/vanderbilt/mac1/OP/rmsf_flexibility_boxplot.png', dpi=300)
plt.close()

# Calculate the average RMSF for residues with increased flexibility
avg_rmsf_high_flex = apo_rmsf_df[apo_rmsf_df['High Flexibility'] == 'Increased Flexibility upon Ligand Binding']['RMSF'].mean()

# Calculate the average RMSF for residues without increased flexibility
avg_rmsf_no_high_flex = apo_rmsf_df[apo_rmsf_df['High Flexibility'] == 'No Increased Flexibility upon Ligand Binding']['RMSF'].mean()

print(f"Average RMSF for residues with increased flexibility: {avg_rmsf_high_flex}")
print(f"Average RMSF for residues without increased flexibility: {avg_rmsf_no_high_flex}")

