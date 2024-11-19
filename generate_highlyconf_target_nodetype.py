import pandas as pd

# File paths
input_file = "NODE_TYPES_highlyconf_targets.txt"

# Genes to add as 'target'
new_genes = ["APC", "ATP1A1", "F2RL1", "HIF1A", "HTR4", "LINC00473", "SCAP"]

# Load the existing file into a DataFrame
df = pd.read_csv(input_file, sep="\t", comment="#", names=["Node", "Node type"], skiprows=1)

# Create a new DataFrame for the genes to add
new_entries = pd.DataFrame({
    "Node": new_genes,
    "Node type": ["target"] * len(new_genes)
})

# Append the new entries
updated_df = pd.concat([df, new_entries], ignore_index=True)

# Save the updated DataFrame to a new file
updated_file_path = "NODE_TYPES_highlyconf_targets.txt"
updated_df.to_csv(updated_file_path, sep="\t", index=False)

updated_file_path
