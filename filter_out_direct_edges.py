import pandas as pd

# Load the network file
network_file = "NETWORK.txt"  # Replace with your network file path
network = pd.read_csv(network_file, sep="\t", names=["#Node1", "Node2", ""])

# Read the file and separate the source and target genes
source_genes = []
target_genes = []

with open("NODE_TYPES_highlyconf_targets.txt", "r") as file:
    for line in file:
        # Skip the header line
        if line.startswith("#Node"):
            continue
        # Parse each line
        gene, gene_type = line.strip().split("\t")
        if gene_type == "source":
            source_genes.append(gene)
        elif gene_type == "target":
            target_genes.append(gene)

# Output the lists
print("Source Genes:", source_genes)
print("Target Genes:", target_genes)

# Remove edges where both nodes are in the gene sets
filtered_network = network[
    ~(
        (network["#Node1"].isin(source_genes) & network["Node2"].isin(target_genes)) |
        (network["#Node1"].isin(target_genes) & network["Node2"].isin(source_genes))
    )
]

# Save the filtered network
output_file = "FILTERED_NETWORK_highlyconf_targets.txt"  # Replace with desired output file path
filtered_network.to_csv(output_file, sep="\t", index=False, header=False)

print(f"Filtered network saved to {output_file}")
