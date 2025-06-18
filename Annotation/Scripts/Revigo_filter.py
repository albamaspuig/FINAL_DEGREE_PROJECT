import pandas as pd
import os

# Folder containing REVIGO output TSVs 
input_folder = './'  

# REVIGO files for each GO domain
revigo_files = {
    'BP': 'Revigo_BP_Table.tsv',
    'CC': 'Revigo_CC_Table.tsv',
    'MF': 'Revigo_MF_Table.tsv'
}

# Filtering threshold for Dispensability
dispensability_threshold = 0.15

# PROCESSING FUNCTION

def filter_revigo_table(file_path, output_path, threshold=dispensability_threshold):
    df = pd.read_csv(file_path, sep="\t")

    # Clean potential quotes from text columns
    df['Name'] = df['Name'].str.replace('"', '')
    df['TermID'] = df['TermID'].str.replace('"', '')

    # Filter by dispensability
    filtered = df[df['Dispensability'] < threshold].copy()

    # Sort by Frequency ascending (more specific GO terms first)
    filtered = filtered.sort_values(by='Frequency', ascending=True)

    # Save filtered table
    filtered[['TermID', 'Name', 'Value', 'Frequency', 'Dispensability']].to_csv(output_path, sep="\t", index=False)

    print(f"Filtered table saved to: {output_path}")
    print(filtered[['TermID', 'Name', 'Value', 'Frequency', 'Dispensability']].head(10))
    return filtered

# Do the filtering for each table

for go_domain, filename in revigo_files.items():
    input_file = os.path.join(input_folder, filename)
    output_file = os.path.join(input_folder, f"Revigo_{go_domain}_filtered.tsv")

    filter_revigo_table(input_file, output_file)

print("Filtering completed for all GO domains.")

