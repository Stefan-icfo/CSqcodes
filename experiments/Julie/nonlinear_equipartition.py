import qcodes as qc
run_ids = [
1262, 1264, 1266, 1298,1300,1302,1304,1308,1310, 1312,1314,1318,1320,1322, 1324, 1328, 1330, 1332, 1334, 1338, 1340, 1342, 1344
]

print(run_ids)
print(len(run_ids))
qc.config["core"]["db_location"] = ".\Data\Raw_data\CD11_D7_C1_part3.db"
def get_metadata(meas_id):
    """
    Retrieve and return the metadata for a single measurement ID.
    """
    # Load the dataset using the given measurement ID
    dataset = qc.load_by_id(meas_id)
    # Return the metadata dictionary
    return dataset.metadata


# Dictionary to store metadata for each run ID
metadata_dict = {}

# Iterate over each run ID and collect metadata
for rid in run_ids:
    try:
        meta = get_metadata(meas_id=rid)
        metadata_dict[rid] = meta
        print(f"Metadata for run ID {rid} retrieved successfully.")
    except Exception as e:
        print(f"Error retrieving metadata for run ID {rid}: {e}")

#max_avg_avg_psd_
# Convert the metadata dictionary into a DataFrame
import pandas as pd
df = pd.DataFrame.from_dict(metadata_dict, orient='index').reset_index()
df.rename(columns={'index': 'run_id'}, inplace=True)

# Define the list of columns you want to KEEP
columns_to_keep = ['run_id', 'max_avg_avg_psd_']  

df = df[columns_to_keep]

# Example temperature list (must match the number of rows in df), in m Kelvin
temperature_list = [35.6766,35.6766,35.6766, 50.0152,50.0152,50.0152, 50.0152, 75.0531,75.0531,75.0531,75.0531,100.38,100.38, 100.38, 100.38, 124.884,124.884,124.884,124.884,149.994,149.994,149.994,149.994]  

# Insert at position 2 (third column)
df.insert(2, 'Temperature', temperature_list)

# Display the DataFrame
print(df)


#####################################################################################
