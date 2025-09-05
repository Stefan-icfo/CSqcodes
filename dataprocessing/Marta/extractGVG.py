import matplotlib.pyplot as plt
from qcodes.dataset import load_by_id, initialise_or_create_database_at

# Percorso al database
db_path = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD13_E3_C2.db"

# Inizializza (o collega) il database
initialise_or_create_database_at(db_path)

# Carica il run con ID 287
dataset = load_by_id(287)

# Controlla che funzioni
print(dataset)

# Estrai in dataframe
data = dataset.to_pandas_dataframe()
print(data.head())  # per vedere le colonne

# Plot
plt.figure(figsize=(6,4))
plt.plot(data["ch1(V)"], data["G_IV (Siemens)"], "o-")
plt.xlabel("ch1 (V)")
plt.ylabel("G_IV (Siemens)")
plt.title("Run 287 - ch1 vs G_IV")
plt.grid(True)
plt.show()





