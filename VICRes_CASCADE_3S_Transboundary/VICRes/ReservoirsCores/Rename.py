import os
import shutil

# Percorso base per le cartelle di destinazione
base_path = "./"

# Cartella sorgente da cui copiare i file
source_folder = os.path.join(base_path, "../ReservoirsCompromizeFlushing")

# Controlla che la cartella sorgente esista
if not os.path.exists(source_folder):
    raise FileNotFoundError(f"La cartella sorgente {source_folder} non esiste.")

# Loop attraverso le cartelle target
for i in range(33001, 33327):
    folder_name = f"ReservoirsOpt{i}"
    folder_path = os.path.join(base_path, folder_name)

    if os.path.exists(folder_path):
        # Rimuove tutti i file nella cartella
        for filename in os.listdir(folder_path):
            file_path = os.path.join(folder_path, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)  # Cancella file o collegamenti
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)  # Cancella directory
            except Exception as e:
                print(f"Errore durante la rimozione di {file_path}: {e}")
    else:
        # Crea la cartella se non esiste
        os.makedirs(folder_path)

    # Copia i file dalla cartella sorgente
    try:
        for filename in os.listdir(source_folder):
            src_file = os.path.join(source_folder, filename)
            dest_file = os.path.join(folder_path, filename)
            if os.path.isfile(src_file):
                shutil.copy2(src_file, dest_file)
    except Exception as e:
        print(f"Errore durante la copia dei file in {folder_path}: {e}")

print("Operazione completata.")