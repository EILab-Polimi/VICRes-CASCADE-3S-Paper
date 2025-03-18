import os
import shutil

# Percorso base per le cartelle di destinazione
base_path = "./"

# Cartella sorgente da cui copiare i file
#source_folder = os.path.join(base_path, "../ReservoirsCompromizeFlushing")


# Loop attraverso le cartelle target
for i in range(33001, 33327):
    folder_name = f"RoutingSetup3S{i}/Results"
    folder_path = os.path.join(base_path, folder_name)

    if os.path.exists(folder_path):
        # Itera solo sui file specificati
        target_files = ["OUTPUT.day", "OUTPUT.day_mm"]
        for filename in target_files:
            file_path = os.path.join(folder_path, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)  # Cancella file o collegamenti
                    print(f"{file_path} rimosso con successo.")
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)  # Se per caso è una directory
                    print(f"{file_path} directory rimossa con successo.")
            except Exception as e:
                print(f"Errore durante la rimozione di {file_path}: {e}")
    

print("Operazione completata.")