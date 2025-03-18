import os

# Specifica il percorso della directory principale
base_directory = "."

# Loop su tutte le cartelle nel range specificato
for i in range(33027, 33340):  # Intervallo di numeri
    folder_name = f"RoutingSetup3S{i}"
    folder_path = os.path.join(base_directory, folder_name)
    
    # Percorso del file configuration.txt
    config_file_path = os.path.join(folder_path, "configuration.txt")
    
    # Controlla se il file esiste nella cartella corrente
    if os.path.exists(config_file_path):
        try:
            # Leggi il contenuto del file
            with open(config_file_path, 'r') as file:
                content = file.read()
            
            # Sostituisci '33026' con il numero della cartella corrente
            updated_content = content.replace("33026", str(i))
            
            # Scrivi il contenuto aggiornato nel file
            with open(config_file_path, 'w') as file:
                file.write(updated_content)
            
            print(f"Modificato con successo: {config_file_path}")
        
        except Exception as e:
            print(f"Errore durante la modifica di {config_file_path}: {e}")
    else:
        print(f"File non trovato: {config_file_path}")
