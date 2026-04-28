import os

# Define the variables based on your folder structure
szas = ['45', '48.2', '60']
albedos = ['0.1', '0.2', '0.3', '0.06', '0.15', '0.25']

# Set the base directory to the current working directory
base_dir = '.'

for sza in szas:
    sza_folder_name = f"SZA_{sza}"
    sza_folder_path = os.path.join(base_dir, sza_folder_name)
    
    # Skip if the SZA folder doesn't exist for some reason
    if not os.path.exists(sza_folder_path):
        print(f"Directory not found: {sza_folder_path}")
        continue

    for albedo in albedos:
        albedo_folder_path = os.path.join(sza_folder_path, albedo)
        
        # Define the original and new file paths
        original_file_name = f"Earth_100pc_{sza}.txt"
        original_file_path = os.path.join(albedo_folder_path, original_file_name)
        
        new_file_name = f"Earth_100pc_{sza}_{albedo}.txt"
        new_file_path = os.path.join(sza_folder_path, new_file_name)
        
        # Check if the target file actually exists
        if os.path.exists(original_file_path):
            # Rename the file and move it up to the SZA folder
            os.rename(original_file_path, new_file_path)
            print(f"Success: Renamed to {new_file_name}")
            
            # Remove the now-empty albedo folder
            try:
                os.rmdir(albedo_folder_path)
                print(f"Cleaned up: Removed empty folder '{albedo_folder_path}'")
            except OSError:
                print(f"Warning: Could not remove '{albedo_folder_path}'. It might contain other hidden files.")
        else:
            # Report the missing file
            print(f"File missing: {original_file_path}")

print("\nScript execution complete.")