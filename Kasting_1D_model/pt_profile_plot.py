import os
import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

states = ['0.1', '1', '10', '100', '150']

colors = {
    '0.1': '#1f77b4', # Blue
    '1':   '#ff7f0e', # Orange
    '10':  '#2ca02c', # Green
    '100': '#d62728', # Red
    '150': '#9467bd'  # Purple
}

temperature_blocks = {
    '0.1': """
      DATA T/279.77, 273.99, 269.26, 263.67, 257.43, 
     &       250.58, 243.33, 235.85, 228.32, 221.16, 
     &       214.77, 209.53, 205.26, 201.37, 197.54, 
     &       193.93, 190.91, 188.69, 187.10, 185.87, 
     &       184.76, 183.67, 182.58, 181.55, 180.65, 
     &       179.98, 179.66, 179.73, 180.10, 180.30, 
     &       179.98, 179.09, 177.74, 176.15, 174.74, 
     &       173.77, 173.21, 172.96, 172.90, 172.79, 
     &       172.59, 172.26, 171.76, 171.06, 170.29, 
     &       169.44, 168.65, 167.89, 167.39, 166.98, 
     &       166.92, 167.12, 167.56, 168.67, 169.79, 
     &       171.94, 174.32, 177.12, 180.90, 184.69, 
     &       189.51, 194.91, 200.32, 206.47, 212.81, 
     &       219.16, 225.48, 231.81, 238.13, 245.89, 
     &       256.50, 267.10, 277.70, 292.00, 310.11, 
     &       328.21, 346.31, 364.42, 383.04, 401.81, 
     &       420.58, 439.34, 458.11, 476.88, 493.28, 
     &       508.15, 523.01, 537.88, 552.75, 567.62, 
     &       582.49, 597.36, 608.66, 618.68, 628.70, 
     &       638.73, 648.75, 658.77, 668.79, 678.81/
    """,
    
    '1': """
      DATA T/282.20, 276.75, 272.15, 266.67, 260.75, 
     &       254.29, 247.39, 240.20, 232.87, 225.80, 
     &       219.49, 214.38, 210.40, 207.06, 204.01, 
     &       201.31, 199.19, 197.62, 196.29, 195.36, 
     &       194.95, 195.13, 195.65, 196.18, 196.49, 
     &       196.46, 196.00, 195.09, 193.71, 191.76, 
     &       189.41, 186.74, 183.91, 181.16, 178.98, 
     &       177.27, 176.01, 175.23, 175.08, 175.31, 
     &       175.53, 175.67, 175.67, 175.58, 175.43, 
     &       175.24, 174.97, 174.67, 174.27, 173.87, 
     &       173.43, 172.99, 172.75, 172.52, 172.54, 
     &       172.78, 173.03, 173.83, 174.64, 175.68, 
     &       176.96, 178.25, 179.74, 181.25, 182.78, 
     &       184.37, 185.95, 187.65, 189.45, 191.25, 
     &       194.13, 197.32, 200.51, 206.22, 212.23, 
     &       218.24, 226.86, 236.02, 245.18, 255.40, 
     &       267.15, 278.91, 290.66, 302.94, 315.92, 
     &       328.90, 341.88, 354.86, 367.33, 379.54, 
     &       391.74, 403.94, 416.15, 428.35, 438.57, 
     &       448.63, 458.69, 468.75, 478.81, 488.87/
    """,
    
    '10': """
      DATA T/282.64, 277.26, 272.71, 267.26, 261.39, 
     &       255.00, 248.17, 241.07, 233.85, 226.94, 
     &       220.88, 216.13, 212.41, 209.30, 206.66, 
     &       204.80, 204.21, 204.79, 205.46, 205.96, 
     &       206.33, 206.75, 207.36, 208.28, 209.54, 
     &       211.11, 213.00, 215.19, 217.62, 220.24, 
     &       222.98, 225.40, 227.41, 228.44, 228.67, 
     &       228.40, 227.03, 225.22, 223.01, 220.53, 
     &       218.03, 215.57, 213.12, 210.67, 208.22, 
     &       205.77, 203.32, 200.92, 198.58, 196.45, 
     &       194.57, 192.94, 191.48, 190.15, 188.81, 
     &       187.43, 186.06, 184.72, 183.37, 182.20, 
     &       181.27, 180.33, 180.00, 179.83, 179.80, 
     &       180.49, 181.19, 182.27, 183.73, 185.19, 
     &       186.74, 188.31, 189.86, 190.81, 191.76, 
     &       192.53, 192.73, 192.93, 193.11, 193.25, 
     &       193.40, 194.04, 195.02, 196.01, 198.36, 
     &       201.20, 204.04, 209.21, 214.94, 220.67, 
     &       228.91, 238.26, 247.61, 257.68, 270.37, 
     &       283.05, 295.73, 308.59, 322.59, 336.58/
    """,
    
    '100': """
      DATA T/282.30, 276.89, 272.30, 266.84, 260.92, 
     &       254.46, 247.57, 240.42, 233.19, 226.33, 
     &       220.39, 215.83, 212.30, 209.36, 206.92, 
     &       205.40, 205.51, 207.28, 209.61, 211.81, 
     &       213.69, 215.42, 217.14, 218.89, 220.72, 
     &       222.58, 224.41, 226.22, 227.99, 229.84, 
     &       231.87, 234.06, 236.59, 239.27, 242.07, 
     &       244.89, 247.68, 250.43, 253.14, 255.82, 
     &       258.33, 260.77, 262.73, 264.58, 265.71, 
     &       266.73, 267.05, 267.25, 266.80, 266.17, 
     &       265.04, 263.67, 262.01, 260.09, 258.06, 
     &       255.69, 253.32, 250.64, 247.82, 244.99, 
     &       242.02, 239.03, 236.04, 233.11, 230.20, 
     &       227.30, 224.43, 221.60, 218.77, 216.05, 
     &       213.53, 211.01, 208.65, 206.84, 205.02, 
     &       203.43, 202.67, 201.92, 201.37, 201.55, 
     &       201.73, 202.04, 202.70, 203.36, 203.90, 
     &       204.13, 204.37, 204.38, 203.86, 203.33, 
     &       202.67, 201.71, 200.75, 199.87, 199.14, 
     &       198.40, 198.23, 198.77, 199.31, 201.30/
    """,
    
    '150': """
      DATA T/282.28, 276.90, 272.33, 266.90, 261.00, 
     &       254.56, 247.70, 240.57, 233.38, 226.59, 
     &       220.75, 216.30, 212.86, 209.97, 207.58, 
     &       206.09, 206.21, 208.01, 210.34, 212.60, 
     &       214.56, 216.41, 218.23, 220.09, 222.07, 
     &       224.09, 226.08, 228.03, 229.92, 231.83, 
     &       233.85, 235.97, 238.36, 240.85, 243.46, 
     &       246.11, 248.77, 251.41, 254.08, 256.80, 
     &       259.55, 262.31, 264.89, 267.36, 269.39, 
     &       271.17, 272.41, 273.29, 273.76, 273.80, 
     &       273.56, 272.80, 271.95, 270.53, 269.10, 
     &       267.28, 265.32, 263.31, 260.79, 258.27, 
     &       255.70, 252.80, 249.90, 246.99, 244.05, 
     &       241.09, 238.14, 235.17, 232.19, 229.22, 
     &       226.26, 223.42, 220.58, 217.74, 215.41, 
     &       213.12, 210.83, 209.32, 208.02, 206.72, 
     &       206.10, 205.76, 205.42, 205.53, 205.84, 
     &       206.15, 206.35, 206.50, 206.65, 206.38, 
     &       205.87, 205.36, 204.52, 203.49, 202.46, 
     &       201.49, 200.55, 199.62, 199.48, 199.67/
    """
}
# =====================================================================

def parse_fortran_temps(text_block):
    """
    Takes a raw Fortran data block string, cleans out the 'DATA T/', '&', 
    and '/', and returns a clean Python list of floats.
    """
    if not text_block.strip():
        return []
        
    text = re.sub(r'DATA\s+T\s*/', '', text_block, flags=re.IGNORECASE)
    text = text.replace('/', '')
    text = text.replace('&', '')
    
    numbers = re.findall(r'[-+]?\d*\.\d+|\d+', text)
    return [float(n) for n in numbers]

def extract_pressures(file_path):
    """
    Safely parses the Kasting output file regardless of line wrapping.
    """
    pressures = []
    expected_z = 50000.0  # Starting altitude in cm
    
    try:
        with open(file_path, 'r') as f:
            for line in f:
                parts = line.split()
                if not parts:
                    continue
                try:
                    val = float(parts[0])
                    if abs(val - expected_z) < 1.0: 
                        p_val_hpa = float(parts[1]) / 1000.0
                        pressures.append(p_val_hpa)
                        expected_z += 100000.0  
                except ValueError:
                    pass 
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        
    return pressures

# Initialize Plot
plt.figure(figsize=(10, 12))

global_p_max = -float('inf')
global_p_min = float('inf')

for state in states:
    color = colors[state]
    raw_temp_text = temperature_blocks[state]
    
    temperatures = parse_fortran_temps(raw_temp_text)
    dat_file = f"{state}pc_PT_profile/Single_SZA/OUTPUT_PLOT.dat"
    
    if temperatures and os.path.exists(dat_file):
        pressures = extract_pressures(dat_file)
        
        if len(pressures) == len(temperatures):
            global_p_max = max(global_p_max, max(pressures))
            global_p_min = min(global_p_min, min(pressures))
            
            plt.plot(temperatures, pressures, linestyle='-', color=color, linewidth=2.5)
        else:
            print(f"Warning: Extracted {len(pressures)} pressure levels for {state}% (Expected {len(temperatures)}).")
    else:
        print(f"Warning: Could not find Kasting output file at {dat_file}")

    csv_file_pal = f"waccm6_{state}%PAL_PT_profile.csv"
    csv_file_alt = f"WACCM6_{state}%_PT_profile.csv"
    
    csv_file = csv_file_pal if os.path.exists(csv_file_pal) else csv_file_alt
    
    if os.path.exists(csv_file):
        df_csv = pd.read_csv(csv_file)
        pressures_csv = df_csv['Press [hPa]']
        
        global_p_max = max(global_p_max, pressures_csv.max())
        global_p_min = min(global_p_min, pressures_csv.min())
        
        plt.plot(df_csv['Temp [K]'], pressures_csv, linestyle=':', color=color, linewidth=2.5)
    else:
        print(f"Warning: Could not find WACCM CSV file for {state}%")

# Format the plot
plt.yscale('log')

if global_p_max != -float('inf') and global_p_min != float('inf'):
    plt.ylim(global_p_max, global_p_min)
else:
    plt.ylim(1000, 0.001)

plt.xlabel('Temperature (K)', fontsize=16)
plt.ylabel('Pressure (hPa)', fontsize=16)
plt.xlim(150, 300) 
plt.ylim(0,0.001)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

# --- Create Custom Legend ---
legend_elements = [
    Line2D([0], [0], color='black', linestyle='-', lw=2.5, label='Kasting (Solid)'),
    Line2D([0], [0], color='black', linestyle=':', lw=2.5, label='WACCM (Dotted)'),
    Line2D([0], [0], color=colors['0.1'], linestyle='-', lw=4, label='0.1% PAL'),
    Line2D([0], [0], color=colors['1'], linestyle='-', lw=4, label='1% PAL'),
    Line2D([0], [0], color=colors['10'], linestyle='-', lw=4, label='10% PAL'),
    Line2D([0], [0], color=colors['100'], linestyle='-', lw=4, label='100% PAL'),
    Line2D([0], [0], color=colors['150'], linestyle='-', lw=4, label='150% PAL')
]

plt.legend(handles=legend_elements, loc='best', fontsize=16)

plt.tight_layout()
plt.savefig('Combined_Kasting_WACCM_PT.png', dpi=500)
plt.show()