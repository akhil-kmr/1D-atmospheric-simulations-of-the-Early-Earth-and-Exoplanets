import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

states = ['0.1', '1', '10', '100', '150']

colors = {
    '0.1': '#1f77b4',
    '1':   '#ff7f0e',
    '10':  '#2ca02c',
    '100': '#d62728',
    '150': '#9467bd'
}

def read_atmos_dist(file_path):
    try:
        df = pd.read_csv(file_path, sep='\s+', skiprows=0)
        press_hpa = df['PRESS'] * 1000.0
        temp_k = df['TEMP']
        return press_hpa, temp_k
    except Exception:
        return None, None

plt.figure(figsize=(10, 12))

global_p_max = -float('inf')
global_p_min = float('inf')

for state in states:
    color = colors[state]
    atmos_file = os.path.join(f"{state}pc", "SZA_48.2", "PTZ_mixingratios_out.dist")
    
    if os.path.exists(atmos_file):
        press_atmos, temp_atmos = read_atmos_dist(atmos_file)
        if press_atmos is not None:
            plt.plot(temp_atmos, press_atmos, linestyle='-', color=color, linewidth=2.5)
            global_p_max = max(global_p_max, press_atmos.max())
            global_p_min = min(global_p_min, press_atmos.min())

    csv_file_pal = f"waccm6_{state}%PAL_PT_profile.csv"
    csv_file_alt = f"WACCM6_{state}%_PT_profile.csv"
    csv_file = csv_file_pal if os.path.exists(csv_file_pal) else csv_file_alt
    
    if os.path.exists(csv_file):
        df_waccm = pd.read_csv(csv_file)
        p_col = 'Press [hPa]' if 'Press [hPa]' in df_waccm.columns else 'Pressure'
        t_col = 'Temp [K]' if 'Temp [K]' in df_waccm.columns else 'Temperature'
        press_waccm = df_waccm[p_col]
        temp_waccm = df_waccm[t_col]
        plt.plot(temp_waccm, press_waccm, linestyle=':', color=color, linewidth=2.5)
        global_p_max = max(global_p_max, press_waccm.max())
        global_p_min = min(global_p_min, press_waccm.min())

plt.yscale('log')

if global_p_max != -float('inf') and global_p_min != float('inf'):
    plt.ylim(global_p_max, global_p_min)
else:
    plt.ylim(1013.25, 1e-5)
    

plt.xlabel('Temperature (K)', fontsize=16)
plt.ylabel('Pressure (hPa)', fontsize=16)
plt.xlim(100, 350) 
plt.ylim(0, 0.001)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(True, which="both", ls="-", alpha=0.2)

legend_elements = [
    Line2D([0], [0], color='black', linestyle='-', lw=2.5, label='Atmos (Solid)'),
    Line2D([0], [0], color='black', linestyle=':', lw=2.5, label='WACCM (Dotted)'),
    Line2D([0], [0], color=colors['0.1'], linestyle='-', lw=4, label='0.1% PAL'),
    Line2D([0], [0], color=colors['1'], linestyle='-', lw=4, label='1% PAL'),
    Line2D([0], [0], color=colors['10'], linestyle='-', lw=4, label='10% PAL'),
    Line2D([0], [0], color=colors['100'], linestyle='-', lw=4, label='100% PAL'),
    Line2D([0], [0], color=colors['150'], linestyle='-', lw=4, label='150% PAL')
]

plt.legend(handles=legend_elements, loc='best', fontsize=12, frameon=True)
plt.title('Pressure-Temperature Profiles: Atmos vs WACCM', fontsize=18)
plt.tight_layout()
plt.savefig('atmos_pt_profile.png', dpi=300)
plt.show()