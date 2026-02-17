import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

# ------------------------
# Define filter
# ------------------------

def filter_molecules_by_wavelength(molecules_dict, wavelength_range):
    wmin, wmax = wavelength_range
    filtered = {}

    for planet, mol_list in molecules_dict.items():
        new_list = []
        for mol in mol_list:
            xmin, xmax = mol['gradient']['xrange']

            # Keep molecule if its band overlaps plotting range
            if xmax >= wmin and xmin <= wmax:
                new_list.append(mol)

        if new_list:
            filtered[planet] = new_list

    return filtered


# ------------------------
# Define thick axes
# ------------------------

def thick_axes(top = False, labelleft = True, direction = 'in'):
    # Accessing the axes object and setting linewidth
    plt.gca().spines['top'].set_linewidth(2)  # Top axis
    plt.gca().spines['bottom'].set_linewidth(2)  # Bottom axis
    plt.gca().spines['left'].set_linewidth(2)  # Left axis
    plt.gca().spines['right'].set_linewidth(2)  # Right axis
    plt.tick_params(which = 'major', axis = 'both', direction = direction, labelsize = 15, length = 6, width = 2, labelleft = labelleft, right = True, top = top)
    plt.tick_params(which = 'minor',axis = 'both', direction = direction, labelsize = 15, length = 3, width = 1, labelleft = labelleft, right = True, top = top)


# ------------------------
# Helper: subscript and superscript
# ------------------------
def sub(x):
    return f"$_{{{x}}}$"
def sup(num):
    return r'$^{'+str(num)+'}$'

# ------------------------
# Gradient helper
# ------------------------
def gradient_image(ax, extent, transform, cmap, cmap_range, alpha):
    vmin, vmax = cmap_range
    gradient = np.linspace(vmin, vmax, 256).reshape(-1, 1)
    gradient = np.flipud(gradient)
    ax.imshow(
        gradient,
        aspect='auto',
        extent=extent,
        transform=transform,
        cmap=cmap,
        alpha=alpha,
        interpolation='bicubic'
    )
    
# ------------------------
# Molecule label + gradient helper
# ------------------------
def add_molecule_labels(ax, molecules, xmin, xmax):
    """
    Add molecular labels and gradient overlays.

    Parameters
    ----------
    ax : matplotlib axis
    molecules : list of dict
        Each dict contains:
        - 'label': str, e.g., "CH4"
        - 'x': float, wavelength in μm
        - 'y': float, vertical position in plot units
        - 'color': str, matplotlib color
        - 'gradient': dict (optional) with keys:
            - 'xrange': tuple (xmin, xmax)
            - 'cmap': colormap
            - 'alpha': float
    xmin, xmax : float
        Wavelength plot range
    """
    x_total = xmax - xmin
    for mol in molecules:
        # Add text
        ax.text(mol['x'], mol['y'], mol['label'], fontsize=15,
                color=mol.get('color','black'),
                ha=mol.get('ha','center'))
        
        # Add gradient if specified
        if 'gradient' in mol:
            g = mol['gradient']
            extent = ((g['xrange'][0]-xmin)/x_total,
                      (g['xrange'][1]-xmin)/x_total,
                      0, 1)
            gradient_image(ax, extent=extent, transform=ax.transAxes,
                           cmap=g['cmap'], cmap_range=(0,1),
                           alpha=g.get('alpha',0.2))

# ------------------------
# Load & combine multiple spectra
# ------------------------
def load_multi_spectrum(file_list):
    arrays = [np.genfromtxt(f, comments="#") for f in file_list]
    Spec = np.vstack(arrays)
    Spec = Spec[np.argsort(Spec[:,0])]
    return Spec

# ------------------------
# Main plotting function
# ------------------------
def plot_spectra(files, figsize = (10,5),
                 transit=True, wavelength_range=(0.2,2.0),
                 molecule_options=None, output_file=None,
                 xscale = 'log'):

    plt.figure(figsize=(figsize))
    gs = gridspec.GridSpec(len(files),1, wspace=0.3)

    for i, (name, filelist) in enumerate(files.items()):
        Spec = load_multi_spectrum(filelist)
        

        S = Spec[:,3]
        P = Spec[:,4]

        contrast = P / S

        ax = plt.subplot(gs[i,0])
        ax.set_ylim(0, 80)
        ax.set_xlim(*wavelength_range)
        ax.set_xscale(xscale)
        ax.set_ylabel('Effective altitude [km]', fontsize=15, weight='bold')
        if i == len(files)-1:
            ax.set_xlabel('Wavelength [μm]', fontsize=15, weight='bold')

        # Plot main curve
        if transit:
            ax.plot(Spec[:,0], Spec[:,5], lw=2, color='#0b71bf', label=name)
            print(Spec[:,0])
            print(Spec[:,5])
        else:
            ax.plot(Spec[:,0], contrast, lw=2, color='#1dccb8', label=name)

        # Molecules / gradients
        if molecule_options and name in molecule_options:
            add_molecule_labels(ax, molecule_options[name],
                                xmin=wavelength_range[0],
                                xmax=wavelength_range[1])
            
        ax.legend(loc=4, fontsize=15)
    
    thick_axes(top = True, labelleft = True, direction = 'in')

    if output_file:
        plt.savefig(output_file, dpi=200, bbox_inches='tight')
    plt.show()

