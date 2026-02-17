# run_plot.py
from config import CONFIG
from plot_spectra import plot_spectra

files = {
    'TRAPPIST-1e': [
        "/Users/gregcooke/HWO_project/spectra/Earth_PI_HWO_UV.txt",
        "/Users/gregcooke/HWO_project/spectra/Earth_PI_HWO_VIS.txt",
        "/Users/gregcooke/HWO_project/spectra/Earth_PI_HWO_NIR.txt"
    ],
    'K2-18 b': [
        "/Users/gregcooke/HWO_project/spectra/K2-18b_HWO_UV.txt",
        "/Users/gregcooke/HWO_project/spectra/K2-18b_HWO_VIS.txt",
        "/Users/gregcooke/HWO_project/spectra/K2-18b_HWO_NIR.txt"
    ]
}

plot_config = CONFIG["plot"]
plot_spectra(
    files,
    transit=plot_config["transit"],
    wavelength_range=plot_config["wavelength_range"],
    output_file=plot_config["output_file"]
)
