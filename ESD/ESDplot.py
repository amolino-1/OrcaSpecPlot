import argparse
import numpy as np
import pandas as pd
from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator)

# Initialize parser
parser = argparse.ArgumentParser()
# Add optional arguments
parser.add_argument("-S", "--spectrumfile",  # Option keyword
                    type=argparse.FileType('r'),  # Open and read the file specified by the argument
                    required=True,  # Indicate whether an argument is required or optional
                    help='ESD output <basename>.spectrum',  # Help message
                    metavar='\b')  # Alternate display name in help message
parser.add_argument("-o", "--output",
                    type=str,
                    nargs='?',  # Number of times argument can be used, '?' = single value, which can be optional
                    default='ESD_spec',  # Default value used when argument is not provided
                    help='plot output name (without file type)')
parser.add_argument("-e", "--energyshift",
                    type=float,
                    nargs='?',
                    default=0.0,
                    help='Energy correction (float) in eV')
parser.add_argument("-xmin",
                    type=float,
                    nargs='?',
                    default=200,
                    help='x-axis MINIMUM value in wavelength')
parser.add_argument("-xmax",
                    type=float,
                    nargs='?',
                    default=1000,
                    help='x-axis MAXIMUM value in wavelength')

# Read arguments from commandline and assign to variables
args = parser.parse_args()
spectrum = args.spectrumfile
energycorrection_ev = args.energyshift
xmin = args.xmin
xmax = args.xmax
output = args.output


# Convert cm to inches as figsize input only accepts Imperial units
def cm_to_inch(value):
    return value / 2.54


# Apply energy correction in eV and then convert back to wavelength
def energy_correction(value):
    electronvolt = 1239.8424122  # Planck constant (eV•s) * Speed of Light (m/s)
    energy_in_ev = (electronvolt / value) + energycorrection_ev
    return electronvolt / energy_in_ev


# Create Pandas dataframe (separated by whitespace)
dataframe = pd.read_csv(spectrum, delim_whitespace=True)
# Remove all empty columns
dataframe.dropna(how='all', axis=1, inplace=True)

# Find λmax (wavelength)
# Get λmax from uncorrected spectrum
indx_max = np.argmax(dataframe['TotalSpectrum'])
x_max = dataframe.iloc[indx_max]
x_max_energy = x_max['Energy']
# Get λmax from energy corrected spectrum
x_corrected_max_energy = energy_correction(x_max_energy)

# Set Subplots, figure size and font size
fig, ax1 = plt.subplots(figsize=(cm_to_inch(12.5), cm_to_inch(6.99)), dpi=cm_to_inch(600))
plt.rcParams["font.size"] = 10

# Plot Data
# Plot NORMALISED ESD Spectrum
ax1.plot(dataframe['Energy'], dataframe['TotalSpectrum'] / max(dataframe['TotalSpectrum']),
         color='Black',
         linestyle='-',
         linewidth=1.0,
         marker='',
         label=f'ESD Spec, λ$_{{max}}$ = {x_max_energy:.1f} nm')
# Plot NORMALISED Energy Corrected ESD Spectrum
ax1.plot(energy_correction(dataframe['Energy']), dataframe['TotalSpectrum'] / max(dataframe['TotalSpectrum']),
         color='Red',
         linestyle='--',
         linewidth=1.0,
         marker='',
         label=f'{energycorrection_ev} eV Correction, λ$_{{max}}$ = {x_corrected_max_energy:.1f} nm')
# Set axis labels
ax1.set_xlabel('Wavelength / nm')
ax1.set_ylabel('Normalised Intensity')
# Set axis range
# xmin and xmax optional arguments, default xmin = 200 nm, xmax = 1000 nm
ax1.set_xlim([xmin, xmax])
# ax1.xaxis.set_major_locator(MultipleLocator(100))
# ax1.xaxis.set_minor_locator(MultipleLocator(25))
ax1.set_ylim([0, 1.2])
ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Plot vertical lines for λmax
# Plot vertical line for λmax of uncorrected spectrum
ax1.vlines(x=[x_max], ymin=[0], ymax=[1],
           colors=['Black'], linestyle='-', linewidth=1)
# Plot vertical line for λmax of energy corrected spectrum
ax1.vlines(x=[x_corrected_max_energy], ymin=[0], ymax=[1],
           colors=['Red'], linestyle='--', linewidth=1)

# Plot legend
plt.legend(loc='upper right',
           frameon=True,
           edgecolor='Black',
           framealpha=0.7,
           handlelength=1.5,
           fontsize=6)

plt.tight_layout()
plt.savefig(f'{output}.png',
            bbox_inches='tight')

print(f"FINISHED - ESD Plot saved to {output}.png")
