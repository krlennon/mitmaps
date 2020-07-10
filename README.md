# MITMAPS

An open-source software package developed for fast analysis of rheological data obtained using the multi-tone MAPS experimental protocol. See the following paper for details about the MAPS experimental protocol:

K. R. Lennon, M. Geri, G. H. McKinley, and J. W. Swan, ["Medium Amplitude Parallel Superposition (MAPS) Rheology, Part 2: Experimental Protocols and Data Analysis"](https://arxiv.org/abs/2006.09465), arXiv preprint (2020).

## Compatibility

MITMAPS has been tested using [Python](https://www.python.org/downloads/) 3.7.6 with [NumPy](https://numpy.org/install/) 1.18.1, [Matplotlib](https://matplotlib.org/downloads.html) 3.1.3, and [SciPy](https://www.scipy.org/install.html) 1.4.1.

The data used to test this software was obtained from the TRIOS v5.0.0 software supplied by TA Instruments. File I/O details are confined to the source file `mapsio.py` and depend on the specific file format exported by the TRIOS software. The remainder of the software functions independently from the specifics of the input data file format.

## Installation

To install the MITMAPS software package, download all of the files in the `src` directory and place these files in the same local directory.

## Usage

### Specifying Input Data Formats

The MITMAPS software requires `.txt` file formats for input data files. Small amplitude oscillatory shear (SAOS) data should be exported from Oscillatory Frequency Sweep experiments in TRIOS, and multi-tone MAPS data should be exported in a separate file from either Arbitrary Wave or Multiwave experiments in TRIOS.

The SAOS data file should contain separate columns for the oscillation frequency (rad/s), storage modulus (Pa), and loss modulus (Pa). The default assumption is that the frequency resides in the fourth column, the storage modulus in the first column, and the loss modulus in the second column. The [example SAOS data file](Example Data/saos.txt) demonstrates this format. Any adjustments to the input file structure can be made by editing the `set_configuration` function in the file `configuration.py`. For example, if the frequency, storage modulus, and loss modulus are in the first, second, and third columns, respectively, the following modifications should be made:

```python
w_col = 1
Gprime_col = 2
Gdoubprime_col = 3
```

Similarly, modification to the input units can be made in the `configuration.py` file. For example, a frequency reported in units of 1/s and moduli in units of kPa are reflected by the unit conversion factors:

```python
w_units = 2*np.pi
G_units = 1000
```

The MAPS data file(s) should contain separate columns for the time (s), stress (Pa), and strain (%). By default, it is assumed that for strain-controlled experiments, these quantities reside in the first, second, and third columns, and for stress-controlled experiments in the fourth, sixth, and fifth columns, respectively. These presets are reflected in the `configuration.py` file and can be adjusted as needed with the following lines:

```python
if control == "strain" or control == None:
	t_col = 1
	stress_col = 2
	strain_col = 3
elif control == "stress":
	t_col = 4
	stress_col = 6
	strain_col = 5
```

Similarly, if the units of these quantities are different than the presents, such as time in ms, stress in kPa, and strain in strain units, unit conversion factors can be specified:

```python
t_units = 0.001
stress_units = 1000
strain_units = 100
```

Different Arbitrary Wave or Multiwave experiments in the same MAPS data file should be separated by the line `[step]` and labelled `Sine Strain - n` or `Arbitrary Wave - n`. These experiments should be organized hierarchically by the input tone set, amplitude, and frequency, in which case the `sort_order` variable in the file `MITMAPS.py` should be set to `"amplitude"`; or by the input tone set, frequency, and amplitude, in which case `sort_order` should be set to `"frequency"`. The [example MAPS data file](Example Data/MAPS Data/maps.txt) demonstrates the appropriate format for a MAPS data file with `sort_order = "amplitude"`.

Multiple MAPS data files may be anaylzed simultaneously as a part of the same MAPS frequency sweeps(s). All MAPS data files to be analyzed simultaneously by the software should be placed in a dedicated directory, even if only one MAPS data file exists.

### Example `MITMAPS.py` File

After configuring the input data file settings, users should modify the beginning lines of the file `MITMAPS.py`. Below is an example of these lines, which when executed will input the example [SAOS](Example Data/saos.txt) and [MAPS](Example Data/MAPS Data/maps.txt) data files, analyze the data, and generate Bode plots and Nyquist diagrams the third order complex viscosity in MAPS frequency sweeps with `[n1,n2,n3] = [5,6,9]` and `[1,4,16]`.

```python
# Select "experimental" mode for experimental data analysis, or "simulation" mode for model predictions/simulations only
mode = "experimental"

# Linear Response
LR_file = "../Example Data/saos.txt"
LR_fit = "Maxwell"

# Experimental mode (either "stress" or "strain")
MAPS_control = "stress"

# MAPS response
MAPS_folder = "../Example Data/MAPS Data"
MAPS_tones = [[5,6,9],[1,4,16]]
MAPS_freqs = [1.28,0.64,0.32,0.16]
sort_order = "amplitude"
plot_var = "eta"

# Constitutive models
full_model = None
maps_models = [crm_eta3]
extra_params = []

# Additional options
plotLR = False
gapLoading = False
outputTable = False
tssComp = False
```

See the [documentation](docs/settings.md) for more details on the settings and options within this block.

## Contributing
Inquiries about the MITMAPS software can be directed to krlennon[at]mit.edu.

## License
[GNU General Public License v3.0](https://choosealicense.com/licenses/gpl-3.0/)
