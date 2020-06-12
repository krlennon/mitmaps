##  MITMAPS - Software for processing MAPS rheological data
##  Copyright (C) 2020 Kyle R. Lennon, Michela Geri, Gareth H. McKinley, James W. Swan
##
##  MITMAPS is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  MITMAPS is distributed in the hope that it will be useful
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with MITMAPS.  If not, see <https://www.gnu.org/licenses/>.

MITMAPS is an open-source software package developed for fast analysis of rheological data obtained using the MAPS protocol (see "Medium Amplitude Parallel Superposition Rheology, Part 2: Experimental Protocols and Data Analysis" by Kyle R. Lennon, Gareth H. McKinley, and James W. Swan). We hope to actively maintain this software at <https://github.com/krlennon/mitmaps>. If issues with the version of this software downloaded from the Supporting Information of the aforementioned paper are encountered, we encourage users to update to the latest version published on GitHub. Inquirites about the software can be directed to krlennon[at]mit.edu.

=============
Compatibility
=============

This software has been tested using Python 3.7.6, Matplotlib 3.1.3, and SciPy 1.4.1. The software uses Python3 specific language. See <https://www.python.org/downloads/> for information on downloading the latest release of Python, <https://matplotlib.org/downloads.html> for information on downloading the latest release of Matplotlib, and <https://www.scipy.org/install.html> for information on downloading the latest release of SciPy.

All data used to test this software was obtained from the TRIOS v5.0.0 software. The details of file I/O for the software package are contained in the file "mapsio.py" and depend on the specific file format exported by the TRIOS software. The remainder of the software functions independently of the specific input file formats.

============
Installation
============

To install the MITMAPS sofware package, download the files "configuration.py", "mapsio.py", "mapsmodels.py", "mapsplot.py", "mapssim.py", "models.py", "process.py", and "MITMAPS.py", and place all files in the same directory. 

=====================
Required Data Formats
=====================

If the format of the input files match those of the examples provided with this software (see "saos.txt" and "maps.txt"), then it should not be necessary to edit any portion of the software outisde of the highlighted section in the "MITMAPS.py" file in order to work the software.

These example files were obtained by exporting the raw time-series data from an Arbitrary Wave experiment in the TRIOS software provided by TA Instruments ("maps.txt") and by exporting the processed data from an Oscillatory Frequency Sweep experiment in the TRIOS software provided by TA Instruments ("saos.txt"). 

It is important that any data used to obtain a fit of a material's LVE behvaior be input in this format, and as a separate file from the MAPS data. The software will read the firt "Frequency Sweep" experiment contained in the specified linear response file (NOTE: The SAOS experiment must be labelled "Frequency sweep - n" in the input file so that it can be detected by the software). The software assumes that the first column in the Frequency Sweep experiment is the "Storage Modulus" in units of Pa, that the second column is the "Loss modulus" in units of Pa, and that the fourth column is the "Angular frequency" in units of rad/s. The content of the remaining columns is unimportant. If these details and ONLY these details are different in a user's file, the software can be adjusted simply by editing the relevant lines of the "configuration.py" file.

The software is able to consolidate multiple files containing MAPS data into a single data set. Therefore, all files containing MAPS data should be collected in a single directory, and this directory will be specified in the "MITMAPS.py" file (implementation details can be found below in the "Implementation" section). Even if only one file containing MAPS data is used, this file should be placed in a directory. The directory containing MAPS data files should contain no other files.

The software will detect different MAPS experiments within the same data file via the line "[step]" which should separate each experiment in the data file. The data file may contain experiments (such as Amplitude Sweeps or Frequency Sweeps) that precede the MAPS experiments. All experiments up until the first experiment labelled "Sine Strain - n" or "Arbitrary Wave - n" will be ignored by the software. After the first experiment labelled "Sine Strain - n" or "Arbitrary Wave - n" is detected, it and all remaning experiments in the data file will be treated by the software as MAPS experiments. If the experiment is a stress-controlled MAPS experiment, the software expects that the fourth column contains the time points ("Step time") in units of s, that the fifth column contains the "Strain" in %, and that the sixth column contains the "Stress" in units of Pa. If the experiment is a strain-controlled MAPS experiment, the software assumes that the first column contains the time points in units of s, that the second column contains the "Stress" in units of Pa, and that the third column contains the "Strain" in %. If these details (the columns in which different quantities appear, and the units of these quantities) and ONLY these details are different in a user's file, the software can be adjusted simply by editing the relevant lines of the "configureation.py" file.

==============
Implementation
==============

To facilitate the use of MAPS rheology and the MITMAPS software, we have endeavored to reduce the amount of modifications necessary to adapt the software to new data to a minimum. The details of the data analysis are encapsulated primarily by the files "process.py" and "mapsplot.py". We use the methods in these files to provide a level of abstraction to many computational details. The file "MITMAPS.py" presents the abstracted workflow of the software in a more straighforward, proceeding through the steps of data input, data processing, and plotting. We hope that this abstraction provides a simplified interface to the MAPS data for future explorations in data analysis, and additional flexibilities to suit the interests of different users. Users interested in such adaptations should see the below section on "Backend Details".

To use the MITMAPS software to produce the set of Bode and Nyquist plots shown in Lennon et. al., users need only modify a highlighted block of lines at the beginning of the "MITMAPS.py" file. The purpose of each line is documented within that file. A more detailed explanation is presented below:

mode - the mode in which the software should operate. "experimental" mode indicates that MAPS data will be provided, and can be plotted against MAPS constitutive model solutions. "simulation" mode should be selected to either plot only model predictions (for any model in "mapsmodels.py", or to simulate MAPS data for a differential constitutive model (for any model in "models.py").

LR_file - The path to the .txt file containing the linear response (SAOS) data -- see the "Required Data Formats" section for more information on the required format of this .txt file. If no linear response data is available, the user may skip features involving LVE data by setting LR_file to the value None. NOTE: interpolation of linear response data is necessary for certain data processing steps. Currently, the software will interpolate by fitting the linear response data to a single-mode Maxwell model plus a Newtonian solvent mode. More options for linear models will be added in future updates. Can be specified in either "experimental" or "simulation" mode.

MAPS_control - The controlled variable (stress or strain) for the MAPS experiments. For the present release, this should be set to either "stress" or "strain", for stress- or strain-controlled experiments, respectively. Only needed for "experimental" mode (all "simulation" mode simulated MAPS experiments are run in strain control).

MAPS_folder - The path to the directory containing all .txt files with MAPS data. This must always be a directory, not a single file, even if only one MAPS data file is available. See "Required Data Formats" for more details on the required formats of the .txt files. NOTE: in reading multiple .txt files contained within this directory, the folders will be ordered by name. Be sure to check that there are no hidden files within this directory, as these will be read by the software causing runtime errors. Only needed for "experimental" mode.

MAPS_tones - The sets of input tones {n1,n2,n3} used for each MAPS experiment (either experimental or simulated). These should be specified as a list of lists, with each input tone set as an inner list. For example, if the data sets contained in MAPS_folder represent only MAPS experiments with {n1,n2,n3} = {5,6,9}, then MAPS_tones should be set to [[5,6,9]]. If the data sets represent experiments with multiple tone sets, for example {5,6,9} and {1,4,6}, MAPS_tones should contain each set as a list, e.g. [[5,6,9],[1,4,6]]. In "experimental" mode, when multiple nsets are present, they should be listed in the order in which they first appear in the .txt file(s) (see MAPS_folder documentation for how multiple files in the same directory are ordered).

MAPS_freqs - The fundamental frequencies w0 used in the MAPS experiments (either experimental or simulated). For "experimental" mode, these frequencies should be listed in the order in which they appear in the .txt file(s) (see MAPS_folder documentation for how multiple files in the same directory are ordered). See sort_order documentation below for more details on the meaning of ordering in these .txt files. This variable should be a list, e.g. [1.28, 0.64, 0.32, 0.16] for the frequencies 1.28 rad/s, 0.64 rad/s, 0.32 rad/s, and 0.16 rad/s. The range over which analytical MAPS solutions to constitutive models (those specified in maps_models) is plotted depends on the maximum and minimum value contained in MAPS_freqs.

sort_order - The outer sorted variable in the .txt files contained in MAPS_folder. The Arbitrary Wave experiments are assumed to be sorted hierarchically, with the input tone set on the highest level always. Within a given input tone set, the next level of sorting determines the value of sort_order. If the next highest level of sorting is the amplitude of the imposed multi-tone signal, then sort_order should be set to "amplitude". If the next highest level is insteady the fundamental frequency, then sort_order should be set to "frequency". For example, the following order of experiments with different amplitudes (g) and frequencies (w0):
	- g = 1, w0 = 1
	- g = 1, w0 = 2
	- g = 2, w0 = 1
	- g = 2, w0 = 2
corresponds to sort_order = "amplitude". The example data file provided with this software also corresponds to an "amplitude" sorting. Only needed for "experimental" mode.

plot_var - The MAPS response function to be plotted on Bode and Nyquist plots. The value can be set to either "G", "eta", "J", or "phi". Note that, if stress-controlled data is provided, then linear response data must also be provided to convert to "G" or "eta", and that if strain-controlled data is provided, then linear response data must also be provided to convert to "J" or "phi". Should be specified in both "experiment" and "simulation" modes.

full_model - The constitutive model that will be simulated in "simulation" mode. See the file "models.py" for a full list of presently available constitutive models for simulation. The simulation will be run in strain control with the input tone sets specified by MAPS_tones and the fundamental frequencies specified in MAPS_freqs. The parameters for this constitutive model may be specified by the LVE data in LR_file plus any extra nonlinear parameters in extra_params, or solely by parameters in extra_params (if LR_file = None). The MAPS response function that will be displayed by the simulation is still set by plot_var. Not needed for "experimental" mode.

models - Any constitutive models to be plotted on the same axes as the experimental MAPS data. See the file "mapsmodels.py" for options. Note that the models specified are specific to the MAPS material function, therefore should be specified consistently with plot_var. These must be specified as a list, even if only one model is selected. To exclude any model predictions from the plots, set to None. Can be specified in either "experimental" or "simulation" mode.

extra_params - Any extra parameters that should be specified to parameterize the constitutive models, in addition to those fit from the linear response data. If no linear response data is provided, these should include parameters that would have otherwise been fit from the linear response data. For example, if LR_file is provided, and the Co-rotational Maxwell model is specified, then extra_params should be set to a blank list: []. If LR_file is provided, and the Giesekus model is specified, then extra_params should be set to a list containing the value of \alpha: [\alpha]. If LR_file is not provided and the Giesekus model is specified, then extra_params should be set to a list containing \eta_0, \lambda_1, and \alpha: [\eta_0, \lamdba_1, \alpha]. NOTE: due to the structure of extra_params, models with incompatible nonlinear parameters (e.g. Phan-Thien and Tanner and Giesekus) cannot be simultaneously specified.

plotLR - Boolean value specifying whether or not to plot the linear response (must be set to False if LR_file is set to None). If True, then plots of the the complex modulus and complex compliance (storage and loss components vs. frequency) are generated and displayed when the script has completed execution.

gapLoading - Boolean value specifying whether or not to plot the gap loading limit anlaysis (must be set to False if LR_file is set to None). If True, then a plot of the wavelength (\lambda_s) and penetration depth (d_s) as a function of frequency are plotted (following the analysis in Appendix C of Lennon et. al.) and the minimum values recorded for both linear rheological measurements and MAPS measurements. These values correspond to the gap loading limit in each regime.

outputTable - Boolean value specifying whether to output all obtained values of each MAPS response function as a .csv file. If True, the .csv file will be saved under the name "output.csv", and will contain as the first three columns the three MAPS frequency coordinates (w1,w2,w3), and as the remaining columns the values of the third order MAPS functions associated with those coordinates, in SI units.

tssComp - Boolean value specifying whether to create a set of parity plots comparing all obtained values of the MAPS response function to predictions obtained by the assumption of time-strain separability. Must be set to false if LR_file is set to None.
