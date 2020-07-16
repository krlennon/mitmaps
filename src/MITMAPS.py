##  Copyright (C) 2020 Kyle R. Lennon, Michela Geri, Gareth H. McKinley, James W. Swan
##
##  This file is part of MITMAPS.
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

from mapsio import *
from process import *
from mapsplot import *
from mapsmodels import *
from mapssim import *
from models import *
import sys

##############################
## EDIT THE FOLLOWING LINES ##
##############################

# Select "experimental" mode for experimental data analysis, or "simulation" mode for model predictions/simulations only
# Note that simulations are only conducted in strain control (MAPS_control will be set to "strain" by default)
mode = "experimental"

# Linear response
LR_file = "../Example Data/saos.txt" # Path to file with linear response data (set to None to skip)
LR_fit = "Maxwell" # Fitting method for LR data ("Maxwell" or "interpn" with n = 1, 2, or 3)

# Experimental mode (either "stress" or "strain")
MAPS_control = "stress"

# MAPS response
MAPS_folder = "../Example Data/MAPS Data" # Path to folder with MAPS data
MAPS_tones = [[5,6,9],[1,4,16]] # Input tone sets for MAPS signals
MAPS_freqs = [1.28,0.64,0.32,0.16] # Fundamental frequencies in the MAPS sweeps
sort_order = "amplitude" # Outer sorted variable ("amplitude" or "frequency")
plot_var = "eta" # MAPS response function to plot ("G", "eta", "J", or "phi")

# Constitutive models
full_model = None # the constitutive model to simulate in "simulation" mode (set to None to plot only analytical MAPS solution)
maps_models = [crm_eta3] # MAPS model to plot (specific to MAPS response function)
extra_params = [] # Parameters in addition to those regressed from LVE response

# Additional options
plotLR = False # Choose whether to plot the linear response (both complex modulus and complex compliance)
gapLoading = False # Choose whether to plot the gap loading limit analysis (only if LVE data + fit provided)
outputTable = False # Choose whether to output the MAPS response functions values obtained by the experiments as a .csv file
tssComp = False # Choose whether to run a TSS comparison

###################################
## EDITS BELOW HERE ARE OPTIONAL ##
###################################

# Check for consistency in inputs
if (plot_var == "J" or plot_var == "phi") and MAPS_control == "strain" and LR_file == None and mode == "experimental":
    sys.exit("Error: MAPS response for " + plot_var + " cannot be determined from strain-controlled data alone. To compute "
          + plot_var + ", include stress-controlled data or linear response data.")

elif (plot_var == "G" or plot_var == "eta") and (MAPS_control == "stress" or MAPS_control == "maostress") and LR_file == None and mode == "experimental":
    sys.exit("Error: MAPS response for " + plot_var + " cannot be determined from stress-controlled data alone. To compute "
          + plot_var + ", include strain-controlled data or linear response data.")
elif (LR_file == None) and plotLR and mode == "experimental":
    sys.exit("Error: No linear response data provided for plotting. Provide linear response data or set plotLR to False")
elif (LR_file == None) and gapLoading and mode == "experimental":
    sys.exit("Error: No linear response data provided for gap loading limit analysis. Provide linear response data or set gapLoading to False")
elif (LR_file == None) and tssComp and mode == "experimental":
    sys.exit("Error: No linear response data provided for TSS comparison. Provide linear response data of set tssComp to False")

# Read LR data and fit to a single Maxwell mode
params = extra_params
if LR_file != None:
    sys.stdout.write("\rProcessing linear response data ...")
    LR_exp_types, LR_exp_data = readfile(LR_file)
    LR_data = getLR(LR_exp_data, LR_exp_types)
    sys.stdout.write("\rProcessing linear response data ... done\n")

    # Fit the LR data to a single-mode Maxwell model
    if LR_fit == "Maxwell":
        sys.stdout.write("\rFitting linear response data ...")
        eta0, lam, etainf = fitLR(LR_data)
        LR_model = lambda w: maxwellLR(w,[eta0,lam,etainf])
        sys.stdout.write("\rFitting linear response data ... done\n")

        # Set the model parameters to include the linear response parameters
        params = [eta0,lam,etainf] + extra_params

    # Fit the LR data to kth order splines
    elif "interp" in LR_fit:
        LR_model = interpolateLR(LR_data,int(LR_fit[6]))

    # If maps_model is tss, set LR_model
    if tss_eta3 in maps_models:
        ind = maps_models.index(tss_eta3)
        maps_models[ind] = lambda w1,w2,w3,params: params[0]*tss_eta3(w1,w2,w3,LR_model)
    elif tss_J3 in maps_models:
        ind = maps_models.index(tss_J3)
        maps_models[ind] = lambda w1,w2,w3,params: params[0]*tss_J3(w1,w2,w3,LR_model)

    # Plot linear response data if specified
    if plotLR:
        plot_linear_response(LR_data,LR_model)

    # Plot gap loading limit analysis if specified
    if gapLoading:
        wmin = np.min(MAPS_freqs)
        wmax = 3*np.max(np.max(np.array(MAPS_tones)))*np.max(MAPS_freqs)
        gap_loading(LR_model,1000,[wmin,wmax])

# Input the MAPS data in "data" mode
if mode == "experimental" and MAPS_folder != None:
    # Read experimental MAPS data
    sys.stdout.write("\rProcessing MAPS data ...")
    exp_types, exp_data = readfolder(MAPS_folder)

    # Sort the data with sortdata(data, types, input sets, input fundamentals)
    tones = ['_'.join([str(elem) for elem in tone_set]) for tone_set in MAPS_tones]
    freqs = [str(elem) for elem in MAPS_freqs]
    data = sortdata(exp_data, exp_types, tones, freqs, sort_order)

    # Build experiment objects (make_experiments is very slow (performance limiting) due to i/o structure)
    # Arguments are sorted data (data) and the control protocol (MAPS_control)
    experiments = make_experiments(data, MAPS_control)

# Simulate a model response in "simulation" mode
if mode == "simulation" and full_model != None:
    sys.stdout.write("\rSimulating MAPS experiments ...")
    # Simulate MAPS experiments with a constitutive model (no simulated noise)
    experiments = []
    for tone_set in MAPS_tones:
        experiments += sim_multitone_eq(full_model,params,0.2,tone_set,MAPS_freqs,[0,0])

    # Set MAPS_control to "strain" - only currently available setting for simulations
    MAPS_control = "strain"

# Process the experiments for "data" mode and "simulation" mode
if (mode == "experimental" and MAPS_folder != None) or (mode == "simulation" and full_model != None):
    # Process the data for each experiment
    for experiment in experiments:
        # Preprocess the data
        experiment.trim(ncycles=5)
        experiment.mean_subtract()

        # Process the data
        experiment.get_MAPS()

        # Convert between stress and strain control (requires LR data)
        if LR_file != None:
            experiment.interconvert(LR_model)

    # If TSS comparison is selected, make the plot
    if tssComp:
        tss_comparison(experiments,LR_model)

    # Set overlap to 1 and display messages
    overlap = 1
    if mode == "experimental":
        sys.stdout.write("\rProcessing MAPS data ... done\n")
    elif mode == "simulation":
        sys.stdout.write("\rSimulating MAPS experiments ... done\n")

    # Make figures from the experiments
    # Arguments are the experiment array, the MAPS function to plot, and the verbosity (0 - don't display or save,
    # 1 - display, don't save, 2 - display and save)
    sys.stdout.write("\rPlotting MAPS data ...")
    figures = make_figures(experiments, plot_var, verbose=0)

# If no data or full model are provided, initialize figures to plot only the analytical MAPS prediction
else:
    figures = {}
    for tone_set in MAPS_tones:
        # Generate the figure name
        name = '_'.join(map(str, tone_set))

        # Get the MAPS coordinates for the input tone set
        sums, mapscoords = get_maps_coords(tone_set)

        # Remove repeated sums
        rem = []
        k = 0
        for nset in mapscoords:
            # Determine if a sum is repeated
            c = list(sums).count(sum(nset))
            if c > 1:
                rem += [k]

            # Correct the signs to reflect the conventions of the subspaces
            if sum(np.sign(nset)) == -1:
                mapscoords[k,:] = -1*mapscoords[k,:]
            k += 1

        # Remove repeated sums
        rem.reverse()
        for k in rem:
            mapscoords = np.delete(mapscoords,k,0)

        # Make the figure object
        figures[name] = MAPSFigure(name,0,mapscoords,[],[])

        # Set overlap to 0
        overlap = 0

# Plot a model prediction
# Arguments are the MAPS model solution, model parameters, frequency range, verbosity
# (0 - don't display or save, 1 - display, don't save, 2 - display and save), and whether to overlap with data
# (1) or make a separate figure (0)

# First set the minimum and maximum frequency to consider
all_tones = [tone for tone_set in MAPS_tones for tone in tone_set]
wmin = np.min(MAPS_freqs)
wmax = np.max(MAPS_freqs)

# Plot model predictions for each figure object
if maps_models != None:
    for figure in figures:
        figures[figure].plot_model(maps_models, params, [wmin,wmax], overlap=overlap, verbose=0)

sys.stdout.write("\rPlotting MAPS data ... done\n")

# Output the data to a CSV file
if outputTable:
    output_MAPS(experiments)

plt.show()
