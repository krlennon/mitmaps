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

import os
import pandas as pd
import numpy as np
from configuration import *
from process import *

def readfolder(foldername):
    ''' Read the contents of a folder containing MAPS experiment files.
    Inputs:
        foldername - name (including path) of folder to be read
    Outputs:
        exp_types - type of each experiment
        exp_data - data from each experiment
    '''

    # Get file names
    files = [f for f in os.listdir(foldername) if os.path.isfile(os.path.join(foldername,f))]
    files.sort()

    # For every file, call readfile and append outputs
    exp_types = []
    exp_data = []
    for f in files:
        if f[0] == '.':
            continue
        elif f[-4:] != ".txt":
            continue
        new_types, new_data = readfile(os.path.join(foldername,f))
        exp_types += new_types
        exp_data += new_data

    return exp_types, exp_data

def readfile(filename):
    """
    Read the contents of an ARES G2 text file for a MAPS frequency sweep.
    Inputs:
        filename - name (including path) of file to be read
    Outputs:
        exp_types - type of each experiment
        exp_data - data from each experiment
    """

    # Open the file and read lines
    with open(filename, 'r') as text:
        lines = text.read().splitlines()

    # List the experiments
    indices = [i for i, x in enumerate(lines) if x == '[step]']
    exp_types = [lines[i+1] for i in indices]

    # Get the data associated with each experiment
    n = len(lines)
    next_indices = indices[1:] + [n]
    exp_data = [lines[indices[i]+4:next_indices[i]-1] for i in range(0,len(indices))]
    return exp_types, exp_data

def writefolder(data, parent):
    """
    Write the separate experiments to a heirarchy of files/folders.
    Assumed three amplitudes per fundamental per set.
    Inputs:
        data - the data to be written (dictionary)
        parent - the name (including path) of the parent folder
    Outputs:
        None
    """

    # Create the parent directory
    try:
        os.mkdir(parent)
    except FileExistsError:
        print('Error: Folder already exists')
        return

    # Get the names of the harmonic sets and fundamentals
    sets = data['sets']
    fundamentals = data['fundamentals']

    # Loop over all entries
    for ii in range(0,len(sets)):
        exp_set = data[sets[ii]]
        child = parent + '/' + sets[ii]
        os.mkdir(child)
        for jj in range(0,len(fundamentals)):
            exp_subset = exp_set[fundamentals[jj]]
            path = child + '/' + fundamentals[jj]
            os.mkdir(path)

            # Write each experiment to a file
            for kk in range(0,3):
                with open(path + '/' + str(kk+1), 'w') as f:
                    for line in exp_subset[kk]:
                        f.write("%s\n" % line)

def sortdata(exp_data, exp_types, sets, fundamentals, sort_order):
    """
    Sort the separated data into a dictionary based on assumed file format.
    Inputs:
        exp_data - the experimental data
        exp_types - the type of each experiment
        sets - the harmonic sets contained in the data set
        fundamentals - the fundamental frequencies in the data set
        sort_order - the outer sorted variable in the dataset
    Outputs:
        data - sorted dictionary
    """

    # Create the dictionary registers
    data = {'sets': sets, 'fundamentals': fundamentals}

    # Find the number of non-MAPS experiments
    for i in range(0,len(exp_types)):
        if 'Sine Strain' in exp_types[i] or 'Arbitrary Wave' in exp_types[i]:
            burn = i
            break

    # Sort through the data
    N = len(exp_types) - burn
    Nset = len(sets)
    Nfun = len(fundamentals)
    Nrep = int(N/(Nset*Nfun))
    for i in range(0,Nset):
        data[sets[i]] = {}
        for j in range(0,Nfun):
            if sort_order == "amplitude":
                n = i*Nfun*Nrep + j + burn
            elif sort_order == "frequency":
                n = i*Nfun*Nrep + j*Nrep + burn
            newdata = []
            for k in range(0,Nrep):
                if sort_order == "amplitude":
                    newdata += [exp_data[n + k*Nfun]]
                elif sort_order == "frequency":
                    newdata += [exp_data[n + k]]
            data[sets[i]][fundamentals[j]] = newdata
    return data

def mergedata(data1, data2):
    """
    Merge two data sets.
    Inputs:
        data1, data2 - data sets (dictionaries) to be merged
    Outputs:
        data - the merged data set
    """

    # Create the new dictionary and assign sets and fundamentals
    data = {'sets': data1['sets'] + data2['sets'], 'fundamentals':
            data1['fundamentals'] + data2['fundamentals']}

    # Add the experiments from the separate data sets to the new data set
    for set1 in data1['sets']:
        for fund1 in data1['fundamentals']:
            data[set1][fund1] = data1[set1][fund1]
    for set2 in data2['sets']:
        for fund2 in data2['fundamentals']:
            data[set2][fund2] = data2[set2][fund2]

    return data

def getLR(exp_data, exp_types):
    """
    Get the linear response data from a data set.
    """

    # Get the configuration
    config = set_configuration(None)
    w_col = config["w_col"]
    Gprime_col = config["Gprime_col"]
    Gdoubprime_col = config["Gdoubprime_col"]
    t_col = config["t_col_LR"]
    w_units = config["w_units"]
    G_units = config["G_units"]

    for i in range(0, len(exp_types)):
        if 'Frequency' in exp_types[i]:
            # Read the data and convert to w, G*
            num_data = [line.split('\t') for line in exp_data[i]]
            w = np.array([w_units*float(line[w_col - 1]) for line in num_data])
            Gr = np.array([G_units*float(line[Gprime_col - 1]) for line in num_data])
            Gi = np.array([G_units*float(line[Gdoubprime_col - 1]) for line in num_data])
            t = np.array([float(line[t_col - 1]) for line in num_data])
            G = Gr + 1j*Gi

            return [w,G,t]

def output_MAPS(experiments,MAPS_folder):
    """
    Output the MAPS data in tabular form
    """

    # Initialize DataFrame
    df = pd.DataFrame([])

    # Append each experiment to DataFrame
    for experiment in experiments:
        subs = []
        RGBs = []
        w1 = []
        w2 = []
        w3 = []
        for ii in range(0,int(experiment.mapscoords.size/3)):
            nset = experiment.mapscoords[ii,:]
            subs += [get_subspace(nset)]
            RGBs += [get_barycentric(nset)]
            nset = list(nset)
            nset.sort()
            w2 += [nset[0]*experiment.base]
            w3 += [nset[1]*experiment.base]
            w1 += [nset[2]*experiment.base]

        G3 = experiment.G3
        eta3 = experiment.eta3
        J3 = experiment.J3
        phi3 = experiment.phi3

        # Create a temporary DataFrame
        temp = pd.DataFrame({"w1 (rad/s)":w1, "w2 (rad/s)":w2, "w3 (rad/s)":w3, "G3 (Pa)":G3, "eta3 (Pa s3)":eta3,
            "J3 (Pa-3)":J3, "phi3 (Pa-3 s-1)":phi3, "Subspace":subs, "[R,G,B]":RGBs})

        # Append to main DataFrame
        df = df.append(temp, ignore_index=True)

    # Write to CSV
    df.to_csv(MAPS_folder + "/output.csv", index=False)

def make_experiments(data, control='strain'):
    """
    From the raw imported data, make an array of Experiment objects. ARES G2 Version.
    Inputs:
        data - data read from file (dictionary)
        control - whether the experiment is 'strain' controlled or 'stress' controlled
    Outputs:
        experiments - array of Experiment objects
    """

    # Get configuration
    configs = set_configuration(control)
    t_col = configs["t_col"]
    stress_col = configs["stress_col"]
    strain_col = configs["strain_col"]
    t_units = configs["t_units"]
    stress_units = configs["stress_units"]
    strain_units = configs["strain_units"]

    # Initialize object list
    experiments = []
    i = 0

    # For all harmonic sets
    for hset in data['sets']:

        # For all fundamental frequencies
        for fund in data['fundamentals']:
            exp_data = data[hset][fund]
            harmonics = [int(n) for n in hset.split('_')]
            fund = float(fund)

            # Create an object
            if control == 'strain':
                experiments += [StrainControlled(fund, harmonics)]
            elif control == 'stress':
                experiments += [StressControlled(fund, harmonics)]
            elif control == 'maostress':
                experiments += [MAOStressControlled(fund, harmonics)]
            elif control == 'micro':
                experiments += [MicroMAPS(fund,harmonics,data["calibration"])]
            else:
                print('Error: Control ' + control + ' is invalid')
                return None

            # For all trials
            n = 0
            for trial in exp_data:
                # This line and the if/else block are performance limiting - fix i/o structure
                if control == 'strain':
                    num_data = [line.split('\t') for line in trial]
                    num_data = [[t_units*float(line[t_col - 1]), stress_units*float(line[stress_col - 1]), strain_units*float(line[strain_col - 1])/100] for line in num_data]

                    # Add the trial to an object
                    experiments[i].addtrials(np.array(num_data).T)
                elif control == 'stress' or control == 'maostress':
                    num_data = [line.split('\t') for line in trial]
                    num_data = [[t_units*float(line[t_col - 1]), strain_units*float(line[strain_col - 1])/100, stress_units*float(line[stress_col - 1])] for line in num_data]

                    # Add the trial to an object
                    experiments[i].addtrials(np.array(num_data).T)
                elif control == 'micro':
                    expeirments[i].addtrials(trial)

                n += 1

            i += 1

    return experiments
