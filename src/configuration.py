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

import numpy as np

def set_configuration(control):
    """
    Sets confiration parameters for file I/O
    """
    ### EDIT THE FOLLOWING LINES TO CHANGE THE CONFIGURATION FOR MAPS FILES ###
    if control == "strain" or control == None:
        ### EDIT THE FOLLOWING LINES FOR STRAIN-CONTROLLED EXPERIMENTS ###
        t_col = 1 # The column (index-1) containing the time value
        stress_col = 2 # The column (index-1) containing the stress time-series
        strain_col = 3 # The column (index-1) containing the strain time-series
    elif control == "stress" or control == "maostress":
        ### EDIT THE FOLLOWING LINES FOR STRESS-CONTROLLED EXPERIMENTS ###
        t_col = 4 # The column (index-1) containing the time value
        stress_col = 6 # The column (index-1) containing the stress time-series
        strain_col = 5 # The column (index-1) containing the strain time-series

    t_units = 1 # Unit correction for time (assuming units of s)
    #t_units = 0.001 <- example if units of time are in ms (0.001 s = 1 ms)
    stress_units = 1 # Unit correction for stress (assuming units of Pa)
    #stress_units = 1000 <- example if units of stress are in kPa (1000 Pa = 1 kPa)
    strain_units = 1 # Unit correction for strain (assuming units of %)
    #strain_units = 100 <- example if units are strain units (100% = 1 strain unit)

    ### EDIT THE FOLLOWING LINES TO CHANGE THE CONFIGURATION FOR SAOS FILES ###
    w_col = 4 # The column (index-1) containing the angular frequency
    Gprime_col = 1 # The column (index-1) containing the storage modulus
    Gdoubprime_col = 2 # The column (index-1) containing the loss modulus
    t_col_LR = 6 # The column (index-1) containing the step time
    w_units = 1 # Unit correction for angular frequency (assuming units of rad/s)
    #w_units = 2*np.pi <- example if units of angular frequency are 1/s (2\pi rad/s = 1/s)
    G_units = 1 # Unit correction for the moduli (assuming units of Pa)
    # G_units = 1000 <- example if units of the moduli are kPa (1000 Pa = 1 kPa)

    configs = {"t_col":t_col, "stress_col":stress_col, "strain_col":strain_col, "t_units":t_units, "stress_units":stress_units, "strain_units":strain_units, "w_col":w_col, "Gprime_col":Gprime_col, "Gdoubprime_col":Gdoubprime_col, "t_col_LR":t_col_LR, "w_units":w_units, "G_units":G_units}

    return configs

