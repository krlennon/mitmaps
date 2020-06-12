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

from scipy.integrate import odeint
import numpy as np
from models import *
from process import *

def sim_multitone_eq(model, params, amp, harmonics, fundamentals, epsilon):
    """
    Simulate a multitone strain-controlled MAPS experiment, in which all tones have equal amplitude.
    Inputs:
        model - the model as a system of ODEs
        params - list of model parameters
        amp - the amplitude of the strain tones
        harmonics - the harmonics of the input tones
        fundamentals - list of fundamental frequencies to sweep
        epsilon - the noise floor (relative)
    Outputs:
        experiments - array of experiment objects
    """
    experiments = []
    for w0 in fundamentals:
        # Run a simulation for three amplitudes
        exp = StrainControlled(w0, harmonics)
        for f in [1, 0.5, 0.5]:
            # Generate the input signal
            s = lambda t: sum([f*amp*n*w0*np.cos(n*w0*t) for n in harmonics])
            gamma = lambda t: sum([f*amp*np.sin(n*w0*t) for n in harmonics])

            # Define the model
            dydt = lambda y, t: model(t, y, params, s)

            # Run a simulation
            t0 = 0
            tf = 20*2*np.pi/w0
            time = np.linspace(t0, tf, 2048*10)
            strain = [gamma(t) for t in time]
            sol = odeint(dydt, [0,0,0,0], time)
            stress = sol[:,0]

            # Add noise
            strain += epsilon[0]*np.random.randn(len(strain))
            stress += epsilon[1]*np.random.randn(len(stress))

            # Add the trial to the experiment object
            exp.addtrials([np.array(time), np.array(stress), np.array(strain)])

        experiments += [exp]

    return experiments

def sim_multitone_free(model, params, amp, harmonics, fundamentals):
    """
    Simulate a multitone strain-controlled MAPS experiment, in which tones have free amplitudes.
    Inputs:
        model - the model as a system of ODEs
        params - list of model parameters
        amp - the amplitude of the strain tones
        harmonics - the harmonics of the input tones
        fundamentals - list of fundamental frequencies to sweep
    Outputs:
        experiments - array of experiment objects
    """
    experiments = []
    for w0 in fundamentals:
        # Run a simulation for three amplitudes
        exp = StrainControlled(w0, harmonics)
        for f in [[1, 1, 1],[0.5, 1, 1],[1, 0.5, 1],[1, 1, 0.5]]:
            # Generate the input signal
            s = lambda t: sum([f[i]*amp*n*w0*np.cos(n*w0*t) for i,n in enumerate(harmonics)])
            gamma = lambda t: sum([f[i]*amp*np.sin(n*w0*t) for i,n in enumerate(harmonics)])

            # Define the model
            dydt = lambda y, t: model(t, y, params, s)

            # Run a simulation
            t0 = 0
            tf = 10*2*np.pi/w0
            time = np.linspace(t0, tf, 2048*10)
            strain = [gamma(t) for t in time]
            sol = odeint(dydt, [0,0,0,0], time)
            stress = sol[:,0]

            # Add the trial to the experiment object
            exp.addtrials([np.array(time), np.array(stress), np.array(strain)])

        experiments += [exp]

    return experiments
