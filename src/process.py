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
from mapsmodels import *
import numpy as np
from itertools import permutations
import matplotlib.pyplot as plt
from scipy.optimize import fmin, curve_fit, least_squares
from scipy.interpolate import interp1d as interp
import sys

class MAPSExperiment:
    """
    Template for an experiment object, consisting of one fundamental frequency at one tone set
    run at multiple characteristic amplitudes.
    """

    def __init__(self, fundamental, harmonics, *args):
        self.base = fundamental
        self.harmonics = harmonics
        self.sums, self.mapscoords = get_maps_coords(harmonics)
        self.times = []
        self.Xs = []
        self.Ys = []
        self.amplitudes = []
        self.freqs = []
        self.XFTs = []
        self.YFTs = []

        for arg in args:
            # Store time vectors
            self.times += [arg[0]]

            # Store a strain signal
            self.Ys += [arg[1]]

            # Store a stress signal
            self.Xs += [arg[2]]

    def addtrials(self, *args):
        """
        Add trials (run at a new amplitude) to the data set.
        """
        for arg in args:
            self.times += [arg[0]]
            self.Ys += [arg[1]]
            self.Xs += [arg[2]]

    def trim(self, **kwargs):
        """
        Trim the time series to remove transient features.
        """
        # Determine the signal period
        T = 2*np.pi/self.base

        for key, value in kwargs.items():
            # If a number of burn-in cycles are specified
            if key == "nburn":
                for i in range(0,len(self.times)):
                    N = int(round(self.times[i][-1]/T))
                    n = len(self.times[i])
                    self.times[i] = self.times[i][int(value*n/N):]
                    self.Xs[i] = self.Xs[i][int(value*n/N):]
                    self.Ys[i] = self.Ys[i][int(value*n/N):]

            # If a total number of cycles is specified (take from the end of the time series)
            elif key == "ncycles":
                t = 2*np.pi*value/self.base
                for i in range(0,len(self.times)):
                    N = int(round(self.times[i][-1]/T))
                    n = len(self.times[i])
                    self.times[i] = self.times[i][int((N-value)*n/N):]
                    self.Xs[i] = self.Xs[i][int((N-value)*n/N):]
                    self.Ys[i] = self.Ys[i][int((N-value)*n/N):]

    def mean_subtract(self):
        """
        Mean subtract the X and Y time series.
        """
        for i in range(0,len(self.times)):
            self.Xs[i] = self.Xs[i] - sum(self.Xs[i])/len(self.Xs[i])
            self.Ys[i] = self.Ys[i] - sum(self.Ys[i])/len(self.Ys[i])

    def fft(self, nburn=1):
        """
        Fourier transform the stress and strain signals. Default burn in time is 2 cycles.
        """
        for i in range(0,len(self.times)):
            n = len(self.times[i])
            samplingtime = (self.times[i][-1] - self.times[i][0])/n
            self.freqs += [2*np.pi*np.fft.rfftfreq(n, d=samplingtime)]
            self.XFTs += [np.fft.rfft(self.Xs[i])/len(self.freqs[i])/2]
            self.YFTs += [np.fft.rfft(self.Ys[i])/len(self.freqs[i])/2]

    def plotft(self):
        """
        Plot the Fourier transform of X and Y.
        """
        # Plot the FT of X and Y
        fig1, ax1 = plt.subplots(1,1)
        fig2, ax2 = plt.subplots(1,1)
        for i in range(0,len(self.freqs)):
            n = len(self.freqs[i])
            ax1.semilogy(self.freqs[i], np.abs(self.XFTs[i][0:n]), marker='o', linestyle='--')
            ax2.semilogy(self.freqs[i], np.abs(self.YFTs[i][0:n]), marker='o', linestyle='--')

        ax1.set_xlim([0, 4*np.max(self.harmonics)*self.base])
        ax2.set_xlim([0, 4*np.max(self.harmonics)*self.base])
        fig1.suptitle("Input Fourier Transforms")
        fig1.suptitle("Output Fourier Transforms")

    def get_Xval(self, f):
        """
        Find the value of the FT X signal at frequency f for each trial.
        """
        idx = (np.abs(self.freqs[0] - f)).argmin()
        return [X[idx] for X in self.XFTs]

    def get_Yval(self, f):
        """
        Find the value of the FT Y signal at frequency f for each trial.
        """
        idx = (np.abs(self.freqs[0] - f)).argmin()
        return [Y[idx] for Y in self.YFTs]

    def get_Xamps(self):
        """
        Find the amplitude of all input tones.
        """
        self.amplitudes = [self.get_Xval(n*self.base) for n in self.harmonics]

    def get_Yamps(self):
        """
        Find the amplitude of all output tones.
        """
        self.y_amplitudes = [self.get_Yval(3*n*self.base) for n in self.harmonics]

    def regress1(self, n):
        """
        Determine the value of G1 at frequency n*w0.
        """
        # Get the stress and strain values at n*w0
        X = np.array(self.amplitudes[self.harmonics.index(n)])
        Y = np.array(self.get_Yval(n*self.base))
        b = np.append(np.real(Y), np.imag(Y))

        # Construct the system of equations
        Vr = np.array([np.real(X), -np.imag(X),
            np.real(X)**3 + np.real(X)*np.imag(X)**2, -np.real(X)**2*np.imag(X) - np.imag(X)**3]).T
        Vi = np.array([np.imag(X), np.real(X),
            np.real(X)**2*np.imag(X) + np.imag(X)**3, np.real(X)**3 + np.real(X)*np.imag(X)**2]).T
        V = np.append(Vr, Vi, axis = 0)

        # Solve
        x = np.linalg.lstsq(V, b, rcond=None)[0]
        r1 = x[0]
        i1 = x[1]

        return r1 + 1j*i1

    def regress3(self, results, variance, w):
        """
        Determine the value of G3 for all coordinates (nsets) located at the same frequenecy sum.
        """
        # Find all MAPS coordinates at the frequency sum
        ind = [i for i, wsum in enumerate(self.sums) if wsum == w]

        # Check if there are too many points at w to separate with the number of trials
        if len(ind) + 1 > len(self.amplitudes[0]):
            for i in ind:
                results[i] = None
            return

        # Otherwise, continue operation
        nsets = [self.mapscoords[i] for i in ind]

        # Initialize the Vandermonde matrix with the linear component
        if w in self.harmonics:
            # If we are at a channel with an active LR, regress against that
            ampsum = self.amplitudes[self.harmonics.index(w)]
        else:
            # Otherwise use the sum of the amplitudes
            ampsum = np.sum(np.array(self.amplitudes), axis=0)

        Vr = np.array([np.real(ampsum), -np.imag(ampsum)])
        Vi = np.array([np.imag(ampsum), np.real(ampsum)])

        # For every MAPS point at the frequency sum
        herm_sym = []
        for nset in nsets:
            # Track the signs for Hermitian symmetry
            signs = np.array([np.sign(n) for n in nset])
            sprod = np.prod(signs)
            herm_sym += [np.sign(sum(signs))]

            # Get the strain amplitudes at each n in nset
            ampset = np.array([self.amplitudes[self.harmonics.index(np.abs(n))] for n in nset])
            ampprod = np.prod(ampset, axis=0)

            # Get the stress amplitude at the frequency sum
            Y = self.get_Yval(self.base*sum(nset))
            b = np.append(np.real(Y), np.imag(Y))

            # Set up the Vandermonde matrix
            nperm = len(set(permutations(nset)))
            Vr = np.append(Vr, np.array([nperm*np.real(ampprod), -nperm*sprod*np.imag(ampprod)]),
                    axis=0)
            Vi = np.append(Vi, np.array([nperm*sprod*np.imag(ampprod), nperm*np.real(ampprod)]),
                    axis=0)

        V = np.append(Vr.T, Vi.T, axis = 0)

        # Solve
        x,res = np.linalg.lstsq(V, b, rcond=None)[0:2]
        if np.size(b) > np.size(x):
            cov = np.linalg.inv(np.dot(V.T,V))*res/(np.size(b) - np.size(x))
        else:
            cov = np.zeros((np.size(x),np.size(x)))

        k = 2
        for i in ind:
            self.mapscoords[i] = herm_sym[int((k-2)/2)]*self.mapscoords[i]
            results[i] = (x[k] + herm_sym[int((k-2)/2)]*1j*x[k+1])
            variance[i] = (cov[k,k] + 1j*cov[k+1,k+1])
            k += 2

    def get_viscosity(self):
        """
        Convert from the complex modulus to the complex viscosity.
        """
        self.eta1 = np.array([G_to_eta(self.G1[i], self.base, self.harmonics[i])
            for i in range(0,len(self.G1))])
        self.eta3 = np.array([G_to_eta(self.G3[i], self.base, self.mapscoords[i])
            for i in range(0,len(self.G3))])
        self.eta3_var = (np.abs(np.imag(self.G3_var)*(np.real(self.eta3)/np.imag(self.G3))**2) +
            1j*np.abs(np.real(self.G3_var)*(np.imag(self.eta3)/np.real(self.G3))**2))

    def get_fluidity(self):
        """
        Convert from the complex compliance to the complex fluidity
        """
        self.phi1 = np.array([J_to_phi(self.J1[i], self.base, self.harmonics[i])
            for i in range(0,len(self.J1))])
        self.phi3 = np.array([J_to_phi(self.J3[i], self.base, self.mapscoords[i])
            for i in range(0,len(self.J3))])
        self.phi3_var = (np.abs(np.imag(self.J3_var)*(np.real(self.phi3)/np.imag(self.J3))**2) +
            1j*np.abs(np.real(self.J3_var)*(np.imag(self.phi3)/np.real(self.J3))**2))

class StrainControlled(MAPSExperiment):
    """
    Class definition for the strain controlled subclass of MAPS experiments
    """
    def __init__(self, fundamental, harmonics, *args):
        super().__init__(fundamental, harmonics, *args)

    def interconvert(self, lrmodel):
        """
        Convert from the third order complex modulus to the third order complex compliance
        """
        self.J3 = np.array([G3_to_J3(self.G3[i], self.base, self.mapscoords[i], lrmodel)
            for i in range(0,len(self.G3))])
        self.J3_var = (np.abs(np.real(self.G3_var)*(np.real(self.J3)/np.real(self.G3))**2) +
                1j*np.abs(np.imag(self.G3_var)*(np.imag(self.J3)/np.imag(self.G3))**2))

        # Get the fluidities too
        self.get_fluidity()

    def get_MAPS(self):
        """
        Fully process the data set to get G1, G3, and eta1, eta3
        """
        # Check if fft has been called, and if not, call it
        if len(self.freqs) < len(self.times):
            self.fft()

        # Check if get_strainamps has been called, and if not, call it
        if len(self.amplitudes) == 0:
            self.get_Xamps()

        # Get the linear response values
        self.G1 = np.array([self.regress1(n) for n in self.harmonics])

        # Get the third order response values
        self.G3 = 1j*np.zeros(len(self.sums))
        self.G3_var = 1j*np.zeros(len(self.sums))
        for w in np.unique(self.sums):
            self.regress3(self.G3, self.G3_var, w)

        # Convert to eta1 and eta3
        self.get_viscosity()

        # Convert to J1
        self.J1 = np.array([1/g1 for g1 in self.G1])

class StressControlled(MAPSExperiment):
    """
    Class definition for the stress controlled subclass of MAPS experiments.
    """
    def __init__(self, fundamental, harmonics, *args):
        super().__init__(fundamental, harmonics, *args)

    def interconvert(self, lrmodel):
        """
        Convert from the third order complex compliance to the third order complex modulus.
        """
        self.G3 = np.array([J3_to_G3(self.J3[i], self.base, self.mapscoords[i], lrmodel)
            for i in range(0,len(self.J3))])
        self.G3_var = (np.abs(np.real(self.J3_var)*(np.real(self.G3)/np.real(self.J3))**2) +
                1j*np.abs(np.imag(self.J3_var)*(np.imag(self.G3)/np.imag(self.J3))**2))

        # Get the complex viscosities too
        self.get_viscosity()

    def get_MAPS(self):
        """
        Fully process the data set to get J1 and J3
        """
        # Check if fft has been called, and if not, call it
        if len(self.freqs) < len(self.times):
            self.fft()

        # Check if get_strainamps has been called, and if not, call it
        if len(self.amplitudes) == 0:
            self.get_Xamps()

        # Get the linear response values
        self.J1 = np.array([self.regress1(n) for n in self.harmonics])

        # Get the third order response values
        self.J3 = 1j*np.zeros(len(self.sums))
        self.J3_var = 1j*np.zeros(len(self.sums))
        for w in np.unique(self.sums):
            self.regress3(self.J3, self.J3_var, w)

        # Convert to phi1 and phi3
        self.get_fluidity()

        # Convert to G1
        self.G1 = np.array([1/j1 for j1 in self.J1])

class MAOStressControlled(StressControlled):
    """
    Class definition for the stress controlled subclass of MAOS experiments.
    """
    def __init__(self, fundamental, harmonics, *args):
        super().__init__(fundamental, harmonics, *args)

    def regress3(self, results, variance, w):
        """
        Determine the value of G3 for all coordinates (nsets) located at the same frequenecy sum.
        """
        # Find all MAPS coordinates at the frequency sum
        ind = [i for i, wsum in enumerate(self.sums) if wsum == w]

        # Check if there are too many points at w to separate with the number of trials
        if len(ind) + 1 > len(self.amplitudes[0]):
            for i in ind:
                results[i] = None
            return

        # Otherwise, continue operation
        nsets = [self.mapscoords[i] for i in ind]

        # Initialize the Vandermonde matrix with the linear component
        if w in self.harmonics:
            # If we are at a channel with an active LR, regress against that
            ampsum = self.amplitudes[self.harmonics.index(w)]
        else:
            # Otherwise include a constant term
            ampsum = (1+1j)*np.ones(np.shape(self.amplitudes[0]))

        Vr = np.array([np.real(ampsum), -np.imag(ampsum)])
        Vi = np.array([np.imag(ampsum), np.real(ampsum)])

        # Construct the third order components to the Vandermonde matrix
        amp = self.amplitudes[self.harmonics.index(1)]
        ampprod3 = np.array(amp)**3
        if w == 1:
            nperm3 = 3
            sprod3 = -1
            nperm5 = 20
            sprod5 = 1
        elif w == 3:
            nperm3 = 1
            sprod3 = 1
            nperm5 = 1
            sprod5 = 1

        Vr = np.append(Vr, np.array([nperm3*np.real(ampprod3), -nperm3*sprod3*np.imag(ampprod3)]), axis=0)
        Vi = np.append(Vi, np.array([nperm3*sprod3*np.imag(ampprod3), nperm3*np.real(ampprod3)]), axis=0)

        # Construct the fifth order components to the Vandermonde matrix
        ampprod5 = np.array(amp)**5
        Vr = np.append(Vr, np.array([nperm5*np.real(ampprod5), -nperm5*sprod5*np.imag(ampprod5)]), axis=0)
        Vi = np.append(Vi, np.array([nperm5*sprod5*np.imag(ampprod5), nperm3*np.real(ampprod3)]), axis=0)

        # Solve
        V = np.append(Vr.T, Vi.T, axis = 0)
        Y = self.get_Yval(self.base*w)
        b = np.append(np.real(Y), np.imag(Y))
        x,res = np.linalg.lstsq(V, b, rcond=None)[0:2]
        if np.size(b) > np.size(x):
            #cov = np.linalg.inv(np.dot(V.T,V))*res/(np.size(b) - np.size(x))
            cov = np.zeros((np.size(x),np.size(x)))
        else:
            cov = np.zeros((np.size(x),np.size(x)))

        k = 2
        for i in ind:
            self.mapscoords[i] = self.mapscoords[i]
            results[i] = (x[k] + 1j*x[k+1])
            variance[i] = (cov[k,k] + 1j*cov[k+1,k+1])
            k += 2

    def get_MAPS(self):
        """
        Fully process the data set to get J1 and J3
        """
        # Check if fft has been called, and if not, call it
        if len(self.freqs) < len(self.times):
            self.fft()

        # Check if get_strainamps has been called, and if not, call it
        if len(self.amplitudes) == 0:
            self.get_Xamps()

        # Get the linear response values
        self.J1 = np.array([self.regress1(n) for n in self.harmonics])

        # Get the third order response values
        self.J3 = 1j*np.zeros(len(self.sums))
        self.J3_var = 1j*np.zeros(len(self.sums))
        for w in np.unique(self.sums):
            self.regress3(self.J3, self.J3_var, w)

class MicroMAPS(MAPSExperiment):
    """
    Class definition for the strain controlled subclass of MAPS experiments.
    """
    def __init__(self, fundamental, harmonics, calibration, *args):
        super().__init__(fundamental, harmonics, *args)

        self.diameter = calibration["diameter"]
        self.samplingRate = calibration["samplingRate"]
        self.temperature = calibration["temperature"]
        self.ktrap = calibration["ktrap"]
        self.responsivity = calibration["responsivity"]

    def get_MAPS(self):
        """
        Fully process the data set to get G1, G3, and eta1, eta3
        """
        # Check if fft has been called, and if not, call it
        if len(self.freqs) < len(self.times):
            self.fft()

        # Check if get_strainamps has been called, and if not, call it
        if len(self.amplitudes) == 0:
            self.get_Xamps()
            self.get_Yamps()

        # Sort the data and include only the low amplitude runs
        ind = np.argsort(self.amplitudes)
        amplitudes = np.sort(self.amplitudes)
        YFTs = np.array([self.YFTs[i] for i in ind[0]])
        self.amplitudes = amplitudes[:,0:3]
        self.YFTs = YFTs[0:3,:]

        # Get the linear response values
        self.R1 = np.array([self.regress1(n) for n in self.harmonics])

        # Get the third order response values
        #self.R3 = 1j*np.zeros(len(self.sums))
        #self.R3_var = 1j*np.zeros(len(self.sums))
        #for w in np.unique(self.sums):
        #    self.regress3(self.R3, self.R3_var, w)

        # Get the first order fluid transfer function
        self.H1 = (self.ktrap/self.R1 - self.ktrap)/(3*np.pi*self.diameter*1j*self.base)

def get_maps_coords(hset):
    """
    Get the measureable MAPS coordanates (n1,n2,n3) from the set of input harmonics.
    """
    hset = [-1*n for n in hset] + hset
    sums = []
    coordinates = []
    for i1 in range(0,len(hset)):
        n1 = hset[i1]
        for i2 in range(i1,len(hset)):
            n2 = hset[i2]
            for i3 in range(i2,len(hset)):
                n3 = hset[i3]
                if n1+n2+n3 > 0:
                    sums += [n1+n2+n3]
                    coordinates += [[n1,n2,n3]]

    # Find the unique points
    return np.array(sums), np.array(coordinates)


def G_to_eta(G, w0, nset):
    """
    Convert from the complex modulus to the complex viscosity.
    """
    wset = 1j*w0*nset
    eta = G/np.prod(wset)
    return eta

def J_to_phi(J, w0, nset):
    """
    Convert from the complex compliance to the complex fluidity.
    """
    wset = 1j*w0*nset
    phi = J*np.sum(wset)
    return phi

def G3_to_J3(G3, w0, nset, lrmodel):
    """
    Convert from the complex modulus to the complex compliance.
    """
    G1 = lambda w: lrmodel(w)
    return -1*G3/(G1(w0*nset[0])*G1(w0*nset[1])*G1(w0*nset[2])*G1(w0*sum(nset)))

def J3_to_G3(J3, w0, nset, lrmodel):
    """
    Convert from the complex compliance to the complex modulus.
    """
    J1 = lambda w: 1/lrmodel(w)
    return -1*J3/(J1(w0*nset[0])*J1(w0*nset[1])*J1(w0*nset[2])*J1(w0*sum(nset)))

def fitLR(data,lrmodel):
    """
    Fit the linear response frequency sweep data to a Maxwell mode.
    Returns [eta0, lambda].
    """
    w = data[0]
    G = data[1]
    Gr = np.real(G)
    Gi = np.imag(G)

    # Fit to the Maxwell mode
    obj = lambda p: np.sum(np.abs(lrmodel(w, p) - (Gr + 1j*Gi))**2)
    popt,fopt,_,_,_ = fmin(func=obj, x0=[1,1,1], disp=False, full_output=1)

    # Calculate the uncertainty in the parameters
    eta0 = popt[0]
    tau = popt[1]
    etainf = popt[2]
    dfdeta0 = np.concatenate((tau*w**2/(1 + (tau*w)**2), w/(1 + (tau*w)**2)))
    dfdtau = np.concatenate((eta0*w**2/(1 + (tau*w)**2) - 2*eta0*tau**2*w**4/((1 + (tau*w)**2)**2), -2*eta0*tau*w**3/((1 + (tau*w)**2)**2)))
    dfdetainf = np.concatenate((np.zeros(np.shape(w)), w))
    A = np.matrix([dfdeta0, dfdtau, dfdetainf]).T
    n = len(w)
    V = (fopt/(n - 3))*np.linalg.inv(np.dot(A.T,A))
    pvar = np.diag(V)/n

    return popt, pvar

def maxwellLR(w, p):
    """
    Define the Maxwell mode linear response function.
    """
    eta0 = p[0]
    lam = p[1]
    etainf = p[2]
    G = (1j*eta0*w/(1 + 1j*lam*w) + 1j*w*etainf)
    return G

def interpolateLR(lrdata,k):
    """
    Define a linear interpolant for linear response data (returns an interpolant function).
    """
    wd = lrdata[0]
    G = lrdata[1]
    Gr = np.real(G)
    Gi = np.imag(G)

    # Get fit type
    if k == 1:
        order = "slinear"
    elif k == 2:
        order = "quadratic"
    elif k == 3:
        order = "cubic"
    else:
        sys.exit("Error: Spline interpolant must be of degree 1 <= k <= 3")

    # Define interpolant funtions for the log of G' and G''
    Gr_int = interp(np.log(wd),np.log(Gr),kind=order)
    Gi_int = interp(np.log(wd),np.log(Gi),kind=order)

    # Combine into composite function
    G_int = lambda w: np.exp(Gr_int(np.log(np.abs(w)))) + 1j*np.sign(w)*np.exp(Gi_int(np.log(np.abs(w))))
    return G_int

def fit11const(experiments,eta0,lambda1,lambda2):
    """
    Fit MAPS experiments to the 11 constant model.
    """
    eta3 = np.array([])
    eta3_11 = np.zeros((2,6))
    for experiment in experiments:
        # Get experimental data
        eta3 = np.append(eta3,experiment.eta3)
        coords = experiment.base*experiment.mapscoords
        # Get value of rational functions
        for coord in coords:
            S2 = eta0*f2(coord[0],coord[1],coord[2],lambda1,lambda2)
            S31 = eta0*f31(coord[0],coord[1],coord[2],lambda1,lambda2)
            S32 = eta0*f32(coord[0],coord[1],coord[2],lambda1,lambda2)
            S41 = eta0*f41(coord[0],coord[1],coord[2],lambda1,lambda2)
            S42 = eta0*f42(coord[0],coord[1],coord[2],lambda1,lambda2)
            S5 = eta0*f5(coord[0],coord[1],coord[2],lambda1,lambda2)
            eta3_vect = np.array([[S2, S31, S32, S41, S42, S5]])
            eta3_11 = np.concatenate((eta3_11,eta3_vect),axis=0)

    eta3_11 = eta3_11[2:,:]

    # Filter out nans
    for i in range(len(eta3)-1,-1,-1):
        if np.isnan(eta3[i]):
            eta3 = np.delete(eta3,i)
            eta3_11 = np.delete(eta3_11,i,0)

    # Regress for A (starting from CRM)
    obj = lambda A: obj_fit11const(A,eta3,eta3_11)
    A0 = [0,-lambda1**2,0,0,0,0]
    A = least_squares(obj,A0)
    return A['x']

def obj_fit11const(A,eta3,eta3_11_mat):
    """
    Define the objective constant for the 11-constant fit.
    """
    eta3_11 = np.zeros(np.size(eta3))
    for i in range(0,6):
        eta3_11 = eta3_11 + A[i]*eta3_11_mat[:,i]

    # Process/scale
    obj_mag = np.log(np.abs(eta3)) - np.log(np.abs(eta3_11))
    obj_arg = np.angle(eta3) - np.angle(eta3_11)
    obj_real = np.log(np.abs(np.real(eta3) - np.real(eta3_11)))
    obj_imag = np.log(np.abs(np.imag(eta3) - np.imag(eta3_11)))
    obj = np.append(obj_real,obj_imag)
    return obj

def gap_loading(lrmodel,rho,wrange):
    """
    Determine the gap loading limit at any time-scale (w) in range of interest
    """
    w = np.logspace(np.log10(wrange[0]),np.log10(wrange[1]),1000)
    G1 = lrmodel(w)
    eta1 = G1/(1j*w)
    etamag = np.abs(eta1)
    etaarg = np.arctan(-np.real(eta1)/np.imag(eta1))

    # Compute the wavelength and penetration depth
    lambdas = 2*np.pi/(np.sqrt(w*rho/etamag)*np.cos(etaarg/2))
    ds = 1/(np.sqrt(w*rho/etamag)*np.sin(etaarg/2))

    # Calculate the minimum values over the range of interest
    lambdamin = np.min(lambdas)
    dmin = np.min(ds)
    lmax_1 = np.min([lambdamin,dmin])
    lmax_3 = lmax_1/3

    # Plot the frequency dependence of these quantities (and one third, for MAPS)
    fig,ax = plt.subplots(1,1)
    ax.loglog(w,lambdas,'b-',label=r"$\lambda_s$")
    ax.loglog(w,ds,'r-',label=r"$d_s$")
    ax.loglog(w,lambdas/3,'b--',label=r"$\lambda_s/3$")
    ax.loglog(w,ds/3,'r--',label=r"$d_s/3$")
    ax.loglog(w,np.ones(np.shape(w))*lmax_1,'k',label=r"$l_{1} = $" + '{0:.2E}'.format(lmax_1))
    ax.loglog(w,np.ones(np.shape(w))*lmax_3,'k--',label=r"$l_{3} = $" + '{0:.2E}'.format(lmax_3))
    ax.legend()
    fig.suptitle("Gap Loading Limit Analysis")

