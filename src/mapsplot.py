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

from process import *
from mapsmodels import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import timeit

class MAPSFigure():
    """
    Class for MAPS figure objects.
    """

    def __init__(self, name, w0, mapscoords, data, variances):
        self.name = name
        self.mapscoords = mapscoords
        barycoords = [get_barycentric(nset) for nset in mapscoords]
        self.barysets = list(np.unique(barycoords, axis=0))
        self.barysets = [list(el) for el in self.barysets]
        self.A = [ [] for n in range(0,len(self.barysets))]
        self.B = [ [] for n in range(0,len(self.barysets))]
        self.C = [ [] for n in range(0,len(self.barysets))]
        self.D = [ [] for n in range(0,len(self.barysets))]

        # Sort data into triangles and by coordinate
        self.triangles = []
        self.baryind = []
        self.norms = []
        for nset in mapscoords:
            nset.sort()
            if nset[1] < 0:
                nset = [-nset[2], -nset[1], -nset[0]]
            if nset[0] > 0:
                self.triangles += "A"
            elif -1*nset[0] >= nset[1]:
                if -1*nset[0] > nset[2]:
                    self.triangles += "D"
                else:
                    self.triangles += "C"
            else:
                self.triangles += "B"

            self.norms += [sum(np.abs(nset))]

        self.barylabels = []
        for bary in barycoords:
            self.barylabels += [self.barysets.index(bary)]

        self.add_data(w0, data, variances)

    def add_data(self, w0, data, variances):
        """
        Add data to each triangle corresponding to certain barycentric coordinates.
        """
        # Sort into triangles and by barycoord
        for i in range(0,len(data)):
            # Check if the data point is empty
            if self.triangles[i] == "A":
                self.A[self.barylabels[i]] += [[data[i], w0*self.norms[i], variances[i]]]
            elif self.triangles[i] == "B":
                self.B[self.barylabels[i]] += [[data[i], w0*self.norms[i], variances[i]]]
            elif self.triangles[i] == "C":
                self.C[self.barylabels[i]] += [[data[i], w0*self.norms[i], variances[i]]]
            else:
                self.D[self.barylabels[i]] += [[data[i], w0*self.norms[i], variances[i]]]

    def plot_data(self, verbose=1):
        """
        Plot the data for the current figure.
        """
        # Determine the number of sampled subspaces
        m = 1
        for set in self.B:
            if len(set) > 0:
                m = 2
                break

        # Create a 2 by m subplot array
        fig1, ax1 = plt.subplots(2,m)
        fig2, ax2 = plt.subplots(2,m)
        fig3, ax3 = plt.subplots(2,m)

        # Build the plots and save
        self.make_plots((fig1,ax1), (fig2,ax2), (fig3,ax3),
                self.A, self.B, self.C, self.D, m)
        if verbose == 1:
            plt.show()
        elif verbose == 2:
            plt.show()
            fig1.savefig(self.name + '_bode_mag.eps')
            fig2.savefig(self.name + '_bode_arg.eps')
            fig3.savefig(self.name + '_nyquist.eps')

        self.bodemag = (fig1,ax1)
        self.bodearg = (fig2,ax2)
        self.nyquist = (fig3,ax3)

    def make_plots(self, bodemag, bodearg, nyquist, A, B, C, D, m, mark=True, line='dashed'):
        """
        Create plots given data.
        """
        # Unpack the figures
        fig1, ax1 = bodemag
        fig1.suptitle("Bode Plot - Magnitude")
        fig2, ax2 = bodearg
        fig2.suptitle("Bode Plot - Argument")
        fig3, ax3 = nyquist
        fig3.suptitle("Nyquist Plot")

        # For each triangle, plot all isobarycentric curves
        markers = ['o','s','^','*','P','X','d','v','<','>','p']
        fills = ["full","none"]
        wmin = 10000.0
        wmax = 0.0
        magmin = 1E10
        magmax = 0.0
        for i in range(0, len(self.barysets)):
            # Determine the marker style and fill
            if mark:
                mstyle = markers[int(np.floor(i/2))]
                fstyle = fills[int(np.mod(i,2))]
            else:
                mstyle = "None"
                fstyle = "none"

            # Make the Bode plots and Nyquist diagram
            for j in range(0,2*m):
                if j == 0:
                    data = np.array(A[i])
                    k,l = [0,0]
                elif j == 1:
                    data = np.array(C[i])
                    k,l = [1,0]
                elif j == 2:
                    data = np.array(B[i])
                    k,l = [0,1]
                else:
                    data = np.array(D[i])
                    k,l = [1,1]

                # Assign the appropriate axes
                if m < 2:
                    ny_ax = ax3[k]
                    bm_ax = ax1[k]
                    ba_ax = ax2[k]
                else:
                    ny_ax = ax3[k,l]
                    bm_ax = ax1[k,l]
                    ba_ax = ax2[k,l]

                if data.any():
                    data = data[data[:,1].argsort()]

                    # Plot the Nyquist diagram
                    ny_ax.plot(np.real(data[:,0]), np.imag(data[:,0]),
                            color=self.barysets[i], marker=mstyle, fillstyle=fstyle, linestyle=line)

                    # Plot error bars on the Bode plots (for experiments with mark = 'o')
                    if mark:
                        # Get the maximum and minimum frequencies
                        if np.real(data[0,1]) < wmin:
                            wmin = np.real(data[0,1])
                        if np.real(data[-1,1]) > wmax:
                            wmax = np.real(data[-1,1])

                        # Calculate the uncertainty of the abs and arg
                        a = np.real(data[:,0])
                        b = np.imag(data[:,0])
                        abs_var = ((a**2)*np.real(data[:,2]) +
                                (b**2)*np.imag(data[:,2]))/(np.abs(data[:,0])**2)
                        arg_var = ((((1/(1 + (b/a)**2))*(1/a))**2)*np.imag(data[:,2]) +
                                (((1/(1 + (b/a)**2))*(b/(a**2)))**2)*np.real(data[:,2]))

                        # Plot on appropriate subplot
                        mag = np.abs(data[:,0])
                        arg = np.angle(data[:,0])
                        bm_ax.set_xscale('log',nonposx='clip')
                        bm_ax.set_yscale('log',nonposy='clip')
                        bm_ax.errorbar(np.real(data[:,1]), mag,
                                yerr=np.sqrt(abs_var), color=self.barysets[i],
                                marker=mstyle, fillstyle=fstyle, linestyle=line)
                        ba_ax.set_xscale('log',nonposx='clip')
                        ba_ax.errorbar(np.real(data[:,1]), arg*2/np.pi,
                                yerr=np.sqrt(arg_var)*2/np.pi, color=self.barysets[i],
                                marker=mstyle, fillstyle=fstyle, linestyle=line)

                        # Get the maximum and minimum magnitudes
                        if np.min(mag) < magmin:
                            magmin = np.min(mag)
                        if np.max(mag) > magmax:
                            magmax = np.max(mag)

                    else:
                        bm_ax.loglog(np.real(data[:,1]), np.abs(data[:,0]),
                                color=self.barysets[i], marker=mstyle, fillstyle=fstyle, linestyle=line)
                        ba_ax.semilogx(np.real(data[:,1]), np.angle(data[:,0])*2/np.pi,
                                color=self.barysets[i], marker=mstyle, fillstyle=fstyle, linestyle=line)



        # Set the appearance of the axes (top + bottom ticks, ticks facing in, enforce
        # ticks every decade)
        for k in [0,1]:
            for l in [0,m-1]:
                if m < 2:
                    bm_ax = ax1[k]
                    ba_ax = ax2[k]
                    ny_ax = ax3[k]
                else:
                    bm_ax = ax1[k,l]
                    ba_ax = ax2[k,l]
                    ny_ax = ax3[k,l]

                # Adjust some plotting parameters
                bm_ax.tick_params(which='both', direction='in', top=True, right=True)
                ba_ax.tick_params(which='both', direction='in', top=True, right=True)
                ny_ax.tick_params(which='both', direction='in', top=True, right=True)
                ba_ax.set_ylim([-2.2,2.2])
                ba_ax.yaxis.set_major_formatter(tck.FuncFormatter(arg_formatter))
                ba_ax.yaxis.set_major_locator(tck.MultipleLocator(base=1.0))
                ba_ax.yaxis.set_minor_locator(tck.MultipleLocator(0.2))

                if mark:
                    ba_ax.set_xlim([wmin/2,wmax*2])
                    bm_ax.set_xlim([wmin/2,wmax*2])
                    bm_ax.set_ylim([magmin/2,magmax*2])

                # More plotting parameters
                #ba_ax.xaxis.set_major_locator(tck.LogLocator(base=10,numticks=6))
                #ba_ax.xaxis.set_minor_locator(tck.LogLocator(base=10,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=8))

                #bm_ax.xaxis.set_major_locator(tck.LogLocator(base=10,numticks=6))
                #bm_ax.xaxis.set_minor_locator(tck.LogLocator(base=10,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=8))

                #bm_ax.yaxis.set_major_locator(tck.LogLocator(base=10,numticks=17))
                #bm_ax.yaxis.set_minor_locator(tck.LogLocator(base=10,subs=(0.33,0.67,1.0),numticks=17))
                #bm_ax.set_yticklabels(["","",r"$10^{-15}$","","",r"$10^{-12}$","","",r"$10^{-9}$","","",r"$10^{-6}$","","",r"$10^{-3}$","","",r"$10^0$"])

    def plot_model(self, models, params, wlim, overlap=0, verbose=0):
        """
        Plot the model prediction over the given given data range.
        """
        # Determine the number of sampled subspaces
        m = 1
        for set in self.B:
            if not self.A[0] or len(set) > 0:
                m = 2
                break

        # Create the subplots
        if overlap == 1:
            fig1, ax1 = self.bodemag
            fig2, ax2 = self.bodearg
            fig3, ax3 = self.nyquist
        else:
            fig1, ax1 = plt.subplots(2,2)
            fig2, ax2 = plt.subplots(2,2)
            fig3, ax3 = plt.subplots(2,2)

        # Plot up to two models on the same axes
        j = 0
        linespecs = ['-','--']
        for model in models:
            # Generate the synthetic data
            A = [[] for n in range(0, len(self.barysets))]
            B = [[] for n in range(0, len(self.barysets))]
            C = [[] for n in range(0, len(self.barysets))]
            D = [[] for n in range(0, len(self.barysets))]
            w = np.logspace(np.log10(wlim[0]/(30*np.max(self.mapscoords))),np.log10(wlim[1]*30*np.max(self.mapscoords)),1000)

            # Make curve for every coordinate
            for i in range(0,len(self.mapscoords)):
                nset = self.mapscoords[i]
                triangle = self.triangles[i]
                barylabel = self.barylabels[i]

                # If in overlap mode, check if data exists at the given barylabel
                start_time = timeit.default_timer()
                if overlap == 1 or overlap == 0:
                    if triangle == 'A':
                        data = self.A[barylabel]
                    elif triangle == 'B':
                        data = self.B[barylabel]
                    elif triangle == 'C':
                        data = self.C[barylabel]
                    elif triangle == 'D':
                        data = self.D[barylabel]

                    if data and np.isnan(data[0][0]):
                        continue

                w1 = nset[0]*w
                w2 = nset[1]*w
                w3 = nset[2]*w
                if triangle == 'A':
                    eta3 = list(model(w1,w2,w3,params))
                    wv = list(np.abs(w1) + np.abs(w2) + np.abs(w3))

                    # Organize in the same way as experimental data
                    A[barylabel] += [[eta3[n],wv[n]] for n in range(0,len(eta3))]
                elif triangle == 'B':
                    eta3 = list(model(w1,w2,w3,params))
                    wv = list(np.abs(w1) + np.abs(w2) + np.abs(w3))
                    B[barylabel] += [[eta3[n],wv[n]] for n in range(0,len(eta3))]
                elif triangle == 'C':
                    eta3 = list(model(w1,w2,w3,params))
                    wv = list(np.abs(w1) + np.abs(w2) + np.abs(w3))
                    C[barylabel] += [[eta3[n],wv[n]] for n in range(0,len(eta3))]
                elif triangle == 'D':
                    eta3 = list(model(w1,w2,w3,params))
                    wv = list(np.abs(w1) + np.abs(w2) + np.abs(w3))
                    D[barylabel] += [[eta3[n],wv[n]] for n in range(0,len(eta3))]

            # Plot
            self.make_plots((fig1,ax1), (fig2,ax2), (fig3,ax3), A, B, C, D, m, False, linespecs[j])
            j += 1

        if verbose == 1:
            plt.show()
        elif verbose == 2:
            plt.show()
            fig1.savefig(self.name + '_bode_mag_model.eps')
            fig2.savefig(self.name + '_bode_arg_model.eps')
            fig3.savefig(self.name + '_nyquist_model.eps')

        if overlap == 1:
            self.bodemag = (fig1,ax1)
            self.bodearg = (fig2,ax2)
            self.nyquist = (fig3,ax3)
        else:
            self.bodemag_model = (fig1,ax1)
            self.bodearg_model = (fig2,ax2)
            self.nyquist_model = (fig3,ax3)

def get_barycentric(nset):
    # Calculate the barycentric coordinates from the associated frequency set
    rv = np.array([0,0,1.0])
    gv = np.array([0,1./2.,1./2])
    bv = np.array([1./3.,1./3.,1./3.])
    A = np.linalg.norm(np.cross(gv-rv,bv-rv))

    # Calculate the barycentric coords for each point via cross products
    coord = np.sort(np.abs(nset))
    coord = coord.astype(float)/np.sum(coord)
    R = np.linalg.norm(np.cross(coord-gv, bv-gv))/A
    G = np.linalg.norm(np.cross(coord-rv, bv-rv))/A
    B = np.linalg.norm(np.cross(coord-rv, gv-rv))/A

    return [R,G,B]

def make_figures(experiments, mapsfun='eta', verbose=1):
    """
    Make figure objects from the experiment objects.
    """
    figures = {}
    names = []

    for experiment in experiments:
        # Make the figure name
        name = '_'.join(map(str, experiment.harmonics))
        names += [name]

        # Either create or append a dictionary entry for the figure object
        if mapsfun == 'eta':
            try:
                figures[name].add_data(experiment.base, experiment.eta3, experiment.eta3_var)
            except KeyError:
                figures[name] = MAPSFigure(name, experiment.base, experiment.mapscoords,
                        experiment.eta3, experiment.eta3_var)
        elif mapsfun == 'G':
            try:
                figures[name].add_data(experiment.base, experiment.G3, experiment.G3_var)
            except KeyError:
                figures[name] = MAPSFigure(name, experiment.base, experiment.mapscoords,
                        experiment.G3, experiment.G3_var)
        elif mapsfun == 'J':
            try:
                figures[name].add_data(experiment.base, experiment.J3, experiment.J3_var)
            except KeyError:
                figures[name] = MAPSFigure(name, experiment.base, experiment.mapscoords,
                        experiment.J3, experiment.J3_var)
        else:
            print('Error: No MAPS function called ' + mapsfun)

    for name in set(names):
        figures[name].plot_data(verbose)

    return figures

def arg_formatter(val,pos):
    """
    Formatter for Bode plots of the argument of a MAPS function (y-axis at multiple of \pi)
    """
    if val == 0:
        return '0'
    elif val == -1:
        return "-$\pi/2$"
    elif val == -2:
        return "-$\pi$"
    elif val == 1:
        return "$\pi/2$"
    elif val == 2:
        return "$\pi$"

def plot_linear_response(lrdata,lrmodel,marker='o',**kwargs):
    """
    Plot the linear response data and fit (G and J)
    """
    w = lrdata[0]
    G = lrdata[1]
    Gr = np.real(G)
    Gi = np.imag(G)
    J = 1/G
    Jr = np.real(J)
    Ji = -np.imag(J)

    # Plot the data
    try:
        figG,axG,figJ,axJ = kwargs["axes"]
    except KeyError:
        figG,axG = plt.subplots(1,1)
        figJ,axJ = plt.subplots(1,1)

    figG.suptitle("Linear Complex Modulus")
    axG.loglog(w,Gr,'r'+marker)
    axG.loglog(w,Gi,'b'+marker)

    figJ.suptitle("Linear Complex Compliance")
    axJ.loglog(w,Jr,'r'+marker)
    axJ.loglog(w,Ji,'b'+marker)

    # Plot the fitted model
    w_v = np.logspace(np.log10(np.min(w)),np.log10(np.max(w)),100)
    Gpred = lrmodel(w_v)
    Jpred = 1/Gpred

    axG.loglog(w_v,np.real(Gpred),'r')
    axG.loglog(w_v,np.imag(Gpred),'b')

    axJ.loglog(w_v,np.real(Jpred),'r')
    axJ.loglog(w_v,-np.imag(Jpred),'b')

    return figG,axG,figJ,axJ

def tss_comparison(experiments,lrmodel):
    """
    Create parity plots to assess if data is TSS-like
    """
    figtssr,axtssr = plt.subplots(1,1)
    figtssr.suptitle(r"TSS Comparison - $|\eta^*_3|$")
    figtssi,axtssi = plt.subplots(1,1)
    figtssi.suptitle(r"TSS Comparison - arg$\eta^*_3$")
    x = []
    y = []

    # Loop over all experiments
    for experiment in experiments:
        w0 = experiment.base
        eta3 = experiment.eta3
        k = 0

        # Loop over all MAPS coordinates
        for nset in experiment.mapscoords:
            # Get the barycentric coordinates and MAPS measurement
            bary = get_barycentric(np.abs(nset))
            eta3_val = eta3[k]

            # Check if value is NaN
            if not np.isnan(eta3_val):
                y += [np.log(np.abs(eta3_val))]

                # Get the TSS prediction
                w1 = nset[0]*w0
                w2 = nset[1]*w0
                w3 = nset[2]*w0
                eta3tss = tss_eta3(w1,w2,w3,lrmodel)
                x += [np.log(np.abs(eta3tss))]

                # Check if the sign matches
                if np.sign(eta3tss) == np.sign(eta3_val):
                    fill = 'full'
                else:
                    fill = 'none'

                # Create parity plot
                tssarg = np.angle(eta3tss)
                valarg = np.angle(-eta3_val)
                if valarg - tssarg > np.pi:
                    valarg = valarg - 2*np.pi
                elif valarg - tssarg < -np.pi:
                    valarg = valarg + 2*np.pi

                axtssr.plot(np.log(np.abs((eta3tss))),np.log(np.abs((eta3_val))),color=bary,marker='o',fillstyle=fill)
                axtssi.plot(tssarg,valarg,color=bary,marker='o',fillstyle=fill)
            k += 1

    # Regress
    A = np.array([np.array(x),np.ones(np.shape(x))]).T
    y = np.array([y]).T
    b = np.linalg.lstsq(A,y,rcond=None)[0]

    # Regress with slope of 1
    A1 = y - np.array([x]).T
    b1 = np.average(A1)

    # Plot the regressed lines
    xv = np.linspace(np.min(x),np.max(x))
    axtssr.plot(xv,xv + b1*np.ones(np.shape(xv)),'k',label="y = x + {:.2f}".format(b1))
    axtssi.plot(np.linspace(-np.pi,np.pi),np.linspace(-np.pi,np.pi),'k')
    axtssr.legend()
