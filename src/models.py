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

def crm(t, y, params, sfun):
    """ Giesekus Model ODEs """
    tau12, tau11, tau22, tau33 = y
    eta0 = params[0]
    lam = params[1]

    # tau12, tau11, tau22, tau33
    s = sfun(t)
    dt12 = eta0*s/lam + 0.5*s*(tau22 - tau11) - tau12/lam
    dt11 = s*tau12 - tau11/lam
    dt22 = -s*tau12 - tau22/lam
    dt33 = 0
    dydt = [dt12, dt11, dt22, dt33]

    return dydt

def giesekus(t, y, params, sfun):
    """ Giesekus Model ODEs """
    tau12, tau11, tau22, tau33 = y
    eta0 = params[0]
    lam = params[1]
    alpha = params[2]

    # tau12, tau11, tau22, tau33
    s = sfun(t)
    dt12 = (eta0*s/lam + s*tau22 - tau12/lam
            - alpha*tau12*(tau11 + tau22)/eta0)
    dt11 = (2*s*tau12 - tau11/lam - alpha*(tau11**2
        + tau12**2)/eta0)
    dt22 = -tau22/lam - alpha*(tau12**2 + tau22**2)/eta0
    dt33 = 0
    dydt = [dt12, dt11, dt22, dt33]

    return dydt

def johnson_segalman(t, y, params, sfun):
    """ Johnson-Segalman ODEs """
    tau12, tau11, tau22, tau33 = y
    eta0 = params[0]
    lam = params[1]
    xi = params[2]

    # tau12, ta11, tau22, tau33
    s = sfun(t)
    dt12 = eta0*s/lam + 0.5*s*(tau22 - tau11) + 0.5*(1 - xi)*s*(tau11 + tau22)  - tau12/lam
    dt11 = -tau11/lam + s*tau12 + (1 - xi)*s*tau12
    dt22 = -tau22/lam - s*tau12 + (1 - xi)*s*tau12
    dt33 = -tau33/lam
    dydt = [dt12, dt11, dt22, dt33]

    return dydt

def phan_thien_tanner(t, y, params, sfun):
    """ Phan-Thien and Tanner ODEs """
    tau12, tau11, tau22, tau33 = y
    eta0 = params[0]
    lam = params[1]
    xi = params[2]
    epsilon = params[3]

    # tau12, ta11, tau22, tau33
    s = sfun(t)
    dt12 = eta0*s/lam + 0.5*s*(tau22 - tau11) + 0.5*(1 - xi)*s*(tau11 + tau22) - (epsilon/eta0)*tau12*(tau11 + tau22 + tau33) - tau12/lam
    dt11 = -tau11/lam + s*tau12 + (1 - xi)*s*tau12 - (epsilon/eta0)*tau11*(tau11 + tau22 + tau33)
    dt22 = -tau22/lam - s*tau12 + (1 - xi)*s*tau12 - (epsilon/eta0)*tau22*(tau11 + tau22 + tau33)
    dt33 = -tau33/lam - (epsilon/eta0)*tau33*(tau11 + tau22 + tau33)
    dydt = [dt12, dt11, dt22, dt33]

    return dydt
