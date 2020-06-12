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

def f2(w1,w2,w3,lambda1,lambda2):
    # Defines the Omega2 rational function
    f = lambda w: 1/(1 + 1j*lambda1*w)
    return f(w1+w2+w3)*(f(w1+w2) + f(w1+w3) + f(w2+w3))/3

def f31(w1,w2,w3,lambda1,lambda2):
    # Defines the Omega31 rational function
    f = lambda w: 1/(1 + 1j*lambda1*w)
    g = lambda w: (1 + 1j*lambda2*w)/(1 + 1j*lambda1*w)
    return f(w1+w2+w3)*(f(w1+w2)*(g(w1)+g(w2)) + f(w1+w3)*(g(w1)+g(w3)) +
            f(w2+w3)*(g(w2)+g(w3)))/6

def f32(w1,w2,w3,lambda1,lambda2):
    # Defines the Omega32 rational function
    f = lambda w: 1/(1 + 1j*lambda1*w)
    g = lambda w: (1 + 1j*lambda2*w)/(1 + 1j*lambda1*w)
    return f(w1+w2+w3)*(f(w1+w2)*g(w3) + f(w1+w3)*g(w2) + f(w2+w3)*g(w1))/3

def f41(w1,w2,w3,lambda1,lambda2):
    # Defines the Omega41 rational function
    f = lambda w: 1/(1 + 1j*lambda1*w)
    g = lambda w: (1 + 1j*lambda2*w)/(1 + 1j*lambda1*w)
    return f(w1+w2+w3)*(g(w3)*f(w1+w2)*(g(w1)+g(w2)) + g(w2)*f(w1+w3)*(g(w1)+g(w3)) +
            g(w1)*f(w2+w3)*(g(w2)+g(w3)))/6

def f42(w1,w2,w3,lambda1,lambda2):
    # Defines the Omega42 rational function
    f = lambda w: 1/(1 + 1j*lambda1*w)
    g = lambda w: (1 + 1j*lambda2*w)/(1 + 1j*lambda1*w)
    return f(w1+w2+w3)*(f(w1+w2)*g(w1)*g(w2) + f(w1+w3)*g(w1)*g(w3) + f(w2+w3)*g(w2)*g(w3))/3

def f5(w1,w2,w3,lambda1,lambda2):
    # Defines the Omega5 rational function
    f = lambda w: 1/(1 + 1j*lambda1*w)
    g = lambda w: (1 + 1j*lambda2*w)/(1 + 1j*lambda1*w)
    return f(w1+w2+w3)*(f(w1+w2) + f(w1+w3) + f(w2+w3))*g(w1)*g(w2)*g(w3)/3

def eleven_const(w1,w2,w3,params):
    # Define parameters
    eta0 = params[0]
    lambda1 = params[1]
    etainf = params[2]
    lambda2 = params[3]
    A = params[4]

    # Defines the '11 constant' model
    g2 = f2(w1,w2,w3,lambda1,lambda2)
    g31 = f31(w1,w2,w3,lambda1,lambda2)
    g32 = f32(w1,w2,w3,lambda1,lambda2)
    g41 = f41(w1,w2,w3,lambda1,lambda2)
    g42 = f42(w1,w2,w3,lambda1,lambda2)
    g5 = f5(w1,w2,w3,lambda1,lambda2)
    return eta0*(A[0]*g2 + A[1]*g31 + A[2]*g32 + A[3]*g41 + A[4]*g42 + A[5]*g5)

def giesekus_eta3(w1,w2,w3,params):
    # Define parameters
    eta0 = params[0]
    lambda1 = params[1]
    etainf = params[2]
    alpha = params[3]

    # Calculate eta3 from the rational functions
    Omega41 = f41(w1,w2,w3,lambda1,0)
    Omega42 = f42(w1,w2,w3,lambda1,0)
    Omega5 = f5(w1,w2,w3,lambda1,0)
    eta3 = eta0*(lambda1**2)*(-2*alpha*Omega41 - alpha*Omega42 + 2*(alpha**2)*Omega5)
    return eta3

def crm_eta3(w1,w2,w3,params):
    # Define parameters
    eta0 = params[0]
    lambda1 = params[1]
    etainf = params[2]

    # Calculate eta3 from the rational functions
    Omega31 = f31(w1,w2,w3,lambda1,0)
    eta3 = -eta0*(lambda1**2)*Omega31
    return eta3

def crm_J3(w1,w2,w3,params):
    # Define parameters
    eta0 = params[0]
    lambda1 = params[1]
    etainf = params[2]

    # Calculate J3 from eta3
    Omega31 = f31(w1,w2,w3,lambda1,0)
    eta3 = -eta0*(lambda1**2)*Omega31
    eta1 = lambda w: eta0/(1 + 1j*lambda1*w) + etainf
    phi3 = -1*eta3/(eta1(w1)*eta1(w2)*eta1(w3)*eta1(w1+w2+w3))
    J3 = phi3/(1j*(w1+w2+w3))
    return J3

def giesekus_J3(w1,w2,w3,params):
    # Define parameters
    eta0 = params[0]
    lambda1 = params[1]
    etainf = params[2]
    alpha = params[3]

    # Calculate J3 from eta3
    eta3 = giesekus_eta3(w1,w2,w3,params)
    eta1 = lambda w: eta0/(1 + 1j*lambda1*w) + etainf
    phi3 = -1*eta3/(eta1(w1)*eta1(w2)*eta1(w3)*eta1(w1+w2+w3))
    J3 = phi3/(1j*(w1+w2+w3))
    return J3

def crj_eta3(w1,w2,w3,params):
    # Define parameters
    eta0 = params[0]
    lambda1 = params[1]
    lambda2 = params[2]
    
    # Calculate eta3 from the rational functions
    Omega2 = f2(w1,w2,w3,lambda1,lambda2)
    Omega31 = f31(w1,w2,w3,lambda1,lambda2)
    #eta3 = -eta0*(lambda1**2)*Omega31 + eta0*lambda1*lambda2*Omega2
    eta3 = -eta0*lambda1*lambda1*Omega31 + eta0*lambda1*lambda2*Omega2
    return eta3

def crj_J3(w1,w2,w3,params):
    # Define parameters
    eta0 = params[0]
    lambda1 = params[1]
    lambda2 = params[2]

    # Calculate J3 from eta3
    eta3 = crj_eta3(w1,w2,w3,params)
    eta1 = lambda w: eta0*(1 + 1j*lambda2*w)/(1 + 1j*lambda1*w)
    phi3 = -1*eta3/(eta1(w1)*eta1(w2)*eta1(w3)*eta1(w1+w2+w3))
    J3 = phi3/(1j*(w1+w2+w3))
    return J3

def ptt_eta3(w1,w2,w3,params):
    # Define parameters
    eta0 = params[0]
    lambda1 = params[1]
    etainf = params[2]
    xi = params[3]
    epsilon = params[4]

    # Calculate eta3 from the rational functions
    Omega31 = f31(w1,w2,w3,lambda1,0)
    Omega41 = f41(w1,w2,w3,lambda1,0)
    eta3 = -eta0*(lambda1**2)*xi*(2-xi)*Omega31 - 2*eta0*(lambda1**2)*epsilon*(1-xi)*Omega41
    return eta3

def ptt_J3(w1,w2,w3,params):
    # Define parameters
    eta0 = params[0]
    lambda1 = params[1]
    etainf = params[2]
    xi = params[3]
    epsilon = params[4]

    # Calculate J3 from eta3
    eta3 = ptt_eta3(w1,w2,w3,params)
    eta1 = lambda w: eta0/(1 + 1j*lambda1*w)
    phi3 = -1*eta3/(eta1(w1)*eta1(w2)*eta1(w3)*eta1(w1+w2+w3))
    J3 = phi3/(1j*(w1+w2+w3))
    return J3

def johnson_segalman_eta3(w1,w2,w3,params):
    # Define parameters
    eta0 = params[0]
    lambda1 = params[1]
    etainf = param[2]
    xi = params[2]

    # Calculate eta3 from the rational functions
    Omega31 = f31(w1,w2,w3,lambda1,0)
    eta3 = -eta0*(lambda1**2)*xi*(2-xi)*Omega31
    return eta3

def johnson_segalman_J3(w1,w2,w3,params):
    # Define parameters
    eta0 = params[0]
    lambda1 = params[1]
    etainf = params[2]
    xi = params[3]

    # Calculate J3 from eta3
    eta3 = ptt_eta3(w1,w2,w3,params)
    eta1 = lambda w: eta0/(1 + 1*lambda1*w) + etainf
    phi3 = -1*eta3/(eta1(w1)*eta1(w2)*eta1(w3)*eta1(w1+w2+w3))
    J3 = phi3/(1j(w1+w2+w3))
    return J3

def wormlike_eta3(w1,w2,w3,params):
    # Define parameters
    eta0 = params[0]
    lambda1 = params[1]
    etainf = params[2]
    d = params[3]
    r = params[4]

    # Define the birth rate nonlinearities
    f_B = lambda w1,w2,w3: (w1*w2*lambda1**2/(1 + 1j*lambda1*w1))*(1/(1 + 1j*lambda1*(w1+w2))
        - 1/(1 + 1j*lambda1*(w1+w2+w3)))
    G3B = (1/12)*(1 - r)*(eta0/lambda1)*(f_B(w1,w2,w3) + f_B(w1,w3,w2) + f_B(w2,w1,w3) + f_B(w2,w3,w1)
            + f_B(w3,w1,w2) + f_B(w3,w2,w1))
    f_D = lambda w1,w2,w3: (1j*lambda1*w1*w2/(w1+w2))*(1/(1 + 1j*lambda1*w1))*(1/(1 + 1j*lambda1*(w1+w2))
            - 1/(1 + 1j*lambda1*(w1+w2+w3)) + 1/(1 + 1j*lambda1*w3) - 1)
    G3D = (1/12)*(1 + r)*(eta0/lambda1)*(f_D(w1,w2,w3) + f_D(w1,w3,w2) + f_D(w2,w1,w3) + f_D(w2,w3,w1)
            + f_D(w3,w1,w2) + f_D(w3,w2,w1))
    G1 = lambda w: (eta0/lambda1)*(1j*w*lambda1/(1 + 1j*w*lambda1))
    G3Q = d*(G1(w1+w2+w3) - G1(w1+w2) - G1(w1+w3) - G1(w2+w3) + G1(w1) + G1(w2) + G1(w3))

    # Define the total nonlinearity
    G3 = G3B + G3D + G3Q
    eta3 = 1j*G3/(w1*w2*w3)
    return eta3

def wormlike_J3(w1,w2,w3,params):
    # Define parameters
    eta0 = params[0]
    lambda1 = params[1]
    etainf = params[2]

    # Calculate eta3 from the rational functions
    eta3 = wormlike_eta3(w1,w2,w3,params)
    eta1 = lambda w: eta0/(1 + 1j*lambda1*w) + etainf
    phi3 = -1*eta3/(eta1(w1)*eta1(w2)*eta1(w3)*eta1(w1+w2+w3))
    J3 = phi3/(1j*(w1+w2+w3))
    return J3

def tss_eta3(w1,w2,w3,lrmodel,params):
    # Compute G3
    G3 = (lrmodel(w1+w2+w3,params) - lrmodel(w1+w2,params) - lrmodel(w1+w3,params) - lrmodel(w2+w3,params)
            + lrmodel(w1,params) + lrmodel(w2,params) + lrmodel(w3,params))
    eta3 = G3/(-1j*w1*w2*w3)
    return eta3
