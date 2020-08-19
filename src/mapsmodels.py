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

def f1(w1,w2,w3,lambda1,lambda2,zeta,r,s):
    # Defines the Omega1 rational function
    f = lambda w,rzeta,s: (1/(1 + rzeta + 1j*lambda1*w))*(1/(1 + s*1j*lambda1*w))
    return f(w1+w2+w3,0,0)*(f(w1+w2,r*zeta,s) + f(w1+w3,r*zeta,s) + f(w2+w3,r*zeta,s))/3

def f2(w1,w2,w3,lambda1,lambda2,zeta,r,s):
    # Defines the Omega2 rational function
    f = lambda w,rzeta,s: 1/(1 + 1j*lambda1*w)
    g = lambda w: (1 + 1j*lambda2*w)/(1 + 1j*lambda1*w)
    return f(w1+w2+w3,0,0)*(f(w1+w2,r*zeta,s)*(g(w1)+g(w2)) + f(w1+w3,r*zeta,s)*(g(w1)+g(w3)) +
            f(w2+w3,r*zeta,s)*(g(w2)+g(w3)))/6

def f3(w1,w2,w3,lambda1,lambda2,zeta,r,s):
    # Defines the Omega3 rational function
    f = lambda w,rzeta,s: 1/(1 + 1j*lambda1*w)
    g = lambda w: (1 + 1j*lambda2*w)/(1 + 1j*lambda1*w)
    return f(w1+w2+w3,0,0)*(f(w1+w2,r*zeta,s)*g(w3) + f(w1+w3,r*zeta,s)*g(w2) + f(w2+w3,r*zeta,s)*g(w1))/3

def f4(w1,w2,w3,lambda1,lambda2,zeta,r,s):
    # Defines the Omega4 rational function
    f = lambda w,rzeta,s: 1/(1 + 1j*lambda1*w)
    g = lambda w: (1 + 1j*lambda2*w)/(1 + 1j*lambda1*w)
    return f(w1+w2+w3,0,0)*(g(w3)*f(w1+w2,r*zeta,s)*(g(w1)+g(w2)) + g(w2)*f(w1+w3,r*zeta,s)*(g(w1)+g(w3)) +
            g(w1)*f(w2+w3,r*zeta,s)*(g(w2)+g(w3)))/6

def f5(w1,w2,w3,lambda1,lambda2,zeta,r,s):
    # Defines the Omega5 rational function
    f = lambda w,rzeta,s: 1/(1 + 1j*lambda1*w)
    g = lambda w: (1 + 1j*lambda2*w)/(1 + 1j*lambda1*w)
    return f(w1+w2+w3,0,0)*(f(w1+w2,r*zeta,s)*g(w1)*g(w2) + f(w1+w3,r*zeta,s)*g(w1)*g(w3) + f(w2+w3,r*zeta,s)*g(w2)*g(w3))/3

def f6(w1,w2,w3,lambda1,lambda2,zeta,r,s):
    # Defines the Omega6 rational function
    f = lambda w,rzeta,s: 1/(1 + 1j*lambda1*w)
    g = lambda w: (1 + 1j*lambda2*w)/(1 + 1j*lambda1*w)
    return f(w1+w2+w3,0,0)*(f(w1+w2,r*zeta,s) + f(w1+w3,r*zeta,s) + f(w2+w3,r*zeta,s))*g(w1)*g(w2)*g(w3)/3

def f7(w1,w2,w3,lambda1,lambda2,zeta,r,s):
    # Defines the Omega7 rational function
    f = lambda w,rzeta,s: 1/(1 + 1j*lambda1*w)
    return f(w1+w2+w3,0,0)

def f8(w1,w2,w3,lambda1,lambda2,zeta,r,s):
    # Defines the Omega8 rational function
    f = lambda w,rzeta,s: 1/(1 + 1j*lambda1*w)
    g = lambda w: (1 + 1j*lambda2*w)/(1 + 1j*lambda1*w)
    return f(w1+w2+w3,0,0)*(g(w1) + g(w2) + g(w3))/3

def f9(w1,w2,w3,lambda1,lambda2,zeta,r,s):
    # Defines the Omega9 rational function
    f = lambda w,rzeta,s: 1/(1 + 1j*lambda1*w)
    g = lambda w: (1 + 1j*lambda2*w)/(1 + 1j*lambda1*w)
    return f(w1+w2+w3,0,0)*(g(w1)*g(w2) + g(w2)*g(w3) + g(w3)*g(w1))/3

def f10(w1,w2,w3,lambda1,lambda2,zeta,r,s):
    # Defines the Omega10 rational function
    f = lambda w,rzeta,s: 1/(1 + 1j*lambda1*w)
    g = lambda w: (1 + 1j*lambda2*w)/(1 + 1j*lambda1*w)
    return f(w1+w2+w3,0,0)*g(w1)*g(w2)*g(w3)

def quadratic_maxwell_eta3(w1,w2,w3,params):
    # Define parameters
    eta0 = params[0]
    lambda1 = params[1]
    etainf = params[2]
    lambda2 = params[3]
    A = params[4]

    # Defines the '11 constant' model
    g2 = f1(w1,w2,w3,lambda1,lambda2,0,0,0)
    g31 = f2(w1,w2,w3,lambda1,lambda2,0,0,0)
    g32 = f3(w1,w2,w3,lambda1,lambda2,0,0,0)
    g41 = f4(w1,w2,w3,lambda1,lambda2,0,0,0)
    g42 = f5(w1,w2,w3,lambda1,lambda2,0,0,0)
    g5 = f6(w1,w2,w3,lambda1,lambda2,0,0,0)
    return eta0*(A[0]*g2 + A[1]*g31 + A[2]*g32 + A[3]*g41 + A[4]*g42 + A[5]*g5)

def quadratic_maxwell_J3(w1,w2,w3,params):
    # Define parameters
    eta0 = params[0]
    lambda1 = params[1]
    etainf = params[2]
    lambda2 = params[3]

    # Calculate J3 from eta3
    eta3 = quadratic_maxwell_eta3(w1,w2,w3,params)
    eta1 = lambda w: eta0*(1 + 1j*lambda2*w)/(1 + 1j*lambda1*w) + etainf
    phi3 = -1*eta3/(eta1(w1)*eta1(w2)*eta1(w3)*eta1(w1+w2+w3))
    J3 = phi3/(1j*(w1+w2+w3))
    return J3

def cubic_maxwell_eta3(w1,w2,w3,params):
    # Define parameters
    eta0 = params[0]
    lambda1 = params[1]
    etainf = params[2]
    lambda2 = params[3]
    zeta = params[4]
    a0 = params[5]
    a3 = params[6]
    b = params[7]

    # Compute the sum over a0 terms
    a0_terms = (a0[0]*f1(w1,w2,w3,lambda1,lambda2,zeta,0,0) + a0[1]*f2(w1,w2,w3,lambda1,lambda2,zeta,0,0) 
            + a0[2]*f3(w1,w2,w3,lambda1,lambda2,zeta,0,0) + a0[3]*f4(w1,w2,w3,lambda1,lambda2,zeta,0,0)
            + a0[4]*f5(w1,w2,w3,lambda1,lambda2,zeta,0,0) + a0[5]*f6(w1,w2,w3,lambda1,lambda2,zeta,0,0)
            + a0[6]*f7(w1,w2,w3,lambda1,lambda2,zeta,0,0) + a0[7]*f8(w1,w2,w3,lambda1,lambda2,zeta,0,0)
            + a0[8]*f9(w1,w2,w3,lambda1,lambda2,zeta,0,0) + a0[9]*f10(w1,w2,w3,lambda1,lambda2,zeta,0,0))

    # Compute the sum over a3 terms
    a3_terms = (a3[0]*f1(w1,w2,w3,lambda1,lambda2,zeta,3,0) + a3[1]*f2(w1,w2,w3,lambda1,lambda2,zeta,3,0)
             + a3[2]*f3(w1,w2,w3,lambda1,lambda2,zeta,3,0) + a3[3]*f4(w1,w2,w3,lambda1,lambda2,zeta,3,0)
             + a3[4]*f5(w1,w2,w3,lambda1,lambda2,zeta,3,0) + a3[5]*f6(w1,w2,w3,lambda1,lambda2,zeta,3,0))

    # Compute the sum over b terms
    b_terms = (b[0]*f1(w1,w2,w3,lambda1,lambda2,zeta,3,1) + b[1]*f2(w1,w2,w3,lambda1,lambda2,zeta,3,1)
             + b[2]*f3(w1,w2,w3,lambda1,lambda2,zeta,3,1) + b[3]*f4(w1,w2,w3,lambda1,lambda2,zeta,3,1)
             + b[4]*f5(w1,w2,w3,lambda1,lambda2,zeta,3,1) + b[5]*f6(w1,w2,w3,lambda1,lambda2,zeta,3,1))

    return eta0*(a0_terms + a3_terms + b_terms)

def cubic_maxwell_J3(w1,w2,w3,params):
    # Define parameters
    eta0 = params[0]
    lambda1 = params[1]
    etainf = params[2]
    lambda2 = params[3]

    # Calculate J3 from eta3
    eta3 = cubic_maxwell_eta3(w1,w2,w3,params)
    eta1 = lambda w: eta0*(1 + 1j*lambda2*w)/(1 + 1j*lambda1*w) + etainf
    phi3 = -1*eta3/(eta1(w1)*eta1(w2)*eta1(w3)*eta1(w1+w2+w3))
    J3 = phi3/(1j*(w1+w2+w3))
    return J3

def giesekus_eta3(w1,w2,w3,params):
    # Define parameters
    eta0 = params[0]
    lambda1 = params[1]
    etainf = params[2]
    alpha = params[3]

    # Calculate eta3 from the rational functions
    Omega4 = f4(w1,w2,w3,lambda1,0,0,0,0)
    Omega5 = f5(w1,w2,w3,lambda1,0,0,0,0)
    Omega6 = f6(w1,w2,w3,lambda1,0,0,0,0)
    eta3 = eta0*(lambda1**2)*(-2*alpha*Omega4 - alpha*Omega5 + 2*(alpha**2)*Omega6)
    return eta3

def crm_eta3(w1,w2,w3,params):
    # Define parameters
    eta0 = params[0]
    lambda1 = params[1]
    etainf = params[2]

    # Calculate eta3 from the rational functions
    Omega2 = f2(w1,w2,w3,lambda1,0,0,0,0)
    eta3 = -eta0*(lambda1**2)*Omega2
    return eta3

def crm_J3(w1,w2,w3,params):
    # Define parameters
    eta0 = params[0]
    lambda1 = params[1]
    etainf = params[2]

    # Calculate J3 from eta3
    Omega2 = f2(w1,w2,w3,lambda1,0,0,0,0)
    eta3 = -eta0*(lambda1**2)*Omega2
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
    Omega1 = f1(w1,w2,w3,lambda1,lambda2,0,0,0)
    Omega2 = f2(w1,w2,w3,lambda1,lambda2,0,0,0)
    eta3 = -eta0*lambda1*lambda1*Omega2 + eta0*lambda1*lambda2*Omega1
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
    Omega2 = f2(w1,w2,w3,lambda1,0,0,0,0)
    Omega4 = f4(w1,w2,w3,lambda1,0,0,0,0)
    eta3 = -eta0*(lambda1**2)*xi*(2-xi)*Omega2 - 2*eta0*(lambda1**2)*epsilon*(1-xi)*Omega4
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
    Omega2 = f2(w1,w2,w3,lambda1,0,0,0,0)
    eta3 = -eta0*(lambda1**2)*xi*(2-xi)*Omega2
    return eta3

def johnson_segalman_J3(w1,w2,w3,params):
    # Define parameters
    eta0 = params[0]
    lambda1 = params[1]
    etainf = params[2]
    xi = params[3]

    # Calculate J3 from eta3
    eta3 = johnson_segalman_eta3(w1,w2,w3,params)
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

def tss_eta3(w1,w2,w3,lrmodel):
    # Compute G3
    G3 = (lrmodel(w1+w2+w3) - lrmodel(w1+w2) - lrmodel(w1+w3) - lrmodel(w2+w3)
            + lrmodel(w1) + lrmodel(w2) + lrmodel(w3))
    eta3 = G3/(-1j*w1*w2*w3)
    return eta3

def tss_J3(w1,w2,w3,lrmodel):
    # Compute eta3
    eta3 = tss_eta3(w1,w2,w3,lrmodel)

    # Compute phi3 and J3
    eta1 = lambda w: lrmodel(w)/(1j*w)
    phi3 = -1*eta3/(eta1(w1)*eta1(w2)*eta1(w3)*eta1(w1+w2+w3))
    J3 = phi3/(1j*(w1+w2+w3))
    return J3
