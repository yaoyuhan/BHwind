#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 15:09:45 2018

@author: Yuhan Yao, yyao@caltech.edu

Reference: 
    David L. Meier 1979-1982, paper I-IV
    David L. Meier (2012): Black Hole Astrophysics: The Engine Paradigm
    
Aim: Supercritical wind stucture
"""
import os, sys
tempdir = os.getcwd()+'/'
sys.path.append(tempdir)
from slim_ss import find_x_trans_
import numpy as np


# constants
PI = 3.141592653589793          # yao   
MSOL = 1.9884754153381438e+33   # yao: Solar mass in grammes!
MPROT = 1.672621898e-24         # yao: Proton mass
C = 29979245800.0               # yao: Speed of light
G = 6.674079999999999e-08       # yao: Gravitational constant
STEFAN_BOLTZMANN = 5.6696e-5    # yao: Stefan boltzmann constant
NA = 6.022140857e+23            # R = BOLTZMANN(k_B) * NA
BOLTZMANN = 1.38064852e-16      # yao: Boltzmann constant
EV2ERGS = 1.602192e-12          # yao: from eV to erg

K0 = 1.55e+24
epsilon = 0.1
kes = 0.342
mu = 0.61 


def wind_structure(alpha = 0.03, m = 10, mdot = 100,
                   estimate = False, Lbb = 1e+39, kTbb = 0.1, isbb = False):
    # constants from slim disk solution
    C2 = 6 * G * MSOL / C**2
    Risco = C2 * m
    M = m * MSOL
    Ledd = 4 * PI * G * M * C / kes
    Medd = 4 * PI * G * M / (C*kes*epsilon)
    Mdot = Medd * mdot
    L_pt = G * M * Mdot / Risco
    
    # k_H    
    C4 = pow(3*C**5/(4*STEFAN_BOLTZMANN*epsilon*kes*G*MSOL), 0.25) * \
            pow(6, -5/8) # k_Tc
    C5 = C**2 / (epsilon * G * MSOL * kes) * pow(6, -3/2) # k_rho
    C6 = C**4 / (epsilon * G * MSOL * kes) * pow(6, -5/2) # k_P
    C7 = pow(6, -0.5) * C #k_Vr
    C8 = 4./3
    '''
    if isbb==True:
        m = Lbb / (4*PI*G*MSOL*C/kes)
        C1 = 2.4
    
        Tbb = kTbb * 1000 * EV2ERGS / BOLTZMANN # Kelvin
        
        temp1 = C1**(-15./28) * C2**(-4./7) * C4**2 * C5**(-6./7) * (kes*K0)**(-2./7)
        mdot1 = pow(Tbb/temp1, -28/25) * alpha**(2/5) * m**(-6./25)
        
        zeta2 = C8**(1./4) * C6**(1./4) * C5**(-1./4) * C7**(-1./2)
        temp2 = C1**(-5./11) * C2**(-12./11) * C4**(32./11) * C5**(-18./11) *\
                (kes*K0)**(-6./11) * zeta2**(10./11) 
        mdot2 = pow(Tbb/temp2, -11./15) * alpha**(-1./3) * m**(-2./15) 
    '''
    x_trans = find_x_trans_(mdot = mdot)
    x_i = x_trans
    Ri = x_i * m * C2
    f_i = 1-pow(x_i, -0.5)
    
    C1 = x_trans/mdot # ~2.4
    
    if estimate==True:
        C1 = 2.4
    
    # Injection Region
    k_rhoi = pow(C1, -1.5) * C5
    k_Pi = pow(C1, -2.5) * C6
    k_Vri = pow(C1, -0.5) * C7
    k_Ti = pow(C1, -5/8) * C4
    
    rhoi = k_rhoi * pow(m, -1) * pow(mdot, -1/2) * pow(alpha, -1) * pow(f_i, -1/2)
    Pi = k_Pi * pow(m, -1) * pow(mdot, -3/2) * pow(alpha, -1) * pow(f_i, 1/2)
    Vri = k_Vri * pow(mdot, -1/2) * alpha
    Ti = k_Ti * pow(m, -1/4) * pow(mdot, -3/8) * pow(alpha, -1) * pow(f_i, 1/8)

    # Rsg: gas sonic point
    k_Rsg =(BOLTZMANN*NA/mu)**(1./3) * C1**(9./8) * C4**(1/3.)* C7**(-2./3)*C2 
    Rsg = k_Rsg * alpha**(-3./4) * m**(11/12.) * mdot**(29./24) * pow(f_i, 1/24)
    
    # Rad: adiabatic radius       
    k_Rad = (C8*C6/C5)**(0.25) * C7**(-0.5) * C1 * C2
    Rad = k_Rad * m * mdot * alpha**(-0.5) * pow(f_i, 1/4)
    zeta = (C8*C6/C5)**(0.25) * C7**(-0.5) * alpha**(-0.5) * pow(f_i, 1/4)
    
    # Rs: sonic point
    k_Rs = C5**(-1./2) * C6**(1./2) * C7**(-1) * C8**(0.5) * C1*C2
    Rs = k_Rs * alpha**(-1) * m * mdot * pow(f_i, 0.25)
    
    # Rsc: scattersphere
    k_Rsc = C1**(1./2) * C2**2 * kes * C5**(5./4) * C8**(-1./4) * C7**0.5 *\
            C6**(-1./4)  
    Rsc = k_Rsc * alpha**(-1/2) * m * mdot**(3./2) * pow(f_i, -3/4)
    
    # Rstar: photosphere
    ## subsonic solution
    k_Rstar1 = C1**(51/56.) * C2**(11./7) * C4**(-1) * C5**(6./7) * \
                (kes*K0)**(2/7) 
    R_star1 = k_Rstar1 * m**(27./28.) * mdot**(85./56.) * \
                alpha**(-17. / 28.) * pow(f_i, -31/56)
    
    k_Tstar1 = pow(C1, -15/28) * (kes*K0)**(-2/7)  * C4**2 * pow(C2, -4/7)*\
                pow(C5, -6/7)
    T_star1 = k_Tstar1 * pow(alpha, 5/14) * pow(m, -3/14) * pow(mdot, -25/28)*\
                pow(f_i, 17/28)
    
    k_tau_es_star1 = pow(kes, 3/7) * pow(K0, -4/7) * pow(C1, -9/28) *\
                        pow(C2, -1/7) * pow(C5, -5/7) * pow(C4, 2) 
    tau_es_star1 = k_tau_es_star1 * pow(alpha, 3/14) * pow(mdot, -15/28) *\
                    pow(m, 1/14) * pow(f_i, 17/28)
    if R_star1 < Rad:
        R_star = R_star1
        T_star = T_star1
        tau_es_star = tau_es_star1
        subCase = 'sub'
    ## supersonic solution
    else:
        R_star = pow(K0*kes, 8/11) * pow(rhoi, 24/11) * pow(Ti, -28/11) *\
                    pow(Ri, 27/11) * pow(zeta, -17/11)
        
        T_star = Ti * zeta**(-1/4) * (R_star/Ri)**(-3/4)        
        tau_es_star = kes * rhoi * pow(zeta, -3) / R_star * Ri**2
        subCase = 'super'
    
    info = {'alpha':        alpha,
            'm':            m,
            'mdot':         mdot,
            'L_pt':         L_pt,
            'Ledd':         Ledd,
            'M':            M,
            'Mdot':         Mdot,
            'Ri':           Ri,
            'rhoi':         rhoi,
            'Risco':        Risco,
            'Rsg':          Rsg,
            'Rad':          Rad,
            'zeta':         zeta,
            'Rs':           Rs,
            'Rsc':          Rsc,
            'R_star':       R_star,
            'T_star':       T_star,
            'tau_es_star':  tau_es_star,
            'subCase':      subCase}
    return info
    
    
    
def find_TOUTSIDE(info1, info2, F_mid, R_mid, H_mid, const_mid, beta,                  
                  n_theta = 100, n_phi = 200, quick = True):
    '''
    Derive outflow irradiationg by Integration.
    
    Parameters:
        info1: standard disk structure, returned from function 
               standard_structure in the slim_ss package
    '''
    # plt.loglog(R_mid, F_mid)
    # find out Hsc
    # R_trans = info1['R_trans']
    R_mi = info1['R_mi']
    R_om = info1['R_om']
    alpha = info2['alpha']
    m = info2['m']
    mdot = info2['mdot']
    Rsc = info2['Rsc']
    Ledd = info2['Ledd']
    
    xsc = Rsc / info2['Risco']
    fsc = 1 - pow(xsc, -0.5)
    if Rsc < R_mi:
        print ('@Yao: error occurs: Rsc inside of R_mi!')
    elif Rsc < R_om:
        print ('@Yao: R_mi < Rsc < R_om')
        k_H3 = pow(C, -2.3) * pow(G*MSOL, 0.9) * pow(8/9*STEFAN_BOLTZMANN, -0.1) *\
            pow(6, 21./20) * pow(epsilon, -0.2) * pow(mu*MPROT/BOLTZMANN, -0.4)*\
            pow(kes,-0.1)
        Hsc = k_H3 * pow(alpha, -0.1) * pow(m, 0.9) * pow(mdot*fsc, 0.2) * pow(xsc, 21./20)
    else:
        print ('@Yao: Rsc > R_om')
        k_H4 = pow(C, -12/5) * pow(G*MSOL, 0.9) * pow(8*STEFAN_BOLTZMANN/9/K0, -1/20)*\
            pow(6, 9/8)*pow(epsilon*kes, -3/20) * pow(mu*MPROT/BOLTZMANN, -3/8)
        Hsc = k_H4 * pow(alpha, -0.1) * pow(m, 0.9) * pow(mdot*fsc, 3/20) * pow(xsc, 9/8)
    
    # integration angle limits
    Theta_0 = np.ones(len(F_mid))*PI/2  
    Theta_1 = np.zeros(len(F_mid))    
    Theta_2 = np.ones(len(F_mid))*PI/2  
    ix = R_mid > 1.01*Rsc
    Theta_0[ix] = PI/2 - np.arcsin(Hsc / Rsc)   
    Theta_2[ix] = PI/2 - np.arctan(H_mid[ix] / R_mid[ix])
    Theta_1[ix] = Theta_2[ix] - np.arccos(Rsc / np.sqrt(R_mid[ix]**2 + H_mid[ix]**2))
    F_mid[ix] = 10
    
    if quick==True:
        id_start = np.where(ix==1)[0][0]
        id_end = np.where(ix==1)[0][-1]
        ids_float_ = np.logspace(np.log10(id_start),np.log10(id_end),num=15)
        ids_float = np.hstack([ids_float_,
                               ids_float_[0]*0.8+ids_float_[1]*0.2,
                               ids_float_[0]*0.6+ids_float_[1]*0.4,
                               ids_float_[0]*0.3+ids_float_[1]*0.7,
                               ids_float_[-2]*0.8+ids_float_[-1]*0.2,
                               ids_float_[-2]*0.6+ids_float_[-1]*0.4,
                               ids_float_[-2]*0.3+ids_float_[-1]*0.7])
        ids_float = ids_float[np.argsort(ids_float)]
        ids_int = np.array(ids_float,dtype=int)
    for i in range(len(F_mid)):
        if (quick==True and i in ids_int) or (quick==False and ix[i]==1):
            Rnow = R_mid[i]
            Hnow = H_mid[i]        
            OP = np.mat([[Rnow,0.,Hnow]])
            anow = np.mat([[1.,0.,-(1./(1+const_mid[i]))*Rnow/Hnow]])    
            thetas = np.linspace(Theta_1[i], Theta_0[i], num=n_theta, endpoint=True)
            d_theta = (Theta_0[i]-Theta_1[i])/len(thetas)
            for j in range(n_theta):
                thetanow = thetas[j]
                phis = np.linspace(-(PI/2-Theta_1[i]), PI/2-Theta_1[i], 
                                   num=n_phi, endpoint=True)
                d_phi = (PI - 2*Theta_1[i])/n_phi
                for k in range(n_phi):
                    phinow = phis[k]
                    OQ = np.mat([[Rsc*np.sin(thetanow)*np.cos(phinow),
                                  Rsc*np.sin(thetanow)*np.sin(phinow),
                                  Rsc*np.cos(thetanow)]])
                    QP = OP-OQ
                    cosPHI = QP*anow.T/np.sqrt(QP*QP.T*anow*anow.T)
                    D2 = QP*QP.T
                    F_mid[i] += Ledd * np.sin(thetanow) * d_theta * \
                                d_phi / D2 * cosPHI    
        if quick==True:
            print (i)
        elif i%50==0:
            print (i)
    if quick==False:             
        F_mid[ix] *= (1-beta) / (4*PI)**2 *2
        return F_mid
    else:
        F_mid[ids_int] *= (1-beta) / (4*PI)**2 *2
        
        
        return F_mid