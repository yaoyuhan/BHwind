#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 23 13:24:22 2018

Author: Yuhan Yao, yyao@caltech.edu

Reference: 
    Frank et al. (2002): Accretion Power in Astrophysics
    Watarai, 2006, ApJ, 648, 523
    Kato et al. (2008): Blakc-Hole Accretion Disks -- Towards a New Paradigm
    David L. Meier (2012): Black Hole Astrophysics: The Engine Paradigm

Aim: Three typical regions of standard steady black hole disk +
     Self-similar solutions of slim disk.
     Derive the coefficients of basic physical quantities.
"""

import numpy as np


# constants
MSOL = 1.9884754153381438e+33   # yao: Solar mass in grammes!
MPROT = 1.672621898e-24         # yao: Proton mass
STEFAN_BOLTZMANN = 5.6696e-5    # yao: Stefan boltzmann constant
C = 29979245800.0               # yao: Speed of light
G = 6.674079999999999e-08       # yao: Gravitational constant
PI = 3.141592653589793          # yao   
BOLTZMANN = 1.38064852e-16      # yao: Boltzmann constant

K0 = 1.55e+24
X = 0.71
mu = 0.61 
kes = 0.342
epsilon=0.1

    
def find_x_trans_(mdot=100):
    x = np.logspace(0.1, 4)
    f = 1-pow(x, -0.5)
    const2 = pow(mdot/4/epsilon , 2)
    y = abs(x**2/(const2*f)-1)
    if min(y)>0.001:
        xbound = x[np.argsort(y)[:2]]
        x = np.linspace(min(xbound), max(xbound))
        f = 1-pow(x, -0.5)
        const2 = pow(mdot/4/epsilon , 2)
        y = abs(x**2/(const2*f)-1)
    result = x[np.argsort(y)[0]]
    return result


def find_x_MI_(mdot=100, alpha=0.03, m=10):
    x = np.logspace(1, 6)
    f = 1-pow(x, -0.5)
    const_RMI = pow(BOLTZMANN/(mu*MPROT), -8/21) * pow(C, 2/7) * \
                pow(epsilon, -16/21) * pow(kes * G*MSOL, 2/21) /6 *\
                pow(8*STEFAN_BOLTZMANN/9, -6/7) * pow(4*STEFAN_BOLTZMANN/3, 20/21)
    temp = const_RMI*pow(m*alpha, 2/21)*pow(mdot*f, 16/21)
    y = abs(x/temp-1)
    if min(y)>0.001:
        xbound = x[np.argsort(y)[:2]]
        x = np.linspace(min(xbound), max(xbound))
        f = 1-pow(x, -0.5)
        temp = const_RMI*pow(m*alpha, 2/21)*pow(mdot*f, 16/21)
        y = abs(x/temp-1)
    result = x[np.argsort(y)[0]]
    return result


def find_x_OM_(mdot=100, m=10):
    x = np.logspace(1, 6)
    f = 1-pow(x, -0.5)
    const_ROM = pow(BOLTZMANN/(mu*MPROT), 1/3) /6. * \
                pow(kes*9*C**2/K0/epsilon/8/STEFAN_BOLTZMANN, 2/3)
    temp = const_ROM*pow(mdot*f, 2/3)
    y = abs(x/temp-1)
    if min(y)>0.001:
        xbound = x[np.argsort(y)[:2]]
        x = np.linspace(min(xbound), max(xbound))
        f = 1-pow(x, -0.5)
        temp = const_ROM*pow(mdot*f, 2/3)
        y = abs(x/temp-1)
    result = x[np.argsort(y)[0]]
    return result

    
def standard_structure(alpha = 0.03, m = 10, mdot = 100, Rout = None):
    x_trans = find_x_trans_(mdot = mdot)
    
    x_MI = find_x_MI_(mdot=mdot, alpha=alpha, m=m)
    
    x_OM = find_x_OM_(mdot=mdot, m=m)
    
    Risco = 6*G*MSOL*m/C**2
    if Rout == None:
        Rout = Risco * 5e+5
    
    xout = Rout/Risco
    ###########################################################################
    # 1. slim innermost disk: kes >> kff/bf, Pr >> Pg, Qadv >~ Qrad
    ###########################################################################
    x1 = np.logspace(0.0001, np.log10(x_trans), num=100, endpoint=False)
    f1 = 1 - pow(x1, -0.5)
    R1 = x1 * Risco
    ff=1
    
    # P
    k_P1 = C**4 / (epsilon * G * MSOL * kes) * pow(6, -5/2)
    P1 = k_P1 * pow(m*alpha, -1) * pow(x1, -2.5) * mdot * pow(f1/ff, 0.5)
    
    # rho
    k_rho1 = C**2 / (epsilon * G * MSOL * kes) * pow(6, -3/2)
    rho1 = k_rho1 * pow(m*alpha, -1) * pow(x1*ff, -1.5) * mdot * pow(f1, -0.5)
    
    # H
    k_H1 = 6 * G* MSOL / C**2
    H1 = k_H1 * m * x1 * pow(f1*ff, 0.5)
    
    # T
    k_Tc1 = pow(3*C**5/(4*STEFAN_BOLTZMANN*epsilon*kes*G*MSOL), 0.25) * pow(6, -5/8)
    Tc1 = k_Tc1 * pow(mdot/alpha/m, 0.25) * pow(f1/ff, 1/8) * pow(x1, -5/8)
    
    # plt.loglog(R1, H1)
    '''
    # Vr
    k_Vr1 = pow(6, -0.5) * C 
    Vr1 = k_Vr1 * alpha * ff * pow(x1, -1/2)
    '''
    
    ###########################################################################
    # 2. standard inner disk: kes >> kff/bf, Pr >> Pg
    ###########################################################################
    x2 = np.logspace(np.log10(x_trans), np.log10(x_MI), num =200, endpoint=False)
    f2 = 1 - pow(x2, -0.5)
    R2 = x2 * Risco
    
    # P
    k_P2 = C**4 / (9 * np.sqrt(6) * G * MSOL * kes)
    P2 = k_P2 * pow(m*alpha, -1) * pow(x2, -1.5)
    
    # rho
    k_rho2 = 16 * np.sqrt(6) / 9. * (epsilon * C)**2 / (G * MSOL * kes)
    rho2 = k_rho2 * pow(mdot*f2, -2) * pow(alpha*m, -1) * pow(x2, 1.5)
    
    # H
    k_H2 = 3 * G * MSOL / (2 * epsilon * C**2)
    H2 = k_H2 * m * mdot * f2
    
    # Tc
    k_Tc2 = pow(C**5 / (12 * np.sqrt(6) * G * MSOL * kes * STEFAN_BOLTZMANN), 1/4)
    Tc2 = k_Tc2 * pow(m*alpha, -0.25) * pow(x2, -3/8)
    
    # plt.loglog(R2, H2)
    
    ###########################################################################
    # 3. standard middle disk: kes >> kff/bf, Pr << Pg
    ###########################################################################
    x3 = np.logspace(np.log10(x_MI), np.log10(x_OM), num=200, endpoint=False)
    f3 = 1 - pow(x3, -0.5)
    R3 = x3 * Risco
    
    # P
    k_P3 = pow(C, 4.3) * pow(G*MSOL*kes, -0.9) * pow(8/9*STEFAN_BOLTZMANN, 0.1) * \
            pow(6, -51./20) * pow(epsilon, -0.8) * pow(mu*MPROT/BOLTZMANN, 0.4)
    P3 = k_P3 * pow(m*alpha, -0.9) * pow(mdot*f3, 0.8) * pow(x3, -51./20)
    
    # rho
    k_rho3 = pow(mu*MPROT/BOLTZMANN, 1.2) * pow(8/9*STEFAN_BOLTZMANN, 0.3)*\
                pow(C, 2.9) * pow(G*MSOL*kes, -0.7) * pow(epsilon, -0.4)*\
                pow(6, -33./20)
    rho3 = k_rho3 * pow(m*alpha, -0.7) * pow(mdot*f3, 0.4) * pow(x3, -33./20)
    
    # H
    k_H3 = pow(C, -2.3) * pow(G*MSOL, 0.9) * pow(8/9*STEFAN_BOLTZMANN, -0.1) *\
            pow(6, 21./20) * pow(epsilon, -0.2) * pow(mu*MPROT/BOLTZMANN, -0.4)*\
            pow(kes,-0.1)
    H3 = k_H3 * pow(alpha, -0.1) * pow(m, 0.9) * pow(mdot*f3, 0.2) * pow(x3, 21./20)
    
    # Tc
    k_Tc3 = pow(G*MSOL*kes, -0.2) * pow(epsilon, -0.4) * pow(C, 1.4) * \
            pow(mu*MPROT/BOLTZMANN, 0.2) * pow(8/9*STEFAN_BOLTZMANN, -0.2) *\
            pow(6, -0.9)
    Tc3 = k_Tc3 * pow(m*alpha, -0.2) * pow(mdot*f3, 0.4) * pow(x3, -0.9)   
    
    # plt.loglog(R3, H3)
    
    ###########################################################################
    # 4. standard outer disk: kes << kff/bf, Pr << Pg
    ###########################################################################
    x4 = np.logspace(np.log10(x_OM), np.log10(xout), num=700, endpoint=True)
    f4 = 1 - pow(x4, -0.5)
    R4 = x4 * Risco
    
    # P
    k_P4 = pow(mu*MPROT/BOLTZMANN, 3/8) * pow(8*STEFAN_BOLTZMANN/9/K0, 1/20)*\
                pow(6, -21./8) * pow(C, 22/5) * pow(epsilon*kes, -17/20)*\
                pow(G*MSOL, -0.9)
    P4 = k_P4 * pow(m*alpha, -0.9) * pow(mdot*f4, 17/20) * pow(x4, -21./8)
    
    # rho
    k_rho4 = pow(mu*MPROT/BOLTZMANN, 9/8) * pow(8*STEFAN_BOLTZMANN/9/K0,3/20)*\
                pow(6, -15./8) * pow(C, 16/5) * pow(epsilon*kes, -11/20)*\
                pow(G*MSOL, -0.7)
    rho4 = k_rho4 * pow(m*alpha, -0.7) * pow(mdot*f4, 11/20) * pow(x4, -15./8)
    
    # H
    k_H4 = pow(C, -12/5) * pow(G*MSOL, 0.9) * pow(8*STEFAN_BOLTZMANN/9/K0, -1/20)*\
            pow(6, 9/8)*pow(epsilon*kes, -3/20) * pow(mu*MPROT/BOLTZMANN, -3/8)
    H4 = k_H4 * pow(alpha, -0.1) * pow(m, 0.9) * pow(mdot*f4, 3/20) * pow(x4, 9/8)
    
    # Tc
    k_Tc4 = pow(epsilon*kes, -0.3)* pow(C,6/5) * pow(G*MSOL,-1/5)*pow(6, -3/4)*\
            pow(mu*MPROT/BOLTZMANN, 0.25) * pow(8*STEFAN_BOLTZMANN/9/K0, -0.1)
    Tc4 = k_Tc4 * pow(m*alpha, -0.2) * pow(mdot*f4, 0.3) * pow(x4, -3/4)
    
    # plt.loglog(R4, H4)
    
    R = np.hstack([R1, R2, R3, R4])
    rho = np.hstack([rho1, rho2, rho3, rho4])
    H = np.hstack([H1, H2, H3, H4])
    P = np.hstack([P1, P2, P3, P4])
    Tc = np.hstack([Tc1, Tc2, Tc3, Tc4])
    
    # define const1 = ( d(lnH)/d(lnR) - 1 )
    const1 = np.hstack([np.zeros(R1.shape),
                        np.zeros(R2.shape),
                        np.ones(R3.shape)*(1/20.),
                        np.ones(R4.shape)*(1/8)])
    
    # plt.loglog(R, P) 
    # plt.loglog(R, Tc)     
    # plt.plot(R, H) 
    # plt.plot(R, rho) 
    
    info = {'R':        R,
            'rho':      rho,
            'H':        H,
            'P':        P,
            'Tc':       Tc,
            'R_trans':  R2[0],
            'R_mi':     R3[0],
            'R_om':     R4[0],
            'const1':   const1}
    return info