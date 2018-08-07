#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 00:04:47 2018

@author: yaoyuhan
"""
import os, sys
tempdir = os.getcwd()+'/'
sys.path.append(tempdir)
from wind import wind_structure, find_TOUTSIDE
from slim_ss import standard_structure
import numpy as np
ptex = tempdir.split('BHwind')[0]+'supertex/'

MSOL = 1.9884754153381438e+33   # yao: Solar mass in grammes!
MPROT = 1.672621898e-24         # yao: Proton mass
MELECTRON = 9.10938356e-28      # yao: Electron mass
YEAR = 3.1556926e7              # yao: Number of seconds in a year (from google)
STEFAN_BOLTZMANN = 5.6696e-5    # yao: Stefan boltzmann constant
G = 6.674079999999999e-08       # yao: Gravitational constant
PI = 3.141592653589793          # yao   
BOLTZMANN = 1.38064852e-16      # yao: Boltzmann constant
H = 6.62607004e-27              # Yao: Planck's constant
C = 29979245800.0               # yao: Speed of light
PC = 3.085677581467192e+18      # yao: pc
THOMPSON = 6.6524587158e-25     # yao: Thompson cross-section for an electron (from wiki)
EV2ERGS = 1.602192e-12          # yao: from eV to erg
ANGSTROM = 1.0e-8               # yao: Definition of an Angstrom in units of this code, e.g. cm 
NA = 6.022140857e+23            # R = BOLTZMANN(k_B) * NA

epsilon=0.1
mu = 0.61                       # Mean molecular weight
kes = 0.342
K0 = 1.55e+24


def planck_lambda(T, lamb):
    # convert to cm for planck equation
    lamb2 = lamb * ANGSTROM
    x = H * C / (BOLTZMANN * T * lamb2)
    x = np.array(x, dtype = np.float128)
    Bnu = (2. * H * C**2 ) /  (lamb2**5. ) / (np.exp(x) - 1. )
    # convert back to ANGSTROM   
    return Bnu*ANGSTROM


class outSpec_comp():
    def __init__(self, wave, specw, specw1, specw2_inside, specw2_outside,
                 R_mid, F_mid, F_intrinsic):
        self.wave = wave
        self.flux = specw
        self.flux_intrinsic = specw1
        self.flux_irrin = specw2_inside
        self.flux_irrout = specw2_outside
        self.R_mid = R_mid
        self.F_mid = F_mid
        self.F_intrinsic = F_intrinsic
        
        
def standirr_spec(alpha=0.03, m = 10, mdot = 100, Rout = None,
                  beta = 0.7, nfreq = 1000, w1 = 1000, w2 = 16000):
    '''
    Compute the SED of disks irradiated by central X-ray accretion luminosity
    '''
    info1 = standard_structure(alpha = alpha, m = m, mdot = mdot, Rout = Rout)
    info2 = wind_structure(alpha = alpha, m = m, mdot = mdot)
    
    R_star = info2['R_star']
    Risco = info2['Risco']
    Mdot = info2['Mdot']
    M = info2['M']
    L_pt = info2['L_pt']
    
    if Rout==None:
        Rout = Risco * 5e+5
    else:
        Rout = Rout
    
    # ---- disk radii from R_star to R_out ----
    R_start = max(R_star, Risco)
    if R_start==Risco:
        print ("@Yao: error -- Rmin = Risco")
        
    R_mid = info1['R']
    ix = R_mid > R_start
    R_mid = R_mid[ix]
    const1_mid = info1['const1'][ix]
    H_mid = info1['H'][ix]
    x_mid = R_mid / Risco    
    f_mid = 1 - x_mid**(-0.5)
    
    R_left = np.hstack([1.5*R_mid[0]-0.5*R_mid[1], (R_mid[:-1]+R_mid[1:])*0.5])
    R_right = np.hstack([(R_mid[:-1]+R_mid[1:])*0.5, 1.5*R_mid[-1]-0.5*R_mid[-2]])
    
    # ---- intrinsic flux ---- 
    Tisco4 = 3 * G * M * Mdot * (8 * PI * STEFAN_BOLTZMANN * Risco**3)**(-1) 
    T_intrinsic = pow(Tisco4 * x_mid**(-3) * f_mid, 0.25)
    # plt.loglog(R_mid, T_intrinsic)
    
    # ---- irradiated flux ----
    F_mid = L_pt / (4*PI*R_mid**2) * (1-beta) * H_mid / R_mid * const1_mid
    T_outflow = pow(F_mid/STEFAN_BOLTZMANN, 0.25)
    # plt.loglog(R_mid, T_outflow)
    
    # ---- add them together ----
    T_disk = (T_intrinsic**4 + T_outflow**4) ** 0.25
    # plt.loglog(R_mid, T_disk)
    
    area = PI * ( R_right**2.0 - R_left**2.0) * 2.0 # in cm^2
    wave = np.logspace(np.log10(w1), np.log10(w2), nfreq)
    specw = np.zeros(len(wave))
    specw1 = np.zeros(len(wave))
    specw2 = np.zeros(len(wave))
    
    for i in range(len(T_disk)):
        specw += planck_lambda(T_disk[i], wave)*area[i]
        specw1 += planck_lambda(T_intrinsic[i], wave)*area[i]
        specw2 += planck_lambda(T_outflow[i], wave)*area[i]
    D = 10*PC
    specw *= 1./D**2
    specw1 *= 1./D**2
    specw2 *= 1./D**2
    newinfo = {'wave':              wave,
               'flux':              specw,
               'flux_intrinsic':    specw1,
               'flux_irr':          specw2}
    info2.update(newinfo)
    return info2     
    
        
def windirr_spec(alpha=0.03, m = 10, mdot = 100, Rout=None,
                 beta = 0.7, nfreq = 1000, w1 = 1000, w2 = 16000,
                 quick=True):
    '''
    Compute the SED of disks irradiated by an outflow wind
    (Assuming all central accretion luminosity is blocked by wind)
    '''
    info1 = standard_structure(alpha = alpha, m = m, mdot = mdot, Rout = Rout)
    info2 = wind_structure(alpha = alpha, m = m, mdot = mdot)
    
    subCase = info2['subCase']
    Rad = info2['Rad']
    Rsc = info2['Rsc']
    R_star = info2['R_star']
    T_star = info2['T_star']
    rhoi = info2['rhoi']
    Ri = info2['Ri']
    zeta = info2['zeta']
    tau_es_star = info2['tau_es_star']
    Risco = info2['Risco']
    Mdot = info2['Mdot']
    M = info2['M']
    
    if Rout==None:
        Rout = Risco * 5e+5
    else:
        Rout = Rout
    
    # ---- disk radii from R_star to R_out ----
    R_start = max(R_star, Risco)
    if R_start==Risco:
        print ("@Yao: error -- Rmin = Risco")
        
    R_mid = info1['R']
    ix = R_mid > R_start
    R_mid = R_mid[ix]
    H_mid = info1['H'][ix]
    const_mid = info1['const1'][ix]
    x_mid = R_mid / Risco    
    f_mid = 1 - x_mid**(-0.5)
    
    R_left = np.zeros(R_mid.shape)
    R_left[0] = R_mid[0]-(R_mid[1]-R_mid[0])*0.5
    for i in range(1, len(R_left)):
        R_left[i] = 2*R_mid[i-1]-R_left[i-1]
    R_right = np.hstack([R_left[1:], 2*R_mid[-1]-R_left[-1]])
    
    if subCase == 'sub':
        rho_es_mid = np.zeros(R_mid.shape)
        q1 = R_mid < Rad
        q2 = R_mid >= Rad
        rho_es_mid[q1] = rhoi * (R_mid[q1]/Ri)**(-3)
        rho_es_mid[q2] = rhoi * zeta**(-1) * (R_mid[q2]/Ri)**(-2)
    elif subCase == 'super':
        rho_es_mid = rhoi * zeta**(-1) * (R_mid/Ri)**(-2)
    
    tau_es_mid = kes * R_mid * rho_es_mid
    tau_es_diff = tau_es_star - tau_es_mid
    
    # ---- intrinsic flux ---- 
    Tisco4 = 3 * G * M * Mdot * (8 * PI * STEFAN_BOLTZMANN * Risco**3)**(-1) 
    T_intrinsic = pow(Tisco4 * x_mid**(-3) * f_mid * (1+0.75*tau_es_mid)**(-1), 0.25)
    F_intrinsic = STEFAN_BOLTZMANN * T_intrinsic**4
    # plt.loglog(R_mid, T_intrinsic)
    
    # ---- irradiated flux ----
    F_mid = 4*PI*R_star**2 * STEFAN_BOLTZMANN * T_star**4 *\
            (1-beta) / (0.75*tau_es_diff+1) / (4*PI*R_mid**2) / (0.75*tau_es_mid+1)
    F_mid = find_TOUTSIDE(info1, info2, F_mid, R_mid, H_mid, const_mid, beta,     
                          n_theta = 100, n_phi = 200, quick=quick)
    T_outflow = pow(F_mid/STEFAN_BOLTZMANN, 0.25)
    # plt.loglog(R_mid, T_outflow)
    
    # ---- add them together ----
    T_disk = pow((F_mid+F_intrinsic)/STEFAN_BOLTZMANN, 0.25)
    # plt.loglog(R_mid, T_disk)
    
    area = PI * ( R_right**2.0 - R_left**2.0) * 2.0 # in cm^2
    
    wave = np.logspace(np.log10(w1), np.log10(w2), nfreq)
    specw = np.zeros(len(wave))
    specw1 = np.zeros(len(wave))
    specw2_inside = np.zeros(len(wave))
    specw2_outside = np.zeros(len(wave))
    
    for i in range(len(T_disk)):
        specw += planck_lambda(T_disk[i], wave)*area[i]
        specw1 += planck_lambda(T_intrinsic[i], wave)*area[i]
        if R_mid[i] <= Rsc:
            specw2_inside += planck_lambda(T_outflow[i], wave)*area[i]
        else:
            specw2_outside += planck_lambda(T_outflow[i], wave)*area[i]
    D = 10*PC
    specw *= 1./D**2
    specw1 *= 1./D**2
    specw2_inside *= 1./D**2
    specw2_outside *= 1./D**2
    spec = outSpec_comp(wave, specw, specw1, specw2_inside, specw2_outside,
                        R_mid, F_mid, F_intrinsic)
    return spec