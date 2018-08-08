#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 22 09:24:18 2018

@author: yaoyuhan
"""
import os, sys
tempdir = os.getcwd()+'/'
sys.path.append(tempdir)

import numpy as np
import matplotlib.pyplot as plt
from specSee import standirr_spec, windirr_spec
from wind import wind_structure
from slim_ss import standard_structure
ptex = tempdir.split('BHwind')[0]+'supertex/'
psave = tempdir.split('BHwind')[0]+'saving/'


STEFAN_BOLTZMANN = 5.6696e-5    # yao: Stefan boltzmann constant
MSOL = 1.9884754153381438e+33   # yao: Solar mass in grammes!
G = 6.674079999999999e-08       # yao: Gravitational constant
C = 29979245800.0               # yao: Speed of light
PC = 3.085677581467192e+18      # yao: pc
C8 = 8176.1081773999285
PI = 3.141592653589793          # yao   
kes = 0.342


def outer_geometry(alpha = 0.03, m = 10, mdot = 100): 
    # m = 10
    # alpha = 0.03
    # mdot = 100
    info1 = standard_structure(alpha = alpha, m = m, mdot = mdot)
    info2 = wind_structure(alpha = alpha, m = m, mdot = mdot)
    Rsc = info2['Rsc']
    
    plt.rcParams['figure.figsize'] = [12, 5]
    plt.rcParams['font.size'] = 28
    plt.figure()
    plt.plot(info1['R'], info1['H'], 'k')
    plt.xlabel('$R$ (linear scale)')
    plt.ylabel('$H$ (linear scale)')
    theta = np.linspace(0,2*PI,100)
    x1 = Rsc*np.cos(theta)
    y1 = Rsc*np.sin(theta)
    plt.plot(x1, y1, 'r-')
    plt.xlim(1e+9,1.2e+12)
    plt.ylim(1e+9,0.4e+12)
    plt.xticks([])
    plt.yticks([])
    plt.text(1e+10, Rsc+0.5e+10, r'$R_{\rm sc}$', color='r')
    plt.tight_layout()
    plt.savefig(ptex+'irr.pdf')
    
    
def compare_spectrum():  
    alpha = 0.03
    Rout = None
    beta = 0.7
    nfreq = 1000
    w1 = 1000
    w2 = 12000
    spout0 = windirr_spec(alpha = alpha, m = 10, mdot = 100, Rout = Rout,
                          beta = beta, nfreq = nfreq, w1 = w1, w2 = w2)
    
    stand0 = standirr_spec(alpha=alpha, m = 10, mdot = 100, Rout = Rout,
                           beta = beta, nfreq = nfreq, w1 = w1, w2 = w2)
    
    # same m, lower mdot
    spout1 = windirr_spec(alpha = alpha, m = 10, mdot = 500, Rout = Rout,
                          beta = beta, nfreq = nfreq, w1 = w1, w2 = w2)
    
    stand1 = standirr_spec(alpha = alpha, m = 10, mdot = 500, Rout = Rout,
                           beta = beta, nfreq = nfreq, w1 = w1, w2 = w2)
    
    # lower m, same mdot
    spout2 = windirr_spec(alpha = alpha, m = 1, mdot = 100, Rout = Rout,
                          beta = beta, nfreq = nfreq, w1 = w1, w2 = w2)
    
    stand2 = standirr_spec(alpha = alpha, m = 1, mdot = 100, Rout = Rout,
                           beta = beta, nfreq = nfreq, w1 = w1, w2 = w2)
    
    
    plt.rcParams['figure.figsize'] = [7.5, 9.]
    plt.rcParams['font.size'] = 15
    plt.figure()
    
    ax1 = plt.subplot(321)
    ax1.loglog(spout0.wave, spout0.flux,'k-',label='All')
    ax1.loglog(spout0.wave, spout0.flux_intrinsic,'r:',label='intrinsic')
    ax1.loglog(spout0.wave, spout0.flux_irrin,'b:',label='inside R_sc')
    ax1.loglog(spout0.wave, spout0.flux_irrout,'g:',label='outside R_sc')
    ax1.set_title('Outflow Irradiation',fontsize=16)
    ax1.set_xticklabels([])
    
    ax2 = plt.subplot(322)
    ax2.loglog(stand0['wave'], stand0['flux'],'k-',label='All')
    ax2.loglog(stand0['wave'], stand0['flux_intrinsic'],'r:',label='intrinsic')
    ax2.loglog(stand0['wave'], stand0['flux_irr'],'c:',label='irradiated')
    ax2.set_title('Hard X-ray Irradiation',fontsize=16)
    ax2.set_xticklabels([])
    
    ax3 = plt.subplot(323)
    ax3.loglog(spout1.wave, spout1.flux,'k-',label='All')
    ax3.loglog(spout1.wave, spout1.flux_intrinsic,'r:')
    ax3.loglog(spout1.wave, spout1.flux_irrin,'b:')
    ax3.loglog(spout1.wave, spout1.flux_irrout,'g:')
    ax3.set_xticklabels([])
    
    ax4 = plt.subplot(324)
    ax4.loglog(stand1['wave'], stand1['flux'],'k-')
    ax4.loglog(stand1['wave'], stand1['flux_intrinsic'],'r:',label='intrinsic')
    ax4.loglog(stand1['wave'], stand1['flux_irr'],'c:',label='irradiated')
    ax4.set_xticklabels([])

    ax5 = plt.subplot(325)
    ax5.loglog(spout2.wave, spout2.flux,'k-')
    ax5.loglog(spout2.wave, spout2.flux_intrinsic,'r:',label='intrinsic')
    ax5.loglog(spout2.wave, spout2.flux_irrin,'b:',label='inside '+r'$R_{\rm sc}$')
    ax5.loglog(spout2.wave, spout2.flux_irrout,'g:',label='outside '+r'$R_{\rm sc}$')
    
    ax6 = plt.subplot(326)
    ax6.loglog(stand2['wave'], stand2['flux'],'k-')
    ax6.loglog(stand2['wave'], stand2['flux_intrinsic'],'r:',label='intrinsic')
    ax6.loglog(stand2['wave'], stand2['flux_irr'],'c:',label='irradiated')
    
    ax2.set_yticklabels([])
    ax4.set_yticklabels([])
    ax6.set_yticklabels([])
    
    ax1.tick_params(direction='in', which='both')
    ax1.yaxis.set_ticks_position('both')
    ax2.tick_params(direction='in', which='both')
    ax2.yaxis.set_ticks_position('both')
    ax3.tick_params(direction='in', which='both')
    ax3.yaxis.set_ticks_position('both')
    ax4.tick_params(direction='in', which='both')
    ax4.yaxis.set_ticks_position('both')
    ax5.tick_params(direction='in', which='both')
    ax5.yaxis.set_ticks_position('both')
    ax6.tick_params(direction='in', which='both')
    ax6.yaxis.set_ticks_position('both')
    
    ymax2 = 7e-4
    ymin2 = 4e-8
    ax1.set_ylim(ymin2, ymax2)
    ax2.set_ylim(ymin2, ymax2)
    
    ymax4 = 4e-3
    ymin4 = 1.1e-7
    ax3.set_ylim(ymin4, ymax4)
    ax4.set_ylim(ymin4, ymax4)
    
    ymax6 = 7e-5
    ymin6 = 1.5e-9
    ax5.set_ylim(ymin6, ymax6)
    ax6.set_ylim(ymin6, ymax6)
    
    ax1.hlines(1e-5, 900, 12100, color='gold',linewidth=2,alpha=0.3)
    ax2.hlines(1e-5, 900, 12100, color='gold',linewidth=2,alpha=0.3)
    ax3.hlines(1e-5, 900, 12100, color='gold',linewidth=2,alpha=0.3)
    ax4.hlines(1e-5, 900, 12100, color='gold',linewidth=2,alpha=0.3)
    ax5.hlines(1e-5, 900, 12100, color='gold',linewidth=2,alpha=0.3)
    ax6.hlines(1e-5, 900, 12100, color='gold',linewidth=2,alpha=0.3)
    
    ax5.legend(loc='upper right',frameon=False)
    ax6.legend(loc='upper right',frameon=False)
    
    
    ax5.set_xlabel(r'$\lambda$'+' ('+r'$\AA$'+')')
    ax6.set_xlabel(r'$\lambda$'+' ('+r'$\AA$'+')')
    ax3.set_ylabel(r'$f_{\lambda}$'+' (erg'+r'$\cdot$'+'cm'+r'$^{-2}$'+
                   r'$\cdot$'+'s'+r'$^{-1}$'+r'$\cdot$'+r'$\AA^{-1}$'+')')
    
    
    plt.tight_layout(h_pad=0.2, w_pad = 0.2)  
    ax2.text(1000, ymin2*7, 'm=%d'%stand0['m'])
    ax2.text(1000, ymin2*2, '$\dot m$=%d'%stand0['mdot'])
    ax4.text(1000, ymin4*7, 'm=%d'%stand1['m'])
    ax4.text(1000, ymin4*2, '$\dot m$=%d'%stand1['mdot'])
    ax6.text(1000, ymin6*7, 'm=%d'%stand2['m'])
    ax6.text(1000, ymin6*2, '$\dot m$=%d'%stand2['mdot'])

    plt.savefig(ptex+'compare.pdf')
    
    
def see_integration():
    alpha = 0.03
    Rout = 3e+13
    beta = 0.7
    nfreq = 1000
    w1 = 1000
    w2 = 12000
    
    m1 = 10
    mdot1 = 100
    spout1 = windirr_spec(alpha = alpha, m = m1, mdot = mdot1, Rout = Rout,
                          beta = beta, nfreq = nfreq, w1 = w1, w2 = w2,
                          quick=False)
    
    np.savetxt(psave+'R1', spout1.R_mid)
    np.savetxt(psave+'F1', spout1.F_mid)
    np.savetxt(psave+'F10', spout1.F_intrinsic)
    np.savetxt(psave+'m1', np.array([m1]))
    np.savetxt(psave+'mdot1', np.array([mdot1]))
    '''
    m2 = 30
    mdot2 = 100
    spout2 = windirr_spec(alpha = alpha, m = m2, mdot = mdot2, Rout = Rout,
                          beta = beta, nfreq = nfreq, w1 = w1, w2 = w2,
                          quick=False)
    
    np.savetxt(psave+'R2', spout2.R_mid)
    np.savetxt(psave+'F2', spout2.F_mid)
    np.savetxt(psave+'F20', spout2.F_intrinsic)
    np.savetxt(psave+'m2', np.array([m2]))
    np.savetxt(psave+'mdot2', np.array([mdot2]))
    
    m3 = 10
    mdot3 = 10
    spout3 = windirr_spec(alpha = alpha, m = m3, mdot = mdot3, Rout = Rout,
                          beta = beta, nfreq = nfreq, w1 = w1, w2 = w2,
                          quick=False)
    np.savetxt(psave+'R3', spout3.R_mid)
    np.savetxt(psave+'F3', spout3.F_mid)
    np.savetxt(psave+'F30', spout3.F_intrinsic)
    np.savetxt(psave+'m3', np.array([m3]))
    np.savetxt(psave+'mdot3', np.array([mdot3]))
    
    m5 = 1
    mdot5 = 100
    spout5 = windirr_spec(alpha = alpha, m = m5, mdot = mdot5, Rout = Rout,
                          beta = beta, nfreq = nfreq, w1 = w1, w2 = w2,
                          quick=False)
    np.savetxt(psave+'R5', spout5.R_mid)
    np.savetxt(psave+'F5', spout5.F_mid)
    np.savetxt(psave+'F50', spout5.F_intrinsic)
    np.savetxt(psave+'m5', np.array([m5]))
    np.savetxt(psave+'mdot5', np.array([mdot5]))
    
    m6 = 10
    mdot6 = 500
    spout6 = windirr_spec(alpha = alpha, m = m6, mdot = mdot6, Rout = Rout,
                          beta = beta, nfreq = nfreq, w1 = w1, w2 = w2,
                          quick=False)
    np.savetxt(psave+'R6', spout6.R_mid)
    np.savetxt(psave+'F6', spout6.F_mid)
    np.savetxt(psave+'F60', spout6.F_intrinsic)
    np.savetxt(psave+'m6', np.array([m6]))
    np.savetxt(psave+'mdot6', np.array([mdot6]))
    '''
    
    R1 = np.loadtxt(psave+'R1')
    F1 = np.loadtxt(psave+'F1')
    F10 = np.loadtxt(psave+'F10')
    m1 = np.loadtxt(psave+'m1')
    mdot1 = np.loadtxt(psave+'mdot1')
    peak_wv1 = 2897.7729*1e+4/pow((F1+F10)/STEFAN_BOLTZMANN, 0.25)
    ix = (peak_wv1>1000)&(peak_wv1<10000)
    r1uv = R1[ix][0]
    r1ir = R1[ix][-1]
    
    R2 = np.loadtxt(psave+'R2')
    F2 = np.loadtxt(psave+'F2')
    F20 = np.loadtxt(psave+'F20')
    m2 = np.loadtxt(psave+'m2')
    mdot2 = np.loadtxt(psave+'mdot2')
    peak_wv2 = 2897.7729*1e+4/pow((F2+F20)/STEFAN_BOLTZMANN, 0.25)
    ix = (peak_wv2>1000)&(peak_wv2<10000)
    r2uv = R2[ix][0]
    r2ir = R2[ix][-1]
    
    R3 = np.loadtxt(psave+'R3')
    F3 = np.loadtxt(psave+'F3')
    F30 = np.loadtxt(psave+'F30')
    m3 = np.loadtxt(psave+'m3')
    mdot3 = np.loadtxt(psave+'mdot3')
    peak_wv3 = 2897.7729*1e+4/pow((F3+F30)/STEFAN_BOLTZMANN, 0.25)
    ix = (peak_wv3>1000)&(peak_wv3<10000)
    r3uv = R3[ix][0]
    r3ir = R3[ix][-1]
    
    R5 = np.loadtxt(psave+'R5')
    F5 = np.loadtxt(psave+'F5')
    F50 = np.loadtxt(psave+'F50')
    m5 = np.loadtxt(psave+'m5')
    mdot5 = np.loadtxt(psave+'mdot5')
    peak_wv5 = 2897.7729*1e+4/pow((F5+F50)/STEFAN_BOLTZMANN, 0.25)
    ix = (peak_wv5>1000)&(peak_wv5<10000)
    r5uv = R5[ix][0]
    r5ir = R5[ix][-1]
    
    R6 = np.loadtxt(psave+'R6')
    F6 = np.loadtxt(psave+'F6')
    F60 = np.loadtxt(psave+'F60')
    m6 = np.loadtxt(psave+'m6')
    mdot6 = np.loadtxt(psave+'mdot6')
    peak_wv6 = 2897.7729*1e+4/pow((F6+F60)/STEFAN_BOLTZMANN, 0.25)
    ix = (peak_wv6>1000)&(peak_wv6<10000)
    r6uv = R6[ix][0]
    r6ir = R6[ix][-1]
    
    plt.rcParams['figure.figsize'] = [7.5, 9.]
    plt.rcParams['font.size'] = 16
    plt.figure()
    yy = 2e+8
    xleft = 1.5e+8
    xright = 4e+13
    yup = 5e+19
    ydown = 9e+7
    yyy = 3e+8
    
    ax1 = plt.subplot(321)
    plt.loglog(R1, F1, 'g-', linewidth=1)    
    plt.loglog(R2, F2, 'b-', linewidth=1)    
    plt.loglog(R3, F3, 'c-', linewidth=1)    
    plt.loglog(R5, F5, 'm-', linewidth=1)
    plt.loglog(R6, F6, 'r-', linewidth=1)
    ax1.set_ylim(ydown, yup)
    ax1.set_xlim(xleft, xright)
    ax1.set_xticklabels([])
    ax1.text(0.5e+9, yyy, r'$F_{\rm irr}$'+' (erg cm'+r'$^{-2}$'+' s'+r'$^{-1}$)')
    
    ax2 = plt.subplot(322)
    ax2.loglog(R3, F3, 'c-', label='$m$=%d, $\dot m$=%d'%(m3, mdot3))    
    ax2.loglog(R3, F30, 'c--')
    ax2.set_ylim(ydown, yup)
    ax2.set_xlim(xleft, xright)
    ax2.set_yticklabels([])
    ax2.set_xticklabels([])
    ax2.text(0.5e+9, yyy, '$m$=%d, $\dot m$=%d'%(m3, mdot3))
    ax2.text(1e+11, 1e+15, r'$F_{\rm irr}$', color='c')
    ax2.text(1e+11, 6e+11, r'$F_{0}$', color='c')
    ax2.plot([r3uv, r3ir], [yy, yy], 'k')
    
    ax3 = plt.subplot(323)
    ax3.loglog(R1, F1, 'g-', label='$m$=%d, $\dot m$=%d'%(m1, mdot1))    
    ax3.loglog(R1, F10, 'g--')
    ax3.set_xticklabels([])
    ax3.set_ylim(ydown, yup)
    ax3.set_xlim(xleft, xright)
    ax3.text(0.5e+9, yyy, '$m$=%d, $\dot m$=%d'%(m1, mdot1))
    ax3.text(1e+10, 1e+17, r'$F_{\rm irr}$',color='g')
    ax3.text(1e+10, 3e+14, r'$F_{0}$',color='g')
    ax3.plot([r1uv, r1ir], [yy, yy], 'k')
    
    ax4 = plt.subplot(324)
    ax4.loglog(R2, F2, 'b-', label='$m$=%d, $\dot m$=%d'%(m2, mdot2))    
    ax4.loglog(R2, F20, 'b--')
    ax4.set_xticklabels([])
    ax4.set_yticklabels([])
    ax4.set_ylim(ydown, yup)
    ax4.set_xlim(xleft, xright)
    ax4.text(0.5e+9, yyy, '$m$=%d, $\dot m$=%d'%(m2, mdot2))
    ax4.text(3e+11, 5e+14, r'$F_{\rm irr}$',color='b')
    ax4.text(3e+11, 1e+12, r'$F_{0}$',color='b')
    ax4.plot([r2uv, r2ir], [yy, yy], 'k')
    
    ax5 = plt.subplot(325)
    ax5.loglog(R5, F5, 'm-', label='$m$=%d, $\dot m$=%d'%(m5, mdot5))    
    ax5.loglog(R5, F50, 'm--')
    ax5.set_ylim(ydown, yup)
    ax5.set_xlim(xleft, xright)
    ax5.text(0.5e+9, yyy, '$m$=%d, $\dot m$=%d'%(m5, mdot5))
    ax5.set_xlabel('$R$ (cm)')
    ax5.text(5e+9, 2e+16, r'$F_{\rm irr}$',color='m')
    ax5.text(5e+9, 1e+14, r'$F_{0}$',color='m')
    ax5.plot([r5uv, r5ir], [yy, yy], 'k')
    
    ax6 = plt.subplot(326)
    ax6.loglog(R6, F6, 'r-', label='$m$=%d, $\dot m$=%d'%(m6, mdot6))    
    ax6.loglog(R6, F60, 'r--')
    ax6.set_ylim(ydown, yup)
    ax6.set_xlim(xleft, xright)
    ax6.set_yticklabels([])
    ax6.text(0.5e+9, yyy, '$m$=%d, $\dot m$=%d'%(m6, mdot6))
    ax6.set_xlabel('$R$ (cm)')
    ax6.text(2e+10, 0.5e+16, r'$F_{\rm irr}$',color='r')
    ax6.text(2e+10, 4e+14, r'$F_{0}$',color='r')
    ax6.plot([r6uv, r6ir], [yy, yy], 'k')
    
    ax1.tick_params(direction='in', which='both')
    ax1.yaxis.set_ticks_position('both')
    ax2.tick_params(direction='in', which='both')
    ax2.yaxis.set_ticks_position('both')
    ax3.tick_params(direction='in', which='both')
    ax3.yaxis.set_ticks_position('both')
    ax4.tick_params(direction='in', which='both')
    ax4.yaxis.set_ticks_position('both')
    ax5.tick_params(direction='in', which='both')
    ax5.yaxis.set_ticks_position('both')
    ax6.tick_params(direction='in', which='both')
    ax6.yaxis.set_ticks_position('both')
    
    plt.tight_layout(h_pad=0.2, w_pad=0.2)
    plt.savefig(ptex+'see_integration.pdf')