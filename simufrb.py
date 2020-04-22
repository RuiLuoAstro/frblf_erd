import numpy as np
import random as rd
from scipy import integrate
from scipy.interpolate import interp1d
import time
import warnings
import sys
import argparse
 
from frb_util import *

dis = AstroDistribution()
cos = Cosmology()
er = EventRate()
tel = Telescope()

def Simu_FRBs(phis, alpha, logls, logl0, mu_w, sigma_w, dnu, ns, fov, npol, g, tsys, bw, sn0):
    res = np.zeros((ns, 10))
    ns0 = ns
    nt = 0
    lamda = er.rate_2d(sn0, bw, npol, g, tsys, dnu, phis, alpha, logls, logl0, mu_w, sigma_w)*fov
    while (ns>0):
        # Sampling luminosities
        vlogL = np.arange(logl0-1., 48., (48.-logl0+1.)/10000)
        vlik = dis.Schechter_log(vlogL, 1, alpha, logls)
        vlogL = Sampling1D(vlogL, vlik, logl0, 47, ns0)
        vlnEps = np.random.uniform(-np.log(2), 0, ns0)
        vEps = np.exp(vlnEps)
        # Sampling galaxy redshifts
        vZg = np.arange(0, 3.1, 3.1/10000)
        vlik = dis.Distribution_volume(vZg)
        vZ = Sampling1D(vZg, vlik, 0, 3.0, ns0)
        # Sampling pulse width
        vlogW0 = np.arange(-0.5, 1.5, 2./10000)
        vlik = dis.dis_logw(vlogW0, mu, sigma)
        vlogW0 = Sampling1D(vlogW0, vlik, -0.4, 1.4, ns0)
        vW = np.power(10., vlogW0)*(1+vZ)
        # Sampling DM of host galaxies
        vDMH0 = np.arange(0, 5001., 5001./10000)
        vlik = dis.Distribution_HostGalaxyDM(vDMH0, fgalaxy_type=fgt)
        vDMH0 = Sampling1D(vDMH0, vlik, 0, 5000, ns0)
        vDMH = vDMH0*np.sqrt(dis.SFR(vZ))/np.sqrt(dis.SFR(0))
        # Sampling DM of local sources
        vDMS = np.random.uniform(0, 50, ns0)
        # Calculate DM of IGM with redshift
        vDMI = cos.DispersionMeasure_IGM(vZ)
        # Sum up to extragalactic DM
        vDME = (vDMH+vDMS)/(1+vZ)+vDMI
        # Set up the selection criteria
        vft = tel.RMEq(sn0, g, tsys, npol, bw, vW)
        vFlux = vEps*cos.Luminosity_to_Flux(vZ, np.power(10., vlogL) , dnu)
        nlen = len(vFlux[vFlux>vft])
        # Sampling event arrival time
        larray = np.repeat(lamda, nlen)
        vT = rd.expovariate(larray)
        if nlen > ns:
            nlen = ns
        res[nt:(nt+nlen),0] = vFlux[vFlux>vft][0:nlen]
        res[nt:(nt+nlen),1] = vW[vFlux>vft][0:nlen]
        res[nt:(nt+nlen),2] = vT[0:nlen]
        res[nt:(nt+nlen),3] = vDME[vFlux>vft][0:nlen]
        res[nt:(nt+nlen),4] = vft[vFlux>vft][0:nlen]
        res[nt:(nt+nlen),5] = vlogL[vFlux>vft][0:nlen]
        res[nt:(nt+nlen),6] = vZ[vFlux>vft][0:nlen]
        res[nt:(nt+nlen),7] = vDMI[vFlux>vft][0:nlen]
        res[nt:(nt+nlen),8] = vDMH[vFlux>vft][0:nlen]
        res[nt:(nt+nlen),9] = vDMS[vFlux>vft][0:nlen]
        nt = nt+nlen
        ns = ns-nlen
        pct = float(nt)/ns0*100
        print ("%0.1f%% mock FRBs have been simulated." %pct)
    return res


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='FRB sample simulator')
    parser.add_argument('-ns', action='store', dest='Ns', type=int, help='FRB number')
    parser.add_argument('-phis', action='store', dest='phis', type=float, help='Charactreristic event rate density in units of Gpc^{-3} yr^{-1}')
    parser.add_argument('-alpha', action='store', dest='alpha', type=float, help='Power-law index of luminosity function')
    parser.add_argument('-logls', action='store', dest='logls', type=float, help='Upper cut-off luminosity')
    parser.add_argument('-logl0', action='store', dest='logl0', type=float, help='Lower cut-off luminosity')
    parser.add_argument('-dnu', action='store', dest='dnu', type=float, help='Reference intrinsic spectral width in units of MHz')
    parser.add_argument('-mu', action='store', dest='mu', type=float, help='Mean of logarithmic intrinsic width distribution')
    parser.add_argument('-sig', action='store', dest='sigma', type=float, help='Standard deviation of logarithmic intrinsic width distribution')
    parser.add_argument('-fgt', action='store', dest='fgt', type=str, help='Host galaxy type')
    parser.add_argument('-ga', action='store', dest='gain', type=float, help='Telescope gain in units of K/Jy')
    parser.add_argument('-npol', action='store', dest='npol', type=int, help='Polarization channel number')
    parser.add_argument('-bw', action='store', dest='bw', type=float, help='Bandwidth in units of MHz')
    parser.add_argument('-ts', action='store', dest='tsys', type=float, help='Telescope system temperature in units of K')
    parser.add_argument('-sn0', action='store', dest='sn0', type=float, help='Detection threshold of signal to noise ratio')
    parser.add_argument('-fov', action='store', dest='fov', type=float, help='Field of view in units of deg^2')
    parser.add_argument('-out', action='store', dest='output', help='Save the output file as txt.')
    
    args = parser.parse_args()
    Ns = args.Ns
    phis = args.phis
    alpha = args.alpha
    logls = args.logls
    logl0 = args.logl0
    dnu = args.dnu
    mu = args.mu
    sigma = args.sigma
    fgt = args.fgt
    gain = args.gain
    npol = args.npol
    bw = args.bw
    Ts = args.tsys
    sn0 = args.sn0
    fov = args.fov
    
    output = args.output
    
    res = Simu_FRBs(phis, alpha, logls, logl0, mu, sigma, dnu, Ns, fov, npol, gain, Ts, bw, sn0)
    np.savetxt(output, res, delimiter=' ', header="#S W T DMe thres logL Z DMi DMh DMs", comments="")
