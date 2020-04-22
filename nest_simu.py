#!/bin/python

import numpy as np
from scipy import interpolate
import time
import pymultinest
import warnings
import sys
import argparse

from frb_util import *
dis = AstroDistribution()
cos = Cosmology()
er = EventRate()
lf = Loadfiles()

dnu = 1000.

sn0 = 10.
npol = 2.

# Two surveys
g1 = 0.7
bw1 = 300 
tsys1 = 30
fov1 = 0.55

g2 = 0.05
bw2 = 300
tsys2 = 100
fov2 = 30

def lnlik(vpar):
    #global vLOGFLUX1, vDME1, vLOGW1, vN1, vT1
    #global vLOGFLUX2, vDME2, vLOGW2, vN2, vT2
    try:
        norm1 = dis.Norm1D(sn0, bw1, npol, g1, tsys1, dnu, vpar[1], vpar[2], vpar[3], vpar[4], vpar[5])
        norm2 = dis.Norm1D(sn0, bw2, npol, g2, tsys2, dnu, vpar[1], vpar[2], vpar[3], vpar[4], vpar[5])
        loglik_fdm1 = np.sum(dis.log_distr_fdmw(dnu, vLOGFLUX1, vDME1, vLOGW1, vpar[1], vpar[2], vpar[3], vpar[4], vpar[5], gtype=fgt)-np.log(norm1))
        loglik_fdm2 = np.sum(dis.log_distr_fdmw(dnu, vLOGFLUX2, vDME2, vLOGW2, vpar[1], vpar[2], vpar[3], vpar[4], vpar[5], gtype=fgt)-np.log(norm2))
        rho1 = er.rate_2d(sn0, bw1, npol, g1, tsys1, dnu, vpar[0], vpar[1], vpar[2], vpar[3],  vpar[4], vpar[5])
        rho2 = er.rate_2d(sn0, bw2, npol, g2, tsys2, dnu, vpar[0], vpar[1], vpar[2], vpar[3],  vpar[4], vpar[5])
        rho = np.array([rho1, rho2])
        loglik_poi = np.sum(er.log_dis_poi(rho, vN, vFOV, vT))
        res = loglik_fdm1 + loglik_fdm2 + loglik_poi
        return res
    except:
        print 'Numerical error: @', vpar
        return -1e99

def myprior(cube, ndim, nparams):
    for i in range(ndim):
        cube[i] = vpar_range[i,0]+cube[i] *(vpar_range[i,1]-vpar_range[i,0])
    if bolupper:
        cube[0] = np.power(10, cube[0])
        cube[3] = np.log10(cube[3])

def myloglike(cube, ndim, nparams):
    cube2 = np.zeros(ndim)
    for i in range(0,ndim):
        cube2[i] = cube[i]
    res = lnlik(cube2)
    return res

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='LF Measurements Program for Simulated FRBs')
    parser.add_argument('-f1', action='store', dest='simu1', type=str, help='Input Survey 1')
    parser.add_argument('-f2', action='store', dest='simu2', type=str, help='Input Survey 2')
    parser.add_argument('-o', action='store', dest='fout', type=str, help='Save the output file')
    parser.add_argument('-g', action='store', dest='fgt', type=str, help='Host galaxy type')
    parser.add_argument('-upper', dest='bolupper', type=bool, help='Bool option: choosing flat prior')
    args = parser.parse_args()
    simu1 = args.simu1
    simu2 = args.simu2
    fout = args.fout
    fgt = args.fgt
    bolupper = args.bolupper
    
    cat1 = lf.LoadSimuData(simu1)
    vLOGFLUX1 = np.log10(cat1['S'])
    vLOGW1 = np.log10(cat1['W'])
    vDME1 = cat1['DMe']
    vdT1 = cat1['T']
    vN1 = len(vLOGFLUX1)
    vT1 = np.sum(vdT1)

    cat2 = lf.LoadSimuData(simu2)
    vLOGFLUX2 = np.log10(cat2['S'])
    vLOGW2 = np.log10(cat2['W'])
    vDME2 = cat2['DMe']
    vdT2 = cat2['T']
    vN2 = len(vLOGFLUX2)
    vT2 = np.sum(vdT2)

    vN = np.array([vN1, vN2])
    vFOV = np.array([fov1, fov2])
    vT = np.array([vT1, vT2])
    
    if bolupper:
        vpara=np.array([1., -3.0, 43.0, 1e37, -1, 0.01])
        vparb=np.array([5.0, 1.1, 47.0, 1e42, 2, 1.0])
    else:
        vpara=np.array([10., -3.0, 43.0, 37.0, -1, 0.01])
        vparb=np.array([1e5, 1.1, 47.0, 42.0, 2, 1.0])
    
    vpar_range=np.dstack((vpara.transpose(),vparb.transpose()))[0,:,:]
    
    print '------------par range-----------'
    print vpar_range
    a1 = time.clock()
    print myloglike(vpara, len(vpara),len(vpara))
    a2 = time.clock()
    print a1,a2
    print "Running Nest Sampling ..."
    # run MultiNest
    pymultinest.run(myloglike, myprior, len(vpara),
                    importance_nested_sampling = False,
                    resume = False,
                    verbose = True,
                    sampling_efficiency = 'model',
                    n_live_points = 1000,
                    outputfiles_basename='nest_out/simu/'+fout)
