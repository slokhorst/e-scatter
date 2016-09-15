#!/usr/bin/env python3
import numpy as np
from numpy import (pi)
from scipy.special import gamma
# from scipy.optimize import curve_fit
from scipy.optimize import minimize
import matplotlib.pyplot as plt


def PSD_Palasantzas1D(k, s, xi, h):
    p0 = 2*s*s*xi*(np.sqrt(np.pi)*gamma(h+0.5)/gamma(h))
    P = p0/np.power(1+np.square(k*xi), h+1/2)
    return P


def calcPSD(X):
    kn = X.shape[-1]
    muX = np.repeat(np.mean(X), kn)
    P = kn*np.square(np.absolute(np.fft.fft(X-muX)/kn))
    K = 2*np.pi/(kn-1)*np.arange(0, kn)
    return P, K


def plotPSD(P, show=True):
    kn = P.shape[-1]
    while len(P.shape) > 1:
        print("WARNING: averaging {} PSDs".format(P.shape[0]))
        P = np.mean(P, 0)
    K = 2*np.pi/(kn-1) * np.arange(0, kn)
    plt.loglog(K[1:-1]/(2*pi), P[1:-1], 'xk')
    plt.title('power spectral density')
    plt.xlabel('wave number [px$^{-1}$]')
    plt.ylabel('power [px$^{3}$]')
    plt.grid(True)
    if show:
        plt.show()


def cfitPSD(P, model, plot=True):
    kn = P.shape[-1]
    while len(P.shape) > 1:
        print("WARNING: averaging {} PSDs".format(P.shape[0]))
        P = np.mean(P, 0)
    K = 2*np.pi/(kn-1)*np.arange(0, kn)
    s0 = np.sqrt(np.sum(P)/(kn-1))
    # dy = (kn-1) / kn
    model = model.lower()
    if model == 'palasantzas':
        def fitFunc(k, s, xi, rx):
            return np.log(
                PSD_Palasantzas1D(k, s, xi, rx) +
                PSD_Palasantzas1D(2*pi - k, s, xi, rx) +
                s0**2 - s**2)
    else:
        raise 'Unknown model "{}"'.format(model)

    def err(p):
        return np.mean((fitFunc(K[1:-1], *p)-np.log(P[1:-1]))**2)

    res = minimize(
        err,
        x0=[s0/2, 1, 0.5],
        bounds=[(0.01*s0, 0.99*s0), (0, None), (0, 1)],
        method='L-BFGS-B',
        tol=1e-10)

    if not res.success:
        raise Exception('failed to fit data: {}'.format(res.message))
    popt = res.x
    noiseLevel = np.sqrt(s0**2-popt[0]**2)
    print('sigma0={:.2f}, noise-level={:.2f}'.format(s0, noiseLevel))
    print('fitted sigma={:.2f} xi={:.2f} H={:.2f}'.format(
        popt[0], popt[1], popt[2]))

    def fitResult(k):
        return fitFunc(k, *popt)
    if plot:
        plotPSD(P, show=False)
        plt.loglog(
            K[1:-1]/(2*pi), np.exp(fitResult(K[1:-1])), '-r',
            np.c_[K[1], K[-1]]/(2*pi),
            np.c_[noiseLevel**2, noiseLevel**2], '--r')

        plt.title(
            'PSD fit\n model: {}, parameters: $\sigma={:.2f}$,'
            ' $\\xi={:.2f}$, $H={:.2f}$\n$\sigma_0={:.2f}$,'
            ' $\sigma_N={:.2f}$'.format(
                model, popt[0], popt[1], popt[2], s0, noiseLevel))
        plt.show()
    return popt
