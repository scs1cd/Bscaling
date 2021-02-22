# -*- coding: utf-8 -*-
"""
Created on Sat Jan  9 16:15:55 2021

@author: Chris Davies
"""

import numpy as np
import pandas as pd
import os, copy
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from collections import Counter

def idxStr(l1, s):
    idx = []
    i = 0
    length = len(l1)
    while i < length:
        if s == l1[i]:
            idx.append(i)
        i += 1
    return idx

def convertDictKeys(indict):
    """
    Convert keys of input dictionary into
    numpy arrays
    """
    for key in indict:
        indict[key] = np.asarray(indict[key])
    return indict

def combineDataDict(datadict, quiet=False):
    """
    This function merges all datasets in datadict
    and returns a merged dictionary.
    """
    alldatadict = {"plot":True, "rmsINT":{}, "rmsCMB":{}, "dipCMB":{}, "plotp":{}}

    firstkey = True
    for key in datadict:
        if not quiet:
            print("\nDataset : %s" %(key))
        if datadict[key]["plot"]:
            if firstkey:
                Rmall    = datadict[key]["Rm"] 
                Eall     = datadict[key]["E"]
                Leall    = datadict[key]["rmsINT"]["Le"]
                Lermscmb = datadict[key]["rmsCMB"]["Le"]
                Ledipcmb = datadict[key]["dipCMB"]["Le"]
                PAall    = datadict[key]["p"]
                fohmall  = datadict[key]["fohm"]
                bdipall  = datadict[key]["bdip"]
                firstkey = False
            else:
                Rmall    = np.concatenate((Rmall,datadict[key]["Rm"]))
                Eall     = np.concatenate((Eall,datadict[key]["E"]))
                Leall    = np.concatenate((Leall,datadict[key]["rmsINT"]["Le"]))
                Lermscmb = np.concatenate((Lermscmb,datadict[key]["rmsCMB"]["Le"]))
                Ledipcmb = np.concatenate((Ledipcmb,datadict[key]["dipCMB"]["Le"]))
                PAall    = np.concatenate((PAall,datadict[key]["p"]))
                fohmall  = np.concatenate((fohmall,datadict[key]["fohm"]))
                bdipall  = np.concatenate((bdipall,datadict[key]["bdip"]))
        else:
            if not quiet:
                print("Plot set to False.")

    alldatadict["Rm"] = Rmall
    alldatadict["E"] = Eall
    alldatadict["rmsINT"]["Le"] = Leall
    alldatadict["rmsCMB"]["Le"] = Lermscmb
    alldatadict["dipCMB"]["Le"] = Ledipcmb 
    alldatadict["p"] = PAall
    alldatadict["fohm"] = fohmall
    alldatadict["bdip"] = bdipall 

    return alldatadict

def getEarthEstimates(quiet=True):
    # Brmax = (a/r)^3 g10
    ri    = 1221e3
    rc    = 3480e3
    rsurf = 6371e3
    omega = 7.272e-5
    dd    = rc-ri
    rho   = 1e4
    mu0   = 4.*np.pi*1.e-7

    # Min and max CMB dipole field strengths for Earth: 20,000 nT and 40,000 nT
    Le_earth_dipmin = 20000e-9 * np.sqrt(2.) * (rsurf/rc)**3 / (np.sqrt(rho*4*np.pi*1e-7)*omega*dd)
    Le_earth_dipmax = 40000e-9 * np.sqrt(2.) * (rsurf/rc)**3 / (np.sqrt(rho*4*np.pi*1e-7)*omega*dd)
    # Min and max core field strengths for Earth: 1 mT and 10 mT
    Le_earth_rmsmin = 1e-3  / (np.sqrt(rho*4*np.pi*1e-7)*omega*dd)
    Le_earth_rmsmax = 10e-3 / (np.sqrt(rho*4*np.pi*1e-7)*omega*dd)
    # Min and max CMB field strengths for Earth: 0.5 mT and 1 mT
    Le_earth_rmscmbmin = 0.5e-3 / (np.sqrt(rho*4*np.pi*1e-7)*omega*dd)
    Le_earth_rmscmbmax = 1e-3 / (np.sqrt(rho*4*np.pi*1e-7)*omega*dd)
    if not quiet:
        print('Le_earth CMB dipole min/max  = ', Le_earth_dipmin, Le_earth_dipmax)
        print('Le_earth Surf dipole min/max = ', Le_earth_dipmin/(rsurf/rc)**3, Le_earth_dipmax/(rsurf/rc)**3)
        print('Le_earth RMS total  min/max  = ', Le_earth_rmsmin, Le_earth_rmsmax)
        print('Le_earth RMS CMB  min/max    = ', Le_earth_rmscmbmin, Le_earth_rmscmbmax)
        # Min buoyancy pow for Earth: 0.1 TW
    P_earth_min = 1e11 * 3. / (4*np.pi*(rc**3-ri**3)) / rho / omega**3 / dd**2
    # Max buoyancy pow for Earth: 5 TW
    P_earth_max = 5e12 * 3. / (4*np.pi*(rc**3-ri**3)) / rho / omega**3 / dd**2
    if not quiet:
        print('Power min/max = ', P_earth_min, P_earth_max)

    return Le_earth_dipmin, Le_earth_dipmax, Le_earth_rmsmin, Le_earth_rmsmax, Le_earth_rmscmbmin, Le_earth_rmscmbmax, P_earth_min, P_earth_max
  
def getPlotProperties(datadict, categorise=False):

    if categorise:
        for key in datadict:
            if datadict[key]["plot"]:
                datadict[key]["plotp"]["Col"] = np.array(datadict[key]["Rm"])/np.array(datadict[key]["Pm"])
                if (datadict[key]["dataset"]=="Mixed"):
                    datadict[key]["plotp"]["marker"] = "*"
                    datadict[key]["plotp"]["size"] = 150
                    datadict[key]["plotp"]["cmap"] = "Blues"
                    datadict[key]["plotp"]["edgecolor"] = "blue"
                    datadict[key]["plotp"]["label"] = "Mixed"
                elif (datadict[key]["dataset"]=="FTFF"):
                    datadict[key]["plotp"]["marker"] = "^"
                    datadict[key]["plotp"]["size"] = 150
                    datadict[key]["plotp"]["cmap"] = "Oranges"
                    datadict[key]["plotp"]["edgecolor"] = "orange"
                    datadict[key]["plotp"]["label"] = "FTFF"
                elif (datadict[key]["dataset"]=="FFFF"):
                    datadict[key]["plotp"]["marker"] = "s"
                    datadict[key]["plotp"]["size"] = 150
                    datadict[key]["plotp"]["cmap"] = "Greens"
                    datadict[key]["plotp"]["edgecolor"] = "green"
                    datadict[key]["plotp"]["label"] = "FFFF"
                elif (datadict[key]["dataset"]=="FF0F"):
                    datadict[key]["plotp"]["marker"] = "v"
                    datadict[key]["plotp"]["size"] = 150
                    datadict[key]["plotp"]["cmap"] = "Purples"
                    datadict[key]["plotp"]["edgecolor"] = "purple"
                    datadict[key]["plotp"]["label"] = "FF0F"
                elif (datadict[key]["dataset"]=="FTFT"):
                    datadict[key]["plotp"]["marker"] = "o"
                    datadict[key]["plotp"]["size"] = 150
                    datadict[key]["plotp"]["cmap"] = "Reds"
                    datadict[key]["plotp"]["edgecolor"] = "red"
                    datadict[key]["plotp"]["label"] = "FTFT"
                elif (datadict[key]["dataset"]=="CE"):
                    datadict[key]["plotp"]["marker"] = "D"
                    datadict[key]["plotp"]["size"] = 150
                    datadict[key]["plotp"]["cmap"] = "Greys"
                    datadict[key]["plotp"]["edgecolor"] = "grey"
                    datadict[key]["plotp"]["label"] = "CE"
                else:
                    raise ValueError("Not valid dataset")
    else:
        for key in datadict:
            if datadict[key]["plot"]:
                datadict[key]["plotp"]["Col"] = np.array(datadict[key]["Rm"])/np.array(datadict[key]["Pm"])
                if (datadict[key]["dataset"]=="L"):
                    datadict[key]["plotp"]["marker"] = "*"
                    datadict[key]["plotp"]["size"] = 150
                    datadict[key]["plotp"]["cmap"] = "Blues"
                    datadict[key]["plotp"]["edgecolor"] = "blue"
                    datadict[key]["plotp"]["label"] = "Leeds"
                elif (datadict[key]["dataset"]=="A"):
                    datadict[key]["plotp"]["marker"] = "^"
                    datadict[key]["plotp"]["size"] = 150
                    datadict[key]["plotp"]["cmap"] = "Oranges"
                    datadict[key]["plotp"]["edgecolor"] = "orange"
                    datadict[key]["plotp"]["label"] = "Aubert et al (2009)"
                elif (datadict[key]["dataset"]=="UC"):
                    datadict[key]["plotp"]["marker"] = "v"
                    datadict[key]["plotp"]["size"] = 150
                    datadict[key]["plotp"]["cmap"] = "Purples"
                    datadict[key]["plotp"]["edgecolor"] = "purple"
                    datadict[key]["plotp"]["label"] = "Christensen FF0F"
                elif (datadict[key]["dataset"]=="UCt"):
                    datadict[key]["plotp"]["marker"] = "^"
                    datadict[key]["plotp"]["size"] = 150
                    datadict[key]["plotp"]["cmap"] = "Purples"
                    datadict[key]["plotp"]["edgecolor"] = "purple"
                    datadict[key]["plotp"]["label"] = "Christensen FTFF"
                elif (datadict[key]["dataset"]=="Y"):
                    datadict[key]["plotp"]["marker"] = "o"
                    datadict[key]["plotp"]["size"] = 150
                    datadict[key]["plotp"]["cmap"] = "Reds"
                    datadict[key]["plotp"]["edgecolor"] = "red"
                    datadict[key]["plotp"]["label"] = "Yadav FTFT"
                elif (datadict[key]["dataset"]=="APath"):
                    datadict[key]["plotp"]["marker"] = "D"
                    datadict[key]["plotp"]["size"] = 150
                    datadict[key]["plotp"]["cmap"] = "Greys"
                    datadict[key]["plotp"]["edgecolor"] = "grey"
                    datadict[key]["plotp"]["label"] = "CE"
                elif (datadict[key]["dataset"]=="S"):
                    datadict[key]["plotp"]["marker"] = "s"
                    datadict[key]["plotp"]["size"] = 150
                    datadict[key]["plotp"]["cmap"] = "Greens"
                    datadict[key]["plotp"]["edgecolor"] = "green"
                    datadict[key]["plotp"]["label"] = "Schwaiger et al (2019)"
                else:
                    raise ValueError("Not valid dataset")

    return datadict

def plotSimulations(ax, datadict=None, alldatadict=None, earthdict=None, field="rmsINT",
                    xrange=[0.,1.], yrange=[0.,1.], colorbar=True, cbarrange=[0.,1.]):

    if field not in ("rmsINT","rmsCMB","dipCMB"):
        raise ValueError("Not valid field provided.")

    xmin = xrange[0]; xmax=xrange[1]
    ymin = yrange[0]; ymax=yrange[1]
    cbarmin = cbarrange[0]; cbarmax = cbarrange[1]

    # Now start the plot
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])

    if earthdict is not None:
        ax.add_patch(Rectangle(xy=(earthdict["p"]["min"],earthdict[field]["min"]),
                               width=(earthdict["p"]["max"]-earthdict["p"]["min"]),
                               height=(earthdict[field]["max"]-earthdict[field]["min"]),
                               linewidth=1, color='black', fill=True))
    for key in datadict:
        if datadict[key]["plot"]:
            plt.scatter(datadict[key]["p"], datadict[key][field]["Le"]/datadict[key]["fohm"]**0.5,
                        s=datadict[key]["plotp"]["size"], marker=datadict[key]["plotp"]["marker"],
                        c=np.log10(datadict[key]["plotp"]["Col"]), vmin=cbarmin, vmax=cbarmax,
                        cmap=plt.get_cmap(datadict[key]["plotp"]["cmap"]), edgecolor=datadict[key]["plotp"]["edgecolor"], 
                        label=datadict[key]["plotp"]["label"])
    
    legend_xpos = 0.01 # x-position of legend (>1, outside of main plot)
    legend_ypos = 0.93; legend_dy   = 0.065 # ypos of legend
    
    iplt = 0
    for key in datadict:
        if datadict[key]["plot"]:
            if (key=="UC"):
                slope_str = "$m$^ = "
            elif (key=="UCt"):
                slope_str = "$m$v = "
            else:
                slope_str = "$m$  = "
     
            plt.text(legend_xpos, legend_ypos-iplt*legend_dy,
                     slope_str+str(np.round(datadict[key][field]["m"],2)) +\
                     "$\pm$"+str(np.round(datadict[key][field]["res"],2)) +\
                     ", SSR="+str(np.round(datadict[key][field]["ssr"],2)),
                     transform=ax.transAxes, color=datadict[key]["plotp"]["edgecolor"])
            iplt += 1

    plt.text(legend_xpos, legend_ypos-iplt*legend_dy, "$m$  = "+str(np.round(alldatadict[field]["m"],2))+\
             "$\pm$"+str(np.round(alldatadict[field]["res"],2))+\
             ", SSR="+str(np.round(alldatadict[field]["ssr"],2)), transform=ax.transAxes, color='black')
    if colorbar:
        # Plot color bar
        cbar = plt.colorbar()
        cbar.set_label("log $Re$")
    
    ax.set_yscale('log')
    ax.set_xscale('log')
    
    return ax, legend_xpos, legend_ypos

def getPlotTitle(myfdip=None, fdip_range= [None]*2, myEr=None, Er_range=[None]*2):
    title_str = ""
    comma = False
    if (myfdip == 1): 
        title_str = "$f_{dip}>0.5$"
        comma = True
    elif (myfdip == 2):
        title_str = "$0.35< f_{dip}<0.80$"
    elif (myfdip == 3):
        title_str = "$0.40< f_{dip}<0.80$"
        comma = True        

    if (myEr == 1):
        if comma:
            title_str += ", "+ "EM/EK $\geq %.1f$" %(Er_range[0])
        else:
            title_str += "EM/EK $\geq %.1f$" %(Er_range[0])
            comma = True

    return title_str

def fits(P, Le, fohm):
    if (len(P) != len(Le) or len(P) != len(fohm) or len(Le) != len(fohm)):
        raise ValueError("Le, P, ohm arrays have different dimensions.")
    else:
        nsims = len(P)
    if nsims>2:
        fp,ssr,rnk,sv,rc = np.polyfit(np.log10(P), np.log10(Le/fohm**0.5), 1, full=True)
        [m, c], cov      = np.polyfit(np.log10(P), np.log10(Le/fohm**0.5), 1, cov=True, full=False)
        res              = np.sqrt(cov[0,0]) # error on slope
    elif nsims==2:
        dy = np.log10(Le[1]/fohm[1]**0.5)-np.log10(Le[0]/fohm[0]**0.5)
        dx = np.log10(P[1])-np.log10(P[0])
        m = dy/dx
        c = 0.5*(np.log10(Le[1]/fohm[1]**0.5)+np.log10(Le[0]/fohm[0]**0.5))-0.5*m*(np.log10(P[1])+np.log10(P[0]))
        ssr = [0.]
        res = 0.
    else:
        raise ValueError("Only 1 data point. No fit is possible.")

    return ssr[0], m, c, res

def getFitBounds(x, c, c_sd, m):
    """
    Get fit error bounds on the slope m.
    c    : prefactor
    c_sd : prefactor error
    m    : slope
    """
    fit_1sdp = 10**(c+c_sd) * x**m
    fit_1sdm = 10**(c-c_sd) * x**m
    fit_2sdp = 10**(c+2.*c_sd) * x**m
    fit_2sdm = 10**(c-2.*c_sd) * x**m
    fit_3sdp = 10**(c+3.*c_sd) * x**m
    fit_3sdm = 10**(c-3.*c_sd) * x**m

    return fit_1sdp, fit_1sdm, fit_2sdp, fit_2sdm, fit_3sdp, fit_3sdm

def shellVolume(ar):
    """
    Calculates spherical shell volume (in units of the
    shell gap D) from aspect ratio ar.
    """
    vol = 4./3.*np.pi*((1.-ar)**(-3)-(1./ar-1.)**(-3))
    return vol

def prefacError(x, y, model=[None,None], plot=False, quiet=False):
    """
    Calculates errors on the power law prefactor.
    This follows Aubert et al. GJI 179, 1414, 2009 (sec. 2.3).
    Specifically, the standard deviation sd_err returned by this function
    is \sigma in Aubert's paper.
    x     : Buoyancy power P of sims data points
    y     : Le/fohm^0.5 of sims data points
    model : list of [10^(prefac), slope] from the fit of the sims
    """
    prefac = model[0]
    slope = model[1]
    err = np.log10(y) - (np.log10(prefac) + slope*np.log10(x))
    mean_err = np.mean(err)
    sd_err = np.std(err)
    
    if not quiet:
        print('\nError estimate on prefactor')
        print('Error (mean, std dev): %s %s' %(mean_err,sd_err))
    if plot:
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)
        xpltmin = np.min(err); xpltmax = np.max(err)
        counts, bins, patches = ax.hist(err, bins=15,
                                        normed=True, align='mid', color='grey', edgecolor='dimgrey', rwidth=0.85)
        ypltmin = 0.; ypltmax = np.max(counts)
        ax.plot([mean_err]*2, [ypltmin,ypltmax], c='r', ls='-', lw=2.)
        ax.plot([mean_err+sd_err]*2, [ypltmin,ypltmax], c='r', ls='--', lw=2.)
        ax.plot([mean_err-sd_err]*2, [ypltmin,ypltmax], c='r', ls='--', lw=2.)
        ax.set_xlim([xpltmin,xpltmax])
        ax.set_ylim([ypltmin,ypltmax])
        ax.set_xlabel('Vertical error e')
        ax.set_ylabel('Probability density')
        plt.show(block=False)
        plt.tight_layout()
        plt.savefig('./eardme/fig/Err_distribution.pdf',format='pdf')
        del ax
    return [mean_err, sd_err]

def saveFitValues(filename=None, datadict=None, alldatadict=None):
    """
    Store fit values in output file
    """
    fout = open(filename, "w")
    delim = "     "
    fout.write("--- B rms ---\n")
    fout.write("Dataset            "+delim+"slope m    "+"std. err. "+" SSR\n")
    for key in datadict:
        if datadict[key]["plot"]:
            fout.write(datadict[key]["plotp"]["label"]+delim+\
                str(np.round(datadict[key]["rmsINT"]["m"],4))+delim+\
                str(np.round(datadict[key]["rmsINT"]["res"],4))+delim+\
                str(np.round(datadict[key]["rmsINT"]["ssr"],4))+"\n")
    fout.write("All                "+delim+str(np.round(alldatadict["rmsINT"]["m"],4))+delim+\
               str(np.round(alldatadict["rmsINT"]["res"],4))+delim+\
               str(np.round(alldatadict["rmsINT"]["ssr"],4))+"\n")
    fout.write("\n")
    fout.write("--- B rms cmb ---\n")
    fout.write("Dataset            "+delim+"slope m    "+"std. err. "+" SSR\n")
    for key in datadict:
        if datadict[key]["plot"]:
            fout.write(datadict[key]["plotp"]["label"]+delim+\
                str(np.round(datadict[key]["rmsCMB"]["m"],4))+delim+\
                str(np.round(datadict[key]["rmsCMB"]["res"],4))+delim+\
                str(np.round(datadict[key]["rmsCMB"]["ssr"],4))+"\n")
    fout.write("All                "+delim+str(np.round(alldatadict["rmsCMB"]["m"],4))+delim+\
               str(np.round(alldatadict["rmsCMB"]["res"],4))+delim+\
               str(np.round(alldatadict["rmsCMB"]["ssr"],4))+"\n")
    fout.write("\n")
    fout.write("--- B dip cmb ---\n")
    fout.write("Dataset            "+delim+"slope m    "+"std. err. "+" SSR\n")
    for key in datadict:
        if datadict[key]["plot"]:
            fout.write(datadict[key]["plotp"]["label"]+delim+\
                str(np.round(datadict[key]["dipCMB"]["m"],4))+delim+\
                str(np.round(datadict[key]["dipCMB"]["res"],4))+delim+\
                str(np.round(datadict[key]["dipCMB"]["ssr"],4))+"\n")
    fout.write("All                "+delim+str(np.round(alldatadict["dipCMB"]["m"],4))+delim+\
               str(np.round(alldatadict["dipCMB"]["res"],4))+delim+\
               str(np.round(alldatadict["dipCMB"]["ssr"],4))+"\n")
    fout.close()

def savePrefacValues(filename=None, indict=None, l_prefac_err=False):
    # Store pre-factors in output file
    fout = open(filename, "w")
    delim = "    "
    fout.write("--- B rms ---\n")
    fout.write("Scaling     "+delim+"Prefactor c1\n")
    fout.write("IMA    1/3  "+delim+str(round(indict["rmsINT"]["cIMA"],10))+"\n")
    fout.write("MAC    1/4  "+delim+str(round(indict["rmsINT"]["cMAC"],10))+"\n")
    fout.write("IMAC   2/5  "+delim+str(round(indict["rmsINT"]["cIMAC"],10))+"\n")
    fout.write("IMACd  1/5  "+delim+str(round(indict["rmsINT"]["cIMACd"],10))+"\n")
    fout.write("IMACi  3/10 "+delim+str(round(indict["rmsINT"]["cIMACi"],10))+"\n")
    fout.write("\n")
    fout.write("--- B dip cmb ---\n")
    if not l_prefac_err:
        fout.write("Scaling     "+delim+"Prefactor c1\n")
        fout.write("IMA    1/3  "+delim+str(round(indict["dipCMB"]["cIMA"],10))+"\n")
        fout.write("MAC    1/4  "+delim+str(round(indict["dipCMB"]["cMAC"],10))+"\n")
        fout.write("IMAC   2/5  "+delim+str(round(indict["dipCMB"]["cIMAC"],10))+"\n")
        fout.write("IMACd  1/5  "+delim+str(round(indict["dipCMB"]["cIMACd"],10))+"\n")
        fout.write("IMACi  3/10 "+delim+str(round(indict["dipCMB"]["cIMACi"],10))+"\n")
    else:
        fout.write("Scaling     "+delim+"Prefactor c1"+delim+" Sigma\n")
        fout.write("IMA    1/3  "+delim+str(round(indict["dipCMB"]["cIMA"],10))+delim+str(round(indict["dipCMB"]["cIMA_sd"],10))+"\n")
        fout.write("MAC    1/4  "+delim+str(round(indict["dipCMB"]["cMAC"],10))+delim+str(round(indict["dipCMB"]["cMAC_sd"],10))+"\n")
        fout.write("IMAC   2/5  "+delim+str(round(indict["dipCMB"]["cIMAC"],10))+delim+str(round(indict["dipCMB"]["cIMAC_sd"],10))+"\n")
        fout.write("IMACd  1/5  "+delim+str(round(indict["dipCMB"]["cIMACd"],10))+delim+str(round(indict["dipCMB"]["cIMACd_sd"],10))+"\n")
        fout.write("IMACi  3/10 "+delim+str(round(indict["dipCMB"]["cIMACi"],10))+delim+str(round(indict["dipCMB"]["cIMACi_sd"],10))+"\n")
    fout.write("\n")
    fout.write("--- Check for VDM script on B dip cmb ---\n")
    # Check for VDM script: Get Le for CMB dipole field at user-defined value of P
    P0 = 1.e-12
    fout.write("P0          "+delim+str(P0)+"\n")
    fout.write("IMA    1/3  "+delim+str((10**indict["dipCMB"]["cIMA"])*P0**(1/3.))+"\n")
    fout.write("MAC    1/4  "+delim+str((10**indict["dipCMB"]["cMAC"])*P0**0.25)+"\n")
    fout.write("IMAC   2/5  "+delim+str((10**indict["dipCMB"]["cIMAC"])*P0**0.40)+"\n")
    fout.write("IMACd  1/5  "+delim+str((10**indict["dipCMB"]["cIMACd"])*P0**0.20)+"\n")
    fout.write("IMACi  3/10 "+delim+str((10**indict["dipCMB"]["cIMACi"])*P0**0.30)+"\n")
    fout.close()

def fitForceScalings(indict, quiet=False):
    """
    Calculates prefactor of all force scalings
    and updates the dictionary
    """
    # - rms core field
    indict["rmsINT"]["mIMA"] = 1/3.
    indict["rmsINT"]["cIMA"] = np.mean(np.log10(indict["rmsINT"]["Le"]/indict["fohm"]**0.5) - \
                                       indict["rmsINT"]["mIMA"]*np.log10(indict["p"]))

    c1_err_IMA = prefacError(indict["p"], indict["rmsINT"]["Le"]/indict["fohm"]**0.5,
                            model=[10**indict["rmsINT"]["cIMA"],indict["rmsINT"]["mIMA"]], plot=False, quiet=False)
    indict["rmsINT"]["cIMA_sd"] = c1_err_IMA[1]
    if not quiet:
        print('Core RMS, IMA   m=1/3 =0.33', indict["rmsINT"]["cIMA"])

    indict["rmsINT"]["mMAC"] = 0.25
    indict["rmsINT"]["cMAC"] = np.mean(np.log10(indict["rmsINT"]["Le"]/indict["fohm"]**0.5) - \
                                       indict["rmsINT"]["mMAC"]*np.log10(indict["p"]))
    c1_err_MAC = prefacError(indict["p"], indict["rmsINT"]["Le"]/indict["fohm"]**0.5,
                            model=[10**indict["rmsINT"]["cMAC"],indict["rmsINT"]["mMAC"]], plot=False, quiet=False)
    indict["rmsINT"]["cMAC_sd"] = c1_err_MAC[1]
    if not quiet:
        print('Core RMS, MAC   m=1/4 =0.25', indict["rmsINT"]["cMAC"])

    indict["rmsINT"]["mIMAC"] = 0.40
    indict["rmsINT"]["cIMAC"] = np.mean(np.log10(indict["rmsINT"]["Le"]/indict["fohm"]**0.5) - \
                                       indict["rmsINT"]["mIMAC"]*np.log10(indict["p"]))
    c1_err_IMAC = prefacError(indict["p"], indict["rmsINT"]["Le"]/indict["fohm"]**0.5,
                              model=[10**indict["rmsINT"]["cIMAC"],indict["rmsINT"]["mIMAC"]], plot=False, quiet=False)
    indict["rmsINT"]["cIMAC_sd"] = c1_err_IMAC[1]
    if not quiet:
        print('Core RMS, IMAC  m=2/5 =0.40', indict["rmsINT"]["cIMAC"])

    indict["rmsINT"]["mIMACd"] = 0.20
    indict["rmsINT"]["cIMACd"] = np.mean(np.log10(indict["rmsINT"]["Le"]/indict["fohm"]**0.5) - \
                                       indict["rmsINT"]["mIMACd"]*np.log10(indict["p"]))
    c1_err_IMACd = prefacError(indict["p"], indict["rmsINT"]["Le"]/indict["fohm"]**0.5,
                              model=[10**indict["rmsINT"]["cIMACd"],indict["rmsINT"]["mIMACd"]], plot=False, quiet=False)
    indict["rmsINT"]["cIMACd_sd"] = c1_err_IMACd[1]
    if not quiet:
        print('Core RMS, IMACd m=1/5 =0.20', indict["rmsINT"]["cIMACd"])

    indict["rmsINT"]["mIMACi"] = 0.30
    indict["rmsINT"]["cIMACi"] = np.mean(np.log10(indict["rmsINT"]["Le"]/indict["fohm"]**0.5) - \
                                       indict["rmsINT"]["mIMACi"]*np.log10(indict["p"]))
    c1_err_IMACi = prefacError(indict["p"], indict["rmsINT"]["Le"]/indict["fohm"]**0.5,
                              model=[10**indict["rmsINT"]["cIMACi"],indict["rmsINT"]["mIMACi"]], plot=False, quiet=False)
    indict["rmsINT"]["cIMACi_sd"] = c1_err_IMACi[1]
    if not quiet:
        print('Core RMS, IMACi m=3/10=0.30', indict["rmsINT"]["cIMACi"])

    # - rms CMB field
    indict["rmsCMB"]["mIMA"] = 1/3.
    indict["rmsCMB"]["cIMA"] = np.mean(np.log10(indict["rmsCMB"]["Le"]/indict["fohm"]**0.5) - \
                                       indict["rmsCMB"]["mIMA"]*np.log10(indict["p"]))
    c1_err_IMA = prefacError(indict["p"], indict["rmsCMB"]["Le"]/indict["fohm"]**0.5,
                            model=[10**indict["rmsCMB"]["cIMA"],indict["rmsCMB"]["mIMA"]], plot=False, quiet=False)
    indict["rmsCMB"]["cIMA_sd"] = c1_err_IMA[1]
    if not quiet:
        print('CMB RMS, IMA   m=1/3 =0.33', indict["rmsCMB"]["cIMA"])

    indict["rmsCMB"]["mMAC"] = 0.25
    indict["rmsCMB"]["cMAC"] = np.mean(np.log10(indict["rmsCMB"]["Le"]/indict["fohm"]**0.5) - \
                                       indict["rmsCMB"]["mMAC"]*np.log10(indict["p"]))
    c1_err_MAC = prefacError(indict["p"], indict["rmsCMB"]["Le"]/indict["fohm"]**0.5,
                            model=[10**indict["rmsCMB"]["cMAC"],indict["rmsCMB"]["mMAC"]], plot=False, quiet=False)
    indict["rmsCMB"]["cMAC_sd"] = c1_err_MAC[1]
    if not quiet:
        print('CMB RMS, MAC   m=1/4 =0.25', indict["rmsCMB"]["cMAC"])

    indict["rmsCMB"]["mIMAC"] = 0.40
    indict["rmsCMB"]["cIMAC"] = np.mean(np.log10(indict["rmsCMB"]["Le"]/indict["fohm"]**0.5) - \
                                       indict["rmsCMB"]["mIMAC"]*np.log10(indict["p"]))
    c1_err_IMAC = prefacError(indict["p"], indict["rmsCMB"]["Le"]/indict["fohm"]**0.5,
                              model=[10**indict["rmsCMB"]["cIMAC"],indict["rmsCMB"]["mIMAC"]], plot=False, quiet=False)
    indict["rmsCMB"]["cIMAC_sd"] = c1_err_IMAC[1]
    if not quiet:
        print('CMB RMS, IMAC  m=2/5 =0.40', indict["rmsCMB"]["cIMAC"])

    indict["rmsCMB"]["mIMACd"] = 0.20
    indict["rmsCMB"]["cIMACd"] = np.mean(np.log10(indict["rmsCMB"]["Le"]/indict["fohm"]**0.5) - \
                                       indict["rmsCMB"]["mIMACd"]*np.log10(indict["p"]))
    c1_err_IMACd = prefacError(indict["p"], indict["rmsCMB"]["Le"]/indict["fohm"]**0.5,
                              model=[10**indict["rmsCMB"]["cIMACd"],indict["rmsCMB"]["mIMACd"]], plot=False, quiet=False)
    indict["rmsCMB"]["cIMACd_sd"] = c1_err_IMACd[1]
    if not quiet:
        print('CMB RMS, IMACd m=1/5 =0.20', indict["rmsCMB"]["cIMACd"])

    indict["rmsCMB"]["mIMACi"] = 0.30
    indict["rmsCMB"]["cIMACi"] = np.mean(np.log10(indict["rmsCMB"]["Le"]/indict["fohm"]**0.5) - \
                                       indict["rmsCMB"]["mIMACi"]*np.log10(indict["p"]))
    c1_err_IMACi = prefacError(indict["p"], indict["rmsCMB"]["Le"]/indict["fohm"]**0.5,
                              model=[10**indict["rmsCMB"]["cIMACi"],indict["rmsCMB"]["mIMACi"]], plot=False, quiet=False)
    indict["rmsCMB"]["cIMACi_sd"] = c1_err_IMACi[1]
    if not quiet:
        print('CMB RMS, IMACi m=3/10=0.30', indict["rmsCMB"]["cIMACi"])

    # - dipole CMB field
    indict["dipCMB"]["mIMA"] = 1/3.
    indict["dipCMB"]["cIMA"] = np.mean(np.log10(indict["dipCMB"]["Le"]/indict["fohm"]**0.5) - \
                                       indict["dipCMB"]["mIMA"]*np.log10(indict["p"]))

    c1_err_IMA = prefacError(indict["p"], indict["dipCMB"]["Le"]/indict["fohm"]**0.5,
                            model=[10**indict["dipCMB"]["cIMA"],indict["dipCMB"]["mIMA"]], plot=False, quiet=False)
    indict["dipCMB"]["cIMA_sd"] = c1_err_IMA[1]
    if not quiet:
        print('CMB dip, IMA   m=1/3 =0.33', indict["dipCMB"]["cIMA"])

    indict["dipCMB"]["mMAC"] = 0.25
    indict["dipCMB"]["cMAC"] = np.mean(np.log10(indict["dipCMB"]["Le"]/indict["fohm"]**0.5) - \
                                       indict["dipCMB"]["mMAC"]*np.log10(indict["p"]))
    c1_err_MAC = prefacError(indict["p"], indict["dipCMB"]["Le"]/indict["fohm"]**0.5,
                            model=[10**indict["dipCMB"]["cMAC"],indict["dipCMB"]["mMAC"]], plot=False, quiet=False)
    indict["dipCMB"]["cMAC_sd"] = c1_err_MAC[1]
    if not quiet:
        print('CMB dip, MAC   m=1/4 =0.25', indict["dipCMB"]["cMAC"])

    indict["dipCMB"]["mIMAC"] = 0.40
    indict["dipCMB"]["cIMAC"] = np.mean(np.log10(indict["dipCMB"]["Le"]/indict["fohm"]**0.5) - \
                                       indict["dipCMB"]["mIMAC"]*np.log10(indict["p"]))
    c1_err_IMAC = prefacError(indict["p"], indict["dipCMB"]["Le"]/indict["fohm"]**0.5,
                              model=[10**indict["dipCMB"]["cIMAC"],indict["dipCMB"]["mIMAC"]], plot=False, quiet=False)
    indict["dipCMB"]["cIMAC_sd"] = c1_err_IMAC[1]
    if not quiet:
        print('CMB dip, IMAC  m=2/5 =0.40', indict["dipCMB"]["cIMAC"])

    indict["dipCMB"]["mIMACd"] = 0.20
    indict["dipCMB"]["cIMACd"] = np.mean(np.log10(indict["dipCMB"]["Le"]/indict["fohm"]**0.5) - \
                                       indict["dipCMB"]["mIMACd"]*np.log10(indict["p"]))
    c1_err_IMACd = prefacError(indict["p"], indict["dipCMB"]["Le"]/indict["fohm"]**0.5,
                              model=[10**indict["dipCMB"]["cIMACd"],indict["dipCMB"]["mIMACd"]], plot=False, quiet=False)
    indict["dipCMB"]["cIMACd_sd"] = c1_err_IMACd[1]
    if not quiet:
        print('CMB dip, IMACd m=1/5 =0.20', indict["dipCMB"]["cIMACd"])

    indict["dipCMB"]["mIMACi"] = 0.30
    indict["dipCMB"]["cIMACi"] = np.mean(np.log10(indict["dipCMB"]["Le"]/indict["fohm"]**0.5) - \
                                       indict["dipCMB"]["mIMACi"]*np.log10(indict["p"]))
    c1_err_IMACi = prefacError(indict["p"], indict["dipCMB"]["Le"]/indict["fohm"]**0.5,
                              model=[10**indict["dipCMB"]["cIMACi"],indict["dipCMB"]["mIMACi"]], plot=False, quiet=False)
    indict["dipCMB"]["cIMACi_sd"] = c1_err_IMACi[1]
    if not quiet:
        print('CMB dip, IMACi m=3/10=0.30', indict["dipCMB"]["cIMACi"])

    return indict

def filter_table(infname=None, outfname=None, dataset="Leeds", fdip_range=None, EkOPm_range=None, EMoEK_range=None, datadict=None,
                 categorise=False, chk=0, myfohm=0):
    """
    This function filters the original tables of simulations by fdip, E/Pm, and EM/EK
    and saves out the filtered table (preserving the original format).
    infname     : Input file name
    outfname    : Output file name of filtered data
    dataset     : String ("Leeds","Aubert","Yadav","Christensen", or "ChristensenT"). Make sure it is consistent with infname.
    fdip_range  : list of min and max values of fdip
    EkOPm_range : set to range filter by E/Pm
    EMoEK_range : set to range to filter by EM/EK
    """

    if fdip_range is not None:
        filter_fdip = True
        print("\nReading file: %s" %(infname))
        fdip_min = fdip_range[0]
        fdip_max = fdip_range[1]
        print("\nFilter fdip (min,max) : %s %s" %(fdip_min,fdip_max))
    else:
        filter_fdip = False
    if EkOPm_range is not None:
        EkOPm_min = EkOPm_range[0]; EkOPm_max = EkOPm_range[1]
        print('\nFiltering by Ek/Pm (min,max) = %s %s' %(EkOPm_min,EkOPm_max))
        if (EkOPm_min > EkOPm_max):
            raise ValueError('Not valid range of Ek/Pm: min > max.')
        filter_EkOPm = True
    else:
        filter_EkOPm = False
    if EMoEK_range is not None:
        filter_Er = True
        Er_min = EMoEK_range[0]; Er_max = EMoEK_range[1]
        print('\nFiltering by EM/EK (min,max) = %s %s' %(Er_min,Er_max))
        if (Er_min > Er_max):
            raise ValueError('Not valid range of EM/EK: min > max.')
    else:
        filter_Er = False
    # -------------------
    # - Leeds simulations
    # -------------------
    if dataset == "Leeds":
        df = pd.read_csv(infname)
        if filter_fdip:
            # Filter fdip
            mask_fdip = (df['cmb_diptyAve'] > fdip_min) & (df['cmb_diptyAve'] < fdip_max)
            new_df = df[mask_fdip]
        else:
            new_df = df
        if filter_Er:
            mask_Er = (new_df['EmagAve'].values/new_df['EkinAve'].values >= Er_min) &\
                      (new_df['EmagAve'].values/new_df['EkinAve'].values <= Er_max)
            new_df = new_df[mask_Er]
        if filter_EkOPm:
            mask_EkOPm = (new_df['Ek'].values*2./new_df['Pm'].values > EkOPm_min) & (new_df['Ek'].values*2./new_df['Pm'].values < EkOPm_max)
            new_df = new_df[mask_EkOPm]
        nsims = len(new_df['Ek'].values)
        datadict["L"]["nsims"] = nsims
        if (nsims<2):
            datadict["L"]["plot"] = False
            print("\nNot plotting L dataset!\n")
        datadict["L"]["d"] = new_df.to_dict("list")
        datadict["L"]["d"] = convertDictKeys(datadict["L"]["d"])
        # Get field strengths and fits
        datadict["L"]["E"]  = datadict["L"]["d"]["Ek"]
        datadict["L"]["Pm"] = datadict["L"]["d"]["Pm"]
        datadict["L"]["Rm"] = datadict["L"]["d"]["RmAve"]
        
        # Calculate shell volume
        ar = datadict["L"]["d"]["ar"] # aspect ratio
        volS = np.asarray([shellVolume(ar[i]) for i in range(len(ar))])

        if myfohm == 1:
            datadict["L"]["fohm"] = np.ones(len(datadict["L"]["d"]["EkinDisAve"]))
        else:
            datadict["L"]["fohm"] = datadict["L"]["d"]["EmagDisAve"] / \
                                    (datadict["L"]["d"]["EkinDisAve"]+datadict["L"]["d"]["EmagDisAve"])
        # get field strengths and power
        datadict["L"]["rmsINT"]["Le"] = np.sqrt(8.0) * (datadict["L"]["E"]/datadict["L"]["Pm"]) * \
                                        np.sqrt(datadict["L"]["d"]["EmagAve"]/volS)       # Manual eq. (3.38)
        datadict["L"]["rmsCMB"]["Le"] = 2. * np.sqrt(datadict["L"]["E"]/datadict["L"]["Pm"] * datadict["L"]["d"]["ElsCMBAve"])
        datadict["L"]["dipCMB"]["Le"] = 2. * np.sqrt(datadict["L"]["E"]/datadict["L"]["Pm"] * datadict["L"]["d"]["ElsDipCMBAve"])
        datadict["L"]["p"] = 8.0 * (datadict["L"]["E"]/datadict["L"]["Pm"])**3 * \
                             (datadict["L"]["d"]["EkinDisAve"]+datadict["L"]["d"]["EmagDisAve"]) / volS
        # get fdip and bdip
        datadict["L"]["bdip"] = datadict["L"]["rmsINT"]["Le"]/datadict["L"]["dipCMB"]["Le"]
        datadict["L"]["fdip"] = datadict["L"]["d"]["cmb_diptyAve"]

        if categorise:
            # store temp BCs and other useful properties for the categorisation of models
            datadict["L"]["TBC"] = datadict["L"]["d"]["TBC"]
            datadict["L"]["Ir"] = datadict["L"]["d"]["Ir"]
            datadict["L"]["ar"] = datadict["L"]["d"]["ar"]
        
        if not categorise:
            # make fit
            datadict["L"]["rmsINT"]["ssr"], datadict["L"]["rmsINT"]["m"], datadict["L"]["rmsINT"]["c"], datadict["L"]["rmsINT"]["res"] = \
                fits(datadict["L"]["p"], datadict["L"]["rmsINT"]["Le"], datadict["L"]["fohm"])
            datadict["L"]["rmsCMB"]["ssr"], datadict["L"]["rmsCMB"]["m"], datadict["L"]["rmsCMB"]["c"], datadict["L"]["rmsCMB"]["res"] = \
                fits(datadict["L"]["p"], datadict["L"]["rmsCMB"]["Le"], datadict["L"]["fohm"])
            datadict["L"]["dipCMB"]["ssr"], datadict["L"]["dipCMB"]["m"], datadict["L"]["dipCMB"]["c"], datadict["L"]["dipCMB"]["res"] = \
                fits(datadict["L"]["p"], datadict["L"]["dipCMB"]["Le"], datadict["L"]["fohm"])

    # ---------------------------------
    # - Aubert et al (2009) simulations
    # ---------------------------------
    # NB - Aubert uses different radii! chi dependence of parameters?
    elif dataset == "Aubert":
        cols_names = ['E','RaQ','Pr','Pm','chi','fi','Ro','Lo','bdip','fdip','lbar','tauDiss','p','fohm']
        E,RaQ,Pr,Pm,chi,fi,Ro,Lo,bdip,fdip,lbar,tauDiss,p,fohm = np.loadtxt(infname,\
            usecols=tuple(np.arange(len(cols_names))), skiprows=1, unpack='true')
        # Create data dictionary
        data_dict = {'E':E,'RaQ':RaQ,'Pr':Pr,'Pm':Pm,'chi':chi,'fi':fi,'Ro':Ro,'Lo':Lo,'bdip':bdip,\
                     'fdip':fdip,'lbar':lbar,'tauDiss':tauDiss,'p':p,'fohm':fohm}
        # Create dataframe
        df = pd.DataFrame(data_dict, columns=cols_names)
        
        if filter_fdip:
            mask_fdip = (df['fdip'] > fdip_min) & (df['fdip'] < fdip_max)
            new_df = df[mask_fdip]
        else:
            new_df = df
        if filter_Er:
            mask_Er = (new_df['Lo'].values**2/new_df['Ro'].values**2 >= Er_min) & (new_df['Lo'].values**2/new_df['Ro'].values**2 <= Er_max)
            new_df = new_df[mask_Er]
        if filter_EkOPm:
            mask_EkOPm = (new_df['E'].values/new_df['Pm'].values > EkOPm_min) & (new_df['E'].values/new_df['Pm'].values < EkOPm_max)
            new_df = new_df[mask_EkOPm]
        nsims = len(new_df['E'].values)
        datadict["A"]["nsims"] = nsims
        if (nsims<2):
            datadict["A"]["plot"] = False
            print("\nNot plotting A dataset!\n")
        datadict["A"]["d"] = new_df.to_dict("list")
        datadict["A"]["d"] = convertDictKeys(datadict["A"]["d"])
        # get field strengths and make fit
        datadict["A"]["E"]    = datadict["A"]["d"]["E"]
        datadict["A"]["Pm"]   = datadict["A"]["d"]["Pm"]
        datadict["A"]["Rm"]   = datadict["A"]["d"]["Ro"] / datadict["A"]["E"] * datadict["A"]["Pm"]
        datadict["A"]["fohm"] = datadict["A"]["d"]["fohm"]

        if myfohm == 1:
            fohmA = np.ones(len(datadict["A"]["d"]["fohm"]))
        else:
            datadict["A"]["fohm"] = datadict["A"]["d"]["fohm"]
        
        datadict["A"]["p"]            = datadict["A"]["d"]["p"]  
        datadict["A"]["bdip"]         = datadict["A"]["d"]["bdip"]
        datadict["A"]["fdip"]         = datadict["A"]["d"]["fdip"]
        datadict["A"]["rmsINT"]["Le"] = datadict["A"]["d"]["Lo"]  
        datadict["A"]["rmsCMB"]["Le"] = datadict["A"]["d"]["Lo"] / (datadict["A"]["d"]["bdip"] * datadict["A"]["d"]["fdip"])
        datadict["A"]["dipCMB"]["Le"] = datadict["A"]["d"]["Lo"] / datadict["A"]["d"]["bdip"]
        if not categorise:
            # fits
            datadict["A"]["rmsINT"]["ssr"], datadict["A"]["rmsINT"]["m"], datadict["A"]["rmsINT"]["c"], datadict["A"]["rmsINT"]["res"] = \
                fits(datadict["A"]["p"], datadict["A"]["rmsINT"]["Le"], datadict["A"]["fohm"])
            datadict["A"]["rmsCMB"]["ssr"], datadict["A"]["rmsCMB"]["m"], datadict["A"]["rmsCMB"]["c"], datadict["A"]["rmsCMB"]["res"] = \
                fits(datadict["A"]["p"], datadict["A"]["rmsCMB"]["Le"], datadict["A"]["fohm"])
            datadict["A"]["dipCMB"]["ssr"], datadict["A"]["dipCMB"]["m"], datadict["A"]["dipCMB"]["c"], datadict["A"]["dipCMB"]["res"] = \
                fits(datadict["A"]["p"], datadict["A"]["dipCMB"]["Le"], datadict["A"]["fohm"])


    # -------------------
    # - Yadav simulations
    # -------------------
    elif dataset == "Yadav":
        cols_names = ['E','Ra','Pm','KE_pol','KE_tor','Re','Rol','Nu','ME_pol','ME_tor','Elsasser','Elsasser_CMB',\
                      'Dip','Dip_CMB','Dip_CMB_l11','Buo_pow','Ohm_diss','l_bar_u','l_bar_B','run_time','l_max','r_max']
        E,Ra,Pm,KE_pol,KE_tor,Re,Rol,Nu,ME_pol,ME_tor,Els,ElsCMB,Dip,DipCMB,DipCMBl11,Buo_pow,Ohm_diss,l_bar_u,l_bar_B,\
            run_time,l_max,r_max = np.loadtxt(infname, usecols=tuple(np.arange(len(cols_names))), skiprows=1, unpack='true')
        # Create data dictionary
        data_dict = {'E':E,'Ra':Ra,'Pm':Pm,'KE_pol':KE_pol,'KE_tor':KE_tor,'Re':Re,'Rol':Rol,'Nu':Nu,\
                     'ME_pol':ME_pol,'ME_tor':ME_tor,'Elsasser':Els,'Elsasser_CMB':ElsCMB,'Dip':Dip,'Dip_CMB':DipCMB,\
                     'Dip_CMB_l11':DipCMBl11,'Buo_pow':Buo_pow,'Ohm_diss':Ohm_diss,'l_bar_u':l_bar_u,'l_bar_B':l_bar_B,\
                     'run_time':run_time,'l_max':l_max,'r_max':r_max}
        # Create data frame
        df = pd.DataFrame(data_dict, columns=cols_names)
        if filter_fdip:
            mask_fdip = (df['Dip_CMB_l11'] > fdip_min) & (df['Dip_CMB_l11'] < fdip_max)
            new_df = df[mask_fdip]
        else:
            new_df = df
        if filter_Er:
            mask_Er = ((new_df['ME_tor'].values+new_df['ME_pol'].values)/(new_df['KE_tor'].values+new_df['KE_pol'].values) >= Er_min) &\
                      ((new_df['ME_tor'].values+new_df['ME_pol'].values)/(new_df['KE_tor'].values+new_df['KE_pol'].values) <= Er_max)
            new_df = new_df[mask_Er]
        if filter_EkOPm:
            mask_EkOPm = (new_df['E'].values/new_df['Pm'].values > EkOPm_min) & (new_df['E'].values/new_df['Pm'].values < EkOPm_max)
            new_df = new_df[mask_EkOPm]
        nsims = len(new_df['E'].values)
        datadict["Y"]["nsims"] = nsims
        if (nsims<2):
            datadict["Y"]["plot"] = False
            print("\nNot plotting Y dataset!\n")
        datadict["Y"]["d"] = new_df.to_dict("list")
        datadict["Y"]["d"] = convertDictKeys(datadict["Y"]["d"])

        # Calculate output measures and perform fits
        datadict["Y"]["E"]  = datadict["Y"]["d"]["E"]
        datadict["Y"]["Pm"] = datadict["Y"]["d"]["Pm"]
        datadict["Y"]["Rm"] = datadict["Y"]["d"]["Re"]*datadict["Y"]["d"]["Pm"]
        if myfohm == 1:
            datadict["Y"]["fohm"] = np.ones(len(datadict["Y"]["d"]["Ohm_diss"]))
        else:
            datadict["Y"]["fohm"] = datadict["Y"]["d"]["Ohm_diss"]/datadict["Y"]["d"]["Buo_pow"]

        datadict["Y"]["rmsINT"]["Le"] = np.sqrt(datadict["Y"]["d"]["Elsasser"]*datadict["Y"]["E"]/datadict["Y"]["Pm"])
        datadict["Y"]["rmsCMB"]["Le"] = np.sqrt(datadict["Y"]["d"]["Elsasser_CMB"]*datadict["Y"]["E"]/datadict["Y"]["Pm"])
        datadict["Y"]["dipCMB"]["Le"] = np.sqrt(datadict["Y"]["d"]["Elsasser_CMB"]*datadict["Y"]["d"]["Dip_CMB"]*datadict["Y"]["E"]/datadict["Y"]["Pm"])
        datadict["Y"]["p"]            = datadict["Y"]["E"]**3 * datadict["Y"]["d"]["Buo_pow"]/shellVolume(0.35)
        datadict["Y"]["bdip"]         = datadict["Y"]["rmsINT"]["Le"]/datadict["Y"]["dipCMB"]["Le"] 
        datadict["Y"]["fdip"]         = datadict["Y"]["d"]["Dip_CMB_l11"]

        if not categorise:
            # fits
            datadict["Y"]["rmsINT"]["ssr"], datadict["Y"]["rmsINT"]["m"], datadict["Y"]["rmsINT"]["c"], datadict["Y"]["rmsINT"]["res"] = \
                fits(datadict["Y"]["p"], datadict["Y"]["rmsINT"]["Le"], datadict["Y"]["fohm"])
            datadict["Y"]["rmsCMB"]["ssr"], datadict["Y"]["rmsCMB"]["m"], datadict["Y"]["rmsCMB"]["c"], datadict["Y"]["rmsCMB"]["res"] = \
                fits(datadict["Y"]["p"], datadict["Y"]["rmsCMB"]["Le"], datadict["Y"]["fohm"])
            datadict["Y"]["dipCMB"]["ssr"], datadict["Y"]["dipCMB"]["m"], datadict["Y"]["dipCMB"]["c"], datadict["Y"]["dipCMB"]["res"] = \
                fits(datadict["Y"]["p"], datadict["Y"]["dipCMB"]["Le"], datadict["Y"]["fohm"])

        if chk == 1:
            # CHECK: Relate Elsasser to EM, total magnetic energy 
            # Em = Els*Els / (2*Pm*E)
            print('*******Yadav***************')
            EfromLe  = datadict["Y"]["rmsINT"]["Le"]**2 * 14.59 / 2.0 / datadict["Y"]["E"]**2
            Efromfile = datadict["Y"]["d"]["ME_pol"]+datadict["Y"]["d"]["ME_tor"]
            print('Energy difference (%) = ', np.abs( (EfromLe-Efromfile) / Efromfile)*100.0)
            print('***************************\n')

    # -------------------------
    # - Christensen simulations
    # -------------------------
    # Ek has no factor 2
    # Fundamental length scale is shell thickness, time scale is viscous diffusion time, magnetic field scale is "Elsasser scale".
    # NOTE that variables names "Els" are actually RMS field strengths!!
    elif dataset == "Christensen":
        cols_names = ['E','Ra','Rac','Pr','Pm','icc','lmax','nr','mc','Re','Rm','Ro','Rol','Ekin','po','ta','lm',
                      'mmd','Nu','Emag','ex','po1','pa','ta1','lm1','mmd1','B','Bsur','B12','Bdip','Bax','pow','jou',
                      'tmag','tvis','trun','tad','Tm','T2','R2R1','Filename']
        str_header0 = '%l%E     Ra  Rac Pr Pm icc lmax nr mc  Re  Rm  Ro  Rol Ekin %po %ta lm mmd  Nu  Emag %ex %po %pa %ta lm mmd B  Bsur  B12  Bdip Bax pow %jou tmag tvis trun tad Tm  T2 R2R1'
        E,Ra,Rac,Pr,Pm,icc,lmax,nr,mc,Re,Rm,Ro,Rol,Ekin,po,ta,lm,mmd,Nu,Emag,ex,po1,pa,ta1,lm1,mmd1,B,Bsur,B12,Bdip,Bax,pow,jou,tmag,tvis,trun,tad,Tm,T2,R2R1 = np.loadtxt(infname, usecols=tuple(np.arange(len(cols_names)-1)), skiprows=2, unpack='true')
        # Create data dictionary
        data_dict = {'E':E,'Ra':Ra,'Rac':Rac,'Pr':Pr,'Pm':Pm,'icc':icc,'lmax':lmax,'nr':nr,'mc':mc,'Re':Re,'Rm':Rm,'Ro':Ro,'Rol':Rol,'Ekin':Ekin,'po':po,'ta':ta,'lm':lm,'mmd':mmd,'Nu':Nu,'Emag':Emag,'ex':ex,'po1':po1,'pa':pa,'ta1':ta1,'lm1':lm1,'mmd1':mmd1,'B':B,'Bsur':Bsur,'B12':B12,'Bdip':Bdip,'Bax':Bax,'pow':pow,'jou':jou,'tmag':tmag,'tvis':tvis,'trun':trun,'tad':tad,'Tm':Tm,'T2':T2,'R2R1':R2R1}
        # Create dataframe
        df = pd.DataFrame(data_dict, columns=cols_names[0:-1])
        if filter_fdip:
            mask_fdip = (df['Bdip'].values/df['B12'].values > fdip_min) & (df['Bdip'].values/df['B12'].values < fdip_max)
            new_df = df[mask_fdip]
        else:
            new_df = df
        if filter_Er:
            mask_Er = ((new_df['Emag'].values/new_df['Ekin'].values) >= Er_min) & ((new_df['Emag'].values/new_df['Ekin'].values) <= Er_max)
            new_df = new_df[mask_Er]
        if filter_EkOPm:
            mask_EkOPm = (new_df['E'].values/new_df['Pm'].values > EkOPm_min) & (new_df['E'].values/new_df['Pm'].values < EkOPm_max)
            new_df = new_df[mask_EkOPm]
        nsims = len(new_df['E'].values)
        datadict["UC"]["nsims"] = nsims
        if (nsims<2):
            datadict["UC"]["plot"] = False
            print("\nNot plotting UC dataset!\n")
        datadict["UC"]["d"] = new_df.to_dict("list")
        datadict["UC"]["d"] = convertDictKeys(datadict["UC"]["d"])
        # Calculate field strengths and do fits
        datadict["UC"]["E"]  = datadict["UC"]["d"]["E"]
        datadict["UC"]["Pm"] = datadict["UC"]["d"]["Pm"]  
        datadict["UC"]["Rm"] = datadict["UC"]["d"]["Rm"]  
        datadict["UC"]["p"]  = 1e7 * datadict["UC"]["d"]["pow"] * datadict["UC"]["E"]**3
        if myfohm == 1:
            datadict["UC"]["fohm"] = np.ones(len(datadict["UC"]["d"]["jou"]))
        else:
            datadict["UC"]["fohm"] = datadict["UC"]["d"]["jou"] / 100

        datadict["UC"]["rmsINT"]["Le"] = np.sqrt(datadict["UC"]["E"]/datadict["UC"]["Pm"]) * datadict["UC"]["d"]["B"] # Eqn 14 of CA06
        datadict["UC"]["rmsCMB"]["Le"] = np.sqrt(datadict["UC"]["E"]/datadict["UC"]["Pm"]) * datadict["UC"]["d"]["Bsur"]
        datadict["UC"]["dipCMB"]["Le"] = np.sqrt(datadict["UC"]["E"]/datadict["UC"]["Pm"]) * datadict["UC"]["d"]["Bdip"]
        datadict["UC"]["bdip"]         = datadict["UC"]["rmsINT"]["Le"]/datadict["UC"]["dipCMB"]["Le"]
        datadict["UC"]["fdip"]         = datadict["UC"]["d"]["Bdip"]/datadict["UC"]["d"]["B12"]
        if not categorise:
            # fits
            datadict["UC"]["rmsINT"]["ssr"], datadict["UC"]["rmsINT"]["m"], datadict["UC"]["rmsINT"]["c"], datadict["UC"]["rmsINT"]["res"] = \
                fits(datadict["UC"]["p"], datadict["UC"]["rmsINT"]["Le"], datadict["UC"]["fohm"])
            datadict["UC"]["rmsCMB"]["ssr"], datadict["UC"]["rmsCMB"]["m"], datadict["UC"]["rmsCMB"]["c"], datadict["UC"]["rmsCMB"]["res"] = \
                fits(datadict["UC"]["p"], datadict["UC"]["rmsCMB"]["Le"], datadict["UC"]["fohm"])
            datadict["UC"]["dipCMB"]["ssr"], datadict["UC"]["dipCMB"]["m"], datadict["UC"]["dipCMB"]["c"], datadict["UC"]["dipCMB"]["res"] = \
                fits(datadict["UC"]["p"], datadict["UC"]["dipCMB"]["Le"], datadict["UC"]["fohm"])

        if chk == 1:
            # CHECK: Relate Elsasser to Em, magnetic energy density (i.e. per unit volume)
            # Em = Els / (2*Pm*E) = Le**2 / (2 E**2) 
            # P is printed to chk against Christensen 10 Fig 2. Note he uses Rc instead of D for 
            #   length and there is a factor 1/4pi in the Ra_Q definition. 
            #   NOTE - I still do not understand how to calculate his F. 
            print('\n*******Christensen FF******')
            EfromLe  = datadict["UC"]["rmsINT"]["Le"]**2 * datadict["UC"]["Pm"] / 2.0 / datadict["UC"]["E"]**2 / datadict["UC"]["Pm"]
            Efromfile = datadict["UC"]["d"]["Emag"]
            print('Energy difference (%) = ', np.abs( (EfromLe-Efromfile) / Efromfile)*100.0)
            print('Power at E=3e-5       = ', datadict["UC"]["p"][26:35] / (1-0.35)**2 *4.0*np.pi)
            print('***************************\n')


    # -------------------------
    # - Christensen simulations
    # -------------------------
    elif dataset == "ChristensenT":
        cols_names = ['E','Ra','Rac','Pr','Pm','icc','qs','lmax','nr','mc','Re','Rm','Ro','Rol','Ekin','po','ta','lm','mmd','Nu','Emag','ex','po1','pa','ta1','lm1','mmd1','B','Bsur','B12','Bdip','Bax','pow','jou','tmag','tvis','trun','tad','Tm','T2','vrat','Filename']
        str_header = '%l%E     Ra  Rac Pr Pm icc  qs  lmax nr mc  Re  Rm  Ro  Rol Ekin %po %ta lm mmd  Nu  Emag %ex %po %pa %ta lm mmd B  Bsur B12  Bdip Bax pow %jou tmag tvis trun tad Tm  T2 vrat'
        E,Ra,Rac,Pr,Pm,icc,qs,lmax,nr,mc,Re,Rm,Ro,Rol,Ekin,po,ta,lm,mmd,Nu,Emag,ex,po1,pa,ta1,lm1,mmd1,B,Bsur,B12,Bdip,Bax,pow,jou,tmag,tvis,trun,tad,Tm,T2,vrat = np.loadtxt(infname, usecols=tuple(np.arange(len(cols_names)-1)), skiprows=2, unpack='true')
        # Create data dictionary
        data_dict = {'E':E,'Ra':Ra,'Rac':Rac,'Pr':Pr,'Pm':Pm,'icc':icc,'qs':qs,'lmax':lmax,'nr':nr,'mc':mc,'Re':Re,'Rm':Rm,'Ro':Ro,'Rol':Rol,'Ekin':Ekin,'po':po,'ta':ta,'lm':lm,'mmd':mmd,'Nu':Nu,'Emag':Emag,'ex':ex,'po1':po1,'pa':pa,'ta1':ta1,'lm1':lm1,'mmd1':mmd1,'B':B,'Bsur':Bsur,'B12':B12,'Bdip':Bdip,'Bax':Bax,'pow':pow,'jou':jou,'tmag':tmag,'tvis':tvis,'trun':trun,'tad':tad,'Tm':Tm,'T2':T2,'vrat':vrat}
        df = pd.DataFrame(data_dict, columns=cols_names[0:-1])
        # Filter fdip
        mask_fdip = (df['Bdip'].values/df['B12'].values > fdip_min) & (df['Bdip'].values/df['B12'].values < fdip_max)
        new_df = df[mask_fdip]
        if filter_Er:
            mask_Er = ((new_df['Emag'].values/new_df['Ekin'].values) >= Er_min) & ((new_df['Emag'].values/new_df['Ekin'].values) <= Er_max)
            new_df = new_df[mask_Er]
        if filter_EkOPm:
            mask_EkOPm = (new_df['E'].values/new_df['Pm'].values > EkOPm_min) & (new_df['E'].values/new_df['Pm'].values < EkOPm_max)
            new_df = new_df[mask_EkOPm]
        nsims = len(new_df['E'].values)
        datadict["UCt"]["nsims"] = nsims
        if (nsims<2):
            datadict["UCt"]["plot"] = False
            print("\nNot plotting UCt dataset!\n")
        datadict["UCt"]["d"] = new_df.to_dict("list")
        datadict["UCt"]["d"] = convertDictKeys(datadict["UCt"]["d"])
        # calculate field strengths and do fits
        datadict["UCt"]["E"] = datadict["UCt"]["d"]["E"]
        datadict["UCt"]["Pm"] = datadict["UCt"]["d"]["Pm"]
        datadict["UCt"]["Rm"] = datadict["UCt"]["d"]["Rm"]
        if myfohm == 1:
            datadict["UCt"]["fohm"] = np.ones(len(datadict["UCt"]["d"]["jou"]))
        else:
            datadict["UCt"]["fohm"] = datadict["UCt"]["d"]["jou"] / 100

        datadict["UCt"]["rmsINT"]["Le"] = np.sqrt(datadict["UCt"]["E"]/datadict["UCt"]["Pm"]) * datadict["UCt"]["d"]["B"]
        datadict["UCt"]["rmsCMB"]["Le"] = np.sqrt(datadict["UCt"]["E"]/datadict["UCt"]["Pm"]) * datadict["UCt"]["d"]["Bsur"]
        datadict["UCt"]["dipCMB"]["Le"] = np.sqrt(datadict["UCt"]["E"]/datadict["UCt"]["Pm"]) * datadict["UCt"]["d"]["Bdip"]
        datadict["UCt"]["p"]            = 1e7 * datadict["UCt"]["d"]["pow"] * datadict["UCt"]["E"]**3
        datadict["UCt"]["bdip"]         = datadict["UCt"]["rmsINT"]["Le"]/datadict["UCt"]["dipCMB"]["Le"]
        datadict["UCt"]["fdip"]         = datadict["UCt"]["d"]["Bdip"]/datadict["UCt"]["d"]["B12"]
        if not categorise:
            # fits
            datadict["UCt"]["rmsINT"]["ssr"], datadict["UCt"]["rmsINT"]["m"], datadict["UCt"]["rmsINT"]["c"], datadict["UCt"]["rmsINT"]["res"] = \
                fits(datadict["UCt"]["p"], datadict["UCt"]["rmsINT"]["Le"], datadict["UCt"]["fohm"])
            datadict["UCt"]["rmsCMB"]["ssr"], datadict["UCt"]["rmsCMB"]["m"], datadict["UCt"]["rmsCMB"]["c"], datadict["UCt"]["rmsCMB"]["res"] = \
                fits(datadict["UCt"]["p"], datadict["UCt"]["rmsCMB"]["Le"], datadict["UCt"]["fohm"])
            datadict["UCt"]["dipCMB"]["ssr"], datadict["UCt"]["dipCMB"]["m"], datadict["UCt"]["dipCMB"]["c"], datadict["UCt"]["dipCMB"]["res"] = \
                fits(datadict["UCt"]["p"], datadict["UCt"]["dipCMB"]["Le"], datadict["UCt"]["fohm"])

        if chk == 1:
            # CHECK: Relate Elsasser to Em, magnetic energy density (i.e. per unit volume)
            # Em = Els / (2*Pm*E) = Le**2 / (2 E**2) 
            print('*******Christensen qT******')
            EfromLe  = datadict["UCt"]["rmsINT"]["Le"]**2 * datadict["UCt"]["Pm"] / 2.0 / datadict["UCt"]["E"]**2 / datadict["UCt"]["Pm"]
            Efromfile = datadict["UCt"]["Emag"]
            print('Energy difference (%) = ', np.abs( (EfromLe-Efromfile) / Efromfile)*100.0)
            print('***************************\n')

    # -------------------------
    # - Aubert 2017 and 2019 simulations
    # -------------------------        
    elif dataset == "APath":
        cols_names = ['eps','E','rmsCMBtotal','RMSCMBl=12','RMSCMBl=1','p','Rm','Pm','fohm','Els']
        
        eps,E,rmsCMBtotal,RMSCMBl12,RMSCMBl1,p,Rm,Pm,fohm,Els = np.loadtxt(infname, usecols=tuple(np.arange(len(cols_names))), skiprows=5, unpack='true')
        # Create data dictionary
        data_dict = {'eps':eps,'E':E,'rmsCMBtotal':rmsCMBtotal,'RMSCMBl=12':RMSCMBl12,'RMSCMBl=1':RMSCMBl1,'p':p,'Rm':Rm,'Pm':Pm,'fohm':fohm,'Els':Els}
        # Create dataframe
        df = pd.DataFrame(data_dict, columns=cols_names)
        # Filter fdip
        mask_fdip = (df['RMSCMBl=1'].values/df['RMSCMBl=12'].values > fdip_min) & (df['RMSCMBl=1'].values/df['RMSCMBl=12'].values < fdip_max)
        new_df = df[mask_fdip]
        #if filter_Er:
        #    mask_Er = ((new_df['Emag'].values/new_df['Ekin'].values) >= Er_min) & ((new_df['Emag'].values/new_df['Ekin'].values) <= Er_max)
        #    new_df = new_df[mask_Er]
        if filter_EkOPm:
            mask_EkOPm = (new_df['E'].values/new_df['Pm'].values > EkOPm_min) & (new_df['E'].values/new_df['Pm'].values < EkOPm_max)
            new_df = new_df[mask_EkOPm]
        nsims = len(new_df['E'].values)
        if (nsims<2):
            datadict["APath"]["plot"] = False
            print("\nNot plotting Apath dataset!\n")
        datadict["APath"]['nsims'] = nsims 
        datadict["APath"]["d"] = new_df.to_dict("list")
        datadict["APath"]["d"] = convertDictKeys(datadict["APath"]["d"])
        
        datadict["APath"]['fohm'] = datadict["APath"]['d']['fohm']
        datadict["APath"]['p']    = datadict["APath"]['d']['p']
        datadict["APath"]['E']    = datadict["APath"]['d']['E']
        datadict["APath"]['Pm']   = datadict["APath"]['d']['Pm']
        datadict["APath"]['Rm']   = datadict["APath"]['d']['Rm']
        datadict["APath"]['fdip'] = datadict["APath"]['d']['RMSCMBl=1']/datadict["APath"]['d']['RMSCMBl=12']
        datadict["APath"]['bdip'] = np.zeros(len(datadict["APath"]['fohm']))
        # get field strengths 
        datadict["APath"]['rmsINT']['Le'] = np.sqrt( (datadict["APath"]['E']/datadict["APath"]['Pm'])   * datadict["APath"]['d']['Els'])
        datadict["APath"]['rmsCMB']['Le'] = np.sqrt( (datadict["APath"]['E']/datadict["APath"]['Pm']) ) * datadict["APath"]['d']['rmsCMBtotal']
        datadict["APath"]['dipCMB']['Le'] = np.sqrt( (datadict["APath"]['E']/datadict["APath"]['Pm']) ) * datadict["APath"]['d']['RMSCMBl=1']
        if not categorise:
            # fit
            datadict["APath"]['rmsINT']['ssr'], datadict["APath"]['rmsINT']['m'], datadict["APath"]['rmsINT']['c'], datadict["APath"]['rmsINT']['res'] = \
                fits(datadict["APath"]['p'], datadict["APath"]['rmsINT']['Le'], datadict["APath"]['fohm'])
            datadict["APath"]['rmsCMB']['ssr'], datadict["APath"]['rmsCMB']['m'], datadict["APath"]['rmsCMB']['c'], datadict["APath"]['rmsCMB']['res'] = \
                fits(datadict["APath"]['p'], datadict["APath"]['rmsCMB']['Le'], datadict["APath"]['fohm'])
            datadict["APath"]['dipCMB']['ssr'], datadict["APath"]['dipCMB']['m'], datadict["APath"]['dipCMB']['c'], datadict["APath"]['dipCMB']['res'] = \
                fits(datadict["APath"]['p'], datadict["APath"]['dipCMB']['Le'], datadict["APath"]['fohm'])

    # -------------------------
    # Schwaiger 2019
    # -------------------------        
    elif dataset == "S":
        # -- DGM WHATS THE DIFFERENCE BETWEEN fdip1 and fdip2. fdip1 FOLLOWS THE DEF IN CA06?
        cols_names = ['E','Pm','rmsCMBtotal','RMSCMBl=1','fdip1','fdip2','p','Rm','fohm', 'Els', 'MEKE']
        
        E,Pm,rmsCMBtotal,RMSCMBl1,fdip1,fdip2,p,Rm,fohm,Els,MEKE = np.loadtxt(infname, usecols=tuple(np.arange(len(cols_names))), skiprows=0, unpack='true')
        # Create data dictionary
        data_dict = {'E':E,'Pm':Pm,'rmsCMBtotal':rmsCMBtotal,'RMSCMBl=1':RMSCMBl1,'fdip1':fdip1,'fdip2':fdip2,'p':p,'Rm':Rm,'fohm':fohm,'Els':Els,'MEKE':MEKE}
        # Create dataframe
        df = pd.DataFrame(data_dict, columns=cols_names)
        # Filter fdip
        mask_fdip = (df['fdip1'].values > fdip_min) & (df['fdip1'].values < fdip_max)
        new_df = df[mask_fdip]
        if filter_Er:
            mask_Er = ((new_df['MEKE'].values) >= Er_min) & ((new_df['MEKE'].values) <= Er_max)
            new_df = new_df[mask_Er]
        if filter_EkOPm:
            mask_EkOPm = (new_df['E'].values/new_df['Pm'].values > EkOPm_min) & (new_df['E'].values/new_df['Pm'].values < EkOPm_max)
            new_df = new_df[mask_EkOPm]
        datadict["S"]["d"] = new_df.to_dict("list")
        datadict["S"]["d"] = convertDictKeys(datadict["S"]["d"])
        ####
        nsims = len(new_df['E'].values)
        if (nsims<2):
            datadict["S"]["plot"] = False
            print("\nNot plotting S dataset!\n")
        datadict["S"]['nsims'] = nsims 
        datadict["S"]["d"] = new_df.to_dict("list")
        datadict["S"]["d"] = convertDictKeys(datadict["S"]["d"])
        
        datadict["S"]['fohm'] = datadict["S"]['d']['fohm']
        datadict["S"]['p']    = datadict["S"]['d']['p']
        datadict["S"]['E']    = datadict["S"]['d']['E']
        datadict["S"]['Pm']   = datadict["S"]['d']['Pm']
        datadict["S"]['Rm']   = datadict["S"]['d']['Rm']
        datadict["S"]['fdip'] = datadict["S"]['d']['fdip1']
        datadict["S"]['bdip'] = np.zeros(len(datadict["S"]['fohm']))
        
        datadict["S"]['rmsINT']['Le'] = np.sqrt( (datadict["S"]['E']/datadict["S"]['Pm'])   * datadict["S"]['d']['Els'])
        datadict["S"]['rmsCMB']['Le'] = np.sqrt( (datadict["S"]['E']/datadict["S"]['Pm']) ) * datadict["S"]['d']['rmsCMBtotal']
        datadict["S"]['dipCMB']['Le'] = np.sqrt( (datadict["S"]['E']/datadict["S"]['Pm']) ) * datadict["S"]['d']['RMSCMBl=1']

        if not categorise:
            # fits
            datadict["S"]['rmsINT']['ssr'], datadict["S"]['rmsINT']['m'], datadict["S"]['rmsINT']['c'], datadict["S"]['rmsINT']['res'] = \
                fits(datadict["S"]['p'], datadict["S"]['rmsINT']['Le'], datadict["S"]['fohm'])
            datadict["S"]['rmsCMB']['ssr'], datadict["S"]['rmsCMB']['m'], datadict["S"]['rmsCMB']['c'], datadict["S"]['rmsCMB']['res'] = \
                fits(datadict["S"]['p'], datadict["S"]['rmsCMB']['Le'], datadict["S"]['fohm'])
            datadict["S"]['dipCMB']['ssr'], datadict["S"]['dipCMB']['m'], datadict["S"]['dipCMB']['c'], datadict["S"]['dipCMB']['res'] = \
                fits(datadict["S"]['p'], datadict["S"]['dipCMB']['Le'], datadict["S"]['fohm'])
    else:
        raise ValueError("Not valid dataset")

    # Save out file
    if (filter_fdip or filter_Er or filter_EkOPm):
        print("\nSaving out file: %s" %(outfname))
        if (dataset=="Leeds"):
            new_df.to_csv(outfname, index=False, sep=",")
        else:
            new_df.to_csv(outfname, index=False, sep=" ")
        if (dataset=="Christensen"):
            # Insert second line header as in original Christensen table
            str_header1 = '                                    10^-3                                                                       10^7    10^-4          10^-4'
            shell_cmd = "sed -i '2i\\"+str_header1+"\\' "+outfname
            os.system(shell_cmd)
        # On screen outputs
        print('Number of models in dataset (original, filtered) = %s %s' %(len(df),len(new_df)))
 
    return new_df, datadict

def plot_bdip(datadict, myfdip):

    Cmax = np.log10(1000.0)#np.log10(np.max(Eall))
    Cmin = np.log10(50.0)#np.log10(np.min(Eall))
    
    # - bdip vs buoyancy power
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)
    for key in datadict:
        if datadict[key]["plot"]:
            plt.scatter(datadict[key]["p"], datadict[key]["bdip"],
                        s=datadict[key]["plotp"]["size"], marker=datadict[key]["plotp"]["marker"],
                        c=np.log10(datadict[key]["plotp"]["Col"]), vmin=Cmin, vmax=Cmax,
                        cmap=plt.get_cmap(datadict[key]["plotp"]["cmap"]), edgecolor=datadict[key]["plotp"]["edgecolor"], label=datadict[key]["plotp"]["label"]) 
    ax.set_xlabel('$P_A$')
    ax.set_ylabel('$b_{dip}$')
    plt.xlim([1.e-10,1.e-3]); plt.ylim([1.e+0,50.])
    plt.xscale("log"); plt.yscale("log")
    if myfdip == 0:
        plt.title("All models")
    elif myfdip == 1: 
        plt.title("$f_{dip}>0.5$")
    elif myfdip == 2:
        #plt.title("$%.2f\leq f_{dip}\leq %.2f$" %(fdip_min,fdip_max))
        plt.title("$%.2f< f_{dip}< %.2f$" %(fdip_min,fdip_max))
    elif myfdip == 3:
        plt.title("$%.2f< f_{dip}< %.2f$" %(fdip_min,fdip_max))

    plt.show(block=False)
    plt.tight_layout()
    plt.savefig('./fig/bdip_vs_P_fdip='+str(myfdip)+'.pdf',format='pdf')
    del ax

def mergeDict(d1, d2):
    if (d1.keys() != d2.keys()):
        print(d1.keys())
        print(d2.keys())
        raise ValueError("d1 and d2 do not have the same keys.")
    
    if (d1["plot"] and d2["plot"]):
        newd = {"rmsINT":{}, "rmsCMB":{}, "dipCMB":{}, "plotp":{}}
        newd["plot"] = True
        newd["dataset"] = d1["dataset"]+"+"+d2["dataset"]
        newd["nsims"] = d1["nsims"]+d2["nsims"]
        for key in d1:
            if (key in ["E","Pm","Rm","fohm","p","bdip","fdip"]):
                newd[key] = np.concatenate((d1[key], d2[key]))
            if (key in ["rmsINT", "rmsCMB", "dipCMB"]):
                newd[key]["Le"] = np.concatenate((d1[key]["Le"], d2[key]["Le"]))
    elif (d1["plot"] and not d2["plot"]):
        newd = copy.deepcopy(d1)
    elif (not d1["plot"] and d2["plot"]):
        newd = copy.deepcopy(d2)
    
    return newd

def filterLeedsDict(d1, fkey="TBC", condition="==FTFF"):
    newd = {"rmsINT":{}, "rmsCMB":{}, "dipCMB":{}}
    # find array elements to filter
    idx = []
    for key in d1:
        if (key in [fkey]):
            for i in range(len(d1[key])):
                if (eval("d1[key][i]"+condition)):
                    idx.append(i)
    if (len(idx)==0):
        print("")
        print("Nothing to filter for "+fkey+condition)
        d1["plot"] = False
        return d1
    else:
        idx = np.asarray(idx)
        # filter all keys
        for key in d1:
            if (key in ["E","Pm","Rm","fohm","p","bdip","fdip","TBC","Ir","ar"]):
                newd[key] = d1[key][idx]
            if (key in ["rmsINT", "rmsCMB", "dipCMB"]):
                newd[key]["Le"] = d1[key]["Le"][idx]
        newd["plot"]    = False
        newd["nsims"]   = len(newd["E"])
        newd["dataset"] = d1["dataset"]
        newd["plotp"]   = d1["plotp"]
        return newd

def redefineDataDict(indict, quiet=False):
    """
    redefines data dictionaries by author to different
    categorisations.
    """
    # clean up author's dictionaries from not useful keys
    for k in indict.keys():
        del indict[k]["d"]

    # define new dictionary
    newdict = {"FTFT":{}, "FTFF":{}, "FF0F":{}, "FFFF":{}, "CE":{}, "Mixed":{}}
    # merge/redefine
    newdict["FTFT"] = mergeDict(indict["S"], indict["Y"])
    newdict["FF0F"] = copy.deepcopy(indict["UC"])
    newdict["CE"]   = copy.deepcopy(indict["APath"])

    # filter Leeds simulations
    dL_Ir0   = filterLeedsDict(indict["L"], fkey="Ir", condition="==0.")
    dL_FTFF  = filterLeedsDict(dL_Ir0, fkey="TBC", condition="=='FTFF'")
    dL_FFFF  = filterLeedsDict(dL_Ir0, fkey="TBC", condition="=='FFFF'")
    dL_mixed = filterLeedsDict(indict["L"], fkey="Ir", condition="!=0.")
   
    # delete unseful stuff 
    del dL_FTFF["TBC"], dL_FTFF["Ir"], dL_FTFF["ar"]
    del dL_FFFF["TBC"], dL_FFFF["Ir"], dL_FFFF["ar"]
    del dL_mixed["TBC"], dL_mixed["Ir"], dL_mixed["ar"]

    # merge with other simulations data
    newdict["FFFF"]  = copy.deepcopy(dL_FFFF)
    newdict["FTFF"]  = mergeDict(indict["UCt"], dL_FTFF)
    newdict["Mixed"] = mergeDict(indict["A"], dL_mixed)
   
    # redefine "dataset" key
    newdict["FTFT"]["dataset"]  = "FTFT" 
    newdict["FFFF"]["dataset"]  = "FFFF"
    newdict["FF0F"]["dataset"]  = "FF0F"
    newdict["FTFF"]["dataset"]  = "FTFF"
    newdict["Mixed"]["dataset"] = "Mixed"
    newdict["CE"]["dataset"]    = "CE"

    # fit field strengths of these merged datasets
    print("")
    for key in newdict:
        print("Fitting new datset: ",key)
        newdict[key]['rmsINT']['ssr'], newdict[key]['rmsINT']['m'], newdict[key]['rmsINT']['c'], newdict[key]['rmsINT']['res'] = \
            fits(newdict[key]['p'], newdict[key]['rmsINT']['Le'], newdict[key]['fohm'])
        newdict[key]['rmsCMB']['ssr'], newdict[key]['rmsCMB']['m'], newdict[key]['rmsCMB']['c'], newdict[key]['rmsCMB']['res'] = \
            fits(newdict[key]['p'], newdict[key]['rmsCMB']['Le'], newdict[key]['fohm'])
        newdict[key]['dipCMB']['ssr'], newdict[key]['dipCMB']['m'], newdict[key]['dipCMB']['c'], newdict[key]['dipCMB']['res'] = \
            fits(newdict[key]['p'], newdict[key]['dipCMB']['Le'], newdict[key]['fohm'])
    return newdict
