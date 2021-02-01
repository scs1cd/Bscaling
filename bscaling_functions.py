# -*- coding: utf-8 -*-
"""
Created on Sat Jan  9 16:15:55 2021

@author: Chris Davies
"""

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def idxStr(l1, s):
    idx = []
    i = 0
    length = len(l1)
    while i < length:
        if s == l1[i]:
            idx.append(i)
        i += 1
    return idx

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
  
def getPlotProperties(datadict):
    for key in datadict:
        if datadict[key]["plot"]:
            # CD - CHK THIS!!!
            datadict[key]["plotp"]["Col"] = np.array(datadict[key]['d']["Rm"])/np.array(datadict[key]['d']["Pm"])
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
                datadict[key]["plotp"]["label"] = "Auber et al (2009)"
            elif (datadict[key]["dataset"]=="UC"):
                datadict[key]["plotp"]["marker"] = "v"
                datadict[key]["plotp"]["size"] = 150
                datadict[key]["plotp"]["cmap"] = "Purples"
                datadict[key]["plotp"]["edgecolor"] = "purple"
                datadict[key]["plotp"]["label"] = "Christensen 0F"
            elif (datadict[key]["dataset"]=="UCt"):
                datadict[key]["plotp"]["marker"] = "^"
                datadict[key]["plotp"]["size"] = 150
                datadict[key]["plotp"]["cmap"] = "Purples"
                datadict[key]["plotp"]["edgecolor"] = "purple"
                datadict[key]["plotp"]["label"] = "Christensen FF"
            elif (datadict[key]["dataset"]=="Y"):
                datadict[key]["plotp"]["marker"] = "o"
                datadict[key]["plotp"]["size"] = 150
                datadict[key]["plotp"]["cmap"] = "Reds"
                datadict[key]["plotp"]["edgecolor"] = "red"
                datadict[key]["plotp"]["label"] = "Yadav et al (2010)"
            elif (datadict[key]["dataset"]=="APath"):
                datadict[key]["plotp"]["marker"] = "D"
                datadict[key]["plotp"]["size"] = 150
                datadict[key]["plotp"]["cmap"] = "Greys"
                datadict[key]["plotp"]["edgecolor"] = "grey"
                datadict[key]["plotp"]["label"] = "Aubert et al (2017); Aubert (2019)"
            elif (datadict[key]["dataset"]=="S"):
                datadict[key]["plotp"]["marker"] = "*"
                datadict[key]["plotp"]["size"] = 150
                datadict[key]["plotp"]["cmap"] = "Greens"
                datadict[key]["plotp"]["edgecolor"] = "green"
                datadict[key]["plotp"]["label"] = "Schwaiger et al (2019)"
            #elif (datadict[key]["dataset"]=="A17"):
            #    print("Symbol properties to be chosen")
            #elif (datadict[key]["dataset"]=="APath"):
            #    print("Symbol properties to be chosen")
            #elif (datadict[key]["dataset"]=="S"):
            #    print("Symbol properties to be chosen")
            else:
                raise ValueError("Not valid dataset")

    return datadict

def plotSimulations(datadict=None, alldatadict=None, earthdict=None, field="rmsINT",
                    xrange=[0.,1.], yrange=[0.,1.], cbarrange=[0.,1.]):

    if field not in ("rmsINT","rmsCMB","dipCMB"):
        raise ValueError("Not valid field provided.")

    xmin = xrange[0]; xmax=xrange[1]
    ymin = yrange[0]; ymax=yrange[1]
    cbarmin = cbarrange[0]; cbarmax = cbarrange[1]

    # Now start the plot
    plt.clf()
    ax  = plt.gca()
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])

    for key in datadict:
        if datadict[key]["plot"]:
            plt.scatter(datadict[key]['d']["p"], datadict[key][field]["Le"]/datadict[key]['d']["fohm"]**0.5,
                        s=datadict[key]["plotp"]["size"], marker=datadict[key]["plotp"]["marker"],
                        c=np.log10(datadict[key]["plotp"]["Col"]), vmin=cbarmin, vmax=cbarmax,
                        cmap=plt.get_cmap(datadict[key]["plotp"]["cmap"]), edgecolor=datadict[key]["plotp"]["edgecolor"], label=datadict[key]["plotp"]["label"])

    ax.add_patch(Rectangle(xy=(earthdict["p"]["min"],earthdict[field]["min"]),
                           width=(earthdict["p"]["max"]-earthdict["p"]["min"]),
                           height=(earthdict[field]["max"]-earthdict[field]["min"]),
                           linewidth=1, color='black', fill=True))
    legend_xpos = 0.01 # x-position of legend (>1, outside of main plot)
    legend_ypos = 0.95; legend_dy   = 0.05
    iplt = 0
    for key in datadict:
        if datadict[key]["plot"]:
            plt.text(legend_xpos, legend_ypos-iplt*legend_dy,
                     "$m$  = "+str(np.round(datadict[key][field]["m"],2)) +\
                     "$\pm$"+str(np.round(datadict[key][field]["res"],2)) +\
                     ", SSR="+str(np.round(datadict[key][field]["ssr"],2)),
                     transform=ax.transAxes, color=datadict[key]["plotp"]["edgecolor"])
            iplt += 1

    plt.text(legend_xpos, legend_ypos-iplt*legend_dy, "$m$  = "+str(np.round(alldatadict[field]["m"],2))+\
             "$\pm$"+str(np.round(alldatadict[field]["res"],2))+\
             ", SSR="+str(np.round(alldatadict[field]["ssr"][0],2)), transform=ax.transAxes, color='black')
    # Plot color bar
    cbar = plt.colorbar()
    cbar.set_label("log $Re$")
    
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.legend(bbox_to_anchor=(legend_xpos-0.1, 1.10), loc=3, ncol=2, borderaxespad=0)
    plt.rcParams["figure.figsize"] = [15,10]
    
    return ax, legend_xpos, legend_ypos

def getPlotTitle(myfdip=None, fdip_range= [None]*2, myEr=None, Er_range=[None]*2):
    title_str = ""
    comma = False
    if (myfdip == 1): 
        title_str = "$f_{dip}>0.5$"
        comma = True
    elif (myfdip == 2):
        title_str = "$0.35< f_{dip}<0.80"
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
        raise ValueError("1 data point. No fit possible")

    return ssr[0], m, c, res

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


def get_fdip(infname=None, dataset="Leeds"):
    """
    Get value of fdip from main data file
    CD - ARE WE DOING loadtxt TWICE??
    """
    if dataset == "Leeds":
        df = pd.read_csv(infname)
        fdip = df['cmb_diptyAve'].values
    elif dataset == "Aubert":
        cols_names = ['E','RaQ','Pr','Pm','chi','fi','Ro','Lo','bdip','fdip','lbar','tauDiss','p','fohm']
        E,RaQ,Pr,Pm,chi,fi,Ro,Lo,bdip,fdip,lbar,tauDiss,p,fohm = np.loadtxt(infname, usecols=tuple(np.arange(len(cols_names))), skiprows=1, unpack='true')
    elif dataset == "Yadav":
        cols_names = ['E','Ra','Pm','KE_pol','KE_tor','Re','Rol','Nu','ME_pol','ME_tor','Elsasser','Elsasser_CMB','Dip','Dip_CMB','Dip_CMB_l11','Buo_pow','Ohm_diss','l_bar_u','l_bar_B','run_time','l_max','r_max']
        E,Ra,Pm,KE_pol,KE_tor,Re,Rol,Nu,ME_pol,ME_tor,Els,ElsCMB,Dip,DipCMB,DipCMBl11,Buo_pow,Ohm_diss,l_bar_u,l_bar_B,run_time,l_max,r_max = np.loadtxt(infname, usecols=tuple(np.arange(len(cols_names))), skiprows=1, unpack='true')
        fdip = DipCMBl11
    elif dataset == "Christensen":
        cols_names = ['E','Ra','Rac','Pr','Pm','icc','lmax','nr','mc','Re','Rm','Ro','Rol','Ekin','po','ta','lm','mmd','Nu','Emag','ex','po1','pa','ta1','lm1','mmd1','B','Bsur','B12','Bdip','Bax','pow','jou','tmag','tvis','trun','tad','Tm','T2','R2R1','Filename']
        E,Ra,Rac,Pr,Pm,icc,lmax,nr,mc,Re,Rm,Ro,Rol,Ekin,po,ta,lm,mmd,Nu,Emag,ex,po1,pa,ta1,lm1,mmd1,B,Bsur,B12,Bdip,Bax,pow,jou,tmag,tvis,trun,tad,Tm,T2,R2R1 = np.loadtxt(infname, usecols=tuple(np.arange(len(cols_names)-1)), skiprows=2, unpack='true')
        fdip = Bdip/B12
    elif dataset == "ChristensenT":
        cols_names = ['E','Ra','Rac','Pr','Pm','icc','qs','lmax','nr','mc','Re','Rm','Ro','Rol','Ekin','po','ta','lm','mmd','Nu','Emag','ex','po1','pa','ta1','lm1','mmd1','B','Bsur','B12','Bdip','Bax','pow','jou','tmag','tvis','trun','tad','Tm','T2','vrat','Filename']
        E,Ra,Rac,Pr,Pm,icc,qs,lmax,nr,mc,Re,Rm,Ro,Rol,Ekin,po,ta,lm,mmd,Nu,Emag,ex,po1,pa,ta1,lm1,mmd1,B,Bsur,B12,Bdip,Bax,pow,jou,tmag,tvis,trun,tad,Tm,T2,vrat = np.loadtxt(infname, usecols=tuple(np.arange(len(cols_names)-1)), skiprows=1, unpack='true')
        fdip = Bdip/B12
    else:
        raise ValueError("Not valid dataset name provided.")
    return fdip

def filter_table(infname=None, outfname=None, dataset="Leeds", fdip_range=None, EkOPm_range=None, EMoEK_range=None, datadict=None):
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
        
    # ---------------------------------
    # - Aubert et al (2009) simulations
    # ---------------------------------
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
        print('Aubert p = ', datadict["A"]["d"]['p'])
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
    # -------------------------
    # - Christensen simulations
    # -------------------------
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

        new_df = df
        nsims = len(new_df['E'].values)
        if (nsims<2):
            datadict["APath"]["plot"] = False
            print("\nNot plotting UC dataset!\n")
        datadict["APath"]['nsims'] = nsims 
        datadict["APath"]['d'] = df        
    # -------------------------
    # Schwaiger 2019
    # -------------------------        
    elif dataset == "S":
        cols_names = ['E','Pm','rmsCMBtotal','RMSCMBl=1','fdip1','fdip2','p','Rm','fohm', 'Els', 'MEKE']
        
        E,Pm,rmsCMBtotal,RMSCMBl1,fdip1,fdip2,p,Rm,fohm,Els,MEKE = np.loadtxt(infname, usecols=tuple(np.arange(len(cols_names))), skiprows=0, unpack='true')
        # Create data dictionary
        data_dict = {'E':E,'Pm':Pm,'rmsCMBtotal':rmsCMBtotal,'RMSCMBl=1':RMSCMBl1,'fdip1':fdip1,'fdip2':fdip2,'p':p,'Rm':Rm,'fohm':fohm,'Els':Els,'MEKE':MEKE}
        # Create dataframe
        df = pd.DataFrame(data_dict, columns=cols_names)

        new_df = df
        nsims = len(new_df['E'].values)
        if (nsims<2):
            datadict["S"]["plot"] = False
            print("\nNot plotting UC dataset!\n")
        datadict["S"]['nsims'] = nsims 
        datadict["S"]['d'] = df        
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