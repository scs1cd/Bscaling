# -*- coding: utf-8 -*-
"""
Created on Sat Jan  9 16:15:55 2021

@author: Chris Davies
"""

import numpy as np
import pandas as pd
import os

def idxStr(l1, s):
    idx = []
    i = 0
    length = len(l1)
    while i < length:
        if s == l1[i]:
            idx.append(i)
        i += 1
    return idx

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

    return ssr, m, c, res

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
            
    # CD - NEED AN ELSE?
    if EMoEK_range is not None:
        filter_Er = True
        Er_min = EMoEK_range[0]; Er_max = EMoEK_range[1]
        print('\nFiltering by EM/EK (min,max) = %s %s' %(Er_min,Er_max))
        if (Er_min > Er_max):
            raise ValueError('Not valid range of EM/EK: min > max.')
    else:
        filter_Er = False
    if EkOPm_range is not None:
        filter_EkOPm = True
    else:
        filter_EkOPm = False
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
        if (nsims<2):
            datadict["L"] = False
            print("\nNot plotting L dataset!\n")
        datadict["nL"] = nsims
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
        if (nsims<2):
            datadict["A"] = False
            print("\nNot plotting A dataset!\n")
        datadict["nA"] = nsims
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
        if (nsims<2):
            datadict["Y"] = False
            print("\nNot plotting Y dataset!\n")
        datadict["nY"] = nsims
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
        if (nsims<2):
            datadict["UC"] = False
            print("\nNot plotting UC dataset!\n")
        datadict["nUC"] = nsims
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
        if (nsims<2):
            datadict["UCt"] = False
            print("\nNot plotting UCt dataset!\n")
        datadict["nUCt"] = nsims
    # -------------------------
    # - Aubert 2017 and 2019 simulations
    # -------------------------        
    elif dataset == "APath":
        cols_names = ['eps','E','rmsCMBtotal','RMSCMBl=12','RMSCMBl=1','p','Rm','Pm','fohm']
        
        eps,E,rmsCMBtotal,RMSCMBl12,RMSCMBl1,p,Rm,Pm,fohm = np.loadtxt(infname, usecols=tuple(np.arange(len(cols_names))), skiprows=5, unpack='true')
        # Create data dictionary
        data_dict = {'eps':eps,'E':E,'rmsCMBtotal':rmsCMBtotal,'RMSCMBl=12':RMSCMBl12,'RMSCMBl=1':RMSCMBl1,'p':p,'Rm':Rm,'Pm':Pm,'fohm':fohm}
        # Create dataframe
        df = pd.DataFrame(data_dict, columns=cols_names)

        new_df = df
        nsims = len(new_df['E'].values)
        if (nsims<2):
            datadict["APath"] = False
            print("\nNot plotting UC dataset!\n")
        datadict["APath"]['YN'] = nsims 
        datadict["APath"]['data'] = df        
        # -------------------------
    # Schwaiger 2019
    # -------------------------        
    elif dataset == "S":
        cols_names = ['E','Pm','rmsCMBtotal','RMSCMBl=1','fdip1','fdip2','p','Rm','fohm']
        
        E,Pm,rmsCMBtotal,RMSCMBl1,fdip1,fdip2,p,Rm,fohm = np.loadtxt(infname, usecols=tuple(np.arange(len(cols_names))), skiprows=0, unpack='true')
        # Create data dictionary
        data_dict = {'E':E,'Pm':Pm,'rmsCMBtotal':rmsCMBtotal,'RMSCMBl=1':RMSCMBl1,'fdip1':fdip1,'fdip2':fdip2,'p':p,'Rm':Rm,'fohm':fohm}
        # Create dataframe
        df = pd.DataFrame(data_dict, columns=cols_names)

        new_df = df
        nsims = len(new_df['E'].values)
        if (nsims<2):
            datadict["APath"] = False
            print("\nNot plotting UC dataset!\n")
        datadict["S"]['YN'] = nsims 
        datadict["S"]['data'] = df        
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

#def plot(x, datadict):
    
#     xmin, xmax, ymin, ymax = x[0], x[1], x[2], x[3]
    
#     plt.clf()
#     ax  = plt.gca()
#     plt.xlim([xmin,xmax])
#     plt.ylim([3e-5,0.2])
#     if datadict["A"]:
#         plt.scatter(PA  ,LeA/fohmA**0.5  ,s=150,marker="^",c=np.log10(ColA),vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Oranges"),edgecolor='orange',label="Aubert et al (2009)")
# if datadict["Y"]:
#     plt.scatter(PY  ,LeY/fohmY**0.5  ,s=150,marker="o",c=np.log10(ColY),vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Reds")   ,edgecolor='red'   ,label="Yadav et al (2016)")
# if datadict["L"]:
#     plt.scatter(PC  ,LeC/fohmC**0.5  ,s=150,marker="*",c=np.log10(ColC),vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Blues")  ,edgecolor='blue'  ,label="Leeds")
# if datadict["UC"]:
#     plt.scatter(PUC ,LeUC /fohmUC**0.5, s=150,marker="v",c=np.log10(ColU),vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Purples"),edgecolor='purple',label="Christensen 0F")
# if datadict["UCt"]:
#     plt.scatter(PUCt,LeUCt/fohmUCt**0.5,s=150,marker="^",c=np.log10(ColUt),vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Purples"),edgecolor='purple',label="Christensen FF")
# plt.loglog(Pfit, fitall, color="black", lw=lw_fit)
# if "IMAC" in plt_extrap_scalings:
#     plt.loglog(Pfit, fitimac, color="dimgrey"  ,linestyle="--", lw=lw_fit, label="$m=2/5$ (IMAC)")
# if "IMACd" in plt_extrap_scalings:
#     plt.loglog(Pfit, fitimacd, color="dimgrey"  ,linestyle="-.", lw=lw_fit, label="$m=1/5$ (IMACd)")
# if "IMACi" in plt_extrap_scalings:
#     plt.loglog(Pfit, fitimaci, color="dimgrey"  ,linestyle=":", lw=lw_fit, label="$m=3/10$ (IMACi)")
# if "IMA" in plt_extrap_scalings:
#     if (calc_prefac_err):
#         plt.loglog(Pfit, fitEn, c=lc_fit[idxIMA], linestyle=ls_fit[idxIMA], lw=lw_fit,
#                    label="$m=1/3$ (QG-MAC), $\sigma=$"+str(np.round(c1_sd_En,3)))
#         plt.fill_between(Pfit, fitEn_1sdm, fitEn_1sdp, color=sc_fit[idxIMA], alpha=sc_alpha, zorder=-1)
#     else:
#         plt.loglog(Pfit, fitEn, c=lc_fit[idxIMA], linestyle=ls_fit[idxIMA], lw=lw_fit, label="$m=1/3$ (QG-MAC)")

# if "MAC" in plt_extrap_scalings:
#     if (calc_prefac_err):
#         plt.loglog(Pfit, fitmac, c=lc_fit[idxMAC], linestyle=ls_fit[idxMAC], lw=lw_fit,
#                    label="$m=1/4$ (MAC), $\sigma=$"+str(np.round(c1_sd_mac,3)))
#         plt.fill_between(Pfit, fitmac_1sdm, fitmac_1sdp, color=sc_fit[idxMAC], alpha=sc_alpha, zorder=-1)
#         plt.loglog(Pfit, fitmac_3sdm, c=lc_fit[idxMAC], linestyle=ls_fit[idxMAC], lw=1., zorder=-1)
#         plt.loglog(Pfit, fitmac_3sdp, c=lc_fit[idxMAC], linestyle=ls_fit[idxMAC], lw=1., zorder=-1)
#     else:
#         plt.loglog(Pfit, fitmac, c=lc_fit[idxMAC], linestyle=ls_fit[idxMAC], lw=lw_fit, label="$m=1/4$ (MAC)")


# ax.add_patch(Rectangle(xy=(P_earth_min,Le_earth_rmsmin) ,width=(P_earth_max-P_earth_min), height=(Le_earth_rmsmax-Le_earth_rmsmin),
#                        linewidth=1, color='black', fill=True))

# legend_xpos = 0.01 # x-position of legend (>1, outside of main plot)
# if datadict["A"]:
#     plt.text(legend_xpos, 0.95, "$m$  = "+str(np.round(mA,2)) +"$\pm$"+str(np.round(resA,2)) +", SSR="+str(np.round(ssrA[0],2)), transform=ax.transAxes, color='orange')
# if datadict["Y"]:
#     plt.text(legend_xpos, 0.90, "$m$  = "+str(np.round(mY,2)) +"$\pm$"+str(np.round(resY,2)) +", SSR="+str(np.round(ssrY[0],2)), transform=ax.transAxes, color='red')
# if datadict["L"]:
#     plt.text(legend_xpos, 0.85, "$m$  = "+str(np.round(mC,2)) +"$\pm$"+str(np.round(resC,2)) +", SSR="+str(np.round(ssrC[0],2)), transform=ax.transAxes, color='blue')
# if datadict["UC"]:
#     plt.text(legend_xpos, 0.80, "$m$v = "+str(np.round(mUC,2))+"$\pm$"+str(np.round(resUC,2))+", SSR="+str(np.round(ssrUC[0],2)), transform=ax.transAxes, color='purple')
# if datadict["UCt"]:
#     plt.text(legend_xpos, 0.75, "$m$^ = "+str(np.round(mUCt,2))+"$\pm$"+str(np.round(resUCt,2))+", SSR="+str(np.round(ssrUCt[0],2)), transform=ax.transAxes, color='purple')
# plt.text(legend_xpos, 0.70, "$m$  = "+str(np.round(mall,2))+"$\pm$"+str(np.round(resall,2))+", SSR="+str(np.round(ssrall[0],2)), transform=ax.transAxes, color='black')
# cbar = plt.colorbar()
# cbar.set_label("log $Re$")
# plt.xlabel('$p_A$')
# if myfohm == 0: 
#     plt.ylabel('$Le^{\\rm rms}_t/f_{ohm}^{1/2}$')
# elif myfohm == 1: 
#     plt.ylabel('$Le^{\\rm rms}_t$')
# ax.set_yscale('log')
# ax.set_xscale('log')
# plt.legend(bbox_to_anchor=(legend_xpos-0.1, 1.10), loc=3, ncol=2, borderaxespad=0)
# plt.rcParams["figure.figsize"] = [15,10]

# title_str = b.getPlotTitle(myfdip=myfdip,
#             myEr=myEr, Er_range=[EMoEK_range[0],EMoEK_range[1]])
# plt.title(title_str)

# file2 = "./fig/Lefohm_PA_Brmsextrap_fdip=" + fdipn + "_fohm=" + fohmn
# file3 = "./fig/Lefohm_PA_Brmsextrap_fdip=" + fdipn + "_fohm=" + fohmn
# if (myEkOPm == 1):
#     file2 += "_" + EkOPm_tag
#     file3 += "_" + EkOPm_tag
# if (myEr == 1):
#     file2 += "_" + EMoEK_tag
#     file3 += "_" + EMoEK_tag
# file2 += ".pdf"
# file3 += ".png"
# plt.savefig(file2, format='pdf',bbox_inches="tight")
# plt.savefig(file3, format='png',bbox_inches="tight")