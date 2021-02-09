import numpy as np
import pandas as pd
import bscaling_functions as b
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size':'12'})
rc('text', usetex=False)

# ---------------------------------------------------------------------------------------
# --- Input parameters ---
# ---------------------------------------------------------------------------------------

calc_prefac_err = True # Calculate and plot prefactor error?
myfdip  = 2 # Use 0 for all fdip values, 1 for fdip > 0.50, 2 for filtering with fdip=(0.35,0.80), 3 for fdip=(0.40,0.80) (see below).
myfohm  = 0 # Use 0 for fohm factor, or 1 for NO fohm factor
myEkOPm = 0 # Use 1 (0) for (not) filtering in a specified range of Ek/Pm values
myEr    = 1 # Use 1 (0) for (not) filtering in specified EM/EK range

if (myEkOPm==1):
    EkOPm_range = [1.e-10,1.e-4] # Range of Ek/Pm to include in the analysis (Ek=\nu/\Omega*D^2 as in Aubert's definition)
else:
    EkOPm_range = None
if (myEr==1):
    EMoEK_range = [2e0,1.e+19]
else:
    EMoEK_range = None
    
# List of scalings to be plotted in extrap figs ("IMA","MAC","IMAC","IMACd","IMACi")
# IMA (or Energy below) corresponds to QG-MAC in the new notation.
plt_extrap_scalings = ["IMA","MAC"]
lc_fit              = ["g",  "darkgrey"] # Line colour of fits
sc_fit              = ["lightgreen", "lightgrey"] # Shading colour for \sigma intervals
sc_alpha            = 0.7
ls_fit              = ["--", ":"] # Line style of fits
lw_fit = 2. # Line-width of fits in plots
chk    = 0 # Use 1 to print on screen checks of energy (quiet otherwise)

# Internal consistency check of Leeds sims. It compares the cmb (total and dipole) field strengths
# calculated from mag_cmb-files and Gauss coeffs-files.
check_gauss_Led = False
# Plot and analyse bdip?
plt_bdip = False

# ---------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------

# Create dictionaries of datasets
datadict = {"L":{  "plot":False, "dataset":"L", "d":{}, "rmsINT":{}, "rmsCMB":{}, "dipCMB":{}, "plotp":{}},
            "A":{  "plot":False, "dataset":"A", "d":{}, "rmsINT":{}, "rmsCMB":{}, "dipCMB":{}, "plotp":{}},
            "Y":{  "plot":False, "dataset":"Y", "d":{}, "rmsINT":{}, "rmsCMB":{}, "dipCMB":{}, "plotp":{}},
            "UC":{ "plot":False, "dataset":"UC", "d":{}, "rmsINT":{}, "rmsCMB":{}, "dipCMB":{}, "plotp":{}},
            "UCt":{"plot":False, "dataset":"UCt", "d":{}, "rmsINT":{}, "rmsCMB":{}, "dipCMB":{}, "plotp":{}},
            "APath":{"plot":False, "dataset":"APath", "d":{}, "rmsINT":{}, "rmsCMB":{}, "dipCMB":{}, "plotp":{}},
            "S":{"plot":True, "dataset":"S", "d":{}, "rmsINT":{}, "rmsCMB":{}, "dipCMB":{}, "plotp":{}}}

earthdict   = {"plot":True, "rmsINT":{}, "rmsCMB":{}, "dipCMB":{}, "p":{}}

# Input datasets
leedsname = "./all_LED_tave_NEW.csv"
yadavname = "./Yadav_etal_pnas_16_data.xls_all"
aubertname= "./aubert2009-all.txt"
christname= "./dynq_mycode.data"
christnamet="./dyntq_mycode.data"
APathname = "./A17_all"
Sname     = "./S19_data"

# fdip filtering
if myfdip == 0:
    fdipn = "0"
    fdip_range = [0.,1.1]
elif myfdip == 1:
    fdipn = "1"
    fdip_range = [0.50,1.1]
elif myfdip == 2:
    fdipn = "2"
    fdip_range = [0.35,0.80]
elif myfdip == 3:
    fdipn = "3"
    fdip_range = [0.40,0.80]
else:
    raise ValueError("Not valid value of myfdip provided.")
filetag = "_fdip"+fdipn

# filtering by Ek/Pm and/or Emag/Ekin
if (myEkOPm == 1 or myEr == 1):
    if (myEkOPm == 1 and myEr != 1):
        filetag += "_EkOPm"
    elif (myEkOPm != 1 and myEr == 1):
        filetag += "_EMoEK"
    elif (myEkOPm == 1 and myEr == 1):
        filetag += "_EkOPm_EMoEK"

if myfohm == 1:
    fohmn = "1"
else: 
    fohmn = "0"

if (myEkOPm == 1):
    EkOPm_min = EkOPm_range[0]; EkOPm_max = EkOPm_range[1]
    EkOPm_tag = "EkOPm="+format(EkOPm_min,'.2e')+"_"+format(EkOPm_max,'.2e')
if (myEr == 1):
    EMoEK_min = EMoEK_range[0]; EMoEK_max = EMoEK_range[1]
    EMoEK_tag = "EMoEK="+format(EMoEK_min,'.1e')+"_"+format(EMoEK_max,'.1e')

# read and filter datasets
if datadict["L"]["plot"]:
    df, datadict = b.filter_table(infname=leedsname  ,outfname=leedsname+filetag, dataset="Leeds", fdip_range=fdip_range, 
                                  EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict,
                                  chk=chk, myfohm=myfohm)
if datadict["Y"]["plot"]:
    df, datadict = b.filter_table(infname=yadavname  ,outfname=yadavname+filetag, dataset="Yadav", fdip_range=fdip_range, 
                                  EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict,
                                  chk=chk, myfohm=myfohm)
if datadict["A"]["plot"]:
    df, datadict = b.filter_table(infname=aubertname ,outfname=aubertname+filetag, dataset="Aubert", fdip_range=fdip_range, 
                                  EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict,
                                  chk=chk, myfohm=myfohm)
if datadict["UC"]["plot"]:
    df, datadict = b.filter_table(infname=christname ,outfname=christname+filetag, dataset="Christensen", fdip_range=fdip_range, 
                                  EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict,
                                  chk=chk, myfohm=myfohm)
if datadict["UCt"]["plot"]:
    df, datadict = b.filter_table(infname=christnamet,outfname= christnamet+filetag, dataset="ChristensenT", fdip_range=fdip_range, 
                                  EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict,
                                  chk=chk, myfohm=myfohm)
if datadict["APath"]["plot"]:
    df, datadict = b.filter_table(infname=APathname,  outfname=APathname+filetag, dataset="APath",fdip_range=fdip_range, 
                                  EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict,
                                  chk=chk, myfohm=myfohm)
if datadict["S"]["plot"]:
    df, datadict = b.filter_table(infname=Sname    ,  outfname=Sname+filetag    , dataset="S"    ,fdip_range=fdip_range, 
                                  EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict,
                                  chk=chk, myfohm=myfohm)

# Output file of fitted slopes and pre-factors
outfpath = "./"
outfname = "scalingSlope_fdip="+fdipn+"_fohm="+fohmn
outfnamepf = "scalingPrefac_fdip="+fdipn+"_fohm="+fohmn
if (myEkOPm == 1):
    outfname += "_"+EkOPm_tag
    outfnamepf += "_"+EkOPm_tag
if (myEr == 1):
    outfname += "_"+EMoEK_tag
    outfnamepf += "_"+EMoEK_tag
outfname += ".txt"
outfnamepf += ".txt"
# ---

# Construct datasets of all simulations
alldatadict = b.combineDataDict(datadict, quiet=False)

# Fit all data
alldatadict["rmsINT"]["ssr"], alldatadict["rmsINT"]["m"], alldatadict["rmsINT"]["c"], alldatadict["rmsINT"]["res"] = \
    b.fits(alldatadict["p"], alldatadict["rmsINT"]["Le"], alldatadict["fohm"])
alldatadict["rmsCMB"]["ssr"], alldatadict["rmsCMB"]["m"], alldatadict["rmsCMB"]["c"], alldatadict["rmsCMB"]["res"] = \
    b.fits(alldatadict["p"], alldatadict["rmsCMB"]["Le"], alldatadict["fohm"])
alldatadict["dipCMB"]["ssr"], alldatadict["dipCMB"]["m"], alldatadict["dipCMB"]["c"], alldatadict["dipCMB"]["res"] = \
    b.fits(alldatadict["p"], alldatadict["dipCMB"]["Le"], alldatadict["fohm"])

# Construct arrays for plotting
xmin = 1e-14
xmax = 1e-3
Pfit = np.linspace(xmin, xmax, 10)
fitall    = 10**alldatadict["rmsINT"]["c"] * Pfit**alldatadict["rmsINT"]["m"] 
fitrmscmb = 10**alldatadict["rmsCMB"]["c"] * Pfit**alldatadict["rmsCMB"]["m"]
fitdipcmb = 10**alldatadict["dipCMB"]["c"] * Pfit**alldatadict["dipCMB"]["m"]

# CD - REMOVE PRINTS?
print('All     (slope, 10^c, std dev c, SSR) =', alldatadict["rmsINT"]["m"], 10**alldatadict["rmsINT"]["c"], \
    alldatadict["rmsINT"]["res"], alldatadict["rmsINT"]["ssr"])
print('RMS CMB (slope, 10^c, std dev c, SSR) =', alldatadict["rmsCMB"]["m"], 10**alldatadict["rmsCMB"]["c"], \
    alldatadict["rmsCMB"]["res"], alldatadict["rmsCMB"]["ssr"])
print('DIP CMB (slope, 10^c, std dev c, SSR) =', alldatadict["dipCMB"]["m"], 10**alldatadict["dipCMB"]["c"], \
    alldatadict["dipCMB"]["res"], alldatadict["dipCMB"]["ssr"])

if calc_prefac_err:
    # - Error estimate on the prefactor for CMB dip field strength
    c_err_dipcmb = b.prefacError(alldatadict["p"], alldatadict["dipCMB"]["Le"]/alldatadict["fohm"]**0.5,
                                 model=[10**alldatadict["dipCMB"]["c"],alldatadict["dipCMB"]["m"]], plot=False, quiet=False)
    alldatadict["dipCMB"]["c_sd"] = c_err_dipcmb[1]
    fitdipcmb_1sdp, fitdipcmb_1sdm, fitdipcmb_3sdp, fitdipcmb_3sdm = b.getFitBounds(Pfit, alldatadict["dipCMB"]["c"],\
                                                                     alldatadict["dipCMB"]["c_sd"], alldatadict["dipCMB"]["m"])

# - Define symbols' properties for plotting
Cmax = np.log10(1000.0)#np.log10(np.max(Eall))
Cmin = np.log10(50.0)#np.log10(np.min(Eall))
datadict = b.getPlotProperties(datadict)

# - Plot bdip stuff
if plt_bdip:
    b.plot_bdip(datadict, myfdip)

# --- Store fitted values in output files
b.saveFitValues(filename=outfpath+outfname, datadict=datadict, alldatadict=alldatadict)

# --- Calculate fits and prefac errors for each scaling on all field strengths
alldatadict = b.fitForceScalings(alldatadict, quiet=False)

# --- save out prefactor file
b.savePrefacValues(filename=outfpath+outfnamepf, indict=alldatadict, l_prefac_err=calc_prefac_err)

earthdict["dipCMB"]["min"], earthdict["dipCMB"]["max"], earthdict["rmsINT"]["min"], earthdict["rmsINT"]["max"],\
    earthdict["rmsCMB"]["min"], earthdict["rmsCMB"]["max"], earthdict["p"]["min"], earthdict["p"]["max"] = b.getEarthEstimates(quiet=False)

##########################
# Brms INTERNAL FIELD
##########################

fitEn    = 10**alldatadict["rmsINT"]["cIMA"]   * Pfit**alldatadict["rmsINT"]["mIMA"]
fitmac   = 10**alldatadict["rmsINT"]["cMAC"]   * Pfit**alldatadict["rmsINT"]["mMAC"]
fitimac  = 10**alldatadict["rmsINT"]["cIMAC"]  * Pfit**alldatadict["rmsINT"]["mIMAC"]
fitimacd = 10**alldatadict["rmsINT"]["cIMACd"] * Pfit**alldatadict["rmsINT"]["mIMACd"]
fitimaci = 10**alldatadict["rmsINT"]["cIMACi"] * Pfit**alldatadict["rmsINT"]["mIMACi"]

# fit scaling bounds
if calc_prefac_err:
    fitEn_1sdp, fitEn_1sdm, fitEn_3sdp, fitEn_3sdm = b.getFitBounds(Pfit, alldatadict["rmsINT"]["cIMA"],\
                                                                    alldatadict["rmsINT"]["cIMA_sd"], alldatadict["rmsINT"]["mIMA"])
    fitmac_1sdp, fitmac_1sdm, fitmac_3sdp, fitmac_3sdm = b.getFitBounds(Pfit, alldatadict["rmsINT"]["cMAC"],\
                                                                         alldatadict["rmsINT"]["cMAC_sd"], alldatadict["rmsINT"]["mMAC"])
# Get plot properties
idxIMA = b.idxStr(plt_extrap_scalings, "IMA")[0]
idxMAC = b.idxStr(plt_extrap_scalings, "MAC")[0]

# Plot simulations
ax, legend_xpos, legend_ypos = b.plotSimulations(datadict=datadict, alldatadict=alldatadict, earthdict=earthdict, field="rmsINT",
                                                 cbarrange=[Cmin,Cmax], xrange=[xmin,xmax], yrange=[3e-5,0.2])

plt.loglog(Pfit ,fitall, color="black")
if "IMAC" in plt_extrap_scalings:
    plt.loglog(Pfit, fitimac, color="dimgrey"  ,linestyle="--", lw=lw_fit, label="$m=2/5$ (IMAC)")
if "IMACd" in plt_extrap_scalings:
    plt.loglog(Pfit, fitimacd, color="dimgrey" ,linestyle="-.", lw=lw_fit, label="$m=1/5$ (IMACd)")
if "IMACi" in plt_extrap_scalings:
    plt.loglog(Pfit, fitimaci, color="dimgrey" ,linestyle=":", lw=lw_fit, label="$m=3/10$ (IMACi)")
if "IMA" in plt_extrap_scalings:
    if (calc_prefac_err):
        plt.loglog(Pfit, fitEn, c=lc_fit[idxIMA], linestyle=ls_fit[idxIMA], lw=lw_fit,
                   label="$m=1/3$ (QG-MAC), $\sigma=$"+str(np.round(alldatadict["rmsINT"]["cIMA_sd"],3)))
        plt.fill_between(Pfit, fitEn_1sdm, fitEn_1sdp, color=sc_fit[idxIMA], alpha=sc_alpha, zorder=-1)
    else:
        plt.loglog(Pfit, fitEn, c=lc_fit[idxIMA], linestyle=ls_fit[idxIMA], lw=lw_fit, label="$m=1/3$ (QG-MAC)")
if "MAC" in plt_extrap_scalings:
    if (calc_prefac_err):
        plt.loglog(Pfit, fitmac, c=lc_fit[idxMAC], linestyle=ls_fit[idxMAC], lw=lw_fit,
                   label="$m=1/4$ (MAC), $\sigma=$"+str(np.round(alldatadict["rmsINT"]["cMAC_sd"],3)))
        plt.fill_between(Pfit, fitmac_1sdm, fitmac_1sdp, color=sc_fit[idxMAC], alpha=sc_alpha, zorder=-1)
        plt.loglog(Pfit, fitmac_3sdm, c=lc_fit[idxMAC], linestyle=ls_fit[idxMAC], lw=1., zorder=-1)
        plt.loglog(Pfit, fitmac_3sdp, c=lc_fit[idxMAC], linestyle=ls_fit[idxMAC], lw=1., zorder=-1)
    else:
        plt.loglog(Pfit, fitmac, c=lc_fit[idxMAC], linestyle=ls_fit[idxMAC], lw=lw_fit, label="$m=1/4$ (MAC)")
# Plot title and axes labels
plt.xlabel('$p_A$')
if myfohm == 0: 
    plt.ylabel('$Le^{\\rm rms}_t/f_{ohm}^{1/2}$')
elif myfohm == 1: 
    plt.ylabel('$Le^{\\rm rms}_t$')

title_str = b.getPlotTitle(myfdip=myfdip, myEr=myEr, Er_range=[EMoEK_range[0],EMoEK_range[1]])
plt.title(title_str)

# save figure
file2 = "./fig/Lefohm_PA_Brmsextrap_fdip=" + fdipn + "_fohm=" + fohmn
file3 = "./fig/Lefohm_PA_Brmsextrap_fdip=" + fdipn + "_fohm=" + fohmn
if (myEkOPm == 1):
    file2 += "_" + EkOPm_tag
    file3 += "_" + EkOPm_tag
if (myEr == 1):
    file2 += "_" + EMoEK_tag
    file3 += "_" + EMoEK_tag
file2 += ".pdf"
file3 += ".png"
plt.savefig(file2, format='pdf',bbox_inches="tight")
plt.savefig(file3, format='png',bbox_inches="tight")

##########################
# Brms CMB FIELD
##########################

# Calculate prefactors of Brms cmb
fitEn    = 10**alldatadict["rmsCMB"]["cIMA"]   * Pfit**alldatadict["rmsCMB"]["mIMA"]
fitmac   = 10**alldatadict["rmsCMB"]["cMAC"]   * Pfit**alldatadict["rmsCMB"]["mMAC"]
fitimac  = 10**alldatadict["rmsCMB"]["cIMAC"]  * Pfit**alldatadict["rmsCMB"]["mIMAC"]
fitimacd = 10**alldatadict["rmsCMB"]["cIMACd"] * Pfit**alldatadict["rmsCMB"]["mIMACd"]
fitimaci = 10**alldatadict["rmsCMB"]["cIMACi"] * Pfit**alldatadict["rmsCMB"]["mIMACi"]

# Plot simulations
ax, legend_xpos, legend_ypos = b.plotSimulations(datadict=datadict, alldatadict=alldatadict, earthdict=earthdict, field="rmsCMB",
                                             cbarrange=[Cmin,Cmax], xrange=[xmin,xmax], yrange=[5e-6,0.2])
plt.loglog(Pfit ,fitrmscmb, color="black")
if "IMAC" in plt_extrap_scalings:
    plt.loglog(Pfit, fitimac,   color="dimgrey"  , linestyle="--", label="$m=2/5$ (IMAC)")
if "IMACd" in plt_extrap_scalings:
    plt.loglog(Pfit, fitimacd,  color="dimgrey"  , linestyle="-.", label="$m=1/5$ (IMACd)")
if "IMACi" in plt_extrap_scalings:
    plt.loglog(Pfit, fitimaci,  color="dimgrey"  , linestyle=":",  label="$m=3/10$ (IMACi)")
if "IMA" in plt_extrap_scalings:
    plt.loglog(Pfit, fitEn, color="darkgrey" , linestyle="--", label="$m=1/3$ (IMA)")
if "MAC" in plt_extrap_scalings:
    plt.loglog(Pfit, fitmac,    color="lightgray", linestyle="--", label="$m=1/4$ (MAC)")

plt.xlabel('$p_A$')
if myfohm == 0: 
    plt.ylabel('$Le^{\\rm rms}_{\\rm cmb}/f_{ohm}^{1/2}$')
elif myfohm == 1: 
    plt.ylabel('$Le^{\\rm rms}_{\\rm cmb}$')

title_str = b.getPlotTitle(myfdip=myfdip, myEr=myEr, Er_range=[EMoEK_range[0],EMoEK_range[1]])
plt.title(title_str)

file2 = "./fig/Lefohm_PA_Brmscmbextrap_fdip=" + fdipn + "_fohm=" + fohmn
file3 = "./fig/Lefohm_PA_Brmscmbextrap_fdip=" + fdipn + "_fohm=" + fohmn
if (myEkOPm == 1):
    file2 += "_" + EkOPm_tag
    file3 += "_" + EkOPm_tag
if (myEr == 1):
    file2 += "_" + EMoEK_tag
    file3 += "_" + EMoEK_tag
file2 += ".pdf"
file3 += ".png"
plt.savefig(file2,format='pdf',bbox_inches="tight")
plt.savefig(file3,format='png',bbox_inches="tight")

##########################
# Bdip cmb plots
##########################

fitEn    = 10**alldatadict["dipCMB"]["cIMA"]   * Pfit**alldatadict["dipCMB"]["mIMA"]
fitmac   = 10**alldatadict["dipCMB"]["cMAC"]   * Pfit**alldatadict["dipCMB"]["mMAC"]
fitimac  = 10**alldatadict["dipCMB"]["cIMAC"]  * Pfit**alldatadict["dipCMB"]["mIMAC"]
fitimacd = 10**alldatadict["dipCMB"]["cIMACd"] * Pfit**alldatadict["dipCMB"]["mIMACd"]
fitimaci = 10**alldatadict["dipCMB"]["cIMACi"] * Pfit**alldatadict["dipCMB"]["mIMACi"]

# Error estimate on the prefactor
if calc_prefac_err:
    fitEn_1sdp, fitEn_1sdm, fitEn_3sdp, fitEn_3sdm = b.getFitBounds(Pfit, alldatadict["dipCMB"]["cIMA"],\
                                                                    alldatadict["dipCMB"]["cIMA_sd"], alldatadict["dipCMB"]["mIMA"])
    fitmac_1sdp, fitmac_1sdm, fitmac_3sdp, fitmac_3sdm = b.getFitBounds(Pfit, alldatadict["dipCMB"]["cMAC"],\
                                                                         alldatadict["dipCMB"]["cMAC_sd"], alldatadict["dipCMB"]["mMAC"])

# Plot simulations
ax, legend_xpos, legend_ypos = b.plotSimulations(datadict=datadict, alldatadict=alldatadict, earthdict=earthdict,
                                                 field="dipCMB",
                                                 cbarrange=[Cmin,Cmax], xrange=[xmin,xmax], yrange=[5e-6,0.2])
plt.loglog(Pfit, fitdipcmb, color="black", lw=lw_fit)
#           label="$m =$"+str(np.round(mdipcmb ,2))+"$\pm$"+str(np.round(resdipcmb ,2))+", SSR="+str(np.round(ssrdipcmb,2)))
#plt.loglog(Pfit, fitdipcmb_3sdp, c="k", ls=':', lw=1.)
#plt.loglog(Pfit, fitdipcmb_3sdm, c="k", ls=':', lw=1.)
if "IMAC" in plt_extrap_scalings:
    plt.loglog(Pfit, fitimac , color="dimgrey", linestyle="--", lw=lw_fit, label="$m=2/5$ (IMAC)")
if "IMACd" in plt_extrap_scalings:
    plt.loglog(Pfit, fitimacd, color="dimgrey", linestyle="-.", lw=lw_fit, label="$m=1/5$ (IMACd)")
if "IMACi" in plt_extrap_scalings:
    plt.loglog(Pfit, fitimaci, color="dimgrey", linestyle=":", lw=lw_fit, label="$m=3/10$ (IMACi)")
if "IMA" in plt_extrap_scalings:
    if (calc_prefac_err):
        plt.loglog(Pfit, fitEn,      c=lc_fit[idxIMA], linestyle=ls_fit[idxIMA], lw=lw_fit,
                   label="$m=1/3$ (QG-MAC), $\sigma=$"+str(np.round(alldatadict["dipCMB"]["cIMA_sd"],3)))
        plt.fill_between(Pfit, fitEn_1sdm, fitEn_1sdp, color=sc_fit[idxIMA], alpha=sc_alpha, zorder=-1)
    else:
        plt.loglog(Pfit, fitEn, color=lc_fit[idxIMA], linestyle=ls_fit[idxIMA], lw=lw_fit,
                   label="$m=1/3$ (QG-MAC)")
if "MAC" in plt_extrap_scalings:
    if (calc_prefac_err):
        plt.loglog(Pfit, fitmac, c=lc_fit[idxMAC], linestyle=ls_fit[idxMAC], lw=lw_fit,
                   label="$m=1/4$ (MAC), $\sigma=$"+str(np.round(alldatadict["dipCMB"]["cMAC_sd"],3)))
        plt.fill_between(Pfit, fitmac_1sdm, fitmac_1sdp, color=sc_fit[idxMAC], alpha=sc_alpha, zorder=-1)
        plt.loglog(Pfit, fitmac_3sdm, c=lc_fit[idxMAC], linestyle=ls_fit[idxMAC], lw=1., zorder=-1)
        plt.loglog(Pfit, fitmac_3sdp, c=lc_fit[idxMAC], linestyle=ls_fit[idxMAC], lw=1., zorder=-1)
    else:
        plt.loglog(Pfit, fitmac, color=lc_fit[idxMAC],linestyle=ls_fit[idxMAC], lw=lw_fit, label="$m=1/4$ (MAC)")

plt.xlabel('$p_A$')
if myfohm == 0: 
    plt.ylabel('$Le_{\\rm{cmb}}^{l=1}/f_{ohm}^{1/2}$')
elif myfohm ==1: 
    plt.ylabel('$Le_{\\rm{cmb}}^{l=1}$')
ax.set_yscale('log')
ax.set_xscale('log')
plt.legend(bbox_to_anchor=(legend_xpos-0.1, 1.1), loc=3, ncol=2, borderaxespad=0)
plt.rcParams["figure.figsize"] = [15,10]

title_str = b.getPlotTitle(myfdip=myfdip, myEr=myEr, Er_range=[EMoEK_range[0],EMoEK_range[1]])
plt.title(title_str)

file2 = "./fig/Lefohm_PA_Bdipcmbextrap_fdip=" + fdipn + "_fohm=" + fohmn
file3 = "./fig/Lefohm_PA_Bdipcmbextrap_fdip=" + fdipn + "_fohm=" + fohmn
if (myEkOPm == 1):
    file2 += "_" + EkOPm_tag
    file3 += "_" + EkOPm_tag
if (myEr == 1):
    file2 += "_" + EMoEK_tag
    file3 += "_" + EMoEK_tag
file2 += ".pdf"
file3 += ".png"
plt.savefig(file2, format='pdf',bbox_inches="tight")
plt.savefig(file3, format='png',bbox_inches="tight")

##########################################################
# Zoomed figs plotted using subfigure
##########################################################

# - B rms core
plt.clf()
plt.figure(figsize=(16,15))
plt.subplot(3,1,1)
ax  = plt.gca()
plt.xlim([1e-10,1e-3])
plt.ylim([1e-4,0.2])

for key in datadict:
    if datadict[key]["plot"]:
        plt.scatter(datadict[key]["p"], datadict[key]["rmsINT"]["Le"]/datadict[key]["fohm"]**0.5,
                        s=datadict[key]["plotp"]["size"], marker=datadict[key]["plotp"]["marker"],
                        c=np.log10(datadict[key]["plotp"]["Col"]), vmin=Cmin, vmax=Cmax,
                        cmap=plt.get_cmap(datadict[key]["plotp"]["cmap"]), edgecolor=datadict[key]["plotp"]["edgecolor"], label=datadict[key]["plotp"]["label"]) 
plt.loglog(Pfit , fitall, color="black")
plt.loglog(5e-12 , 2e-4, marker="s", markersize=30)

legend_xpos = 0.05
legend_ypos = 0.95; legend_dy   = 0.06
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
                 slope_str+str(np.round(datadict[key]["rmsINT"]["m"],2)) +\
                 "$\pm$"+str(np.round(datadict[key]["rmsINT"]["res"],2)) +\
                 ", SSR="+str(np.round(datadict[key]["rmsINT"]["ssr"],2)),
                 transform=ax.transAxes, color=datadict[key]["plotp"]["edgecolor"])
        iplt += 1
plt.text(legend_xpos, legend_ypos-iplt*legend_dy, "$m$  = "+str(np.round(alldatadict["rmsINT"]["m"],2))+\
                                                  "$\pm$"+str(np.round(alldatadict["rmsINT"]["res"],2))+\
                                                  ", SSR="+str(np.round(alldatadict["rmsINT"]["ssr"],2)),
                                                  transform=ax.transAxes, color='black')
if myfohm == 0: 
    plt.ylabel('$Le_t^{rms}/f_{ohm}^{1/2}$')
elif myfohm == 1: 
    plt.ylabel('$Le_t^{rms}$')
ax.set_yscale('log')
ax.set_xscale('log')
if myfdip == 0:
    plt.title("Field strength ($Le$) vs Buoyant Power ($P$) for all models")
elif myfdip == 1:
    plt.title("Field strength ($Le$) vs Buoyant Power ($P$) for $f_{dip}>0.5$")
elif myfdip == 2:
    plt.title("Field strength ($Le$) vs Buoyant Power ($P$) for $%.2f<f_{dip}<%.2f$" %(fdip_range[0],fdip_range[1]))
elif myfdip == 3:
    plt.title("Field strength ($Le$) vs Buoyant Power ($P$) for $%.2f<f_{dip}<%.2f$" %(fdip_range[0],fdip_range[1]))
plt.legend(bbox_to_anchor=(0.73, 0.05), loc=3, borderaxespad=0)
#file2="./fig/Lefohm_PA_Brmsallzoom_fdip=" + fdipn + "_fohm=" + fohmn + ".pdf"
#plt.savefig(file2, format='pdf',bbox_inches="tight")

# Brms - CMB
plt.subplot(3,1,2)
ax  = plt.gca()
plt.xlim([1e-10,1e-3])
plt.ylim([1e-4,0.2])

for key in datadict:
    if datadict[key]["plot"]:
        plt.scatter(datadict[key]["p"], datadict[key]["rmsCMB"]["Le"]/datadict[key]["fohm"]**0.5,
                        s=datadict[key]["plotp"]["size"], marker=datadict[key]["plotp"]["marker"],
                        c=np.log10(datadict[key]["plotp"]["Col"]), vmin=Cmin, vmax=Cmax,
                        cmap=plt.get_cmap(datadict[key]["plotp"]["cmap"]), edgecolor=datadict[key]["plotp"]["edgecolor"], label=datadict[key]["plotp"]["label"]) 
plt.loglog(Pfit,fitrmscmb,color="black")
plt.loglog(5e-12 , 2e-4, marker="s", markersize=30)

legend_ypos = 0.95; legend_dy   = 0.06
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
                 slope_str+str(np.round(datadict[key]["rmsCMB"]["m"],2)) +\
                 "$\pm$"+str(np.round(datadict[key]["rmsCMB"]["res"],2)) +\
                 ", SSR="+str(np.round(datadict[key]["rmsCMB"]["ssr"],2)),
                 transform=ax.transAxes, color=datadict[key]["plotp"]["edgecolor"])
        iplt += 1
plt.text(legend_xpos, legend_ypos-iplt*legend_dy, "$m$  = "+str(np.round(alldatadict["rmsCMB"]["m"],2))+\
                                                  "$\pm$"+str(np.round(alldatadict["rmsCMB"]["res"],2))+\
                                                  ", SSR="+str(np.round(alldatadict["rmsCMB"]["ssr"],2)), transform=ax.transAxes, color='black')
if myfohm == 0: 
    plt.ylabel('$Le_{cmb}^{rms}/f_{ohm}^{1/2}$')
elif myfohm == 1: 
    plt.ylabel('$Le_{cmb}^{rms}/$')
ax.set_yscale('log')
ax.set_xscale('log')
#plt.legend(bbox_to_anchor=(0.05, 0.70), loc=3, borderaxespad=0)
#file2="./fig/Lefohm_PA_Brmscmb_fdip=" + fdipn + "_fohm=" + fohmn + ".pdf"
#plt.savefig(file2, format='pdf',bbox_inches="tight")

# Bdip CMB
plt.subplot(3,1,3)
ax  = plt.gca()
plt.xlim([1e-10,1e-3])
plt.ylim([1e-4,0.2])

for key in datadict:
    if datadict[key]["plot"]:
        plt.scatter(datadict[key]["p"], datadict[key]["dipCMB"]["Le"]/datadict[key]["fohm"]**0.5,
                        s=datadict[key]["plotp"]["size"], marker=datadict[key]["plotp"]["marker"],
                        c=np.log10(datadict[key]["plotp"]["Col"]), vmin=Cmin, vmax=Cmax,
                        cmap=plt.get_cmap(datadict[key]["plotp"]["cmap"]), edgecolor=datadict[key]["plotp"]["edgecolor"], label=datadict[key]["plotp"]["label"]) 
plt.loglog(Pfit,fitdipcmb,color="black")
if (calc_prefac_err):
    plt.loglog(Pfit,fitdipcmb_1sdp,c='k',ls=':',lw=1.)
    plt.loglog(Pfit,fitdipcmb_1sdm,c='k',ls=':',lw=1.)
    #plt.loglog(Pfit,fitdipcmb_3sdp,c='k',ls='--',lw=1.)
    #plt.loglog(Pfit,fitdipcmb_3sdm,c='k',ls='--',lw=1.)

legend_ypos = 0.95; legend_dy   = 0.06
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
                 slope_str+str(np.round(datadict[key]["dipCMB"]["m"],2))+\
                 "$\pm$"+str(np.round(datadict[key]["dipCMB"]["res"],2))+\
                 ", SSR="+str(np.round(datadict[key]["dipCMB"]["ssr"],2)),
                 transform=ax.transAxes, color=datadict[key]["plotp"]["edgecolor"])
        iplt += 1
plt.text(legend_xpos, legend_ypos-iplt*legend_dy, "$m$ = "+str(np.round(alldatadict["dipCMB"]["m"],2))+\
                              "$\pm$"+str(np.round(alldatadict["dipCMB"]["res"],2))+\
                              ", SSR="+str(np.round(alldatadict["dipCMB"]["ssr"],2)), transform=ax.transAxes, color='black')
plt.text(legend_xpos, legend_ypos-(iplt+1)*legend_dy,"$\sigma$ = "+str(np.round(alldatadict["dipCMB"]["c_sd"],4)), transform=ax.transAxes, color='black')
plt.xlabel('$P_A$')
if myfohm == 0: 
    plt.ylabel('$Le_{cmb}^{l=1}/f_{ohm}^{1/2}$')
elif myfohm == 1: 
    plt.ylabel('$Le_{cmb}^{l=1}$')

file2="./fig/Lefohm_PA_Bdipcmb_fdip=" + fdipn + "_fohm=" + fohmn
file3="./fig/Lefohm_PA_Bdipcmb_fdip=" + fdipn + "_fohm=" + fohmn
if (myEkOPm == 1):
    file2 += "_" + EkOPm_tag
    file3 += "_" + EkOPm_tag
if (myEr == 1):
    file2 += "_" + EMoEK_tag
    file3 += "_" + EMoEK_tag
file2 += ".pdf"
file3 += ".png"
plt.savefig(file2, format='pdf',bbox_inches="tight")
plt.savefig(file3, format='png',bbox_inches="tight")
