import numpy as np
import os
import bscaling_functions as b
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size':'18'})
rc('text', usetex=False)

# ---------------------------------------------------------------------------------------
# --- Input parameters ---
# ---------------------------------------------------------------------------------------

calc_prefac_err = True # Calculate and plot prefactor error?
myfdip  = 0 # Use 0 for all fdip values, 1 for fdip > 0.50, 2 for filtering with fdip=(0.35,0.80), 3 for fdip=(0.40,0.80) (see below).
myfohm  = 0 # Use 0 for fohm factor, or 1 for NO fohm factor
myEkOPm = 0 # Use 1 (0) for (not) filtering in a specified range of Ek/Pm values
myEr    = 1 # Use 1 (0) for (not) filtering in specified EM/EK range

if (myEkOPm==1):
    EkOPm_range = [1.e-10,1.e-4] # Range of Ek/Pm to include in the analysis (Ek=\nu/\Omega*D^2 as in Aubert's definition)
else:
    EkOPm_range = None
if (myEr==1):
    EMoEK_range = [2.,1.e+19]
else:
    EMoEK_range = None
    
# -- List of scalings to be plotted in extrap figs ("IMA","MAC","IMAC","IMACd","IMACi")
#    IMA (or Energy below) corresponds to QG-MAC in the new notation.
plt_extrap_scalings = ["IMA", "MAC"]
# = [""]
lc_fit              = ["g",  "g"]                 # Line colour of fits
sc_fit              = ["lightgreen", "lightgrey"] # Shading colour for \sigma intervals
sc_alpha            = 0.5
ls_fit              = ["--", ":"]                 # Line style of fits
lw_fit = 2.                                       # Line-width of fits in plots
chk    = 0 # Use 1 to print on screen checks of energy (quiet otherwise)

# -- Internal consistency check of Leeds sims. It compares the cmb (total and dipole) field strengths
#    calculated from mag_cmb-files and Gauss coeffs-files.
check_gauss_Led = False
# -- Plot and analyse bdip?
plt_bdip = True

# -- Categorise simulations by driving? (plots are by authors otherwise)
#    Categories are: FTFT, FFFF, FF0F, FTFF, Mixed, CE
#    Set all plots to true if categorise = True. 
categorise = True
plt_categ  = ["FTFT", "FF0F", "FTFF", "Mixed", "CE"]
#plt_categ  = ["FTFT", "FF0F", "FTFF", "Mixed"]
# -- write out file of databases for check?
write_check = True
# -- To compare unfiltered and filtered datasets it is easiest to set all filters above to 0
#    and change the name here; then set the required filter, change name here and rerun. 
outf_check  = './data/Dataset_UNFILTERED'

# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------

# Create dictionaries of datasets
datadict = {"L":{    "plot":True, "dataset":"L",     "d":{}, "rmsINT":{}, "rmsCMB":{}, "dipCMB":{}, "plotp":{}},
            "A":{    "plot":True, "dataset":"A",     "d":{}, "rmsINT":{}, "rmsCMB":{}, "dipCMB":{}, "plotp":{}},
            "Y":{    "plot":True, "dataset":"Y",     "d":{}, "rmsINT":{}, "rmsCMB":{}, "dipCMB":{}, "plotp":{}},
            "UC":{   "plot":True, "dataset":"UC",    "d":{}, "rmsINT":{}, "rmsCMB":{}, "dipCMB":{}, "plotp":{}},
            "UCt":{  "plot":True, "dataset":"UCt",   "d":{}, "rmsINT":{}, "rmsCMB":{}, "dipCMB":{}, "plotp":{}},
            "APath":{"plot":True, "dataset":"APath", "d":{}, "rmsINT":{}, "rmsCMB":{}, "dipCMB":{}, "plotp":{}},
            "S":{    "plot":True, "dataset":"S",     "d":{}, "rmsINT":{}, "rmsCMB":{}, "dipCMB":{}, "plotp":{}}}

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
    df, datadict = b.filter_table(infname=leedsname  ,outfname="./data/"+leedsname+filetag, dataset="Leeds", fdip_range=fdip_range, 
                                  EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict,
                                  categorise=categorise, chk=chk, myfohm=myfohm)
if datadict["Y"]["plot"]:
    df, datadict = b.filter_table(infname=yadavname  ,outfname="./data/"+yadavname+filetag, dataset="Yadav", fdip_range=fdip_range, 
                                  EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict,
                                  categorise=categorise, chk=chk, myfohm=myfohm)
if datadict["A"]["plot"]:
    df, datadict = b.filter_table(infname=aubertname ,outfname="./data/"+aubertname+filetag, dataset="Aubert", fdip_range=fdip_range, 
                                  EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict,
                                  categorise=categorise, chk=chk, myfohm=myfohm)
if datadict["UC"]["plot"]:
    df, datadict = b.filter_table(infname=christname ,outfname="./data/"+christname+filetag, dataset="Christensen", fdip_range=fdip_range, 
                                  EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict,
                                  categorise=categorise, chk=chk, myfohm=myfohm)
if datadict["UCt"]["plot"]:
    df, datadict = b.filter_table(infname=christnamet,outfname= "./data/"+christnamet+filetag, dataset="ChristensenT", fdip_range=fdip_range, 
                                  EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict,
                                  categorise=categorise, chk=chk, myfohm=myfohm)
if datadict["APath"]["plot"]:
    df, datadict = b.filter_table(infname=APathname,  outfname="./data/"+APathname+filetag, dataset="APath",fdip_range=fdip_range, 
                                  EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict,
                                  categorise=categorise, chk=chk, myfohm=myfohm)
if datadict["S"]["plot"]:
    df, datadict = b.filter_table(infname=Sname    ,  outfname="./data/"+Sname+filetag    , dataset="S"    ,fdip_range=fdip_range, 
                                  EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict,
                                  categorise=categorise, chk=chk, myfohm=myfohm)

# Categorise simulations by driving and select which to plot
if categorise:
    datadict = b.redefineDataDict(datadict, quiet=True)
    if ("FTFT" in plt_categ) : datadict["FTFT"]["plot"] = True 
    if ("FFFF" in plt_categ) : datadict["FFFF"]["plot"] = True 
    if ("FF0F" in plt_categ) : datadict["FF0F"]["plot"] = True 
    if ("FTFF" in plt_categ) : datadict["FTFF"]["plot"] = True 
    if ("Mixed" in plt_categ): datadict["Mixed"]["plot"] = True 
    if ("CE" in plt_categ)   : datadict["CE"]["plot"] = True

if write_check:
    if (not os.path.isdir('./data')):
        os.mkdir('./data')
    # write filtered datasets in output file:
    b.writefilecheck(datadict, outfiletag=outf_check)

# check by saving out filtered datasets
# ----------

# get output file name of fitted slopes and pre-factors
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

if calc_prefac_err:
    # - Error estimate on the prefactor for RMS internal field strength
    c_err_all = b.prefacError(alldatadict["p"], alldatadict["rmsINT"]["Le"]/alldatadict["fohm"]**0.5,
                                 model=[10**alldatadict["rmsINT"]["c"],alldatadict["rmsINT"]["m"]], plot=False, quiet=False)
    alldatadict["rmsINT"]["c_sd"] = c_err_all[1]
    fitall_1sdp, fitall_1sdm, fitall_2sdp, fitall_2sdm, fitall_3sdp, fitall_3sdm = b.getFitBounds(Pfit, alldatadict["rmsINT"]["c"],\
                                                                     alldatadict["rmsINT"]["c_sd"], alldatadict["rmsINT"]["m"])
    # - Error estimate on the prefactor for RMS CMB field strength
    c_err_rmscmb = b.prefacError(alldatadict["p"], alldatadict["rmsCMB"]["Le"]/alldatadict["fohm"]**0.5,
                                 model=[10**alldatadict["rmsCMB"]["c"],alldatadict["rmsCMB"]["m"]], plot=False, quiet=False)
    alldatadict["rmsCMB"]["c_sd"] = c_err_rmscmb[1]
    fitrmscmb_1sdp, fitrmscmb_1sdm, fitrmscmb_2sdp, fitrmscmb_2sdm, fitrmscmb_3sdp, fitrmscmb_3sdm = b.getFitBounds(Pfit, alldatadict["rmsCMB"]["c"],\
                                                                     alldatadict["rmsCMB"]["c_sd"], alldatadict["rmsCMB"]["m"])

    # - Error estimate on the prefactor for CMB dip field strength
    c_err_dipcmb = b.prefacError(alldatadict["p"], alldatadict["dipCMB"]["Le"]/alldatadict["fohm"]**0.5,
                                 model=[10**alldatadict["dipCMB"]["c"],alldatadict["dipCMB"]["m"]], plot=False, quiet=False)
    alldatadict["dipCMB"]["c_sd"] = c_err_dipcmb[1]
    fitdipcmb_1sdp, fitdipcmb_1sdm, fitdipcmb_2sdp, fitdipcmb_2sdm, fitdipcmb_3sdp, fitdipcmb_3sdm = b.getFitBounds(Pfit, alldatadict["dipCMB"]["c"],\
                                                                     alldatadict["dipCMB"]["c_sd"], alldatadict["dipCMB"]["m"])

# - Define symbols' properties for plotting
Cmax = np.log10(1000.0)#np.log10(np.max(Eall))
Cmin = np.log10(50.0)#np.log10(np.min(Eall))
datadict = b.getPlotProperties(datadict, categorise=categorise)

# --- Store fitted values in output files
b.saveFitValues(filename=outfpath+outfname, datadict=datadict, alldatadict=alldatadict)

# --- Calculate fits and prefac errors for each scaling on all field strengths
alldatadict = b.fitForceScalings(alldatadict, quiet=False)

# --- save out prefactor file
b.savePrefacValues(filename=outfpath+outfnamepf, indict=alldatadict, l_prefac_err=calc_prefac_err)

earthdict["dipCMB"]["min"], earthdict["dipCMB"]["max"], earthdict["rmsINT"]["min"], earthdict["rmsINT"]["max"],\
    earthdict["rmsCMB"]["min"], earthdict["rmsCMB"]["max"], earthdict["p"]["min"], earthdict["p"]["max"] = b.getEarthEstimates(quiet=False)

###############################################################################################################################
# Brms INTERNAL FIELD
###############################################################################################################################

fitEn    = 10**alldatadict["rmsINT"]["cIMA"]   * Pfit**alldatadict["rmsINT"]["mIMA"]
fitmac   = 10**alldatadict["rmsINT"]["cMAC"]   * Pfit**alldatadict["rmsINT"]["mMAC"]
fitimac  = 10**alldatadict["rmsINT"]["cIMAC"]  * Pfit**alldatadict["rmsINT"]["mIMAC"]
fitimacd = 10**alldatadict["rmsINT"]["cIMACd"] * Pfit**alldatadict["rmsINT"]["mIMACd"]
fitimaci = 10**alldatadict["rmsINT"]["cIMACi"] * Pfit**alldatadict["rmsINT"]["mIMACi"]

# fit scaling bounds
if calc_prefac_err:
    fitEn_1sdp, fitEn_1sdm, fitEn_2sdp, fitEn_2sdm, fitEn_3sdp, fitEn_3sdm = b.getFitBounds(Pfit, alldatadict["rmsINT"]["cIMA"],\
                                                                    alldatadict["rmsINT"]["cIMA_sd"], alldatadict["rmsINT"]["mIMA"])
    fitmac_1sdp, fitmac_1sdm, fitmac_2sdp, fitmac_2sdm, fitmac_3sdp, fitmac_3sdm = b.getFitBounds(Pfit, alldatadict["rmsINT"]["cMAC"],\
                                                                         alldatadict["rmsINT"]["cMAC_sd"], alldatadict["rmsINT"]["mMAC"])

# Plot simulations
plt.clf()
plt.rcParams["figure.figsize"] = [16,6]
ax = plt.gca()
ax, dd_xpos, legend_ypos = b.plotSimulations(ax, datadict=datadict, alldatadict=alldatadict, earthdict=earthdict, field="rmsINT",
                                                 cbarrange=[Cmin,Cmax], xrange=[xmin,xmax], yrange=[3e-5,0.2])
plt.loglog(Pfit, fitall, color="k", lw=lw_fit, zorder=2)
if (calc_prefac_err):
    plt.fill_between(Pfit, fitall_2sdm, fitall_2sdp, color="k", alpha=0.1, zorder=-1)
    plt.loglog(Pfit, fitall_1sdm, color="k", lw=0.5, ls='--', zorder=0) 
    plt.loglog(Pfit, fitall_1sdp, color="k", lw=0.5, ls='--', zorder=0) 

if "IMAC" in plt_extrap_scalings:
    plt.loglog(Pfit, fitimac, color="dimgrey"  ,linestyle="--", lw=lw_fit, label="$m=2/5$ (IMAC)")
if "IMACd" in plt_extrap_scalings:
    plt.loglog(Pfit, fitimacd, color="dimgrey" ,linestyle="-.", lw=lw_fit, label="$m=1/5$ (IMACd)")
if "IMACi" in plt_extrap_scalings:
    plt.loglog(Pfit, fitimaci, color="dimgrey" ,linestyle=":", lw=lw_fit, label="$m=3/10$ (IMACi)")
if "IMA" in plt_extrap_scalings:
    idxIMA = b.idxStr(plt_extrap_scalings, "IMA")[0]
    if (calc_prefac_err):
        plt.loglog(Pfit, fitEn, c=lc_fit[idxIMA], linestyle=ls_fit[idxIMA], lw=lw_fit,
                   label="$m=1/3$, $\sigma=$"+str(np.round(alldatadict["rmsINT"]["cIMA_sd"],3)))
        plt.fill_between(Pfit, fitEn_1sdm, fitEn_1sdp, color=sc_fit[idxIMA], alpha=sc_alpha, zorder=-1)
    else:
        plt.loglog(Pfit, fitEn, c=lc_fit[idxIMA], linestyle=ls_fit[idxIMA], lw=lw_fit, label="$m=1/3$ (QG-MAC)")
if "MAC" in plt_extrap_scalings:
    idxMAC = b.idxStr(plt_extrap_scalings, "MAC")[0]
    if (calc_prefac_err):
        plt.loglog(Pfit, fitmac, c=lc_fit[idxMAC], linestyle=ls_fit[idxMAC], lw=lw_fit,
                   label="$m=1/4$, $\sigma=$"+str(np.round(alldatadict["rmsINT"]["cMAC_sd"],3)))
        #plt.fill_between(Pfit, fitmac_1sdm, fitmac_1sdp, color=sc_fit[idxMAC], alpha=sc_alpha, zorder=-1)
        #plt.loglog(Pfit, fitmac_3sdm, c=lc_fit[idxMAC], linestyle=ls_fit[idxMAC], lw=1., zorder=-1)  #3sigma bounds on MAC
        #plt.loglog(Pfit, fitmac_3sdp, c=lc_fit[idxMAC], linestyle=ls_fit[idxMAC], lw=1., zorder=-1)
    else:
        plt.loglog(Pfit, fitmac, c=lc_fit[idxMAC], linestyle=ls_fit[idxMAC], lw=lw_fit, label="$m=1/4$ (MAC)")
# Plot title and axes labels
plt.xlabel('$p_A$')
if myfohm == 0: 
    plt.ylabel('$Le^{\\rm rms}_t/f_{ohm}^{1/2}$')
elif myfohm == 1: 
    plt.ylabel('$Le^{\\rm rms}_t$')
if myEr==1: 
    title_str = b.getPlotTitle(myfdip=myfdip, myEr=myEr, Er_range=[EMoEK_range[0],EMoEK_range[1]])
else:
    title_str = b.getPlotTitle(myfdip=myfdip, myEr=myEr, Er_range=[None]*2)
plt.title(title_str)
leg = plt.legend(bbox_to_anchor=(0.55, 0.01), loc=3, ncol=2, borderaxespad=0, 
                 labelspacing=0.1,handlelength=0.8, columnspacing=0.5)
for lh in leg.legendHandles:
    lh.set_linewidth(3.)

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
#Check ./fig directory exists, otherwise create it
if (not os.path.isdir('./fig')):
    os.mkdir('./fig')
plt.savefig(file2, format='pdf',bbox_inches="tight")
plt.savefig(file3, format='png',bbox_inches="tight")
###############################################################################################################################
# Brms CMB FIELD
###############################################################################################################################

# Calculate prefactors of Brms cmb
fitEn    = 10**alldatadict["rmsCMB"]["cIMA"]   * Pfit**alldatadict["rmsCMB"]["mIMA"]
fitmac   = 10**alldatadict["rmsCMB"]["cMAC"]   * Pfit**alldatadict["rmsCMB"]["mMAC"]
fitimac  = 10**alldatadict["rmsCMB"]["cIMAC"]  * Pfit**alldatadict["rmsCMB"]["mIMAC"]
fitimacd = 10**alldatadict["rmsCMB"]["cIMACd"] * Pfit**alldatadict["rmsCMB"]["mIMACd"]
fitimaci = 10**alldatadict["rmsCMB"]["cIMACi"] * Pfit**alldatadict["rmsCMB"]["mIMACi"]

# Plot simulations
del ax
plt.clf()
plt.rcParams["figure.figsize"] = [16,6]
ax = plt.gca()
ax, legend_xpos, legend_ypos = b.plotSimulations(ax, datadict=datadict, alldatadict=alldatadict, earthdict=earthdict, field="rmsCMB",
                                                 cbarrange=[Cmin,Cmax], xrange=[xmin,xmax], yrange=[5e-6,0.2])
plt.loglog(Pfit ,fitrmscmb, color="k", lw=lw_fit)
if (calc_prefac_err):
    plt.fill_between(Pfit, fitrmscmb_2sdm, fitrmscmb_2sdp, color="k", alpha=0.1, zorder=-1)
    plt.loglog(Pfit, fitrmscmb_1sdm, color="k", lw=0.5, ls='--', zorder=0) 
    plt.loglog(Pfit, fitrmscmb_1sdp, color="k", lw=0.5, ls='--', zorder=0) 
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

if myEr==1:
    title_str = b.getPlotTitle(myfdip=myfdip, myEr=myEr, Er_range=[EMoEK_range[0],EMoEK_range[1]])
else:
    title_str = b.getPlotTitle(myfdip=myfdip, myEr=myEr, Er_range=[None]*2)
plt.title(title_str)
leg = plt.legend(bbox_to_anchor=(0.55, 0.01), loc=3, ncol=2, borderaxespad=0, 
                 labelspacing=0.1,handlelength=0.8, columnspacing=0.5)
for lh in leg.legendHandles:
    lh.set_linewidth(3.)

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

###############################################################################################################################
# Bdip cmb plots
###############################################################################################################################

fitEn    = 10**alldatadict["dipCMB"]["cIMA"]   * Pfit**alldatadict["dipCMB"]["mIMA"]
fitmac   = 10**alldatadict["dipCMB"]["cMAC"]   * Pfit**alldatadict["dipCMB"]["mMAC"]
fitimac  = 10**alldatadict["dipCMB"]["cIMAC"]  * Pfit**alldatadict["dipCMB"]["mIMAC"]
fitimacd = 10**alldatadict["dipCMB"]["cIMACd"] * Pfit**alldatadict["dipCMB"]["mIMACd"]
fitimaci = 10**alldatadict["dipCMB"]["cIMACi"] * Pfit**alldatadict["dipCMB"]["mIMACi"]

# Error estimate on the prefactor
if calc_prefac_err:
    fitEn_1sdp, fitEn_1sdm, fitEn_2sdp, fitEn_2sdm, fitEn_3sdp, fitEn_3sdm = b.getFitBounds(Pfit, alldatadict["dipCMB"]["cIMA"],\
                                                                    alldatadict["dipCMB"]["cIMA_sd"], alldatadict["dipCMB"]["mIMA"])
    fitmac_1sdp, fitmac_1sdm, fitmac_2sdp, fitmac_2sdm, fitmac_3sdp, fitmac_3sdm = b.getFitBounds(Pfit, alldatadict["dipCMB"]["cMAC"],\
                                                                         alldatadict["dipCMB"]["cMAC_sd"], alldatadict["dipCMB"]["mMAC"])

# Plot simulations
del ax
plt.clf()
plt.rcParams["figure.figsize"] = [16,6]
ax = plt.gca()
ax, legend_xpos, legend_ypos = b.plotSimulations(ax, datadict=datadict, alldatadict=alldatadict, earthdict=earthdict,
                                                 field="dipCMB",
                                                 cbarrange=[Cmin,Cmax], xrange=[xmin,xmax], yrange=[5e-6,0.2])
plt.loglog(Pfit, fitdipcmb, color="k", lw=lw_fit)
if (calc_prefac_err):
    plt.fill_between(Pfit, fitdipcmb_2sdm, fitdipcmb_2sdp, color="k", alpha=0.1, zorder=-1)
    plt.loglog(Pfit, fitdipcmb_1sdm, color="k", lw=0.5, ls='--', zorder=0) 
    plt.loglog(Pfit, fitdipcmb_1sdp, color="k", lw=0.5, ls='--', zorder=0) 
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
                   label="$m=1/3$, $\sigma=$"+str(np.round(alldatadict["dipCMB"]["cIMA_sd"],3)))
        plt.fill_between(Pfit, fitEn_1sdm, fitEn_1sdp, color=sc_fit[idxIMA], alpha=sc_alpha, zorder=-1)
    else:
        plt.loglog(Pfit, fitEn, color=lc_fit[idxIMA], linestyle=ls_fit[idxIMA], lw=lw_fit,
                   label="$m=1/3$")
if "MAC" in plt_extrap_scalings:
    if (calc_prefac_err):
        plt.loglog(Pfit, fitmac, c=lc_fit[idxMAC], linestyle=ls_fit[idxMAC], lw=lw_fit,
                   label="$m=1/4$, $\sigma=$"+str(np.round(alldatadict["dipCMB"]["cMAC_sd"],3)))
    else:
        plt.loglog(Pfit, fitmac, color=lc_fit[idxMAC],linestyle=ls_fit[idxMAC], lw=lw_fit, label="$m=1/4$ (MAC)")

plt.xlabel('$p_A$')
if myfohm == 0: 
    plt.ylabel('$Le_{\\rm{cmb}}^{l=1}/f_{ohm}^{1/2}$')
elif myfohm ==1: 
    plt.ylabel('$Le_{\\rm{cmb}}^{l=1}$')
ax.set_yscale('log')
ax.set_xscale('log')

leg = plt.legend(bbox_to_anchor=(0.55, 0.01), loc=3, ncol=2, borderaxespad=0, 
                 labelspacing=0.1,handlelength=0.8, columnspacing=0.5)
for lh in leg.legendHandles:
    lh.set_linewidth(3.)

if myEr==1:
    title_str = b.getPlotTitle(myfdip=myfdip, myEr=myEr, Er_range=[EMoEK_range[0],EMoEK_range[1]])
else:
    title_str = b.getPlotTitle(myfdip=myfdip, myEr=myEr, Er_range=[None]*2)
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

###############################################################################################################################
# Zoomed figs plotted using subfigure
###############################################################################################################################

# - B rms core
del ax
plt.clf()
plt.figure(figsize=(16,15))
plt.subplot(3,1,1)
ax = plt.gca()
ax, legend_xpos, legend_ypos = b.plotSimulations(ax, datadict=datadict, alldatadict=alldatadict, earthdict=None, field="rmsINT",
                                                 colorbar=True, cbarrange=[Cmin,Cmax], xrange=[1e-10,1e-3], yrange=[1e-4,0.2])
plt.loglog(Pfit , fitall, color="black")
plt.loglog(5e-12 , 2e-4, marker="s", markersize=30)
if (calc_prefac_err):
    plt.loglog(Pfit,fitall_1sdp,c='k',ls=':',lw=1.)
    plt.loglog(Pfit,fitall_1sdm,c='k',ls=':',lw=1.)
    plt.fill_between(Pfit, fitall_1sdp, fitall_1sdm, color="k", alpha=0.1, zorder=-1)
    
if myfohm == 0: 
    plt.ylabel('$Le_t^{rms}/f_{ohm}^{1/2}$')
elif myfohm == 1: 
    plt.ylabel('$Le_t^{rms}$')
if myEr==1:
    title_str = b.getPlotTitle(myfdip=myfdip, myEr=myEr, Er_range=[EMoEK_range[0],EMoEK_range[1]])
else:
    title_str = b.getPlotTitle(myfdip=myfdip, myEr=myEr, Er_range=[None]*2)
plt.title("Field strength vs Buoyant power for "+title_str)
leg = plt.legend(bbox_to_anchor=(0.8, 0.05), loc=3, borderaxespad=0)
for lh in leg.legendHandles:
    lh.set_linewidth(3.)
"""
# --- DGM for improving the legend, try something like:
line = Line2D([], [], label='abc', color='red', linewidth=1.5, marker=ur'$\u25CC$',
              markeredgecolor='indigo', markeredgewidth=0.5, markersize=16)
plt.legend(handles=[line], numpoints=1)
"""
# - Brms at CMB
del ax
plt.subplot(3,1,2)
ax  = plt.gca()
ax, legend_xpos, legend_ypos = b.plotSimulations(ax, datadict=datadict, alldatadict=alldatadict, earthdict=None, field="rmsCMB",
                                                 colorbar=True, cbarrange=[Cmin,Cmax], xrange=[1e-10,1e-3], yrange=[1e-4,0.2])
plt.loglog(Pfit,fitrmscmb,color="black")
plt.loglog(5e-12 , 2e-4, marker="s", markersize=30)
if (calc_prefac_err):
    plt.loglog(Pfit,fitrmscmb_1sdp,c='k',ls=':',lw=1.)
    plt.loglog(Pfit,fitrmscmb_1sdm,c='k',ls=':',lw=1.)
    plt.fill_between(Pfit, fitrmscmb_1sdp, fitrmscmb_1sdm, color="k", alpha=0.1, zorder=-1)
    
if myfohm == 0: 
    plt.ylabel('$Le_{cmb}^{rms}/f_{ohm}^{1/2}$')
elif myfohm == 1: 
    plt.ylabel('$Le_{cmb}^{rms}$')

# - Bdip at CMB
del ax
plt.subplot(3,1,3)
ax = plt.gca()
ax, legend_xpos, legend_ypos = b.plotSimulations(ax, datadict=datadict, alldatadict=alldatadict, earthdict=None, field="dipCMB",
                                                 colorbar=True, cbarrange=[Cmin,Cmax], xrange=[1e-10,1e-3], yrange=[1e-4,0.2])
plt.loglog(Pfit,fitdipcmb,color="black")
if (calc_prefac_err):
    plt.loglog(Pfit,fitdipcmb_1sdp,c='k',ls=':',lw=1.)
    plt.loglog(Pfit,fitdipcmb_1sdm,c='k',ls=':',lw=1.)
    plt.fill_between(Pfit, fitdipcmb_1sdp, fitdipcmb_1sdm, color="k", alpha=0.1, zorder=-1)
    #plt.loglog(Pfit,fitdipcmb_3sdp,c='k',ls='--',lw=1.)
    #plt.loglog(Pfit,fitdipcmb_3sdm,c='k',ls='--',lw=1.)
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

# - Plot bdip stuff
if plt_bdip:
    b.plot_bdip(datadict, title_str)
