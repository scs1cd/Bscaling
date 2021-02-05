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
myfdip  = 1 # Use 0 for all fdip values, 1 for fdip > 0.50, 2 for filtering with fdip=(0.35,0.80), 3 for fdip=(0.40,0.80) (see below).
myfohm  = 0 # Use 0 for fohm factor, or 1 for NO fohm factor
myEkOPm = 0 # Use 1 (0) for (not) filtering in a specified range of Ek/Pm values
myEr    = 1 # Use 1 (0) for (not) filtering in specified EM/EK range

if (myEkOPm==1):
    EkOPm_range = [1.e-10,1.e-4] # Range of Ek/Pm to include in the analysis (Ek=\nu/\Omega*D^2 as in Aubert's definition)
else:
    EkOPm_range = None
if (myEr==1):
    EMoEK_range = [2e-15,1.e+19]
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
datadict = {"L":{  "plot":True, "dataset":"L", "d":{}, "rmsINT":{}, "rmsCMB":{}, "dipCMB":{}, "plotp":{}},
            "A":{  "plot":True, "dataset":"A", "d":{}, "rmsINT":{}, "rmsCMB":{}, "dipCMB":{}, "plotp":{}},
            "Y":{  "plot":True, "dataset":"Y", "d":{}, "rmsINT":{}, "rmsCMB":{}, "dipCMB":{}, "plotp":{}},
            "UC":{ "plot":True, "dataset":"UC", "d":{}, "rmsINT":{}, "rmsCMB":{}, "dipCMB":{}, "plotp":{}},
            "UCt":{"plot":True, "dataset":"UCt", "d":{}, "rmsINT":{}, "rmsCMB":{}, "dipCMB":{}, "plotp":{}},
            "APath":{"plot":False, "dataset":"APath", "d":{}, "rmsINT":{}, "rmsCMB":{}, "dipCMB":{}, "plotp":{}},
            "S":{"plot":False, "dataset":"S", "d":{}, "rmsINT":{}, "rmsCMB":{}, "dipCMB":{}, "plotp":{}}}

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

leedsOutName   = leedsname+filetag
yadavOutName   = yadavname+filetag
aubertOutName  = aubertname+filetag
christOutName  = christname+filetag
christOutNameT = christnamet+filetag
APathOutName   = APathname+filetag
SOutName       = APathname+filetag

# read and filter datasets
if datadict["L"]["plot"]:
    df, datadict = b.filter_table(infname=leedsname  ,outfname=leedsOutName, dataset="Leeds", fdip_range=fdip_range, 
                                  EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
if datadict["Y"]["plot"]:
    df, datadict = b.filter_table(infname=yadavname  ,outfname=yadavOutName, dataset="Yadav", fdip_range=fdip_range, 
                                  EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
if datadict["A"]["plot"]:
    df, datadict = b.filter_table(infname=aubertname ,outfname=aubertOutName, dataset="Aubert", fdip_range=fdip_range, 
                                  EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
if datadict["UC"]["plot"]:
    df, datadict = b.filter_table(infname=christname ,outfname=christOutName, dataset="Christensen", fdip_range=fdip_range, 
                                  EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
if datadict["UCt"]["plot"]:
    df, datadict = b.filter_table(infname=christnamet,outfname=christOutNameT, dataset="ChristensenT", fdip_range=fdip_range, 
                                  EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
if datadict["APath"]["plot"]:
    df, datadict = b.filter_table(infname=APathname,  outfname=APathOutName, dataset="APath",fdip_range=fdip_range, 
                                  EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
if datadict["S"]["plot"]:
    df, datadict = b.filter_table(infname=Sname    ,  outfname=SOutName    , dataset="S"    ,fdip_range=fdip_range, 
                                  EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)

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

# Leeds dataset
if datadict["L"]["plot"]:
    datadict["L"]["E"]  = datadict["L"]["d"]["Ek"]
    datadict["L"]["Pm"] = datadict["L"]["d"]["Pm"]
    datadict["L"]["Rm"] = datadict["L"]["d"]["RmAve"]
    
    # Calculate shell volume
    ar = datadict["L"]["d"]["ar"] # aspect ratio
    volS = np.asarray([b.shellVolume(ar[i]) for i in range(len(ar))])

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
    # make fit
    datadict["L"]["rmsINT"]["ssr"], datadict["L"]["rmsINT"]["m"], datadict["L"]["rmsINT"]["c"], datadict["L"]["rmsINT"]["res"] = \
        b.fits(datadict["L"]["p"], datadict["L"]["rmsINT"]["Le"], datadict["L"]["fohm"])
    datadict["L"]["rmsCMB"]["ssr"], datadict["L"]["rmsCMB"]["m"], datadict["L"]["rmsCMB"]["c"], datadict["L"]["rmsCMB"]["res"] = \
        b.fits(datadict["L"]["p"], datadict["L"]["rmsCMB"]["Le"], datadict["L"]["fohm"])
    datadict["L"]["dipCMB"]["ssr"], datadict["L"]["dipCMB"]["m"], datadict["L"]["dipCMB"]["c"], datadict["L"]["dipCMB"]["res"] = \
        b.fits(datadict["L"]["p"], datadict["L"]["dipCMB"]["Le"], datadict["L"]["fohm"])

# Get Aubert data
# NB - Aubert uses different radii! chi dependence of parameters?
if datadict["A"]:
    datadict["A"]["E"]    = datadict["A"]["d"]["E"]
    datadict["A"]["Pm"]   = datadict["A"]["d"]["Pm"]
    datadict["A"]["Rm"]   = datadict["A"]["d"]["Ro"] / datadict["A"]["E"] * datadict["A"]["Pm"]
    datadict["A"]["fohm"] = datadict["A"]["d"]["fohm"]
    EA   , LeA   ,bdipA,fdipA,PA   ,fohmA, dAall,RoA,PmA   = np.loadtxt(aubertname, usecols=(0,7,8,9,12,13,10,6,3), skiprows=1, unpack='true')
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
    # fits
    datadict["A"]["rmsINT"]["ssr"], datadict["A"]["rmsINT"]["m"], datadict["A"]["rmsINT"]["c"], datadict["A"]["rmsINT"]["res"] = \
        b.fits(datadict["A"]["p"], datadict["A"]["rmsINT"]["Le"], datadict["A"]["fohm"])
    datadict["A"]["rmsCMB"]["ssr"], datadict["A"]["rmsCMB"]["m"], datadict["A"]["rmsCMB"]["c"], datadict["A"]["rmsCMB"]["res"] = \
        b.fits(datadict["A"]["p"], datadict["A"]["rmsCMB"]["Le"], datadict["A"]["fohm"])
    datadict["A"]["dipCMB"]["ssr"], datadict["A"]["dipCMB"]["m"], datadict["A"]["dipCMB"]["c"], datadict["A"]["dipCMB"]["res"] = \
        b.fits(datadict["A"]["p"], datadict["A"]["dipCMB"]["Le"], datadict["A"]["fohm"])

# Christensen dataset
# Ek has no factor 2
# Fundamental length scale is shell thickness, time scale is viscous diffusion time, magnetic field scale is "Elsasser scale".
# NOTE that variables names "Els" are actually RMS field strengths!!
if datadict["UC"]["plot"]:
    datadict["UC"]["E"]  = datadict["UC"]["d"]["E"]
    datadict["UC"]["Pm"] = datadict["UC"]["d"]["Pm"]  
    datadict["UC"]["Rm"] = datadict["UC"]["d"]["Rm"]  
    datadict["UC"]["p"]  = 1e7 * datadict["UC"]["d"]["pow"] * datadict["UC"]["E"]**3
    if myfohm == 1:
        datadict["UC"]["fohm"] = np.ones(len(datadict["UC"]["d"]["jou"]))
    else:
        datadict["UC"]["fohm"] = datadict["UC"]["d"]["jou"] / 100

    datadict["UC"]["rmsINT"]["Le"] = np.sqrt(datadict["UC"]["E"]/datadict["UC"]["Pm"]) * datadict["UC"]["d"]["B"]        # Eqn 14 of CA06
    datadict["UC"]["rmsCMB"]["Le"] = np.sqrt(datadict["UC"]["E"]/datadict["UC"]["Pm"]) * datadict["UC"]["d"]["Bsur"]
    datadict["UC"]["dipCMB"]["Le"] = np.sqrt(datadict["UC"]["E"]/datadict["UC"]["Pm"]) * datadict["UC"]["d"]["Bdip"]
    datadict["UC"]["bdip"]         = datadict["UC"]["rmsINT"]["Le"]/datadict["UC"]["dipCMB"]["Le"]
    datadict["UC"]["fdip"]         = datadict["UC"]["d"]["Bdip"]/datadict["UC"]["d"]["B12"]
    # fits
    datadict["UC"]["rmsINT"]["ssr"], datadict["UC"]["rmsINT"]["m"], datadict["UC"]["rmsINT"]["c"], datadict["UC"]["rmsINT"]["res"] = \
        b.fits(datadict["UC"]["p"], datadict["UC"]["rmsINT"]["Le"], datadict["UC"]["fohm"])
    datadict["UC"]["rmsCMB"]["ssr"], datadict["UC"]["rmsCMB"]["m"], datadict["UC"]["rmsCMB"]["c"], datadict["UC"]["rmsCMB"]["res"] = \
        b.fits(datadict["UC"]["p"], datadict["UC"]["rmsCMB"]["Le"], datadict["UC"]["fohm"])
    datadict["UC"]["dipCMB"]["ssr"], datadict["UC"]["dipCMB"]["m"], datadict["UC"]["dipCMB"]["c"], datadict["UC"]["dipCMB"]["res"] = \
        b.fits(datadict["UC"]["p"], datadict["UC"]["dipCMB"]["Le"], datadict["UC"]["fohm"])

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

if datadict["UCt"]["plot"]:
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
    # fits
    datadict["UCt"]["rmsINT"]["ssr"], datadict["UCt"]["rmsINT"]["m"], datadict["UCt"]["rmsINT"]["c"], datadict["UCt"]["rmsINT"]["res"] = \
        b.fits(datadict["UCt"]["p"], datadict["UCt"]["rmsINT"]["Le"], datadict["UCt"]["fohm"])
    datadict["UCt"]["rmsCMB"]["ssr"], datadict["UCt"]["rmsCMB"]["m"], datadict["UCt"]["rmsCMB"]["c"], datadict["UCt"]["rmsCMB"]["res"] = \
        b.fits(datadict["UCt"]["p"], datadict["UCt"]["rmsCMB"]["Le"], datadict["UCt"]["fohm"])
    datadict["UCt"]["dipCMB"]["ssr"], datadict["UCt"]["dipCMB"]["m"], datadict["UCt"]["dipCMB"]["c"], datadict["UCt"]["dipCMB"]["res"] = \
        b.fits(datadict["UCt"]["p"], datadict["UCt"]["dipCMB"]["Le"], datadict["UCt"]["fohm"])

    if chk == 1:
        # CHECK: Relate Elsasser to Em, magnetic energy density (i.e. per unit volume)
        # Em = Els / (2*Pm*E) = Le**2 / (2 E**2) 
        print('*******Christensen qT******')
        EfromLe  = datadict["UCt"]["rmsINT"]["Le"]**2 * datadict["UCt"]["Pm"] / 2.0 / datadict["UCt"]["E"]**2 / datadict["UCt"]["Pm"]
        Efromfile = datadict["UCt"]["Emag"]
        print('Energy difference (%) = ', np.abs( (EfromLe-Efromfile) / Efromfile)*100.0)
        print('***************************\n')

# ---
# Yadav - NOTE - he uses fixed T. Extracts E, Pm, Elasser, Elsasser_CMB, Dip_CMB, Buo_pow, Ohm_diss
if datadict["Y"]["plot"]:
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
    datadict["Y"]["p"]            = datadict["Y"]["E"]**3 * datadict["Y"]["d"]["Buo_pow"]/b.shellVolume(0.35)
    datadict["Y"]["bdip"]         = datadict["Y"]["rmsINT"]["Le"]/datadict["Y"]["dipCMB"]["Le"] 
    datadict["Y"]["fdip"]         = datadict["Y"]["d"]["Dip_CMB_l11"]

    datadict["Y"]["rmsINT"]["ssr"], datadict["Y"]["rmsINT"]["m"], datadict["Y"]["rmsINT"]["c"], datadict["Y"]["rmsINT"]["res"] = \
        b.fits(datadict["Y"]["p"], datadict["Y"]["rmsINT"]["Le"], datadict["Y"]["fohm"])
    datadict["Y"]["rmsCMB"]["ssr"], datadict["Y"]["rmsCMB"]["m"], datadict["Y"]["rmsCMB"]["c"], datadict["Y"]["rmsCMB"]["res"] = \
        b.fits(datadict["Y"]["p"], datadict["Y"]["rmsCMB"]["Le"], datadict["Y"]["fohm"])
    datadict["Y"]["dipCMB"]["ssr"], datadict["Y"]["dipCMB"]["m"], datadict["Y"]["dipCMB"]["c"], datadict["Y"]["dipCMB"]["res"] = \
        b.fits(datadict["Y"]["p"], datadict["Y"]["dipCMB"]["Le"], datadict["Y"]["fohm"])

    if chk == 1:
        # CHECK: Relate Elsasser to EM, total magnetic energy 
        # Em = Els*Els / (2*Pm*E)
        print('*******Yadav***************')
        EfromLe  = datadict["Y"]["rmsINT"]["Le"]**2 * 14.59 / 2.0 / datadict["Y"]["E"]**2
        Efromfile = datadict["Y"]["d"]["ME_pol"]+datadict["Y"]["d"]["ME_tor"]
        print('Energy difference (%) = ', np.abs( (EfromLe-Efromfile) / Efromfile)*100.0)
        print('***************************\n')
# ---
if datadict['APath']["plot"]:
    fohm_dum = datadict["APath"]['d']['fohm']
    P_dum    = datadict["APath"]['d']['p']
    E_dum    = datadict["APath"]['d']['E']
    Pm_dum   = datadict["APath"]['d']['Pm']
    Els_dum  = datadict["APath"]['d']['Els']
    Els12_dum  = datadict["APath"]['d']['rmsCMBtotal']
    Els1_dum   = datadict["APath"]['d']['RMSCMBl=1']
    
    datadict["APath"]['d']['bdip'] = 0.0

    datadict["APath"]['rmsINT']['Le'] = np.sqrt( (E_dum/Pm_dum)   * Els_dum)
    datadict["APath"]['rmsCMB']['Le'] = np.sqrt( (E_dum/Pm_dum) ) * Els12_dum
    datadict["APath"]['dipCMB']['Le'] = np.sqrt( (E_dum/Pm_dum) ) * Els1_dum

    datadict["APath"]['rmsINT']['ssr'], datadict["APath"]['rmsINT']['m'], datadict["APath"]['rmsINT']['c'], datadict["APath"]['rmsINT']['res'] = b.fits(P_dum, datadict["APath"]['rmsINT']['Le'], fohm_dum)
    datadict["APath"]['rmsCMB']['ssr'], datadict["APath"]['rmsCMB']['m'], datadict["APath"]['rmsCMB']['c'], datadict["APath"]['rmsCMB']['res'] = b.fits(P_dum, datadict["APath"]['rmsCMB']['Le'], fohm_dum)
    datadict["APath"]['dipCMB']['ssr'], datadict["APath"]['dipCMB']['m'], datadict["APath"]['dipCMB']['c'], datadict["APath"]['dipCMB']['res'] = b.fits(P_dum, datadict["APath"]['dipCMB']['Le'], fohm_dum)

# ---
if datadict['S']["plot"]:
    fohm_dum   = datadict["S"]['d']['fohm']
    P_dum   = datadict["S"]['d']['p']
    E_dum   = datadict["S"]['d']['E']
    Pm_dum  = datadict["S"]['d']['Pm']
    Els_dum = datadict["S"]['d']['Els']
    Els_cmb = datadict["S"]['d']['rmsCMBtotal']
    Els1_dum  = datadict["S"]['d']['RMSCMBl=1']

    datadict["S"]['d']['bdip'] = 0.0

    datadict["S"]['rmsINT']['Le'] = np.sqrt( (E_dum/Pm_dum) * Els_dum)
    datadict["S"]['rmsCMB']['Le'] = np.sqrt( (E_dum/Pm_dum) ) * Els_cmb
    datadict["S"]['dipCMB']['Le'] = np.sqrt( (E_dum/Pm_dum) ) * Els1_dum


    datadict["S"]['rmsINT']['ssr'], datadict["S"]['rmsINT']['m'], datadict["S"]['rmsINT']['c'], datadict["S"]['rmsINT']['res'] = b.fits(P_dum, datadict["S"]['rmsINT']['Le'], fohm_dum)
    datadict["S"]['rmsCMB']['ssr'], datadict["S"]['rmsCMB']['m'], datadict["S"]['rmsCMB']['c'], datadict["S"]['rmsCMB']['res'] = b.fits(P_dum, datadict["S"]['rmsCMB']['Le'], fohm_dum)
    datadict["S"]['dipCMB']['ssr'], datadict["S"]['dipCMB']['m'], datadict["S"]['dipCMB']['c'], datadict["S"]['dipCMB']['res'] = b.fits(P_dum, datadict["S"]['dipCMB']['Le'], fohm_dum)

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


##########################################################
# Extrapolated figures
##########################################################

earthdict["dipCMB"]["min"], earthdict["dipCMB"]["max"], earthdict["rmsINT"]["min"], earthdict["rmsINT"]["max"],\
    earthdict["rmsCMB"]["min"], earthdict["rmsCMB"]["max"], earthdict["p"]["min"], earthdict["p"]["max"] = b.getEarthEstimates(quiet=False)

##########################
# Brms plots
##########################

fitEn    = 10**alldatadict["rmsINT"]["cIMA"] * Pfit**alldatadict["rmsINT"]["mIMA"]
fitmac   = 10**alldatadict["rmsINT"]["cMAC"] * Pfit**alldatadict["rmsINT"]["mMAC"]
fitimac  = 10**alldatadict["rmsINT"]["cIMAC"] * Pfit**alldatadict["rmsINT"]["mIMAC"]
fitimacd = 10**alldatadict["rmsINT"]["cIMACd"] * Pfit**alldatadict["rmsINT"]["mIMACd"]
fitimaci = 10**alldatadict["rmsINT"]["cIMACi"]* Pfit**alldatadict["rmsINT"]["mIMACi"]

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
    plt.loglog(Pfit, fitimacd, color="dimgrey"  ,linestyle="-.", lw=lw_fit, label="$m=1/5$ (IMACd)")
if "IMACi" in plt_extrap_scalings:
    plt.loglog(Pfit, fitimaci, color="dimgrey"  ,linestyle=":", lw=lw_fit, label="$m=3/10$ (IMACi)")
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
# Brms cmb plots
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
    plt.title("Field strength ($Le$) vs Buoyant Power ($P$) for $%.2f<f_{dip}<%.2f$" %(fdip_min,fdip_max))
elif myfdip == 3:
    plt.title("Field strength ($Le$) vs Buoyant Power ($P$) for $%.2f<f_{dip}<%.2f$" %(fdip_min,fdip_max))
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
plt.text(0.05, 0.61, "$m$ = "+str(np.round(alldatadict["dipCMB"]["m"],2))+\
                              "$\pm$"+str(np.round(alldatadict["dipCMB"]["res"],2))+\
                              ", SSR="+str(np.round(alldatadict["dipCMB"]["ssr"],2)), transform=ax.transAxes, color='black')
plt.text(0.05, 0.56, "$\sigma$ = "+str(np.round(alldatadict["dipCMB"]["c_sd"],4)), transform=ax.transAxes, color='black')
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
