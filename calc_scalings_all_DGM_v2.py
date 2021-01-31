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
myfdip  = 0 # Use 0 for all fdip values, 1 for fdip > 0.50, 2 for filtering with fdip=(0.35,0.80), 3 for fdip=(0.40,0.80) (see below).
myfohm  = 0 # Use 0 for fohm factor, or 1 for NO fohm factor
myEkOPm = 0 # Use 1 (0) for (not) filtering in a specified range of Ek/Pm values
myEr    = 1 # Use 1 (0) for (not) filtering in specified EM/EK range

if (myEkOPm==1):
    EkOPm_range = [1.e-10,1.e-4] # Range of Ek/Pm to include in the analysis (Ek=\nu/\Omega*D^2 as in Aubert's definition)
else:
    EkOPm_range = None
if (myEr==1):
    EMoEK_range = [2.,1.e+9]
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
            "APath":{"plot":True, "dataset":"APath", "d":{}, "rmsINT":{}, "rmsCMB":{}, "dipCMB":{}, "plotp":{}},
            "S":{"plot":True, "dataset":"S", "d":{}, "rmsINT":{}, "rmsCMB":{}, "dipCMB":{}, "plotp":{}}}

alldatadict = {"plot":True, "rmsINT":{}, "rmsCMB":{}, "dipCMB":{}, "plotp":{}}

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

# read datasets, filter and calculate fdip
if datadict["L"]["plot"]:
    df, datadict = b.filter_table(infname=leedsname  ,outfname=leedsOutName, dataset="Leeds", fdip_range=fdip_range, 
                                  EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
    datadict["L"]["fdip"] = b.get_fdip(infname=leedsname, dataset="Leeds")
if datadict["Y"]["plot"]:
    df, datadict = b.filter_table(infname=yadavname  ,outfname=yadavOutName, dataset="Yadav", fdip_range=fdip_range, 
                                  EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
    datadict["Y"]["fdip"] = b.get_fdip(infname=yadavname, dataset="Yadav")
if datadict["A"]["plot"]:
    df, datadict = b.filter_table(infname=aubertname ,outfname=aubertOutName, dataset="Aubert", fdip_range=fdip_range, 
                                  EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
    datadict["A"]["fdip"] = b.get_fdip(infname=aubertname, dataset="Aubert")
if datadict["UC"]["plot"]:
    df, datadict = b.filter_table(infname=christname ,outfname=christOutName, dataset="Christensen", fdip_range=fdip_range, 
                                  EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
    datadict["UC"]["fdip"] = b.get_fdip(infname=christname, dataset="Christensen")
if datadict["UCt"]["plot"]:
    df, datadict = b.filter_table(infname=christnamet,outfname=christOutNameT, dataset="ChristensenT", fdip_range=fdip_range, 
                                  EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
    datadict["UCt"]["fdip"] = b.get_fdip(infname=christnamet, dataset="ChristensenT")
if datadict["APath"]["plot"]:
    df, datadict = b.filter_table(infname=APathname,  outfname=APathOutName, dataset="APath",fdip_range=fdip_range, 
                                  EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
if datadict["S"]["plot"]:
    df, datadict = b.filter_table(infname=Sname    ,  outfname=SOutName    , dataset="S"    ,fdip_range=fdip_range, 
                                  EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
leedsname = leedsOutName
yadavname = yadavOutName
aubertname = aubertOutName
christname = christOutName
christnamet = christOutNameT
APathname = APathOutName

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
    ar,EC,PmC,EmC,VDC,ODC,ElsCmbC,ElsDipCmbC,gauss_ElsCmbC,gauss_ElsDipCmbC,RmC = np.loadtxt(leedsname,
                                                                                  usecols=(76,52,70,63,54,64,58,60,3,5,73),
                                                                                  skiprows=1, unpack='true', delimiter=",")
    if (check_gauss_Led):
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(211)
        #xpltmin = 0.; xpltmax = np.max(ElsCmbC)
        xpltmin = 0.; xpltmax = 5.
        ypltmin = 0.; ypltmax = np.max(gauss_ElsCmbC)
        ax.plot([0.,5.],[0.,5.], ls='-', lw=1., c='k')
        ax.scatter(ElsCmbC,gauss_ElsCmbC)
        for ii in range(len(leeds_df['id'].values)):
            ax.text(ElsCmbC[ii],gauss_ElsCmbC[ii],leeds_df['id'].values[ii])
        ax.set_xlim([xpltmin,xpltmax])
        ax.set_ylim([ypltmin,ypltmax])
        ax.set_xlabel('$\Lambda_{\mathrm{CMB}}$')
        ax.set_ylabel('$\sum_{\ell=1}^{32} W_\ell$')
        ax = fig.add_subplot(212)
        xpltmin = 0.; xpltmax = np.max(ElsDipCmbC)
        ypltmin = 0.; ypltmax = np.max(gauss_ElsDipCmbC)
        ax.plot([0.,5.],[0.,5.], ls='-', lw=1., c='k')
        ax.scatter(ElsDipCmbC,gauss_ElsDipCmbC)
        ax.set_xlim([xpltmin,xpltmax])
        ax.set_ylim([ypltmin,ypltmax])
        ax.set_xlabel('$\Lambda_{\mathrm{CMB}}^{\ell=1}$')
        ax.set_ylabel('$W_1$')
        plt.show(block=False)
        plt.tight_layout()
        plt.savefig('./fig/Check_Emag_Gauss.pdf',format='pdf')
        del ax

    # Calculate shell volume
    volS = np.asarray([b.shellVolume(ar[i]) for i in range(len(ar))])

    if myfohm == 1:
        fohmC = np.ones(len(ODC))
    else:
        fohmC = ODC/(VDC+ODC)
    LeC      = np.sqrt(8.0)* (EC/PmC)  * np.sqrt(EmC/volS)       # Manual eq. (3.38)
    LeCrmscmb= 2. * np.sqrt(EC/PmC * ElsCmbC)
    LeCdipcmb= 2. * np.sqrt(EC/PmC * ElsDipCmbC)
    PC       = 8.0 * (EC/PmC)**3 * (VDC+ODC) / volS
    bdipC    = LeC/LeCdipcmb
    # Store in dictionary and make fit
    datadict["L"]['d']["E"] = EC
    datadict["L"]['d']["Pm"] = PmC
    datadict["L"]['d']["Rm"] = RmC 
    datadict["L"]['d']["fohm"] = fohmC 
    datadict["L"]['d']["p"] = PC
    datadict["L"]['d']["bdip"] = bdipC
    datadict["L"]["rmsINT"]["Le"] = LeC; datadict["L"]["rmsCMB"]["Le"] = LeCrmscmb; datadict["L"]["dipCMB"]["Le"] = LeCdipcmb

    datadict["L"]["rmsINT"]["ssr"], datadict["L"]["rmsINT"]["m"], datadict["L"]["rmsINT"]["c"], datadict["L"]["rmsINT"]["res"] = b.fits(PC, LeC, fohmC)
    datadict["L"]["rmsCMB"]["ssr"], datadict["L"]["rmsCMB"]["m"], datadict["L"]["rmsCMB"]["c"], datadict["L"]["rmsCMB"]["res"] = b.fits(PC, LeCrmscmb, fohmC)
    datadict["L"]["dipCMB"]["ssr"], datadict["L"]["dipCMB"]["m"], datadict["L"]["dipCMB"]["c"], datadict["L"]["dipCMB"]["res"] = b.fits(PC, LeCdipcmb, fohmC)

# Get Aubert data
# NB - Aubert uses different radii! chi dependence of parameters?
if datadict["A"]:
    EA   , LeA   ,bdipA,fdipA,PA   ,fohmA, dAall,RoA,PmA   = np.loadtxt(aubertname, usecols=(0,7,8,9,12,13,10,6,3), skiprows=1, unpack='true')
    RmA = RoA / EA * PmA
    if myfohm == 1:
        fohmA = np.ones(len(fohmA))
    LeAdipcmb    = LeA   /bdipA
    LeArmscmb    = LeA   /(bdipA*fdipA)
    lAall = 2.0*np.pi / (dAall + 0.5) # Jeans' formula.
    # Store in dictionary and make fit
    datadict["A"]["E"]       = EA
    datadict["A"]["Pm"]      = PmA
    datadict["A"]['d']["Rm"] = RmA
    datadict["A"]['d']["fohm"]    = np.array(fohmA)
    datadict["A"]["p"]       = PA
    datadict["A"]["bdip"]    = bdipA
    datadict["A"]["rmsINT"]["Le"] = np.array(LeA); datadict["A"]["rmsCMB"]["Le"] = LeArmscmb; datadict["A"]["dipCMB"]["Le"] = LeAdipcmb

    datadict["A"]["rmsINT"]["ssr"], datadict["A"]["rmsINT"]["m"], datadict["A"]["rmsINT"]["c"], datadict["A"]["rmsINT"]["res"] = b.fits(PA, LeA, fohmA)
    datadict["A"]["rmsCMB"]["ssr"], datadict["A"]["rmsCMB"]["m"], datadict["A"]["rmsCMB"]["c"], datadict["A"]["rmsCMB"]["res"] = b.fits(PA, LeArmscmb, fohmA)
    datadict["A"]["dipCMB"]["ssr"], datadict["A"]["dipCMB"]["m"], datadict["A"]["dipCMB"]["c"], datadict["A"]["dipCMB"]["res"] = b.fits(PA, LeAdipcmb, fohmA)

# Christensen dataset
# Ek has no factor 2
# Fundamental length scale is shell thickness, time scale is viscous diffusion time, magnetic field scale is "Elsasser scale".
# NOTE that variables names "Els" are actually RMS field strengths!!
if datadict["UC"]["plot"]:
    EUC,PmUC,ElsUC,ElsrmscmbUC,ElsdipcmbUC,PUC,fohmUC,RmUC,EmagUC = np.loadtxt(christname, usecols=(0,4,26,27,29,31,32,10,19), skiprows=1, unpack='true')
    PUC       = 1e7 * PUC * EUC**3
    fohmUC    = fohmUC / 100
    if myfohm == 1:
        fohmUC = np.ones(len(fohmUC))
    LeUC      = np.sqrt( (EUC/PmUC) ) * ElsUC        # Eqn 14 of CA06
    LeUCrmscmb= np.sqrt( (EUC/PmUC) ) * ElsrmscmbUC
    LeUCdipcmb= np.sqrt( (EUC/PmUC) ) * ElsdipcmbUC
    bdipUC    = LeUC/LeUCdipcmb
    # Store in dictionary and make fit
    #datadict["UC"]["E"] = EUC
    #datadict["UC"]["Pm"] = PmUC
    #datadict["UC"]["Rm"] = RmUC
    
    datadict["UC"]['d']["fohm"] = fohmUC
    datadict["UC"]['d']["p"]    = PUC
    datadict["UC"]['d']["bdip"] = bdipUC
    datadict["UC"]["rmsINT"]["Le"] = LeUC; datadict["UC"]["rmsCMB"]["Le"] = LeUCrmscmb; datadict["UC"]["dipCMB"]["Le"] = LeUCdipcmb

    datadict["UC"]["rmsINT"]["ssr"], datadict["UC"]["rmsINT"]["m"], datadict["UC"]["rmsINT"]["c"], datadict["UC"]["rmsINT"]["res"] = b.fits(PUC, LeUC, fohmUC)
    datadict["UC"]["rmsCMB"]["ssr"], datadict["UC"]["rmsCMB"]["m"], datadict["UC"]["rmsCMB"]["c"], datadict["UC"]["rmsCMB"]["res"] = b.fits(PUC, LeUCrmscmb, fohmUC)
    datadict["UC"]["dipCMB"]["ssr"], datadict["UC"]["dipCMB"]["m"], datadict["UC"]["dipCMB"]["c"], datadict["UC"]["dipCMB"]["res"] = b.fits(PUC, LeUCdipcmb, fohmUC)

    if chk == 1:
        # CHECK: Relate Elsasser to Em, magnetic energy density (i.e. per unit volume)
        # Em = Els / (2*Pm*E) = Le**2 / (2 E**2) 
        # P is printed to chk against Christensen 10 Fig 2. Note he uses Rc instead of D for 
        #   length and there is a factor 1/4pi in the Ra_Q definition. 
        #   NOTE - I still do not understand how to calculate his F. 
        print('\n*******Christensen FF******')
        EfromLe  = LeUC**2 * PmUC / 2.0 / EUC**2 / PmUC
        Efromfile = EmagUC
        print('Energy difference (%) = ', np.abs( (EfromLe-Efromfile) / Efromfile)*100.0)
        print('Power at E=3e-5       = ', PUC[26:35] / (1-0.35)**2 *4.0*np.pi)
        print('***************************\n')

if datadict["UCt"]["plot"]:
    EUCt,PmUCt,ElsUCt,ElsrmscmbUCt,ElsdipcmbUCt,PUCt,fohmUCt,RmUCt,EmagUCt = np.loadtxt(christnamet, usecols=(0,4,27,28,30,32,33,11,20), skiprows=1, unpack='true')
    LeUCt      = np.sqrt((EUCt/PmUCt)) * ElsUCt
    LeUCtrmscmb= np.sqrt((EUCt/PmUCt)) * ElsrmscmbUCt
    LeUCtdipcmb= np.sqrt((EUCt/PmUCt)) * ElsdipcmbUCt
    PUCt       = 1e7 * PUCt * EUCt**3
    fohmUCt    = fohmUCt / 100
    if myfohm == 1:
        fohmUCt = np.ones(len(fohmUCt))
    bdipUCt    = LeUCt/LeUCtdipcmb
    # Store in dictionary and make fit
    #datadict["UCt"]["E"] = EUCt
    #datadict["UCt"]["Pm"] = PmUCt
    #datadict["UCt"]["Rm"] = RmUCt
    datadict["UCt"]['d']["fohm"] = fohmUCt
    datadict["UCt"]['d']["p"] = PUCt
    datadict["UCt"]['d']["bdip"] = bdipUCt
    datadict["UCt"]["rmsINT"]["Le"] = LeUCt; datadict["UCt"]["rmsCMB"]["Le"] = LeUCtrmscmb; datadict["UCt"]["dipCMB"]["Le"] = LeUCtdipcmb

    datadict["UCt"]["rmsINT"]["ssr"], datadict["UCt"]["rmsINT"]["m"], datadict["UCt"]["rmsINT"]["c"], datadict["UCt"]["rmsINT"]["res"] = b.fits(PUCt, LeUCt, fohmUCt)
    datadict["UCt"]["rmsCMB"]["ssr"], datadict["UCt"]["rmsCMB"]["m"], datadict["UCt"]["rmsCMB"]["c"], datadict["UCt"]["rmsCMB"]["res"] = b.fits(PUCt, LeUCtrmscmb, fohmUCt)
    datadict["UCt"]["dipCMB"]["ssr"], datadict["UCt"]["dipCMB"]["m"], datadict["UCt"]["dipCMB"]["c"], datadict["UCt"]["dipCMB"]["res"] = b.fits(PUCt, LeUCtdipcmb, fohmUCt)

    if chk == 1:
        # CHECK: Relate Elsasser to Em, magnetic energy density (i.e. per unit volume)
        # Em = Els / (2*Pm*E) = Le**2 / (2 E**2) 
        print('*******Christensen qT******')
        EfromLe  = LeUCt**2 * PmUCt / 2.0 / EUCt**2 / PmUCt
        Efromfile = EmagUCt
        print('Energy difference (%) = ', np.abs( (EfromLe-Efromfile) / Efromfile)*100.0)
        print('***************************\n')

# Yadav - NOTE - he uses fixed T. Extracts E, Pm, Elasser, Elsasser_CMB, Dip_CMB, Buo_pow, Ohm_diss
if datadict["Y"]["plot"]:
    EY,PmY,ElY,ElCMBY,dipCMBY,BWY,ODY,ReY,EmPY,EmTY  = np.loadtxt(yadavname, usecols=(0,2,10,11,13,15,16,5,8,9), skiprows=1, unpack='true')
    RmY = ReY/PmY
    if myfohm == 1:
        fohmY = np.ones(len(ODY))
    else:
        fohmY    = ODY/BWY
    LeY      = np.sqrt(ElY   *EY/PmY)
    LeYrmscmb= np.sqrt(ElCMBY*EY/PmY)
    LeYdipcmb= np.sqrt(ElCMBY*dipCMBY*EY/PmY)
    PY       = (EY)**3 * BWY   / b.shellVolume(0.35)
    bdipY    = LeY/LeYdipcmb
    # Store in dictionary and make fit
    #datadict["Y"]["E"] = EY
    #datadict["Y"]["Pm"] = PmY
    datadict["Y"]['d']["Rm"] = RmY
    datadict["Y"]['d']["fohm"] = fohmY
    datadict["Y"]['d']["p"] = PY
    datadict["Y"]['d']["bdip"] = bdipY
    datadict["Y"]["rmsINT"]["Le"] = LeY; datadict["Y"]["rmsCMB"]["Le"] = LeYrmscmb; datadict["Y"]["dipCMB"]["Le"] = LeYdipcmb

    datadict["Y"]["rmsINT"]["ssr"], datadict["Y"]["rmsINT"]["m"], datadict["Y"]["rmsINT"]["c"], datadict["Y"]["rmsINT"]["res"] = b.fits(PY, LeY, fohmY)
    datadict["Y"]["rmsCMB"]["ssr"], datadict["Y"]["rmsCMB"]["m"], datadict["Y"]["rmsCMB"]["c"], datadict["Y"]["rmsCMB"]["res"] = b.fits(PY, LeYrmscmb, fohmY)
    datadict["Y"]["dipCMB"]["ssr"], datadict["Y"]["dipCMB"]["m"], datadict["Y"]["dipCMB"]["c"], datadict["Y"]["dipCMB"]["res"] = b.fits(PY, LeYdipcmb, fohmY)

    if chk == 1:
        # CHECK: Relate Elsasser to EM, total magnetic energy 
        # Em = Els*Els / (2*Pm*E)
        print('*******Yadav***************')
        EfromLe  = LeY**2 * 14.59 / 2.0 / EY**2
        Efromfile = EmPY + EmTY
        print('Energy difference (%) = ', np.abs( (EfromLe-Efromfile) / Efromfile)*100.0)
        print('***************************\n')

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

    print(datadict["S"]['rmsINT']['Le'])

    datadict["S"]['rmsINT']['ssr'], datadict["S"]['rmsINT']['m'], datadict["S"]['rmsINT']['c'], datadict["S"]['rmsINT']['res'] = b.fits(P_dum, datadict["S"]['rmsINT']['Le'], fohm_dum)
    datadict["S"]['rmsCMB']['ssr'], datadict["S"]['rmsCMB']['m'], datadict["S"]['rmsCMB']['c'], datadict["S"]['rmsCMB']['res'] = b.fits(P_dum, datadict["S"]['rmsCMB']['Le'], fohm_dum)
    datadict["S"]['dipCMB']['ssr'], datadict["S"]['dipCMB']['m'], datadict["S"]['dipCMB']['c'], datadict["S"]['dipCMB']['res'] = b.fits(P_dum, datadict["S"]['dipCMB']['Le'], fohm_dum)


# Construct datasets of all simulations
firstkey = True
for key in datadict:
    print("\nDataset : %s" %(key))
    if datadict[key]["plot"]:
        if firstkey:
            Rmall    = datadict[key]['d']["Rm"] 
            Eall     = datadict[key]['d']["E"]
            Leall    = datadict[key]["rmsINT"]["Le"]
            Lermscmb = datadict[key]["rmsCMB"]["Le"]
            Ledipcmb = datadict[key]["dipCMB"]["Le"]
            PAall    = datadict[key]['d']["p"]
            fohmall  = datadict[key]['d']["fohm"]
            bdipall  = datadict[key]['d']["bdip"]
            firstkey = False
        else:
            Rmall    = np.concatenate((Rmall,datadict[key]['d']["Rm"]))
            Eall     = np.concatenate((Eall,datadict[key]['d']["E"]))
            Leall    = np.concatenate((Leall,datadict[key]["rmsINT"]["Le"]))
            Lermscmb = np.concatenate((Lermscmb,datadict[key]["rmsCMB"]["Le"]))
            Ledipcmb = np.concatenate((Ledipcmb,datadict[key]["dipCMB"]["Le"]))
            PAall    = np.concatenate((PAall,datadict[key]['d']["p"]))
            fohmall  = np.concatenate((fohmall,datadict[key]['d']["fohm"]))
            bdipall  = np.concatenate((bdipall,datadict[key]['d']["bdip"]))
    else:
        print("Plot set to False.")

alldatadict["Rm"] = Rmall
alldatadict["E"] = Eall
alldatadict["rmsINT"]["Le"] = Leall
alldatadict["rmsCMB"]["Le"] = Lermscmb
alldatadict["dipCMB"]["Le"] = Ledipcmb 
alldatadict["p"] = PAall
alldatadict["fohm"] = fohmall
alldatadict["bdip"] = bdipall 

# Fit all data
xmin = 1e-14
xmax = 1e-3
Pfit = np.linspace(xmin, xmax, 10)
[mall , call]     ,  covall      = np.polyfit(np.log10(PAall), np.log10(Leall/fohmall**0.5)   , 1, cov=True, full=False)
[mrmscmb, crmscmb], covrmscmb    = np.polyfit(np.log10(PAall), np.log10(Lermscmb/fohmall**0.5), 1, cov=True, full=False)
[mdipcmb, cdipcmb], covdipcmb    = np.polyfit(np.log10(PAall), np.log10(Ledipcmb/fohmall**0.5), 1, cov=True, full=False)
fp,ssrall   ,rnk,sv,rc = np.polyfit(np.log10(PAall), np.log10(Leall/fohmall**0.5)   , 1, full=True)
fp,ssrrmscmb,rnk,sv,rc = np.polyfit(np.log10(PAall), np.log10(Lermscmb/fohmall**0.5), 1, full=True)
fp,ssrdipcmb,rnk,sv,rc = np.polyfit(np.log10(PAall), np.log10(Ledipcmb/fohmall**0.5), 1, full=True)
fitall    = 10**call    * Pfit**mall 
fitrmscmb = 10**crmscmb * Pfit**mrmscmb
fitdipcmb = 10**cdipcmb * Pfit**mdipcmb
resall    = np.sqrt(covall[0,0])
resrmscmb = np.sqrt(covrmscmb[0,0])
resdipcmb = np.sqrt(covdipcmb[0,0])
# Now store in dictionary
alldatadict["rmsINT"]["m"] = mall; alldatadict["rmsINT"]["res"] = resall; alldatadict["rmsINT"]["ssr"] = ssrall
alldatadict["rmsCMB"]["m"] = mrmscmb; alldatadict["rmsCMB"]["res"] = resrmscmb; alldatadict["rmsCMB"]["ssr"] = ssrrmscmb
alldatadict["dipCMB"]["m"] = mdipcmb; alldatadict["dipCMB"]["res"] = resdipcmb; alldatadict["dipCMB"]["ssr"] = ssrdipcmb

print('All     (slope, 10^c, std dev c, SSR) =', mall   , 10**call   , resall, ssrall)
print('RMS CMB (slope, 10^c, std dev c, SSR) =', mrmscmb, 10**crmscmb, resrmscmb, ssrrmscmb)
print('DIP CMB (slope, 10^c, std dev c, SSR) =', mdipcmb, 10**cdipcmb, resdipcmb, ssrdipcmb)

if calc_prefac_err:
    # - Error estimate on the prefactor for CMB dip field strength
    c_err_dipcmb = b.prefacError(PAall, Ledipcmb/fohmall**0.5, model=[10**cdipcmb,mdipcmb], plot=False, quiet=False)
    c_sd = c_err_dipcmb[1]
    fitdipcmb_1sdp = 10**(cdipcmb+c_sd) * Pfit**mdipcmb
    fitdipcmb_1sdm = 10**(cdipcmb-c_sd) * Pfit**mdipcmb
    fitdipcmb_3sdp = 10**(cdipcmb+3.*c_sd) * Pfit**mdipcmb
    fitdipcmb_3sdm = 10**(cdipcmb-3.*c_sd) * Pfit**mdipcmb
    fitdipcmb_5sdp = 10**(cdipcmb+5.*c_sd) * Pfit**mdipcmb
    fitdipcmb_5sdm = 10**(cdipcmb-5.*c_sd) * Pfit**mdipcmb

# - Define symbols' properties for plotting
Cmax = np.log10(1000.0)#np.log10(np.max(Eall))
Cmin = np.log10(50.0)#np.log10(np.min(Eall))

datadict = b.getPlotProperties(datadict)

# - Plot bdip stuff
if plt_bdip: b.plot_bdip(datadict, myfdip)


# --- Store fitted values in output files
b.saveFitValues(filename=outfpath+outfname, datadict=datadict, alldatadict=alldatadict)

##########################################################
# Extrapolated figures
##########################################################

earthdict["dipCMB"]["min"], earthdict["dipCMB"]["max"], earthdict["rmsINT"]["min"], earthdict["rmsINT"]["max"],\
    earthdict["rmsCMB"]["min"], earthdict["rmsCMB"]["max"], earthdict["p"]["min"], earthdict["p"]["max"] = b.getEarthEstimates(quiet=False)

# Calculate prefactors of Brms
c1_En = np.mean(np.log10(Leall/fohmall**0.5) - (1/3.)*np.log10(PAall))
fitEn = 10**c1_En * Pfit**(1/3.)
print('Core RMS, IMA   m=1/3 =0.33', c1_En)
c1_mac = np.mean(np.log10(Leall/fohmall**0.5) - 0.25*np.log10(PAall))
fitmac = 10**c1_mac * Pfit**0.25
print('Core RMS, MAC   m=1/4 =0.25', c1_mac)
c1_imac = np.mean(np.log10(Leall/fohmall**0.5) - 0.40*np.log10(PAall))
print('Core RMS, IMAC  m=2/5 =0.40', c1_imac)
fitimac = 10**c1_imac * Pfit**0.40
c1_imacd = np.mean(np.log10(Leall/fohmall**0.5) - 0.20*np.log10(PAall))
print('Core RMS, IMACd m=1/5 =0.20', c1_imacd)
fitimacd = 10**c1_imacd * Pfit**0.20
c1_imaci = np.mean(np.log10(Leall/fohmall**0.5) - 0.30*np.log10(PAall))
print('Core RMS, IMACi m=3/10=0.30', c1_imaci)
fitimaci = 10**c1_imaci * Pfit**0.30

# Error estimate on the prefactor
if calc_prefac_err:
    c1_err_En = b.prefacError(PAall, Leall/fohmall**0.5, model=[10**c1_En,1/3.], plot=False, quiet=False)
    c1_sd_En = c1_err_En[1]
    fitEn_1sdp = 10**(c1_En+c1_sd_En) * Pfit**(1/3.)
    fitEn_1sdm = 10**(c1_En-c1_sd_En) * Pfit**(1/3.)
    fitEn_3sdp = 10**(c1_En+3.*c1_sd_En) * Pfit**(1/3.)
    fitEn_3sdm = 10**(c1_En-3.*c1_sd_En) * Pfit**(1/3.)
    c1_err_mac = b.prefacError(PAall, Leall/fohmall**0.5, model=[10**c1_mac,0.25], plot=False, quiet=False)
    c1_sd_mac = c1_err_mac[1]
    fitmac_1sdp = 10**(c1_mac+c1_sd_mac) * Pfit**0.25
    fitmac_1sdm = 10**(c1_mac-c1_sd_mac) * Pfit**0.25
    fitmac_3sdp = 10**(c1_mac+3.*c1_sd_mac) * Pfit**0.25
    fitmac_3sdm = 10**(c1_mac-3.*c1_sd_mac) * Pfit**0.25

# Get plot properties
idxIMA = b.idxStr(plt_extrap_scalings, "IMA")[0]
idxMAC = b.idxStr(plt_extrap_scalings, "MAC")[0]

##########################
# Brms plots
##########################

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
                   label="$m=1/3$ (QG-MAC), $\sigma=$"+str(np.round(c1_sd_En,3)))
        plt.fill_between(Pfit, fitEn_1sdm, fitEn_1sdp, color=sc_fit[idxIMA], alpha=sc_alpha, zorder=-1)
    else:
        plt.loglog(Pfit, fitEn, c=lc_fit[idxIMA], linestyle=ls_fit[idxIMA], lw=lw_fit, label="$m=1/3$ (QG-MAC)")
if "MAC" in plt_extrap_scalings:
    if (calc_prefac_err):
        plt.loglog(Pfit, fitmac, c=lc_fit[idxMAC], linestyle=ls_fit[idxMAC], lw=lw_fit,
                   label="$m=1/4$ (MAC), $\sigma=$"+str(np.round(c1_sd_mac,3)))
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

title_str = b.getPlotTitle(myfdip=myfdip,
            myEr=myEr, Er_range=[EMoEK_range[0],EMoEK_range[1]])
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
c1_En_rmscmb = np.mean(np.log10(Lermscmb/fohmall**0.5) - (1/3.)*np.log10(PAall))
fitEn = 10**c1_En_rmscmb * Pfit**(1/3.)
print('CMB RMS, IMA   m=1/3 =0.33', c1_En_rmscmb)
c1_mac_rmscmb = np.mean(np.log10(Lermscmb/fohmall**0.5) - 0.25*np.log10(PAall))
fitmac = 10**c1_mac_rmscmb * Pfit**0.25
print('CMB RMS, MAC   m=1/4 =0.25', c1_mac_rmscmb)
c1_imac_rmscmb = np.mean(np.log10(Lermscmb/fohmall**0.5) - 0.40*np.log10(PAall))
print('CMB RMS, IMAC  m=2/5 =0.40', c1_imac_rmscmb)
fitimac = 10**c1_imac_rmscmb * Pfit**0.40
c1_imacd_rmscmb = np.mean(np.log10(Lermscmb/fohmall**0.5) - 0.20*np.log10(PAall))
print('CMB RMS, IMACd m=1/5 =0.20', c1_imacd_rmscmb)
fitimacd = 10**c1_imacd_rmscmb * Pfit**0.20
c1_imaci_rmscmb = np.mean(np.log10(Lermscmb/fohmall**0.5) - 0.30*np.log10(PAall))
print('CMB RMS, IMACi m=3/10=0.30', c1_imaci_rmscmb)
fitimaci = 10**c1_imaci_rmscmb * Pfit**0.30

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

if myfdip == 0:
    plt.title("All $f_{dip}$")
elif myfdip == 1: 
    plt.title("$f_{dip}>0.5$")
elif myfdip == 2:
    plt.title("$%.2f< f_{dip}< %.2f$" %(fdip_min,fdip_max))
elif myfdip == 3:
    plt.title("$%.2f< f_{dip}< %.2f$" %(fdip_min,fdip_max))

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
plt.savefig(file2, format='pdf',bbox_inches="tight")
plt.savefig(file3, format='png',bbox_inches="tight")

# - Bdip CMB extrapolated
c1_En_dipcmb = np.mean(np.log10(Ledipcmb/fohmall**0.5) - (1/3.)*np.log10(PAall))
fitEn = 10**c1_En_dipcmb * Pfit**(1/3.)
print('DIP CMB, IMA   m=1/3 =0.33', c1_En_dipcmb)
c1_mac_dipcmb = np.mean(np.log10(Ledipcmb/fohmall**0.5) - 0.25*np.log10(PAall))
fitmac = 10**c1_mac_dipcmb * Pfit**0.25
print('DIP CMB, MAC   m=1/4 =0.25', c1_mac_dipcmb)
c1_imac_dipcmb = np.mean(np.log10(Ledipcmb/fohmall**0.5) - 0.40*np.log10(PAall))
fitimac = 10**c1_imac_dipcmb * Pfit**0.40
print('DIP CMB, IMAC  m=2/5 =0.40', c1_imac_dipcmb)
c1_imacd_dipcmb = np.mean(np.log10(Ledipcmb/fohmall**0.5) - 0.20*np.log10(PAall))
fitimacd = 10**c1_imacd_dipcmb * Pfit**0.20
print('DIP CMB, IMACd m=1/5 =0.20', c1_imacd_dipcmb)
c1_imaci_dipcmb = np.mean(np.log10(Ledipcmb/fohmall**0.5) - 0.30*np.log10(PAall))
fitimaci = 10**c1_imaci_dipcmb * Pfit**0.30
print('DIP CMB, IMACi m=3/10=0.30', c1_imaci_dipcmb)

# Error estimate on the prefactor
if calc_prefac_err:
    c1_err_En_dipcmb = b.prefacError(PAall, Ledipcmb/fohmall**0.5, model=[10**c1_En_dipcmb,1/3.], plot=False, quiet=False)
    c1_sd_En_dipcmb = c1_err_En_dipcmb[1]
    fitEn_1sdp = 10**(c1_En_dipcmb+c1_sd_En_dipcmb) * Pfit**(1/3.)
    fitEn_1sdm = 10**(c1_En_dipcmb-c1_sd_En_dipcmb) * Pfit**(1/3.)
    fitEn_3sdp = 10**(c1_En_dipcmb+3.*c1_sd_En_dipcmb) * Pfit**(1/3.)
    fitEn_3sdm = 10**(c1_En_dipcmb-3.*c1_sd_En_dipcmb) * Pfit**(1/3.)
    fitEn_5sdp = 10**(c1_En_dipcmb+5.*c1_sd_En_dipcmb) * Pfit**(1/3.)
    fitEn_5sdm = 10**(c1_En_dipcmb-5.*c1_sd_En_dipcmb) * Pfit**(1/3.)
    c1_err_mac_dipcmb = b.prefacError(PAall, Ledipcmb/fohmall**0.5, model=[10**c1_mac_dipcmb,0.25], plot=False, quiet=False)
    c1_sd_mac_dipcmb = c1_err_mac_dipcmb[1]
    fitmac_1sdp = 10**(c1_mac_dipcmb+c1_sd_mac_dipcmb) * Pfit**0.25
    fitmac_1sdm = 10**(c1_mac_dipcmb-c1_sd_mac_dipcmb) * Pfit**0.25
    fitmac_3sdp = 10**(c1_mac_dipcmb+3.*c1_sd_mac_dipcmb) * Pfit**0.25
    fitmac_3sdm = 10**(c1_mac_dipcmb-3.*c1_sd_mac_dipcmb) * Pfit**0.25
    c1_err_imac_dipcmb = b.prefacError(PAall, Ledipcmb/fohmall**0.5, model=[10**c1_imac_dipcmb,0.40], plot=False, quiet=False)
    c1_sd_imac_dipcmb = c1_err_imac_dipcmb[1]
    c1_err_imacd_dipcmb = b.prefacError(PAall, Ledipcmb/fohmall**0.5, model=[10**c1_imacd_dipcmb,0.20], plot=False, quiet=False)
    c1_sd_imacd_dipcmb = c1_err_imacd_dipcmb[1]
    c1_err_imaci_dipcmb = b.prefacError(PAall, Ledipcmb/fohmall**0.5, model=[10**c1_imaci_dipcmb,0.30], plot=False, quiet=False)
    c1_sd_imaci_dipcmb = c1_err_imaci_dipcmb[1]

# Store pre-factors in output file
fout = open(outfpath+outfnamepf, "w")
delim = "    "
fout.write("--- B rms ---\n")
fout.write("Scaling     "+delim+"Prefactor c1\n")
fout.write("IMA    1/3  "+delim+str(round(c1_En,10))+"\n")
fout.write("MAC    1/4  "+delim+str(round(c1_mac,10))+"\n")
fout.write("IMAC   2/5  "+delim+str(round(c1_imac,10))+"\n")
fout.write("IMACd  1/5  "+delim+str(round(c1_imacd,10))+"\n")
fout.write("IMACi  3/10 "+delim+str(round(c1_imaci,10))+"\n")
fout.write("\n")
fout.write("--- B dip cmb ---\n")
if not calc_prefac_err:
    fout.write("Scaling     "+delim+"Prefactor c1\n")
    fout.write("IMA    1/3  "+delim+str(round(c1_En_dipcmb,10))+"\n")
    fout.write("MAC    1/4  "+delim+str(round(c1_mac_dipcmb,10))+"\n")
    fout.write("IMAC   2/5  "+delim+str(round(c1_imac_dipcmb,10))+"\n")
    fout.write("IMACd  1/5  "+delim+str(round(c1_imacd_dipcmb,10))+"\n")
    fout.write("IMACi  3/10 "+delim+str(round(c1_imaci_dipcmb,10))+"\n")
else:
    fout.write("Scaling     "+delim+"Prefactor c1"+delim+" Sigma\n")
    fout.write("IMA    1/3  "+delim+str(round(c1_En_dipcmb,10))+delim+str(round(c1_sd_En_dipcmb,10))+"\n")
    fout.write("MAC    1/4  "+delim+str(round(c1_mac_dipcmb,10))+delim+str(round(c1_sd_mac_dipcmb,10))+"\n")
    fout.write("IMAC   2/5  "+delim+str(round(c1_imac_dipcmb,10))+delim+str(round(c1_sd_imac_dipcmb,10))+"\n")
    fout.write("IMACd  1/5  "+delim+str(round(c1_imacd_dipcmb,10))+delim+str(round(c1_sd_imacd_dipcmb,10))+"\n")
    fout.write("IMACi  3/10 "+delim+str(round(c1_imaci_dipcmb,10))+delim+str(round(c1_sd_imaci_dipcmb,10))+"\n")
fout.write("\n")
fout.write("--- Check for VDM script on B dip cmb ---\n")
# Check for VDM script: Get Le for CMB dipole field at user-defined value of P
P0 = 1.e-12
fout.write("P0          "+delim+str(P0)+"\n")
fout.write("IMA    1/3  "+delim+str((10**c1_En_dipcmb)*P0**(1/3.))+"\n")
fout.write("MAC    1/4  "+delim+str((10**c1_mac_dipcmb)*P0**0.25)+"\n")
fout.write("IMAC   2/5  "+delim+str((10**c1_imac_dipcmb)*P0**0.40)+"\n")
fout.write("IMACd  1/5  "+delim+str((10**c1_imacd_dipcmb)*P0**0.20)+"\n")
fout.write("IMACi  3/10 "+delim+str((10**c1_imaci_dipcmb)*P0**0.30)+"\n")
fout.close()

# Plot simulations
ax, legend_xpos, legend_ypos = b.plotSimulations(datadict=datadict, alldatadict=alldatadict, earthdict=earthdict,
                                                 field="dipCMB",
                                                 cbarrange=[Cmin,Cmax], xrange=[xmin,xmax], yrange=[5e-6,0.2])
plt.loglog(Pfit, fitdipcmb, color="black", lw=lw_fit)
#           label="$m =$"+str(np.round(mdipcmb ,2))+"$\pm$"+str(np.round(resdipcmb ,2))+", SSR="+str(np.round(ssrdipcmb,2)))
#plt.loglog(Pfit, fitdipcmb_3sdp, c="k", ls=':', lw=1.)
#plt.loglog(Pfit, fitdipcmb_3sdm, c="k", ls=':', lw=1.)
#plt.loglog(Pfit, fitdipcmb_5sdp, c="k", ls=':', lw=1.)
#plt.loglog(Pfit, fitdipcmb_5sdm, c="k", ls=':', lw=1.)
if "IMAC" in plt_extrap_scalings:
    plt.loglog(Pfit, fitimac , color="dimgrey", linestyle="--", lw=lw_fit, label="$m=2/5$ (IMAC)")
if "IMACd" in plt_extrap_scalings:
    plt.loglog(Pfit, fitimacd, color="dimgrey", linestyle="-.", lw=lw_fit, label="$m=1/5$ (IMACd)")
if "IMACi" in plt_extrap_scalings:
    plt.loglog(Pfit, fitimaci, color="dimgrey", linestyle=":", lw=lw_fit, label="$m=3/10$ (IMACi)")
if "IMA" in plt_extrap_scalings:
    if (calc_prefac_err):
        plt.loglog(Pfit, fitEn,      c=lc_fit[idxIMA], linestyle=ls_fit[idxIMA], lw=lw_fit,
                   label="$m=1/3$ (QG-MAC), $\sigma=$"+str(np.round(c1_sd_En_dipcmb,3)))
        plt.fill_between(Pfit, fitEn_1sdm, fitEn_1sdp, color=sc_fit[idxIMA], alpha=sc_alpha, zorder=-1)
    else:
        plt.loglog(Pfit, fitEn, color=lc_fit[idxIMA], linestyle=ls_fit[idxIMA], lw=lw_fit,
                   label="$m=1/3$ (QG-MAC)")
if "MAC" in plt_extrap_scalings:
    if (calc_prefac_err):
        plt.loglog(Pfit, fitmac, c=lc_fit[idxMAC], linestyle=ls_fit[idxMAC], lw=lw_fit,
                   label="$m=1/4$ (MAC), $\sigma=$"+str(np.round(c1_sd_mac_dipcmb,3)))
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

title_str = b.getPlotTitle(myfdip=myfdip,
            myEr=myEr, Er_range=[EMoEK_range[0],EMoEK_range[1]])
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

plt.clf()
plt.figure(figsize=(16,15))
plt.subplot(3,1,1)
ax  = plt.gca()
plt.xlim([1e-10,1e-3])
plt.ylim([1e-4,0.2])

for key in datadict:
    if datadict[key]["plot"]:
        plt.scatter(datadict[key]['d']["p"], datadict[key]["rmsINT"]["Le"]/datadict[key]['d']["fohm"]**0.5,
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
        plt.text(legend_xpos, legend_ypos-iplt*legend_dy,
                 "$m$  = "+str(np.round(datadict[key]["rmsINT"]["m"],2)) +\
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
        plt.scatter(datadict[key]['d']["p"], datadict[key]["rmsCMB"]["Le"]/datadict[key]['d']["fohm"]**0.5,
                        s=datadict[key]["plotp"]["size"], marker=datadict[key]["plotp"]["marker"],
                        c=np.log10(datadict[key]["plotp"]["Col"]), vmin=Cmin, vmax=Cmax,
                        cmap=plt.get_cmap(datadict[key]["plotp"]["cmap"]), edgecolor=datadict[key]["plotp"]["edgecolor"], label=datadict[key]["plotp"]["label"]) 
plt.loglog(Pfit,fitrmscmb,color="black")
plt.loglog(5e-12 , 2e-4, marker="s", markersize=30)

legend_ypos = 0.95; legend_dy   = 0.06
iplt = 0
for key in datadict:
    if datadict[key]["plot"]:
        plt.text(legend_xpos, legend_ypos-iplt*legend_dy,
                 "$m$  = "+str(np.round(datadict[key]["rmsCMB"]["m"],2)) +\
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
        plt.scatter(datadict[key]['d']["p"], datadict[key]["dipCMB"]["Le"]/datadict[key]['d']["fohm"]**0.5,
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
        plt.text(legend_xpos, legend_ypos-iplt*legend_dy,
                 "$m$  = "+str(np.round(datadict[key]["dipCMB"]["m"],2))+\
                 "$\pm$"+str(np.round(datadict[key]["dipCMB"]["res"],2))+\
                 ", SSR="+str(np.round(datadict[key]["dipCMB"]["ssr"],2)),
                 transform=ax.transAxes, color=datadict[key]["plotp"]["edgecolor"])
        iplt += 1
plt.text(0.05, 0.61, "$m$ = "+str(np.round(alldatadict["dipCMB"]["m"],2))+\
                              "$\pm$"+str(np.round(alldatadict["dipCMB"]["res"],2))+\
                              ", SSR="+str(np.round(alldatadict["dipCMB"]["ssr"],2)), transform=ax.transAxes, color='black')
plt.text(0.05, 0.56, "$\sigma$ = "+str(np.round(c_sd,4)), transform=ax.transAxes, color='black')
plt.xlabel('$P_A$')
if myfohm == 0: 
    plt.ylabel('$Le_{cmb}^{l=1}/f_{ohm}^{1/2}$')
elif myfohm == 1: 
    plt.ylabel('$Le_{cmb}^{l=1}$')
plt.legend(bbox_to_anchor=(0.05, 0.70), loc=3, borderaxespad=0)

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
