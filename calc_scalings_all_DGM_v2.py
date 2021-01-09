import numpy as np
import pandas as pd
import bscaling_functions as b
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size':'12'})
rc('text', usetex=False)


# ------------------------
# --- Input parameters ---
# ------------------------
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

# ------------------------
# ------------------------

# Logical values for plotting datasets
datadict = {"L":True, "A":True, "Y":True, "UC":True, "UCt":True, "A17":True}

leedsname = "./all_LED_tave_NEW.csv"
yadavname = "./Yadav_etal_pnas_16_data.xls_all"
aubertname= "./aubert2009-all.txt"
christname= "./dynq_mycode.data"
christnamet="./dyntq_mycode.data"

if myfdip == 0:
    fdipn = "0"

    if (myEkOPm == 1 or myEr == 1):
        filetag = "_"
        if (myEkOPm == 1 and myEr != 1):
            filetag += "EkOPm"
        elif (myEkOPm != 1 and myEr == 1):
            filetag += "EMoEK"
        elif (myEkOPm == 1 and myEr == 1):
            filetag += "EkOPm_EMoEK"

        leedsOutName   = "./all_LED_tave_NEW.csv"+filetag
        yadavOutName   = "./Yadav_etal_pnas_16_data.xls"+filetag
        aubertOutName  = "./aubert2009.txt"+filetag
        christOutName  = "./dynq_mycode.data"+filetag
        christOutNameT = "./dyntq_mycode.data"+filetag
        df, datadict = b.filter_table(infname="./all_LED_tave_NEW.csv",
                                    outfname=leedsOutName, dataset="Leeds",fdip_range=[0.,1.1], EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
        df, datadict = b.filter_table(infname="./Yadav_etal_pnas_16_data.xls_all",
                     outfname=yadavOutName, dataset="Yadav", fdip_range=[0.,1.1], EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
        df, datadict = b.filter_table(infname="./aubert2009-all.txt",
                     outfname=aubertOutName, dataset="Aubert", fdip_range=[0.,1.1], EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
        df, datadict = b.filter_table(infname="./dynq_mycode.data",
                     outfname=christOutName, dataset="Christensen", fdip_range=[0.,1.1], EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
        df, datadict = b.filter_table(infname="./dyntq_mycode.data",
                     outfname=christOutNameT, dataset="ChristensenT", fdip_range=[0.,1.1], EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
        leedsname = leedsOutName
        yadavname = yadavOutName
        aubertname = aubertOutName
        christname = christOutName
        christnamet = christOutNameT

elif myfdip == 1: 
    fdipn = "1"
    leedsname = "./all_LED_tave_NEW_fdip50.csv"
    yadavname = "./Yadav_etal_pnas_16_data.xls_dipolar"
    aubertname= "./aubert2009050.txt"
    christname= "./dynq_mycode.data_fdip50"
    christnamet="./dyntq_mycode.data_fdip50"
    if (myEkOPm == 1 or myEr == 1):
        filetag = "_"
        if (myEkOPm == 1 and myEr != 1):
            filetag += "EkOPm"
        elif (myEkOPm != 1 and myEr == 1):
            filetag += "EMoEK"
        elif (myEkOPm == 1 and myEr == 1):
            filetag += "EkOPm_EMoEK"
        leedsOutName   = "./all_LED_tave_NEW_fdip50.csv"+filetag
        yadavOutName   = "./Yadav_etal_pnas_16_data.xls_dipolar"+filetag
        aubertOutName  = "./aubert2009050.txt"+filetag
        christOutName  = "./dynq_mycode.data_fdip50"+filetag
        christOutNameT = "./dyntq_mycode.data_fdip50"+filetag
        df, datadict = b.filter_table(infname="./all_LED_tave_NEW_fdip50.csv",         
                     outfname=leedsOutName, dataset="Leeds", fdip_range=[0.,1.1], EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
        df, datadict = b.filter_table(infname="./Yadav_etal_pnas_16_data.xls_dipolar",
                     outfname=yadavOutName, dataset="Yadav", fdip_range=[0.,1.1], EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
        df, datadict = b.filter_table(infname="./aubert2009050.txt",
                     outfname=aubertOutName, dataset="Aubert", fdip_range=[0.,1.1], EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
        df, datadict = b.filter_table(infname="./dynq_mycode.data_fdip50",
                     outfname=christOutName, dataset="Christensen", fdip_range=[0.,1.1], EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
        df, datadict = b.filter_table(infname="./dyntq_mycode.data_fdip50", 
                     outfname=christOutNameT, dataset="ChristensenT", fdip_range=[0.,1.1], EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
        leedsname = leedsOutName
        yadavname = yadavOutName
        aubertname = aubertOutName
        christname = christOutName
        christnamet = christOutNameT

elif myfdip == 2:
    fdipn = "2"
    fdip_min = 0.35; fdip_max = 0.80
    leedsOutName   = "./all_LED_tave_NEW.csv_fdip035-080"
    yadavOutName   = "./Yadav_etal_pnas_16_data.xls_fdip035-080"
    aubertOutName  = "./aubert2009.txt_fdip035-080"
    christOutName  = "./dynq_mycode.data_fdip035-080"
    christOutNameT = "./dyntq_mycode.data_fdip035-080"
    if (myEkOPm == 1 or myEr == 1):
        filetag = "_"
        if (myEkOPm == 1 and myEr != 1):
            filetag += "EkOPm"
        elif (myEkOPm != 1 and myEr == 1):
            filetag += "EMoEK"
        elif (myEkOPm == 1 and myEr == 1):
            filetag += "EkOPm_EMoEK"   
        leedsOutName   += filetag
        yadavOutName   += filetag
        aubertOutName  += filetag
        christOutName  += filetag
        christOutNameT += filetag
        df, datadict = b.filter_table(infname="./all_LED_tave_NEW.csv",
                     outfname=leedsOutName, dataset="Leeds", fdip_range=[fdip_min,fdip_max], EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
        df, datadict = b.filter_table(infname="./Yadav_etal_pnas_16_data.xls_all",
                     outfname=yadavOutName, dataset="Yadav", fdip_range=[fdip_min,fdip_max], EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
        df, datadict = b.filter_table(infname="./aubert2009-all.txt",
                     outfname=aubertOutName, dataset="Aubert", fdip_range=[fdip_min,fdip_max], EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
        df, datadict = b.filter_table(infname="./dynq_mycode.data",
                     outfname=christOutName, dataset="Christensen", fdip_range=[fdip_min,fdip_max], EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
        df, datadict = b.filter_table(infname="./dyntq_mycode.data",
                     outfname=christOutNameT, dataset="ChristensenT", fdip_range=[fdip_min,fdip_max], EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
    else:
        df, datadict = b.filter_table(infname="./all_LED_tave_NEW.csv",            outfname=leedsOutName,   dataset="Leeds",        fdip_range=[fdip_min,fdip_max], datadict=datadict)
        df, datadict = b.filter_table(infname="./Yadav_etal_pnas_16_data.xls_all", outfname=yadavOutName,   dataset="Yadav",        fdip_range=[fdip_min,fdip_max], datadict=datadict)
        df, datadict = b.filter_table(infname="./aubert2009-all.txt",              outfname=aubertOutName,  dataset="Aubert",       fdip_range=[fdip_min,fdip_max], datadict=datadict)
        df, datadict = b.filter_table(infname="./dynq_mycode.data",                outfname=christOutName,  dataset="Christensen",  fdip_range=[fdip_min,fdip_max], datadict=datadict)
        df, datadict = b.filter_table(infname="./dyntq_mycode.data",               outfname=christOutNameT, dataset="ChristensenT", fdip_range=[fdip_min,fdip_max], datadict=datadict)
        EkOPm_range = None
        EMoEK_range = None
    leedsname = leedsOutName
    yadavname = yadavOutName
    aubertname = aubertOutName
    christname = christOutName
    christnamet = christOutNameT

elif myfdip == 3:
    fdipn = "3"
    fdip_min = 0.40; fdip_max = 0.80
    leedsOutName   = "./all_LED_tave_NEW.csv_fdip040-080"
    yadavOutName   = "./Yadav_etal_pnas_16_data.xls_fdip040-080"
    aubertOutName  = "./aubert2009.txt_fdip040-080"
    christOutName  = "./dynq_mycode.data_fdip040-080"
    christOutNameT = "./dyntq_mycode.data_fdip040-080"
    if (myEkOPm == 1 or myEr == 1):
        filetag = "_"
        if (myEkOPm == 1 and myEr != 1):
            filetag += "EkOPm"
        elif (myEkOPm != 1 and myEr == 1):
            filetag += "EMoEK"
        elif (myEkOPm == 1 and myEr == 1):
            filetag += "EkOPm_EMoEK"   
        leedsOutName   += filetag
        yadavOutName   += filetag
        aubertOutName  += filetag
        christOutName  += filetag
        christOutNameT += filetag
        df, datadict = b.filter_table(infname="./all_LED_tave_NEW.csv",
                     outfname=leedsOutName, dataset="Leeds", fdip_range=[fdip_min,fdip_max], EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
        df, datadict = b.filter_table(infname="./Yadav_etal_pnas_16_data.xls_all",
                     outfname=yadavOutName, dataset="Yadav", fdip_range=[fdip_min,fdip_max], EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
        df, datadict = b.filter_table(infname="./aubert2009-all.txt",
                     outfname=aubertOutName, dataset="Aubert", fdip_range=[fdip_min,fdip_max], EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
        df, datadict = b.filter_table(infname="./dynq_mycode.data", 
                     outfname=christOutName, dataset="Christensen", fdip_range=[fdip_min,fdip_max], EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
        df, datadict = b.filter_table(infname="./dyntq_mycode.data",
                     outfname=christOutNameT, dataset="ChristensenT", fdip_range=[fdip_min,fdip_max], EkOPm_range=EkOPm_range, EMoEK_range=EMoEK_range, datadict=datadict)
    else:
        df, datadict = b.filter_table(infname="./all_LED_tave_NEW.csv", outfname=leedsOutName, dataset="Leeds", fdip_range=[fdip_min,fdip_max], datadict=datadict)
        df, datadict = b.filter_table(infname="./Yadav_etal_pnas_16_data.xls_all", outfname=yadavOutName, dataset="Yadav", fdip_range=[fdip_min,fdip_max], datadict=datadict)
        df, datadict = b.filter_table(infname="./aubert2009-all.txt", outfname=aubertOutName, dataset="Aubert", fdip_range=[fdip_min,fdip_max], datadict=datadict)
        df, datadict = b.filter_table(infname="./dynq_mycode.data", outfname=christOutName, dataset="Christensen", fdip_range=[fdip_min,fdip_max], datadict=datadict)
        df, datadict = b.filter_table(infname="./dyntq_mycode.data", outfname=christOutNameT, dataset="ChristensenT", fdip_range=[fdip_min,fdip_max], datadict=datadict)
        EkOPm_range = None
        EMoEK_range = None
    leedsname = leedsOutName
    yadavname = yadavOutName
    aubertname = aubertOutName
    christname = christOutName
    christnamet = christOutNameT

else:
    raise ValueError('Not valid value for myfdip')

# Get fdip values and store in separate variables
if datadict["L"]:
    fdipC   = b.get_fdip(infname=leedsname, dataset="Leeds")
if datadict["Y"]:
    fdipY   = b.get_fdip(infname=yadavname, dataset="Yadav")
if datadict["A"]:
    fdipA   = b.get_fdip(infname=aubertname, dataset="Aubert")
if datadict["UC"]:
    fdipUC  = b.get_fdip(infname=christname, dataset="Christensen")
if datadict["UCt"]:
    fdipUCt = b.get_fdip(infname=christnamet, dataset="ChristensenT")

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

# Get Lhuillier values
#t0, bdip035, fdip035, PA035, fohmA035 = np.loadtxt("lhuillier2019", usecols=(6,8,9,12,13), skiprows=0, unpack='true')
#LeA   , bdipall, fdipall, PA   , fohmA    = np.loadtxt("aubert2009-all.txt", usecols=(7,8,9,12,13), skiprows=0, unpack='true')
#LeLh = 0.01*t0

# Christensen
# Ek has no factor 2
# Fundamental length scale is shell thickness, time scale is viscous diffusion time, magnetic field scale is "Elsasser scale".
# NOTE that variables names "Els" are actually RMS field strengths!!
if datadict["UC"]:
    EUC,PmUC,ElsUC,ElsrmscmbUC,ElsdipcmbUC,PUC,fohmUC,RmUC,EmagUC = np.loadtxt(christname, usecols=(0,4,26,27,29,31,32,10,19), skiprows=2, unpack='true')
    LeUC      = np.sqrt( (EUC/PmUC) ) * ElsUC        # Eqn 14 of CA06
    LeUCrmscmb= np.sqrt( (EUC/PmUC) ) * ElsrmscmbUC
    LeUCdipcmb= np.sqrt( (EUC/PmUC) ) * ElsdipcmbUC
    PUC       = 1e7 * PUC * EUC**3
    fohmUC    = fohmUC / 100
    if myfohm == 1:
        fohmUC = np.ones(len(fohmUC)) 
    bdipUC    = LeUC/LeUCdipcmb
    nsimsUC   = len(PUC)
    ssrUC      , mUC      , cUC      , resUC       = b.fits(PUC, LeUC      , fohmUC)
    ssrUCrmscmb, mUCrmscmb, cUCrmscmb, resUCrmscmb = b.fits(PUC, LeUCrmscmb, fohmUC)
    ssrUCdipcmb, mUCdipcmb, cUCdipcmb, resUCdipcmb = b.fits(PUC, LeUCdipcmb, fohmUC)

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

if datadict["UCt"]:
    EUCt,PmUCt,ElsUCt,ElsrmscmbUCt,ElsdipcmbUCt,PUCt,fohmUCt,RmUCt,EmagUCt = np.loadtxt(christnamet, usecols=(0,4,27,28,30,32,33,11,20), skiprows=1, unpack='true')
    LeUCt      = np.sqrt((EUCt/PmUCt)) * ElsUCt
    LeUCtrmscmb= np.sqrt((EUCt/PmUCt)) * ElsrmscmbUCt
    LeUCtdipcmb= np.sqrt((EUCt/PmUCt)) * ElsdipcmbUCt
    PUCt       = 1e7 * PUCt * EUCt**3
    fohmUCt    = fohmUCt / 100
    if myfohm == 1:
        fohmUCt = np.ones(len(fohmUCt))
    bdipUCt    = LeUCt/LeUCtdipcmb
    nsimsUCt   = len(PUCt)
    ssrUCt      , mUCt      , cUCt      , resUCt       = b.fits(PUCt, LeUCt      , fohmUCt)
    ssrUCtrmscmb, mUCtrmscmb, cUCtrmscmb, resUCtrmscmb = b.fits(PUCt, LeUCtrmscmb, fohmUCt)
    ssrUCtdipcmb, mUCtdipcmb, cUCtdipcmb, resUCtdipcmb = b.fits(PUCt, LeUCtdipcmb, fohmUCt)

    if chk == 1:
        # CHECK: Relate Elsasser to Em, magnetic energy density (i.e. per unit volume)
        # Em = Els / (2*Pm*E) = Le**2 / (2 E**2) 
        print('*******Christensen qT******')
        EfromLe  = LeUCt**2 * PmUCt / 2.0 / EUCt**2 / PmUCt
        Efromfile = EmagUCt
        print('Energy difference (%) = ', np.abs( (EfromLe-Efromfile) / Efromfile)*100.0)
        print('***************************\n')

# Leeds dataset
if datadict["L"]:
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
        plt.savefig('./eardme/fig/Check_Emag_Gauss.pdf',format='pdf')
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
    nsimsC   = len(PC)
    ssrC      , mC      , cC      , resC       = b.fits(PC, LeC      , fohmC)
    ssrCrmscmb, mCrmscmb, cCrmscmb, resCrmscmb = b.fits(PC, LeCrmscmb, fohmC)
    ssrCdipcmb, mCdipcmb, cCdipcmb, resCdipcmb = b.fits(PC, LeCdipcmb, fohmC)

# Yadav - NOTE - he uses fixed T. Extracts E, Pm, Elasser, Elsasser_CMB, Dip_CMB, Buo_pow, Ohm_diss
if datadict["Y"]:
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
    nsimsY   = len(PY)
    ssrY      , mY      , cY      , resY       = b.fits(PY, LeY      , fohmY)
    ssrYrmscmb, mYrmscmb, cYrmscmb, resYrmscmb = b.fits(PY, LeYrmscmb, fohmY)
    ssrYdipcmb, mYdipcmb, cYdipcmb, resYdipcmb = b.fits(PY, LeYdipcmb, fohmY)

    if chk == 1:
        # CHECK: Relate Elsasser to EM, total magnetic energy 
        # Em = Els*Els / (2*Pm*E)
        print('*******Yadav***************')
        EfromLe  = LeY**2 * 14.59 / 2.0 / EY**2
        Efromfile = EmPY + EmTY
        print('Energy difference (%) = ', np.abs( (EfromLe-Efromfile) / Efromfile)*100.0)
        print('***************************\n')

# Get Aubert values
# NB - Aubert uses different radii! chi dependence of parameters?
if datadict["A"]:
    EA   , LeA   ,bdipA,fdipA,PA   ,fohmA, dAall,RoA,PmA   = np.loadtxt(aubertname, usecols=(0,7,8,9,12,13,10,6,3), skiprows=1, unpack='true')
    RmA = RoA / EA * PmA
    if myfohm == 1:
        fohmA = np.ones(len(fohmA))
    LeAdipcmb    = LeA   /bdipA
    LeArmscmb    = LeA   /(bdipA*fdipA)
    lAall = 2.0*np.pi / (dAall + 0.5) # Jeans' formula.
    nsimsA = len(PA)
    ssrA      , mA      , cA      , resA       = b.fits(PA, LeA      , fohmA)
    ssrArmscmb, mArmscmb, cArmscmb, resArmscmb = b.fits(PA, LeArmscmb, fohmA)
    ssrAdipcmb, mAdipcmb, cAdipcmb, resAdipcmb = b.fits(PA, LeAdipcmb, fohmA)
    print('cA = ', cA, cArmscmb, cAdipcmb)

# Construct datasets
if not datadict["Y"]:
    raise ValueError('Y dataset is emply. Choose other initial dataset.')
Rmall   = RmY
Eall    = EY
Leall   = LeY
Lermscmb= LeYrmscmb
Ledipcmb= LeYdipcmb
PAall   = PY
fohmall = fohmY
bdipall = bdipY
if datadict["A"]:
    Rmall   = np.concatenate((Rmall,RmA))
    Eall    = np.concatenate((Eall,EA))
    Leall   = np.concatenate((Leall,LeA))
    Lermscmb= np.concatenate((Lermscmb,LeArmscmb))
    Ledipcmb= np.concatenate((Ledipcmb,LeAdipcmb))
    PAall   = np.concatenate((PAall,PA))
    fohmall = np.concatenate((fohmall,fohmA))
    bdipall = np.concatenate((bdipall,bdipA))
if datadict["L"]:
    Rmall   = np.concatenate((Rmall,RmC))
    Eall    = np.concatenate((Eall,EC))
    Leall   = np.concatenate((Leall,LeC))
    Lermscmb= np.concatenate((Lermscmb,LeCrmscmb))
    Ledipcmb= np.concatenate((Ledipcmb,LeCdipcmb))
    PAall   = np.concatenate((PAall,PC))
    fohmall = np.concatenate((fohmall,fohmC))
    bdipall = np.concatenate((bdipall,bdipC))
if datadict["UC"]:
    Rmall   = np.concatenate((Rmall,RmUC))
    Eall    = np.concatenate((Eall,EUC))
    Leall   = np.concatenate((Leall,LeUC))
    Lermscmb= np.concatenate((Lermscmb,LeUCrmscmb))
    Ledipcmb= np.concatenate((Ledipcmb,LeUCdipcmb))
    PAall   = np.concatenate((PAall,PUC))
    fohmall = np.concatenate((fohmall,fohmUC))
    bdipall = np.concatenate((bdipall,bdipUC))
if datadict["UCt"]:
    Rmall   = np.concatenate((Rmall,RmUCt))
    Eall    = np.concatenate((Eall,EUCt))
    Leall   = np.concatenate((Leall,LeUCt))
    Lermscmb= np.concatenate((Lermscmb,LeUCtrmscmb))
    Ledipcmb= np.concatenate((Ledipcmb,LeUCtdipcmb))
    PAall   = np.concatenate((PAall,PUCt))
    fohmall = np.concatenate((fohmall,fohmUCt))
    bdipall = np.concatenate((bdipall,bdipUCt))
#Rmall   = np.concatenate((RmA      , RmC      , RmY      , RmUC      , RmUCt))
#Eall    = np.concatenate((EA       , EC       , EY       , EUC       , EUCt))
#Leall   = np.concatenate((LeA      , LeC      , LeY      , LeUC      , LeUCt))
#Lermscmb= np.concatenate((LeArmscmb, LeCrmscmb, LeYrmscmb, LeUCrmscmb, LeUCtrmscmb))
#Ledipcmb= np.concatenate((LeAdipcmb, LeCdipcmb, LeYdipcmb, LeUCdipcmb, LeUCtdipcmb))
#PAall   = np.concatenate((PA       , PC       , PY       , PUC       , PUCt))
#fohmall = np.concatenate((fohmA    , fohmC    , fohmY    , fohmUC    , fohmUCt))
#bdipall = np.concatenate((bdipA    , bdipC    , bdipY    , bdipUC    , bdipUCt))
# Fits
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

# - Define symbols' colours for plotting
Cmax = np.log10(1000.0)#np.log10(np.max(Eall))
Cmin = np.log10(50.0)#np.log10(np.min(Eall))
if datadict["A"]:
    ColA = RmA / PmA
if datadict["L"]:
    ColC = RmC / PmC
if datadict["Y"]:
    ColY = RmY / PmY
if datadict["UC"]:
    ColU = RmUC/ PmUC
if datadict["UCt"]:
    ColUt= RmUCt/ PmUCt

# - Plot bdip stuff
if plt_bdip:
    """
    # - bdip vs fdip
    fig = plt.figure(figsize=(8,10))
    ax = fig.add_subplot(311)
    plt.scatter(fdipA, bdipA, s=150, marker="^", c=np.log10(ColA), vmin=Cmin, vmax=Cmax,
                cmap=plt.get_cmap("Oranges"), edgecolor='orange',label="Aubert et al (2009)")
    plt.scatter(fdipY, bdipY, s=150, marker="o", c=np.log10(ColY), vmin=Cmin, vmax=Cmax,
                cmap=plt.get_cmap("Reds"),edgecolor='red'   ,label="Yadav et al (2016)")
    plt.scatter(fdipC, bdipC, s=150, marker="*", c=np.log10(ColC), vmin=Cmin, vmax=Cmax,
                cmap=plt.get_cmap("Blues")  ,edgecolor='blue'  ,label="Leeds")
    plt.scatter(fdipUC, bdipUC, s=150, marker="v", c=np.log10(ColU), vmin=Cmin, vmax=Cmax,
                cmap=plt.get_cmap("Purples"),edgecolor='purple',label="Christensen 0F")
    plt.scatter(fdipUCt, bdipUCt, s=150, marker="^", c=np.log10(ColUt), vmin=Cmin, vmax=Cmax,
                cmap=plt.get_cmap("Purples"),edgecolor='purple',label="Christensen FF")
    ax.set_xlabel('$f_{dip}$')
    ax.set_ylabel('$b_{dip}$')
    plt.xlim([0.,1.]); plt.ylim([0.,200.])

    ax = fig.add_subplot(312)
    plt.scatter(fdipA, bdipA, s=150, marker="^", c=np.log10(ColA), vmin=Cmin, vmax=Cmax,
                cmap=plt.get_cmap("Oranges"), edgecolor='orange',label="Aubert et al (2009)")
    plt.scatter(fdipY, bdipY, s=150, marker="o", c=np.log10(ColY), vmin=Cmin, vmax=Cmax,
                cmap=plt.get_cmap("Reds"),edgecolor='red'   ,label="Yadav et al (2016)")
    plt.scatter(fdipC, bdipC, s=150, marker="*", c=np.log10(ColC), vmin=Cmin, vmax=Cmax,
                cmap=plt.get_cmap("Blues")  ,edgecolor='blue'  ,label="Leeds")
    plt.scatter(fdipUC, bdipUC, s=150, marker="v", c=np.log10(ColU), vmin=Cmin, vmax=Cmax,
                cmap=plt.get_cmap("Purples"),edgecolor='purple',label="Christensen 0F")
    plt.scatter(fdipUCt, bdipUCt, s=150, marker="^", c=np.log10(ColUt), vmin=Cmin, vmax=Cmax,
                cmap=plt.get_cmap("Purples"),edgecolor='purple',label="Christensen FF")
    ax.set_xlabel('$f_{dip}$')
    ax.set_ylabel('$b_{dip}$')
    plt.xlim([0.5,1.]); plt.ylim([0.,30.])
    ax = fig.add_subplot(313)
    plt.scatter(fdipA, bdipA, s=150, marker="^", c=np.log10(ColA), vmin=Cmin, vmax=Cmax,
                cmap=plt.get_cmap("Oranges"), edgecolor='orange',label="Aubert et al (2009)")
    plt.scatter(fdipY, bdipY, s=150, marker="o", c=np.log10(ColY), vmin=Cmin, vmax=Cmax,
                cmap=plt.get_cmap("Reds"),edgecolor='red'   ,label="Yadav et al (2016)")
    plt.scatter(fdipC, bdipC, s=150, marker="*", c=np.log10(ColC), vmin=Cmin, vmax=Cmax,
                cmap=plt.get_cmap("Blues")  ,edgecolor='blue'  ,label="Leeds")
    plt.scatter(fdipUC, bdipUC, s=150, marker="v", c=np.log10(ColU), vmin=Cmin, vmax=Cmax,
                cmap=plt.get_cmap("Purples"),edgecolor='purple',label="Christensen 0F")
    plt.scatter(fdipUCt, bdipUCt, s=150, marker="^", c=np.log10(ColUt), vmin=Cmin, vmax=Cmax,
                cmap=plt.get_cmap("Purples"),edgecolor='purple',label="Christensen FF")
    ax.set_xlabel('$f_{dip}$')
    ax.set_ylabel('$b_{dip}$')
    plt.xlim([0.4,0.8]); plt.ylim([0.,70.])
 
    plt.show(block=False)
    plt.tight_layout()
    plt.savefig('./eardme/fig/bdip_vs_fdip_fdipN='+fdipn+'.pdf',format='pdf')
    del ax
    """
    # - bdip vs buoyancy power
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)
    plt.scatter(PA, bdipA, s=150, marker="^", c=np.log10(ColA), vmin=Cmin, vmax=Cmax,
                cmap=plt.get_cmap("Oranges"), edgecolor='orange',label="Aubert et al (2009)")
    plt.scatter(PY, bdipY, s=150, marker="o", c=np.log10(ColY), vmin=Cmin, vmax=Cmax,
                cmap=plt.get_cmap("Reds"),edgecolor='red'   ,label="Yadav et al (2016)")
    plt.scatter(PC, bdipC, s=150, marker="*", c=np.log10(ColC), vmin=Cmin, vmax=Cmax,
                cmap=plt.get_cmap("Blues")  ,edgecolor='blue'  ,label="Leeds")
    plt.scatter(PUC, bdipUC, s=150, marker="v", c=np.log10(ColU), vmin=Cmin, vmax=Cmax,
                cmap=plt.get_cmap("Purples"),edgecolor='purple',label="Christensen 0F")
    plt.scatter(PUCt, bdipUCt, s=150, marker="^", c=np.log10(ColUt), vmin=Cmin, vmax=Cmax,
                cmap=plt.get_cmap("Purples"),edgecolor='purple',label="Christensen FF")
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
    plt.savefig('./eardme/fig/bdip_vs_P_fdip='+fdipn+'.pdf',format='pdf')
    del ax

# --- Store fitted slopes and pre-factors in output files
# Slopes file
fout = open(outfpath+outfname, "w")
delim = "     "
fout.write("--- B rms ---\n")
fout.write("Dataset            "+delim+"slope m    "+"std. err. "+" SSR\n")
if datadict["A"]:
    fout.write("Aubert et al (2009)"+delim+str(np.round(mA,4))+delim+str(np.round(resA,4))+delim+str(np.round(ssrA[0],4))+"\n")
if datadict["Y"]:
    fout.write("Yadav et al  (2016)"+delim+str(np.round(mY,4))+delim+str(np.round(resY,4))+delim+str(np.round(ssrY[0],4))+"\n")
if datadict["L"]:
    fout.write("Leeds              "+delim+str(np.round(mC,4))+delim+str(np.round(resC,4))+delim+str(np.round(ssrC[0],4))+"\n")
if datadict["UC"]:
    fout.write("Christensen        "+delim+str(np.round(mUC,4))+delim+str(np.round(resUC,4))+delim+str(np.round(ssrUC[0],4))+"\n")
if datadict["UCt"]:
    fout.write("Christensen T      "+delim+str(np.round(mUCt,4))+delim+str(np.round(resUCt,4))+delim+str(np.round(ssrUCt[0],4))+"\n")
fout.write("All                "+delim+str(np.round(mall,4))+delim+str(np.round(resall,4))+delim+str(np.round(ssrall[0],4))+"\n")
fout.write("\n")
fout.write("--- B rms cmb ---\n")
fout.write("Dataset            "+delim+"slope m    "+"std. err. "+" SSR\n")
if datadict["A"]:
    fout.write("Aubert et al (2009)"+delim+str(np.round(mArmscmb,4))+delim+str(np.round(resArmscmb,4))+delim+str(np.round(ssrArmscmb[0],4))+"\n")
if datadict["Y"]:
    fout.write("Yadav et al  (2016)"+delim+str(np.round(mYrmscmb,4))+delim+str(np.round(resYrmscmb,4))+delim+str(np.round(ssrYrmscmb[0],4))+"\n")
if datadict["L"]:
    fout.write("Leeds              "+delim+str(np.round(mCrmscmb,4))+delim+str(np.round(resCrmscmb,4))+delim+str(np.round(ssrCrmscmb[0],4))+"\n")
if datadict["UC"]:
    fout.write("Christensen        "+delim+str(np.round(mUCrmscmb,4))+delim+str(np.round(resUCrmscmb,4))+delim+str(np.round(ssrUCrmscmb[0],4))+"\n")
if datadict["UCt"]:
    fout.write("Christensen T      "+delim+str(np.round(mUCtrmscmb,4))+delim+str(np.round(resUCtrmscmb,4))+delim+str(np.round(ssrUCtrmscmb[0],4))+"\n")
fout.write("All                "+delim+str(np.round(mrmscmb,4))+delim+str(np.round(resrmscmb,4))+delim+str(np.round(ssrrmscmb[0],4))+"\n")
fout.write("\n")
fout.write("--- B dip cmb ---\n")
fout.write("Dataset            "+delim+"slope m    "+"std. err. "+" SSR\n")
if datadict["A"]:
    fout.write("Aubert et al (2009)"+delim+str(np.round(mAdipcmb,4))+delim+str(np.round(resAdipcmb,4))+delim+str(np.round(ssrAdipcmb[0],4))+"\n")
if datadict["Y"]:
    fout.write("Yadav et al  (2016)"+delim+str(np.round(mYdipcmb,4))+delim+str(np.round(resYdipcmb,4))+delim+str(np.round(ssrYdipcmb[0],4))+"\n")
if datadict["L"]:
    fout.write("Leeds              "+delim+str(np.round(mCdipcmb,4))+delim+str(np.round(resCdipcmb,4))+delim+str(np.round(ssrCdipcmb[0],4))+"\n")
if datadict["UC"]:
    fout.write("Christensen        "+delim+str(np.round(mUCdipcmb,4))+delim+str(np.round(resUCdipcmb,4))+delim+str(np.round(ssrUCdipcmb[0],4))+"\n")
if datadict["UCt"]:
    fout.write("Christensen T      "+delim+str(np.round(mUCtdipcmb,4))+delim+str(np.round(resUCtdipcmb,4))+delim+str(np.round(ssrUCtdipcmb[0],4))+"\n")
fout.write("All                "+delim+str(np.round(mdipcmb,4))+delim+str(np.round(resdipcmb,4))+delim+str(np.round(ssrdipcmb[0],4))+"\n")
fout.close()


##########################################################
# Extrapolated figures
##########################################################

# Brmax = (a/r)^3 g10
ri    = 1221e3
rc    = 3480e3
rsurf = 6371e3
omega = 7.272e-5
dd    = rc-ri
rho   = 1e4
mu0   = 4.*np.pi*1.e-7 
Bfac  = np.sqrt(mu0*rho)*omega*dd # Denominator of Lenhert number

# Min and max CMB dipole field strengths for Earth: 20,000 nT and 40,000 nT
Le_earth_dipmin = 20000e-9 * np.sqrt(2.) * (rsurf/rc)**3 / (np.sqrt(rho*4*np.pi*1e-7)*omega*dd)
Le_earth_dipmax = 40000e-9 * np.sqrt(2.) * (rsurf/rc)**3 / (np.sqrt(rho*4*np.pi*1e-7)*omega*dd)
# Min and max core field strengths for Earth: 1 mT and 10 mT
Le_earth_rmsmin = 1e-3  / (np.sqrt(rho*4*np.pi*1e-7)*omega*dd)
Le_earth_rmsmax = 10e-3 / (np.sqrt(rho*4*np.pi*1e-7)*omega*dd)
# Min and max CMB field strengths for Earth: 0.5 mT and 1 mT
Le_earth_rmscmbmin = 0.5e-3 / (np.sqrt(rho*4*np.pi*1e-7)*omega*dd)
Le_earth_rmscmbmax = 1e-3 / (np.sqrt(rho*4*np.pi*1e-7)*omega*dd)
print('Le_earth CMB dipole min/max  = ', Le_earth_dipmin, Le_earth_dipmax)
print('Le_earth Surf dipole min/max = ', Le_earth_dipmin/(rsurf/rc)**3, Le_earth_dipmax/(rsurf/rc)**3)
print('Le_earth RMS total  min/max  = ', Le_earth_rmsmin, Le_earth_rmsmax)
print('Le_earth RMS CMB  min/max    = ', Le_earth_rmscmbmin, Le_earth_rmscmbmax)
# Min buoyancy pow for Earth: 0.1 TW
P_earth_min = 1e11 * 3. / (4*np.pi*(rc**3-ri**3)) / rho / omega**3 / dd**2
# Max buoyancy pow for Earth: 5 TW
P_earth_max = 5e12 * 3. / (4*np.pi*(rc**3-ri**3)) / rho / omega**3 / dd**2
print('Power min/max = ', P_earth_min, P_earth_max)

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

# Brms plots
plt.clf()
ax  = plt.gca()
plt.xlim([xmin,xmax])
plt.ylim([3e-5,0.2])
if datadict["A"]:
    plt.scatter(PA  ,LeA/fohmA**0.5  ,s=150,marker="^",c=np.log10(ColA),vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Oranges"),edgecolor='orange',label="Aubert et al (2009)")
if datadict["Y"]:
    plt.scatter(PY  ,LeY/fohmY**0.5  ,s=150,marker="o",c=np.log10(ColY),vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Reds")   ,edgecolor='red'   ,label="Yadav et al (2016)")
if datadict["L"]:
    plt.scatter(PC  ,LeC/fohmC**0.5  ,s=150,marker="*",c=np.log10(ColC),vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Blues")  ,edgecolor='blue'  ,label="Leeds")
if datadict["UC"]:
    plt.scatter(PUC ,LeUC /fohmUC**0.5, s=150,marker="v",c=np.log10(ColU),vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Purples"),edgecolor='purple',label="Christensen 0F")
if datadict["UCt"]:
    plt.scatter(PUCt,LeUCt/fohmUCt**0.5,s=150,marker="^",c=np.log10(ColUt),vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Purples"),edgecolor='purple',label="Christensen FF")
plt.loglog(Pfit, fitall, color="black", lw=lw_fit)
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


ax.add_patch(Rectangle(xy=(P_earth_min,Le_earth_rmsmin) ,width=(P_earth_max-P_earth_min), height=(Le_earth_rmsmax-Le_earth_rmsmin),
                       linewidth=1, color='black', fill=True))

legend_xpos = 0.01 # x-position of legend (>1, outside of main plot)
if datadict["A"]:
    plt.text(legend_xpos, 0.95, "$m$  = "+str(np.round(mA,2)) +"$\pm$"+str(np.round(resA,2)) +", SSR="+str(np.round(ssrA[0],2)), transform=ax.transAxes, color='orange')
if datadict["Y"]:
    plt.text(legend_xpos, 0.90, "$m$  = "+str(np.round(mY,2)) +"$\pm$"+str(np.round(resY,2)) +", SSR="+str(np.round(ssrY[0],2)), transform=ax.transAxes, color='red')
if datadict["L"]:
    plt.text(legend_xpos, 0.85, "$m$  = "+str(np.round(mC,2)) +"$\pm$"+str(np.round(resC,2)) +", SSR="+str(np.round(ssrC[0],2)), transform=ax.transAxes, color='blue')
if datadict["UC"]:
    plt.text(legend_xpos, 0.80, "$m$v = "+str(np.round(mUC,2))+"$\pm$"+str(np.round(resUC,2))+", SSR="+str(np.round(ssrUC[0],2)), transform=ax.transAxes, color='purple')
if datadict["UCt"]:
    plt.text(legend_xpos, 0.75, "$m$^ = "+str(np.round(mUCt,2))+"$\pm$"+str(np.round(resUCt,2))+", SSR="+str(np.round(ssrUCt[0],2)), transform=ax.transAxes, color='purple')
plt.text(legend_xpos, 0.70, "$m$  = "+str(np.round(mall,2))+"$\pm$"+str(np.round(resall,2))+", SSR="+str(np.round(ssrall[0],2)), transform=ax.transAxes, color='black')
cbar = plt.colorbar()
cbar.set_label("log $Re$")
plt.xlabel('$p_A$')
if myfohm == 0: 
    plt.ylabel('$Le^{\\rm rms}_t/f_{ohm}^{1/2}$')
elif myfohm == 1: 
    plt.ylabel('$Le^{\\rm rms}_t$')
ax.set_yscale('log')
ax.set_xscale('log')
plt.legend(bbox_to_anchor=(legend_xpos-0.1, 1.10), loc=3, ncol=2, borderaxespad=0)
plt.rcParams["figure.figsize"] = [15,10]

title_str = b.getPlotTitle(myfdip=myfdip,
            myEr=myEr, Er_range=[EMoEK_range[0],EMoEK_range[1]])
plt.title(title_str)

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

# Brms cmb plots
plt.clf()
ax  = plt.gca()
plt.xlim([xmin,xmax])
plt.ylim([5e-6,0.2])
if datadict["A"]:
    plt.scatter(PA  ,LeArmscmb/fohmA**0.5  ,s=150,marker="^",c=np.log10(ColA),vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Oranges"),
                edgecolor='orange',label="Aubert et al (2009)")
if datadict["Y"]:
    plt.scatter(PY  ,LeYrmscmb/fohmY**0.5  ,s=150,marker="o",c=np.log10(ColY),vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Reds"),
                edgecolor='red'   ,label="Yadav et al (2016)")
if datadict["L"]:
    plt.scatter(PC  ,LeCrmscmb/fohmC**0.5  ,s=150,marker="*",c=np.log10(ColC),vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Blues"),
                edgecolor='blue'  ,label="Leeds")
if datadict["UC"]:
    plt.scatter(PUC ,LeUCrmscmb /fohmUC**0.5, s=150,marker="v",c=np.log10(ColU),vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Purples"),
                edgecolor='purple',label="Christensen 0F")
if datadict["UCt"]:
    plt.scatter(PUCt,LeUCtrmscmb/fohmUCt**0.5,s=150,marker="^",c=np.log10(ColUt),vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Purples"),
                edgecolor='purple',label="Christensen FF")
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
ax.add_patch(Rectangle(xy=(P_earth_min,Le_earth_rmscmbmin),
             width=(P_earth_max-P_earth_min), height=(Le_earth_rmscmbmax-Le_earth_rmscmbmin),
             linewidth=1, color='black', fill=True))

legend_xpos = 0.01 # x-position of legend (>1, outside of main plot)
if datadict["A"]:
    plt.text(legend_xpos, 0.95, "$m$  = "+str(np.round(mArmscmb,2)) +"$\pm$"+str(np.round(resArmscmb,2)) +", SSR="+str(np.round(ssrArmscmb[0],2)),
             transform=ax.transAxes, color='orange')
if datadict["Y"]:
    plt.text(legend_xpos, 0.90, "$m$  = "+str(np.round(mYrmscmb,2)) +"$\pm$"+str(np.round(resYrmscmb,2)) +", SSR="+str(np.round(ssrYrmscmb[0],2)),
             transform=ax.transAxes, color='red')
if datadict["L"]:
    plt.text(legend_xpos, 0.85, "$m$  = "+str(np.round(mCrmscmb,2)) +"$\pm$"+str(np.round(resCrmscmb,2)) +", SSR="+str(np.round(ssrCrmscmb[0],2)),
             transform=ax.transAxes, color='blue')
if datadict["UC"]:
    plt.text(legend_xpos, 0.80, "$m$v = "+str(np.round(mUCrmscmb,2))+"$\pm$"+str(np.round(resUCrmscmb,2))+", SSR="+str(np.round(ssrUCrmscmb[0],2)),
             transform=ax.transAxes, color='purple')
if datadict["UCt"]:
    plt.text(legend_xpos, 0.75, "$m$^ = "+str(np.round(mUCtrmscmb,2))+"$\pm$"+str(np.round(resUCtrmscmb,2))+", SSR="+str(np.round(ssrUCtrmscmb[0],2)),
             transform=ax.transAxes, color='purple')
plt.text(legend_xpos, 0.70, "$m$  = "+str(np.round(mrmscmb,2))+"$\pm$"+str(np.round(resrmscmb,2))+", SSR="+str(np.round(ssrrmscmb[0],2)),
         transform=ax.transAxes, color='black')
cbar = plt.colorbar()
cbar.set_label("log $Re$")
plt.xlabel('$p_A$')
if myfohm == 0: 
    plt.ylabel('$Le^{\\rm rms}_{\\rm cmb}/f_{ohm}^{1/2}$')
elif myfohm == 1: 
    plt.ylabel('$Le^{\\rm rms}_{\\rm cmb}$')
ax.set_yscale('log')
ax.set_xscale('log')
plt.legend(bbox_to_anchor=(legend_xpos-0.1, 1.10), loc=3, ncol=2, borderaxespad=0)
plt.rcParams["figure.figsize"] = [15,10]
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

plt.clf()
ax  = plt.gca()
plt.xlim([xmin,xmax])
plt.ylim([5e-6,0.05])
if datadict["A"]:
    plt.scatter(PA  ,LeAdipcmb/fohmA**0.5  ,s=150,marker="^",c=np.log10(ColA) ,vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Oranges"),
                edgecolor='orange',label="Aubert et al (2009)")
if datadict["Y"]:
    plt.scatter(PY  ,LeYdipcmb/fohmY**0.5  ,s=150,marker="o",c=np.log10(ColY) ,vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Reds"),
                edgecolor='red'   ,label="Yadav et al (2016)")
if datadict["L"]:
    plt.scatter(PC  ,LeCdipcmb/fohmC**0.5  ,s=150,marker="*",c=np.log10(ColC) ,vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Blues"),
                edgecolor='blue'  ,label="Leeds")
if datadict["UC"]:
    plt.scatter(PUC ,LeUCdipcmb/fohmUC**0.5,s=150,marker="v",c=np.log10(ColU),vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Purples"),
                edgecolor='purple',label="Christensen 0F")
if datadict["UCt"]:
    plt.scatter(PUCt,LeUCtdipcmb/fohmUCt**0.5,s=150,marker="^",c=np.log10(ColUt),vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Purples"),
                edgecolor='purple',label="Christensen FF")
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

ax.add_patch(Rectangle(xy=(P_earth_min,Le_earth_dipmin) ,width=(P_earth_max-P_earth_min), height=(Le_earth_dipmax-Le_earth_dipmin),
             linewidth=1, color='black', fill=True))

legend_xpos = 0.01 # x-position of legend (>1, outside of main plot)
if datadict["A"]:
    plt.text(legend_xpos, 0.95, "$m$  = "+str(np.round(mAdipcmb,2)) +"$\pm$"+str(np.round(resAdipcmb,2)) +", SSR="+str(np.round(ssrAdipcmb[0],2)),
             transform=ax.transAxes, color='orange')
if datadict["Y"]:
    plt.text(legend_xpos, 0.90, "$m$  = "+str(np.round(mYdipcmb,2)) +"$\pm$"+str(np.round(resYdipcmb,2)) +", SSR="+str(np.round(ssrYdipcmb[0],2)),
             transform=ax.transAxes, color='red')
if datadict["L"]:
    plt.text(legend_xpos, 0.85, "$m$  = "+str(np.round(mCdipcmb,2)) +"$\pm$"+str(np.round(resCdipcmb,2)) +", SSR="+str(np.round(ssrCdipcmb[0],2)),
             transform=ax.transAxes, color='blue')
if datadict["UC"]:
    plt.text(legend_xpos, 0.80, "$m$v = "+str(np.round(mUCdipcmb,2))+"$\pm$"+str(np.round(resUCdipcmb,2))+", SSR="+str(np.round(ssrUCdipcmb[0],2)),
             transform=ax.transAxes, color='purple')
if datadict["UCt"]:
    plt.text(legend_xpos, 0.75, "$m$^ = "+str(np.round(mUCtdipcmb,2))+"$\pm$"+str(np.round(resUCtdipcmb,2))+", SSR="+str(np.round(ssrUCtdipcmb[0],2)),
             transform=ax.transAxes, color='purple')
plt.text(legend_xpos, 0.70, "$m$  = "+str(np.round(mdipcmb,2))  +"$\pm$"+str(np.round(resdipcmb,2))  +", SSR="+str(np.round(ssrdipcmb[0],2)),
         transform=ax.transAxes, color='black')
cbar = plt.colorbar()
cbar.set_label("log $Re$")
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
if datadict["A"]:
    ax.scatter(PA  ,LeA/fohmA**0.5  ,s=150,marker="^",c=np.log10(ColA) ,vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Oranges"),edgecolor='orange',label="Aubert et al (2009)")
if datadict["Y"]:
    ax.scatter(PY  ,LeY/fohmY**0.5  ,s=150,marker="o",c=np.log10(ColY) ,vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Reds")   ,edgecolor='red'   ,label="Yadav et al (2016)")
if datadict["L"]:
    ax.scatter(PC  ,LeC/fohmC**0.5  ,s=150,marker="*",c=np.log10(ColC) ,vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Blues")  ,edgecolor='blue'  ,label="Leeds")
if datadict["UC"]:
    ax.scatter(PUC ,LeUC/fohmUC**0.5,s=150,marker="v",c=np.log10(ColU),vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Purples"),edgecolor='purple',label="Christensen 0F")
if datadict["UCt"]:
    ax.scatter(PUCt,LeUCt/fohmUCt**0.5,s=150,marker="^",c=np.log10(ColUt),vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Purples"),edgecolor='purple',label="Christensen FF")
plt.loglog(Pfit , fitall, color="black")
plt.loglog(5e-12 , 2e-4, marker="s", markersize=30)
if datadict["A"]:
    plt.text(0.05, 0.90, "$m$  = "+str(np.round(mA,2)) +"$\pm$"+str(np.round(resA,2)) +", SSR="+str(np.round(ssrA[0],2)), transform=ax.transAxes, color='orange')
if datadict["Y"]:
    plt.text(0.05, 0.84, "$m$  = "+str(np.round(mY,2)) +"$\pm$"+str(np.round(resY,2)) +", SSR="+str(np.round(ssrY[0],2)), transform=ax.transAxes, color='red')
if datadict["L"]:
    plt.text(0.05, 0.78, "$m$  = "+str(np.round(mC,2)) +"$\pm$"+str(np.round(resC,2)) +", SSR="+str(np.round(ssrC[0],2)), transform=ax.transAxes, color='blue')
if datadict["UC"]:
    plt.text(0.05, 0.72, "$m$v = "+str(np.round(mUC,2))+"$\pm$"+str(np.round(resUC,2))+", SSR="+str(np.round(ssrUC[0],2)), transform=ax.transAxes, color='purple')
if datadict["UCt"]:
    plt.text(0.05, 0.66, "$m$^ = "+str(np.round(mUCt,2))+"$\pm$"+str(np.round(resUCt,2))+", SSR="+str(np.round(ssrUCt[0],2)), transform=ax.transAxes, color='purple')
plt.text(0.05, 0.61, "$m$  = "+str(np.round(mall,2))+"$\pm$"+str(np.round(resall,2))+", SSR="+str(np.round(ssrall[0],2)), transform=ax.transAxes, color='black')
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
if datadict["A"]:
    ax.scatter(PA  ,LeArmscmb/fohmA**0.5  ,s=150,marker="^",c=np.log10(ColA) ,vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Oranges"),edgecolor='orange')
if datadict["Y"]:
    ax.scatter(PY  ,LeYrmscmb/fohmY**0.5  ,s=150,marker="o",c=np.log10(ColY) ,vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Reds")   ,edgecolor='red'   )
if datadict["L"]:
    ax.scatter(PC  ,LeCrmscmb/fohmC**0.5  ,s=150,marker="*",c=np.log10(ColC) ,vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Blues")  ,edgecolor='blue'  )
if datadict["UC"]:
    ax.scatter(PUC ,LeUCrmscmb/fohmUC**0.5,s=150,marker="v",c=np.log10(ColU),vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Purples"),edgecolor='purple')
if datadict["UCt"]:
    ax.scatter(PUCt,LeUCtrmscmb/fohmUCt**0.5,s=150,marker="^",c=np.log10(ColUt),vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Purples"),edgecolor='purple')
plt.loglog(Pfit,fitrmscmb,color="black")
plt.loglog(5e-12 , 2e-4, marker="s", markersize=30)
if datadict["A"]:
    plt.text(0.05, 0.90, "$m$  = "+str(np.round(mArmscmb,2)) +"$\pm$"+str(np.round(resArmscmb,2)) +", SSR="+str(np.round(ssrArmscmb[0],2)) , transform=ax.transAxes, color='orange')
if datadict["Y"]:
    plt.text(0.05, 0.84, "$m$  = "+str(np.round(mYrmscmb,2)) +"$\pm$"+str(np.round(resYrmscmb,2)) +", SSR="+str(np.round(ssrYrmscmb[0],2)) , transform=ax.transAxes, color='red')
if datadict["L"]:
    plt.text(0.05, 0.78, "$m$  = "+str(np.round(mCrmscmb,2)) +"$\pm$"+str(np.round(resCrmscmb,2)) +", SSR="+str(np.round(ssrCrmscmb[0],2)) , transform=ax.transAxes, color='blue')
if datadict["UC"]:
    plt.text(0.05, 0.72, "$m$v = "+str(np.round(mUCrmscmb,2))+"$\pm$"+str(np.round(resUCrmscmb,2))+", SSR="+str(np.round(ssrUCrmscmb[0],2)), transform=ax.transAxes, color='purple')
if datadict["UCt"]:
    plt.text(0.05, 0.66, "$m$^ = "+str(np.round(mUCtrmscmb,2))+"$\pm$"+str(np.round(resUCtrmscmb,2))+", SSR="+str(np.round(ssrUCtrmscmb[0],2)), transform=ax.transAxes, color='purple')
plt.text(0.05, 0.62, "$m$  = "+str(np.round(mrmscmb,2))  +"$\pm$"+str(np.round(resrmscmb,2))  +", SSR="+str(np.round(ssrrmscmb[0],2))  , transform=ax.transAxes, color='black')
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
if datadict["A"]:
    ax.scatter(PA  ,LeAdipcmb/fohmA**0.5  ,s=150,marker="^",c=np.log10(ColA) ,vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Oranges"),edgecolor='orange')
if datadict["Y"]:
    ax.scatter(PY  ,LeYdipcmb/fohmY**0.5  ,s=150,marker="o",c=np.log10(ColY) ,vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Reds")   ,edgecolor='red'   )
if datadict["L"]:
    ax.scatter(PC  ,LeCdipcmb/fohmC**0.5  ,s=150,marker="*",c=np.log10(ColC) ,vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Blues")  ,edgecolor='blue'  )
if datadict["UC"]:
    ax.scatter(PUC ,LeUCdipcmb/fohmUC**0.5,s=150,marker="v",c=np.log10(ColU),vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Purples"),edgecolor='purple')
if datadict["UCt"]:
    ax.scatter(PUCt,LeUCtdipcmb/fohmUCt**0.5,s=150,marker="^",c=np.log10(ColUt),vmin=Cmin,vmax=Cmax,cmap=plt.get_cmap("Purples"),edgecolor='purple')
plt.loglog(Pfit,fitdipcmb,color="black")
if (calc_prefac_err):
    plt.loglog(Pfit,fitdipcmb_1sdp,c='k',ls=':',lw=1.)
    plt.loglog(Pfit,fitdipcmb_1sdm,c='k',ls=':',lw=1.)
    #plt.loglog(Pfit,fitdipcmb_3sdp,c='k',ls='--',lw=1.)
    #plt.loglog(Pfit,fitdipcmb_3sdm,c='k',ls='--',lw=1.)
if datadict["A"]:
    plt.text(0.05, 0.90, "$m$ = "+str(np.round(mAdipcmb,2)) +"$\pm$"+str(np.round(resAdipcmb,2)) +", SSR="+str(np.round(ssrAdipcmb[0],2)),
             transform=ax.transAxes, color='orange')
if datadict["Y"]:
    plt.text(0.05, 0.84, "$m$ = "+str(np.round(mYdipcmb,2)) +"$\pm$"+str(np.round(resYdipcmb,2)) +", SSR="+str(np.round(ssrYdipcmb[0],2)),
             transform=ax.transAxes, color='red')
if datadict["L"]:
    plt.text(0.05, 0.78, "$m$ = "+str(np.round(mCdipcmb,2)) +"$\pm$"+str(np.round(resCdipcmb,2)) +", SSR="+str(np.round(ssrCdipcmb[0],2)),
             transform=ax.transAxes, color='blue')
if datadict["UC"]:
    plt.text(0.05, 0.72, "$m$v = "+str(np.round(mUCdipcmb,2))+"$\pm$"+str(np.round(resUCdipcmb,2))+", SSR="+str(np.round(ssrUCdipcmb[0],2)),
             transform=ax.transAxes, color='purple')
if datadict["UCt"]:
    plt.text(0.05, 0.66, "$m$^ = "+str(np.round(mUCtdipcmb,2))+"$\pm$"+str(np.round(resUCtdipcmb,2))+", SSR="+str(np.round(ssrUCtdipcmb[0],2)),
             transform=ax.transAxes, color='purple')
plt.text(0.05, 0.61, "$m$ = "+str(np.round(mdipcmb,2))  +"$\pm$"+str(np.round(resdipcmb,2))  +", SSR="+str(np.round(ssrdipcmb[0],2)),
         transform=ax.transAxes, color='black')
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
