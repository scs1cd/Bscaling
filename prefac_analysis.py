import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size':'22'})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

def readPrefac(fname):
    # --- This function reads file for prefactors (output of dynamo scaling script)
    f = open(fname, "r")
    i = 0
    for line in f:
        if (i==10):
            row = line.split()
            prefac_IMA       = float(row[2])
            prefac_IMA_sigma = float(row[3])
        elif (i==11):
            row = line.split()
            prefac_MAC = float(row[2])
            prefac_MAC_sigma = float(row[3])
        i += 1
    return prefac_IMA, prefac_IMA_sigma, prefac_MAC, prefac_MAC_sigma

def readSlope(fname):
    # --- This function reads file for slopes (output of dynamo scaling script)
    f = open(fname, "r")
    i = 0
    for line in f:
        if (i==7):
            row = line.split()
            mrmsall    = float(row[1])
            merrrsmall = float(row[2])
            SSRrmsall  = float(row[3])
        if (i==16):
            row = line.split()
            mrmsrmscmb = float(row[1])
            merrrmscmb = float(row[2])
            SSRrmscmb  = float(row[3])
        if (i==25):
            row = line.split()
            mdipcmb    = float(row[1])
            merrdipcmb = float(row[2])
            SSRdipcmb  = float(row[3])
        i += 1
    return mrmsall, merrrsmall, SSRrmsall, mrmsrmscmb, merrrmscmb, SSRrmscmb, mdipcmb, merrdipcmb, SSRdipcmb


fdipr = ["0","1","2"]
EmoKEr= ["1.0e-19", "1.0e+00","2.0e+00","3.0e+00","4.0e+00","5.0e+00"]
color = ["black", "blue","purple", "orange", "gold", "red"]
sym   = ["s", 'o', '^']
names=["Dataset","m", "merr", "SSR"]

fig1, (ax1a, ax1b) = plt.subplots(1, 2, figsize=(20,6))
fig2, (ax2a, ax2b) = plt.subplots(1, 2, figsize=(20,6))
fig3,         ax3  = plt.subplots(figsize=(12,6))
fig4, (ax4a, ax4b) = plt.subplots(1, 2, figsize=(20,6))

prefac_array1 = np.zeros(len(EmoKEr)*len(fdipr))
prefac_array2 = np.zeros(len(EmoKEr)*len(fdipr))
prefac_sig1 = np.zeros(len(EmoKEr)*len(fdipr))

i = 0 
for fi, fdip in enumerate(fdipr):
    for e, EMoKE in enumerate(EmoKEr):
        slopefile ="scalingSlope_fdip="+fdip+"_fohm=0_EMoEK="+EMoKE+"_1.0e+19.txt"
        prefacfile="scalingPrefac_fdip="+fdip+"_fohm=0_EMoEK="+EMoKE+"_1.0e+19.txt"

        slope  = readSlope(slopefile)
        prefac = readPrefac(prefacfile)

        prefac_array1[i] = 10**prefac[0]
        prefac_array2[i] = 10**prefac[2]
        prefac_sig1[i]   = 10**(prefac[0]+prefac[1])

        #prefac_array1[i] = prefac[0]
        #prefac_array2[i] = prefac[2]
        #print(slopefile,slope)
        print(prefacfile, prefac)

        if e == 0:
            x = ["0", "0"]
        else:
            x = EMoKE.split("e")

        # mSSR
        if fi == 0:
            ax1a.scatter(slope[0], slope[2], color=color[e], label="$E_M/E_K > $"+str(x[0]), s=200, marker=sym[fi])
            ax1b.scatter(slope[6], slope[8], color=color[e], label="$E_M/E_K > $"+str(x[0]), s=200, marker=sym[fi])
        else:
            ax1a.scatter(slope[0], slope[2], color=color[e], s=200, marker=sym[fi])
            ax1b.scatter(slope[6], slope[8], color=color[e], s=200, marker=sym[fi])       

        ax3.scatter(prefac_array1[i], prefac_array2[i], color=color[e], label="ME/KE="+str(EMoKE))

        # cSSR
        ax4a.scatter(prefac_array1[i], slope[2], color=color[e], label="$E_M/E_K > $"+str(EMoKE), s=200, marker=sym[fi])
        ax4b.scatter(prefac_array1[i], slope[8], color=color[e], label="$E_M/E_K > $"+str(EMoKE), s=200, marker=sym[fi])

        i = i + 1

#        df = pd.read_csv('scalingSlope_test.txt',sep='\s+',skiprows=[0,1,8,9,10,17,18,19],header=0,
#            names=('c1','c2','c3','c4'),error_bad_lines=False)


print("Mean c1 m=1/3 = ", np.mean(prefac_array1))
print("Max = ", np.max(prefac_array1), " Min = ", np.min(prefac_array1), 
    " Max-Min = ",(np.max(np.abs(prefac_array1)) - np.min(np.abs(prefac_array1))))
print("Mean c1 m=1/4 = ", np.mean(prefac_array2))
print( (np.max(np.abs(prefac_array1)) - np.min(np.abs(prefac_array1)))/np.mean(np.abs(prefac_array1)) )

np.savetxt("prefac.txt", [prefac_array1, prefac_array2])

ax1a.set_title("RMS Internal Field")
ax1a.set_xlim([0.28, 0.335])
ax1a.set_ylim([0.00, 12.0])
ax1a.set_xlabel("$m$")
ax1a.set_ylabel("SSR")
ax1b.set_title("CMB Dipole Field")
ax1b.set_ylim([0.00, 12.0])
ax1b.set_xlim([0.28, 0.335])
ax1b.set_xlabel("$m$")
ax1a.legend(loc=2, handlelen=0.05, ncol=2)
fig1.savefig('mSSR.pdf',format='pdf')

ax2a.set_title("$m=1/3$")
ax2a.set_xlabel("$c$")
ax2a.set_ylabel("Occurrence")
ax2b.set_xlabel("$c_1$")
ax2b.set_title("$m=1/4$")
ax2a.hist(prefac_array1)
ax2b.hist(prefac_array2)
fig2.savefig('hist.pdf',format='pdf')

ax3.set_xlabel("$c_1, m=1/3$")
ax3.set_ylabel("$c_1, m=1/4$")
fig3.savefig('c1.pdf',format='pdf')

ax4a.set_title("RMS Internal Field")
#ax4a.set_xlim([0.29, 0.33])
ax4a.set_ylim([0.00, 12.0])
ax4a.set_xlabel("$c_1$")
ax4a.set_ylabel("SSR")
ax4b.set_title("CMB Dipole Field")
ax4b.set_ylim([0.00, 12.0])
#ax4b.set_xlim([0.29, 0.33])
ax4b.set_xlabel("$c_1$")
fig4.savefig('cSSR.pdf',format='pdf')

