# Bscaling
This code analyses scaling laws for the RMS internal field, RMS CMB field and CMB dipole field. It consists of a main script (`calc_scalings_all_DGM_v2.py`) which calls various functions defined in `bscaling_functions.py`.
## Input parameters
In the main script, the following parameters can be set:
* calc_prefac_err     : (logical) If True, calculate and plot the scaling prefactor error (\sigma) as a shaded region. The prefactor error \sigma is defined as in Aubert et al GJI 2009, sec. 2.3.
* myfdip              : (integer) Use 0 for all fdip values, 1 for fdip > 0.50, 2 for filtering with fdip=(0.35,0.80), 3 for fdip=(0.40,0.80)
* myfohm              : (integer) Use 0 to apply fohm correction factor, or 1 for NO fohm factor
* myEkOPm             : Use 1 (0) for (not) filtering in a specified range of Ek/Pm values. The range of values is defined by the variable `EkOPm_range`.
* myEr                : Use 1 (0) for (not) filtering in a specified range of magnetic to kinetic energy ratio (EM/EK). The range of values is defined by the variable `EMoEK_range`.

* plt_extrap_scalings : (list of strings) List of theoretical scalings to be plotted in extrapolated figures. Possible values are "IMA","MAC","IMAC","IMACd","IMACi".
* lc_fit              : (list of strings) Line colour for each theoretical scaling.
* sc_fit              : (list of strings) Shading colour for the \sigma intervals.
* sc_alpha            : (float) alpha value for the shading (from 0 to 1. 1 is opaque.)
* ls_fit              : (list of strings) Line style for each theoretical scaling.
* lw_fit              : (float) Line width of the theoretical scaling.
* plt_bdip            : (logical) Plot bdip (ratio of CMb field strength to internal field strength) as function of buoyancy power. 

* chk                 : (integer) Use 1 to print on screen checks of energy (quiet otherwise).
* check_gauss_Leeds   : (logical) If True, perform internal consistency check of Leeds simulations. It makes sure that the CMB (total and dipole) field strengths calculated from the mag_cmb files and from the Gauss coeffs files are the same.

* categorise          : (logical) If True, categorise and plots the simulations by driving. If False, plots are by authors. The categories are: FTFT, FFFF, FF0F, FTFF, Mixed, CE. ***If categorise=True, make sure all the keys `plot` in the simulations dictionary `datadict` are set to True. If categorise=False, the simulations are plotted by author (L=Leeds, Y=Yadav, UCt=Uli Christensen FF, UC=Uli Christensen, chemical (0F at CMB), A=Aubert et al GJI 2009, APath=Aubert et al JFM 2017, S= Schwaiger et al GJI 2019). To select which simulation database to plot, change the keys `plot` in the simulations dictionary `datadict`.***
* write_check         : (logical) If True, write out the simulation data (fdip, fohm, p, Le, etc...) for each dataset
outf_check            : (string) file name of the simulation data check. The file name extension is set automatically to the simulation database. ***To compare unfiltered and filtered datasets it is easiest to set all filters above (myfdip, myEkOPm, myEr) to 0 and change the name here; then set the required filter, change name here and rerun.***
