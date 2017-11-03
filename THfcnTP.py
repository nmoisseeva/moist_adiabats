#=======================================================
# This script is an updated version of the original tephigram work.
# Created by: nmoissee@eoas.ubc.ca May 2017
#=======================================================
#INPUT
Tmin =-100
Tmax =40
THmin =-100.    
THmax = 100.
Pbot = 105
Ptop = 1    #kPa - upper atmosphere limit surface
Plim = 1
degree =10      #degree of polinomial to model the curves
datafile = '%s-%s_T.npy' %(THmin,THmax)
#=======================================================

# #supress warnings
# import warnings
# def fxn():
#     warnings.warn("deprecated", DeprecationWarning)
# with warnings.catch_warnings():
#     warnings.simplefilter("ignore")
#     fxn()

#import necessary packages 
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, minimize
from scipy import spatial
from scipy.interpolate import griddata

#constants
T0 = 273.15     #standard temperature
P0 = 100.      #kPa

#create temp and pressure axes (from Part 1) - DO NOT CHANGE UNLESS PART 1 IS RECALCULATED
PrangeList = np.concatenate((np.arange(Pbot,10,-0.0001),np.arange(10,2,-0.00001),np.arange(2,Ptop,-0.000001)))
Prange = PrangeList[:-1]

#set up temperature ranges
nTw = np.arange(Tmin,Tmax,0.5)
nTwk = [i + T0 for i in nTw]

#load pre-calculated THw
arrayTHw = np.load(datafile)
nThetaW = np.arange(THmin,THmax,0.5)

#create a custom pressure range to aid fitting - DO NOT CHANGE UNLESS PART 1 IS RECALCULATED
segments = [10,2,Plim]  
segIdx = [np.argmin(abs(PrangeList - i)) for i in segments]
PrangeIdx = np.concatenate((np.arange(0,segIdx[0],1000,dtype=int),np.arange(segIdx[0],segIdx[1],10000,dtype=int)))
PrangeFit = Prange[PrangeIdx]




#extract contours of constant temperature, plot
cs = plt.contour(nThetaW,PrangeFit,arrayTHw[:,PrangeIdx].T,nTwk)
plt.gca().invert_yaxis()
plt.show()
plt.close()

arrayTw = []
arrayTax = []
arrayP = []
for nT, T in enumerate(nTwk):
    cs_paths = cs.collections[nT].get_paths()
    if len(cs_paths) > 0:
        path_vals =cs_paths[0]
        vtx = path_vals.vertices
        tnorm = vtx[:,0]

        # tnorm = (vtx[:,0] - nTw[nT])
        arrayTax.append(vtx[:,0])
        arrayTw.append(tnorm)
        arrayP.append(vtx[:,1])

C0axis = []
for nT, T in enumerate(nTwk):
    P0idx = np.argmin(abs(np.array(arrayP[nT][:]) - P0))
    C0axis.append(arrayTw[nT][P0idx])
# P0idx = np.argmin(abs(PrangeFit - P0))
# Tidx = [np.argmin(abs(nTw - i)) for i in nTw]
# P0axis = arrayTHw[Tidx,P0idx]
# C0axis = P0axis - T0


#normailzing by one of the curves removes the non-linearity from the data
refIdx = np.argmax([len(item) for item in arrayTw])     # make sure to pick the longest curve
Tref = arrayTw[refIdx]
Pref = arrayP[refIdx]
Tref_fit = np.poly1d(np.polyfit(Pref,Tref,20))
MAE_Tref = np.mean(abs(Tref-Tref_fit(Pref)))
print('MAE for polynomial fit of Trefcurve: %.2E' %MAE_Tref)
np.savetxt('Trefcoeffs.txt',Tref_fit.coeffs)
np.savetxt("Trefcoeffs_latex.txt", Tref_fit.coeffs, delimiter=' & ', fmt='%2.2e')



# Now model,store coeffs and plot (for specified degree polynomial)
print('Fitting polynomials to normalized curves')
numterms = degree+1
store_args = np.zeros((numterms,len(nTwk)))
tags = ['$\\kappa_{%s}$' %i for i in range(numterms)]
for nT in range(len(nTwk)):
    main_pfit = np.poly1d(np.polyfit(Tref_fit(arrayP[nT]),arrayTw[nT],degree))
    store_args[:,nT] = main_pfit.coeffs
    plt.plot(main_pfit(Tref_fit(arrayP[nT])),arrayP[nT],'r')
    plt.plot(arrayTw[nT],arrayP[nT],'b')
plt.gca().invert_yaxis()
plt.show()
plt.close()

#now do fits for individual parameters
print('Fitting polynomials to curve parameters')
fitFCNs = []
store_coeffs = []
for iDeg in range(numterms):
    pfit = np.poly1d(np.polyfit(C0axis,store_args[iDeg,:],20))
    MAE = np.mean(abs(store_args[iDeg,:] - pfit(C0axis)))
    print('%s MAE = %0.2E' %(tags[iDeg],MAE))
    fitFCNs.append(pfit)
    store_coeffs.append(pfit.coeffs)
np.savetxt('kcoeffsT.txt', store_coeffs)
np.savetxt("kcoeffsT_latex.txt", store_coeffs, delimiter=' & ', fmt='%2.2e', newline=' \\\\\n')


#TESTING THE METHOD======================================
arrayTHfit = []
for nT, T in enumerate(C0axis):
    k = []
    #calculate parameters 
    for iDeg in range(numterms):
        k.append(fitFCNs[iDeg](T))
    #fit the moist adiabats
    storeFit = []
    for nP,P in enumerate(arrayP[nT]):
        #sum up the polynomial terms
        THfit = 0.
        for iDeg in range(numterms):
            THfit = THfit + k[iDeg]*Tref_fit(P)**(degree-iDeg)
            # THfit = THfit + k[iDeg]*P**(degree-iDeg)
        storeFit.append(THfit)
    arrayTHfit.append(storeFit)

arrayDiff = []
MAE = []
for nT, T in enumerate(nTwk):
    diff = arrayTw[nT] - arrayTHfit[nT]
    arrayDiff.append(diff)
    plt.plot(arrayTw[nT],arrayP[nT],'b')
    plt.plot(arrayTHfit[nT],arrayP[nT],'r' )
    MAE.extend(abs(diff))
plt.gca().invert_yaxis()
plt.show()
plt.close()
MAE = np.mean(MAE)
print('FULL DOMAIN MAE: %s' %MAE) 

#regrid error to rectilinear for plotting
Taxis, Paxis = np.arange(Tmin,Tmax,0.5), np.arange(Pbot,Plim,-0.1)
Pidx = [np.argmin(abs(Paxis - i)) for i in np.arange(100,1,-10)]
grid_x, grid_y = np.meshgrid(Taxis,Paxis)
flatten = lambda l: [item for sublist in l for item in sublist]
flatTH, flatP, flatError = np.array(flatten(arrayTax)),np.array(flatten(arrayP)), np.array(flatten(arrayDiff))
gridError = griddata(zip(flatTH,flatP),flatError,(grid_x, grid_y), method='cubic')

#=====================PLOTTING===========================
#plot fit of single curve Tref
# plt.title('$T_{ref}$ POLYNOMIAL FIT')
plt.plot(Tref,Pref,'g')
plt.plot(Tref_fit(Pref),Pref,'r')
plt.gca().invert_yaxis()
plt.xlabel('$\\theta_w$ [$^\\circ$C]')
plt.ylabel('pressure [kPa]')
plt.savefig('./figs/Tref_May.pdf')
plt.show()
plt.close()


#subplot of fits for individual parameters
fig = plt.figure(figsize=(12, 10)) 
# plt.suptitle('FIT PARAMETERS')
for iDeg in range(degree+1):
    plt.subplot(3,4,iDeg+1)
    plt.title(tags[iDeg])
    plt.xlabel('T [$^\\circ$C]',fontsize=8)
    plt.plot(C0axis,store_args[iDeg,:],'g')
    plt.plot(C0axis,fitFCNs[iDeg](C0axis),'r')
    plt.gca().tick_params(labelsize=8.5)
plt.subplots_adjust(top = .92, hspace=0.4, wspace=0.4, left=0.05, right=0.97, bottom=0.05)
plt.savefig('./figs/fit_params_T_May.pdf')
plt.show()
plt.close()


plt.figure(figsize=(8,6))
# plt.title('ERROR CONTOURS')
plt.imshow(gridError,aspect='auto',origin='lower',vmin=-0.01, vmax=0.01,cmap='RdBu_r')
# plt.contourf(grid_x,grid_y,gridError,vmin=-0.05, vmax=0.05,cmap='RdBu_r')
ax = plt.gca()
ax.set_xticks(np.arange(0,len(Taxis),len(Taxis)/7))
ax.set_yticks(Pidx)
ax.set_xticklabels(Taxis[::len(Taxis)/7].astype(int))
ax.set_yticklabels(np.arange(100,1,-10, dtype=int))
cbar = plt.colorbar(format='%.3f')
cbar.set_label('temperature difference [K]')
plt.xlabel('$\\theta_w$ [$^\circ$C]')
plt.ylabel('pressure [kPa]')
plt.savefig('./figs/ErrorT_May.pdf')
plt.show()
plt.close()
