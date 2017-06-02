#=======================================================
# This script is an updated version of the original tephigram work.
# Created by: nmoissee@eoas.ubc.ca May 2017
#=======================================================
#INPUT
Tmin =-80
Tmax =40
THmin =-80    
THmax = 90

Ptop = 1    #kPa - upper atmosphere limit surface (from Part 1)
degree =10      #degree of polinomial to model the curves
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
T0 = 273.16     #standard temperature
P0 = 101.3      #kPa

#create temp and pressure axes (from Part 1) - DO NOT CHANGE UNLESS PART 1 IS RECALCULATED
PrangeList = np.concatenate((np.arange(P0,10,-0.0001),np.arange(10,2,-0.00001),np.arange(2,Ptop,-0.000001)))
Prange = PrangeList[:-1]

#create a custom pressure range to aid fitting
segments = [10,2,Ptop]
segIdx = [np.argmin(abs(PrangeList - i)) for i in segments]
# PrangeIdx = np.concatenate((np.arange(P0,segIdx[0],1,dtype=int),np.arange(segIdx[0],segIdx[1],1000,dtype=int),np.arange(segIdx[1],segIdx[2],100,dtype=int)))

PrangeIdx = np.concatenate((np.arange(0,segIdx[0],1000,dtype=int),np.arange(segIdx[0],segIdx[1],10000,dtype=int)))
PrangeFit = Prange[PrangeIdx]

# PrangeIdx = np.arange(0,1413000,10)
# PrangeFit = Prange[PrangeIdx]

#load pre-calculated THw
arrayTHw = np.load('%s.0-%s.0.npy' %(THmin,THmax))
nThetaW = np.arange(THmin,THmax)

#set up temperature ranges
nTw = np.arange(Tmin,Tmax,0.5)
nTwk = [i + T0 for i in nTw]

cs = plt.contour(nThetaW,PrangeFit,arrayTHw[:,PrangeIdx].T,nTwk)
# ax = plt.gca()
# ax.set_yticks(np.arange(0,len(PrangeFit),100))
# ax.set_yticklabels(PrangeFit[::100])


plt.gca().invert_yaxis()
plt.show()

arrayTw = []
arrayTax = []
arrayP = []
for nT, T in enumerate(nTwk):
    cs_paths = cs.collections[nT].get_paths()
    if len(cs_paths) > 0:
        path_vals =cs_paths[0]
        vtx = path_vals.vertices
        tnorm = (vtx[:,0] - nTw[nT])
        arrayTax.append(vtx[:,0])
        arrayTw.append(tnorm)
        arrayP.append(vtx[:,1])


#normailzing by one of the adiabats removes the non-linearity from the data
refIdx = np.argmax([len(item) for item in arrayTw])
Tref = arrayTw[refIdx]
Pref = arrayP[refIdx]
Tref_fit = np.poly1d(np.polyfit(Pref,Tref,20))
MAE_Tref = np.mean(abs(Tref-Tref_fit(Pref)))
print('MAE for polynomial fit of Trefcurve: %.2E' %MAE_Tref)

# np.savetxt('THrefcoeffs.txt',THref_fit.coeffs)

#plot fit of single adiabat Tref
# plt.title('$\\T_{ref} = \\theta_{%s}$ POLYNOMIAL FIT' %nTwk[refIdx])
plt.plot(Tref,Pref,'g')
plt.plot(Tref_fit(Pref),Pref,'r')
plt.gca().invert_yaxis()
# plt.xlim([0,100])
# plt.ylim([6,5])
plt.xlabel('normalized temperature')
plt.ylabel('pressure [kPa]')
# plt.savefig('./figs/Tref.pdf')
plt.show()
plt.close()


# Now model,store coeffs and plot (for specified degree polynomial)
print('Fitting polynomials to normalized curves')
numterms = degree+1
store_args = np.zeros((numterms,len(nTwk)))
tags = ['k%s' %i for i in range(numterms)]
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
    pfit = np.poly1d(np.polyfit(nTwk,store_args[iDeg,:],25))
    MAE = np.mean(abs(store_args[iDeg,:] - pfit(nTwk[:])))
    print('%s MAE = %0.2E' %(tags[iDeg],MAE))
    fitFCNs.append(pfit)
    store_coeffs.append(pfit.coeffs)
# np.savetxt('kcoeffs.txt', store_coeffs)

#subplot of fits for individual parameters
fig = plt.figure(figsize=(10, 12)) 
plt.suptitle('FIT PARAMETERS')
for iDeg in range(degree):
    plt.subplot(4,4,iDeg+1)
    plt.title(tags[iDeg])
    plt.xlabel('temperature [K]',fontsize=8)
    plt.plot(nTwk,store_args[iDeg,:],'g')
    plt.plot(nTwk,fitFCNs[iDeg](nTwk),'r')
plt.subplots_adjust(top = .92, hspace=0.4, wspace=0.4, left=0.05, right=0.97, bottom=0.05)
plt.savefig('./figs/fit_params_T.pdf')
plt.show()
plt.close()


#TESTING THE METHOD======================================
arrayTHfit = []
for nT, T in enumerate(nTwk):
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
for nT, T in enumerate(nTwk):
    diff = arrayTw[nT] - arrayTHfit[nT]
    arrayDiff.append(diff)
    plt.plot(arrayTw[nT],arrayP[nT],'b')
    plt.plot(arrayTHfit[nT],arrayP[nT],'r' )
plt.gca().invert_yaxis()
plt.show()
plt.close()


# grid_x, grid_y = np.mgrid[-40:Tmax:0.5,P0:2:-0.01]
grid_x, grid_y = np.mgrid[-40:Tmax:0.05,P0:2:-0.05]


flatten = lambda l: [item for sublist in l for item in sublist]

# maxLen = len(sorted(arrayDiff,key=len, reverse=True)[0])
# error2D = np.zeros((len(arrayDiff),maxLen))*np.nan
# p = np.zeros((len(arrayDiff),maxLen))*np.nan

# for nT, T in enumerate(nTwk):
#     vals = arrayDiff[nT][::-1]
#     p[nT,:len(vals)] = arrayP[nT][::-1]
#     error2D[nT,:len(vals)] = vals

flatTH, flatP, flatError = np.array(flatten(arrayTax)),np.array(flatten(arrayP)), np.array(flatten(arrayDiff))

gridError = griddata(zip(flatTH,flatP),flatError,(grid_x, grid_y), method='cubic')

plt.figure(figsize=(8,6))
# plt.title('ABSOLUTE ERROR (C)')

plt.imshow(gridError.T,aspect='auto',origin='lower',vmin=-0.05, vmax=0.05,cmap='RdBu_r')
# plt.gca().invert_yaxis()
plt.colorbar()
plt.show()
plt.close()


plt.figure(figsize=(8,6))
plt.title('ABSOLUTE ERROR (C)')
# for nT,T in enumerate(nTwk):
#     plt.scatter(arrayTw[nT][::10]+nTw[nT],arrayP[nT][::10],s=1,c=arrayDiff[nT],vmin=-0.05,vmax=0.05, cmap='RdBu_r')
# plt.contourf(error2D,cmap='RdBu_r',vmin=-0.05, vmax=0.05)

plt.imshow(error2D.T,aspect='auto',origin='lower',cmap='RdBu_r',vmin=-0.05, vmax=0.05)
# plt.gca().invert_yaxis()
# plt.xlim([-40,40])
plt.colorbar()
plt.savefig('./figs/DiffContoursTHw.pdf')
plt.show()
plt.close()
