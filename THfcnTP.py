#=======================================================
# This script is an updated version of the original tephigram work.
# Created by: nmoissee@eoas.ubc.ca May 2017
#=======================================================
#INPUT
Tmin =-60.
Tmax =40.
Ptop = 1    #kPa - upper atmosphere limit surface
degree =10      #degree of polinomial to model the curves

THmin =-60    
THmax = 60

#=======================================================
#supress warnings
import warnings
def fxn():
    warnings.warn("deprecated", DeprecationWarning)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    fxn()

#import necessary packages 
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, minimize
from scipy import spatial

# #thermodynamic constants
# Rd = 287.058    #[J K^-1 kg^-1] gas constant for dry air
# Rv = 461.5  #[J K^-1 kg^-1] gas constant for water vapour
# Cp = 1006.  #[J K^-1 kg^-1] specific heat of dry air at constant pressure
T0 = 273.16     #standard temperature
P0 = 101.3      #kPa
# e0 = 0.611657   #kPa: adjusted Clausius-Clayperon constant (Koutsoyiannis 2011)

# #derived constants
# Eps = Rd/Rv     #dimensionless
# c1 = Rd/Cp      #dimensionless

#create temp and pressure axes
PrangeList = np.concatenate((np.arange(P0,10,-0.0001),np.arange(10,2,-0.00001),np.arange(2,Ptop,-0.000001)))
Prange = PrangeList[:-1]

#monotonically select points every 0.1 kPa for fitting
segments = [95,11,10]
segIdx = [np.argmin(abs(PrangeList - i)) for i in segments]
PrangeIdx = np.concatenate((np.arange(P0,segIdx[0],1,dtype=int),np.arange(segIdx[0],segIdx[1],1000,dtype=int),np.arange(segIdx[1],segIdx[2],100,dtype=int)))
PrangeFit = Prange[PrangeIdx]


# Psubrange = np.concatenate((np.arange(P0,90,-0.0001),np.arange(90,6,-0.0001),np.arange(6,5,-0.0001)))
# PrangeIdx = [np.argmin(abs(PrangeList - i)) for i in Psubrange]
# PrangeFit = Prange[PrangeIdx]

# PrangeIdx = np.arange(0,1413000,10)
# PrangeFit = Prange[PrangeIdx]

#load pre-calculated THw
arrayTHw = np.load('%s.0-%s.0.npy' %(THmin,THmax))
nThetaW = np.arange(THmin,THmax)

#set up temperature ranges
nTw = np.arange(Tmin,Tmax,0.5)
nTwk = [i + T0 for i in nTw]
# nTwlist = [round(i,2) for i in nTw]   

cs = plt.contour(nThetaW,PrangeFit,arrayTHw[:,PrangeIdx].T,nTwk)
# cs = plt.contour(nThetaW,Prange,arrayTHw[:,:].T,nTwk[:2])

plt.gca().invert_yaxis()
plt.show()

arrayTw = []
arrayP = []
for nT, T in enumerate(nTwk):
    cs_paths = cs.collections[nT].get_paths()
    if len(cs_paths) > 0:
        path_vals =cs_paths[0]
        vtx = path_vals.vertices
        tnorm = (vtx[:,0] - nTw[nT])
        # print tnorm[-1]
        arrayTw.append(tnorm)
        # arrayTw.append(vtx[:,0])
        arrayP.append(vtx[:,1])


#normailzing by one of the adiabats removes the non-linearity from the data
refIdx = np.argmax([len(item) for item in arrayTw])
Tref = arrayTw[refIdx]
Pref = arrayP[refIdx]
Tref_fit = np.poly1d(np.polyfit(Pref,Tref,20))
MAE_Tref = np.mean(abs(Tref-Tref_fit(Pref)))
print('MAE for polynomial fit of Trefcurve: %.2E' %MAE_Tref)

# np.savetxt('THrefcoeffs.txt',THref_fit.coeffs)

# arrayTnorm = []
# for nT, T in enumerate(nTwk[-3:]):
#     arrayTnorm.append(Tref_fit(arrayP[nT])/arrayTw[nT])
#     # plt.plot(arrayTnorm[nT],arrayP[nT],'g')
#     plt.plot(Tref_fit(arrayP[nT]),'r')
#     plt.plot(arrayTw[nT],'b')
# plt.gca().invert_yaxis()
# # plt.xlim([0,5])
# plt.show()

#plot fit of single adiabat THref
# plt.title('$\\theta_{ref} = \\theta_{40C}$ POLYNOMIAL FIT')
plt.plot(Tref,Pref,'g')
plt.plot(Tref_fit(Pref),Pref,'r')
plt.gca().invert_yaxis()
# plt.xlim([0,100])
# plt.ylim([6,5])
plt.xlabel('normalized temperature')
plt.ylabel('pressure [kPa]')
# plt.savefig('./figs/THref.pdf')
plt.show()
plt.close()


# Now model,store coeffs and plot (for specified degree polynomial)
print('Fitting polynomials to normalized curves')
numterms = degree+1
store_args = np.zeros((numterms,len(nTwk)))
tags = ['k%s' %i for i in range(numterms)]
for nT in range(len(nTwk)):
    # main_pfit = np.poly1d(np.polyfit(arrayP[nT],arrayTw[nT],degree))
    main_pfit = np.poly1d(np.polyfit(Tref_fit(arrayP[nT]),arrayTw[nT],degree))
    store_args[:,nT] = main_pfit.coeffs
    # plt.plot(main_pfit(arrayP[nT]),arrayP[nT],'r')
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
    pfit = np.poly1d(np.polyfit(nTwk,store_args[iDeg,:],20))
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
# plt.savefig('./figs/fit_params.pdf')
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
    for nP,P in enumerate(arrayP[nT][::10]):
        #sum up the polynomial terms
        THfit = 0.
        for iDeg in range(numterms):
            THfit = THfit + k[iDeg]*Tref_fit(P)**(degree-iDeg)
            # THfit = THfit + k[iDeg]*P**(degree-iDeg)
        storeFit.append(THfit)
    arrayTHfit.append(storeFit)

arrayDiff = []
for nT, T in enumerate(nTwk):
    diff = arrayTw[nT][::10] - arrayTHfit[nT]
    arrayDiff.append(diff)
    MAE = np.mean(abs(diff),0)
    print MAE


for nT in range(len(nTwk)):
    plt.plot(arrayTw[nT],arrayP[nT],'b')
    plt.plot(arrayTHfit[nT],arrayP[nT][::10],'r' )
plt.gca().invert_yaxis()
plt.show()
plt.close()


plt.figure(figsize=(8,6))
plt.title('ABSOLUTE ERROR (C)')
for nT,T in enumerate(nTwk):
    plt.scatter(arrayTw[nT][::10]+nTw[nT],arrayP[nT][::10],s=1,c=arrayDiff[nT],vmin=-0.1,vmax=0.1, cmap='RdBu_r')
plt.gca().invert_yaxis()
plt.xlim([-40,40])
# plt.colorbar()
# plt.savefig('./figs/DiffContoursTemp.pdf')
# plt.show()


plt.show()
plt.close()


# plt.figure(figsize=(8,6))
# plt.title('ABSOLUTE ERROR (C)')
# plt.imshow(arrayDiffTemp.T,aspect='auto',origin='lower',cmap='RdBu_r',vmin=-0.05, vmax=0.05)
# plt.xlabel("temperature [C]")
# plt.ylabel("pressure [kPa]")
# ax = plt.gca()
# ax.set_xticks(np.arange(0,len(nThetaW),10))
# ax.set_yticks([1300,51300,100300])
# ax.set_xticklabels(nThetaW[::10])
# ax.set_yticklabels([100,50,1])
# plt.colorbar()
# plt.savefig('./figs/DiffContoursTemp.pdf')
# plt.show()

# #=====================PLOTTING===========================


# #plot true and fitted normalized saturated adiabats
# plt.figure(figsize=(8,6))
# plt.title('TRUE AND MODELLED $\\theta_{norm}$')
# plt.plot(Prange,arrayTHnorm[1::10,:].T,color='g', label='directly computed $\\theta_{norm}$')
# plt.plot(Prange,arrayTHfit[1::10,:].T,'r',label='modelled $\\theta_{norm}$')
# plt.gca().invert_xaxis()
# plt.xlim([101,1])
# plt.ylim([0,1.1])
# plt.grid()
# plt.ylabel("normalized moist adiabats")
# plt.xlabel("pressure [kPa]")
# plt.savefig('./figs/THfit.pdf')
# plt.show()

# #plot error distribution contours
# plt.figure(figsize=(8,6))
# plt.title('ERROR CONTOURS OF $\\theta_{norm}$')
# plt.imshow(arrayDiff.T,aspect='auto',origin='lower',cmap='RdBu_r',vmin=-0.005,vmax=0.005)
# plt.xlabel("temperature [C]")
# plt.ylabel("pressure [kPa]")
# ax = plt.gca()
# ax.set_xticks(np.arange(0,len(nThetaW),10))
# ax.set_yticks([1300,51300,100300])
# ax.set_xticklabels(nThetaW[::10])
# ax.set_yticklabels([100,50,1])
# plt.colorbar()
# plt.savefig('./figs/DiffContours.pdf')
# plt.show()

# #plot error distribution contours in degrees
# plt.figure(figsize=(8,6))
# plt.title('ABSOLUTE ERROR (C)')
# plt.imshow(arrayDiffTemp.T,aspect='auto',origin='lower',cmap='RdBu_r',vmin=-0.05, vmax=0.05)
# plt.xlabel("temperature [C]")
# plt.ylabel("pressure [kPa]")
# ax = plt.gca()
# ax.set_xticks(np.arange(0,len(nThetaW),10))
# ax.set_yticks([1300,51300,100300])
# ax.set_xticklabels(nThetaW[::10])
# ax.set_yticklabels([100,50,1])
# plt.colorbar()
# plt.savefig('./figs/DiffContoursTemp.pdf')
# plt.show()

# #plot stuve's diagram with modelled adiabats
# plt.figure(figsize=(9,6))
# plt.title('TESTING')
# plt.plot(arrayTHw[-1,:]-T0,Prange, 'b',label='moist adiabat $\\theta_w$')
# plt.plot(arrayTHd[-1,:]-T0,Prange, 'r--',label='dry adiabat $\\theta_d$')
# plt.plot(arrayTHe[-1,:]-T0,Prange, 'g:',label='equivalent potential temperature $\\theta_e}$')
# plt.plot(arrayTHwm[-1,:]-T0,Prange, 'y',label='modelled moist adiabat $\\theta_{mod}$')
# plt.plot(arrayTHw[0::5,:].T-T0,Prange, 'b')
# plt.plot(arrayTHd[0::10,:].T-T0,Prange, 'r--')
# plt.plot(arrayTHe[0::10,:].T-T0,Prange, 'g:')
# plt.plot(arrayTHwm[0::5,:].T-T0,Prange, 'y')
# plt.ylim([40,101])
# plt.gca().invert_yaxis()
# plt.xlim([-40,40])
# plt.grid()
# plt.xlabel("temperature [C]")
# plt.ylabel("pressure [kPa]")
# plt.legend(loc='upper right',fontsize=12)
# # plt.savefig('./figs/stuve.pdf')
# plt.show()
# plt.close()


# # #plot THw(TP)
# # plt.figure(figsize=(8,6))
# # plt.title('$\\theta_{w}(T,P)$')
# # plt.imshow(WofTP_MAE.T,aspect='auto',origin='lower',cmap='RdBu_r')
# # plt.xlabel("temperature [C]")
# # plt.ylabel("pressure [kPa]")
# # ax = plt.gca()
# # ax.set_xticks(np.arange(0,len(nThetaW),5))
# # # ax.set_yticks([1300,51300,100300])
# # ax.set_xticklabels(nThetaW[::5])
# # # ax.set_yticklabels([100,50,1])
# # plt.colorbar()
# # # plt.savefig('./figs/DiffContoursTemp.pdf')
# # plt.show()