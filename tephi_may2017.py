#=======================================================
# This script is an updated version of the original tephigram work.
# Created by: nmoissee@eoas.ubc.ca April 2017
#=======================================================
#INPUT
Tmin=-100.
Tmax=100.
THmin =-68.
THmax= 42.
Pbot = 105
Ptop = 1    #kPa - upper atmosphere limit surface
Plim = 1
degree =10      #degree of polinomial to model the curves
datafile = '%s-%s.npy' %(Tmin,Tmax)
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
import os.path


#thermodynamic constants
Rd = 287.058    #[J K^-1 kg^-1] gas constant for dry air
Rv = 461.5      #[J K^-1 kg^-1] gas constant for water vapour
Cp = 1006.      #[J K^-1 kg^-1] specific heat of dry air at constant pressure
T0 = 273.16     #standard temperature
P0 = 100.      #kPa
e0 = 0.611657   #kPa: adjusted Clausius-Clayperon constant (Koutsoyiannis 2011)

#derived constants
Eps = Rd/Rv     #dimensionless
c1 = Rd/Cp      #dimensionless

#create temp and pressure axes
Pbrange = np.arange(Pbot,10,-0.0001)
Pmrange = np.arange(10,2,-0.00001)
Ptrange = np.arange(2,Ptop,-0.000001)
PrangeList = np.concatenate((Pbrange,Pmrange,Ptrange))

# PrangeList = np.concatenate((Pbrange,Pmrange))
Prange = PrangeList[:-1]

#create temp and pressure axes
nTw = np.arange(Tmin,Tmax,0.5)
nThetaW = np.arange(THmin,THmax,0.5)

if os.path.isfile(datafile):
    arrayTHw= np.load(datafile)
else:
    #create storage arrays
    arrayTHw = np.empty((len(nTw),len(Prange)))     #theta moist
    arrayTHnorm = np.empty_like(arrayTHw)     #normalized array

    def f_es(T):
        #saturation vapour pressure
        return  e0*np.exp(24.921*(1.-(T0/T)))*((T0/T)**5.06)
    def f_rs(P,es):
        #saturated mixing ratio of water at temperature
        return  Eps*es / (P - es)
    def dTdP(P,T):
        return (c1*T + c3*rs)/(P*(1.+(c2*rs/T**2.)))

    for nT, THw in enumerate(nTw):
        T = THw + T0
        Tz = np.copy(T)         
        print('Current adiabat: %s' %THw)
        for nP,P in enumerate(Prange):
            #update the 'constants'
            Lv = 3.139e6 - 2336 * Tz
            c2 = (Lv**2)/(Rv*Cp)    #[K^2]
            c3 = Lv/Cp      #[K]

            #get moist adiabat
            es = f_es(Tz)
            rs = f_rs(P,es)
            grad = dTdP(P,Tz)
            Pstep = P - PrangeList[nP+1]
            Tz = Tz - grad*Pstep
            arrayTHw[nT,nP] = Tz

    np.save(datafile, arrayTHw)

#monotonically select points every 0.1 kPa for fitting
PrangeIdx = [np.argmin(abs(PrangeList - i)) for i in np.arange(Pbot,Plim,-0.1)]
PrangeFit = Prange[PrangeIdx]
P0idx = np.argmin(abs(Prange - P0))
Tidx = [np.argmin(abs(nTw - i)) for i in nThetaW]


# arrayTHnorm = np.copy(arrayTHw[Tidx,:])
arrayTHnorm = np.empty((len(nThetaW),len(PrangeFit)))*np.nan
subarrayTHw = arrayTHw[Tidx,:]
P0axis = subarrayTHw[:,P0idx]
C0axis = P0axis - T0

arrayTHnorm[:,:] = subarrayTHw[:,PrangeIdx]     #redundant, fix above


#normailzing by one of the adiabats removes the non-linearity from the data
THref = arrayTHnorm[0,:]
THref_fit = np.poly1d(np.polyfit(PrangeFit,THref,20))
MAE_THref = np.mean(abs(THref-THref_fit(PrangeFit)))
print('MAE for polynomial fit of Tmax reference curve: %.2E' %MAE_THref)
np.savetxt('THrefcoeffs.txt',THref_fit.coeffs)
np.savetxt('THrefcoeffs_latex.txt', THref_fit.coeffs.T, delimiter=' & ', fmt='%2.2e')

# for n in range(30): 
#     #normailzing by one of the adiabats removes the non-linearity from the data
#     THref_fit = np.poly1d(np.polyfit(PrangeFit,THref,n))
#     MAE_THref = np.mean(abs(THref-THref_fit(PrangeFit)))
#     print('MAE for polynomial fit of Tmax reference curve with n=%s: %.2E' %(n,MAE_THref))
#     # np.savetxt('THrefcoeffs.txt',THref_fit.coeffs)

# Now model,store coeffs and plot (for specified degree polynomial)
print('Fitting polynomials to normalized curves')
numterms = degree+1
store_args = np.zeros((numterms,len(nThetaW)))
tags = ['k$_{%s}$' %i for i in range(numterms)]
for i in range(len(P0axis)):
    main_pfit = np.poly1d(np.polyfit(THref_fit(PrangeFit),arrayTHnorm[i,:],degree))
    store_args[:,i] = main_pfit.coeffs
    plt.plot(THref_fit(PrangeFit),main_pfit(THref_fit(PrangeFit)),'r')
    plt.plot(THref_fit(PrangeFit),arrayTHnorm[i,:],'b')
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
np.savetxt('kcoeffs.txt', store_coeffs)
np.savetxt("kcoeffs_latex.txt", store_coeffs, delimiter=' & ', fmt='%2.2e', newline=' \\\\\n')


#TESTING THE METHOD======================================
print('Evaluating polynomial fit method....')
arrayTHfit = np.zeros((len(C0axis),len(PrangeFit)))

for nT, T in enumerate(C0axis):
    k = []
    #calculate parameters 
    for iDeg in range(numterms):
        k.append(fitFCNs[iDeg](T))
    #fit the moist adiabats
    for nP,P in enumerate(PrangeFit):
        THrefm = THref_fit(P)        
        THfit = 0.
        #sum up the polynomial terms
        for iDeg in range(numterms):
            THfit = THfit + k[iDeg]*THrefm**(degree-iDeg)
        arrayTHfit[nT,nP] = THfit

arrayDiff = arrayTHnorm-arrayTHfit
MAE = np.mean(abs(arrayDiff.ravel()))
print('FULL DOMAIN MAE: %s' %MAE) 

# #convert back to true adiabats
# arrayTHwm = np.empty_like(arrayTHfit)

# for nT, THw in enumerate(C0axis):
#     Ttop = Top_fit(THw)
#     Tbot = P0axis[nT]
#     # arrayTHwm[nT,:] = arrayTHfit[nT,:]*Ttop + Ttop
#     arrayTHwm[nT,:] = arrayTHfit[nT,:]*(Ttop - Tbot) + Tbot

# arrayDiffTemp = subarrayTHw[:,PrangeIdx]-arrayTHwm
# MAE = np.mean(abs(arrayDiffTemp.ravel()))
# print('FULL DOMAIN MAE: %s' %MAE) 

#=====================PLOTTING===========================

PaxisIdx = [np.argmin(abs(PrangeFit - i)) for i in np.arange(P0,Plim,-10)]
#plot emegram diagram
TrueTHwPlt = subarrayTHw[0::10,PrangeIdx].T - T0
TrueTHwLgnd = subarrayTHw[-1,PrangeIdx]- T0
plt.figure(figsize=(9,6))
# plt.title('EMAGRAM')
plt.plot(TrueTHwPlt,PrangeFit, 'k')
plt.plot(TrueTHwLgnd,PrangeFit, 'k',label='"true" $\\theta_w$')
plt.plot(arrayTHfit[0::10,:].T-T0,PrangeFit, 'r--')
plt.plot(arrayTHfit[-1,:]-T0,PrangeFit, 'r--',label='modelled $\\theta_w$')
ax = plt.gca()
ax.invert_yaxis()
plt.ylim([105,10])
plt.xlim([-100,40])
plt.yscale('log')
ax.set_yticks(PrangeFit[PaxisIdx])
ax.set_yticklabels(PrangeFit[PaxisIdx].round())
plt.grid()
plt.xlabel("temperature [$^\circ$C]")
plt.ylabel("pressure [kPa]")
plt.legend(loc='upper right',fontsize=12)
plt.savefig('./figs/emagram.pdf')
plt.show()
plt.close()

#plot fit of single adiabat THref
# plt.title('$\\theta_{ref} = \\theta_{-70}$ POLYNOMIAL FIT')
plt.plot(THref,PrangeFit,'g')
plt.plot(THref_fit(PrangeFit),PrangeFit,'r')
ax = plt.gca()
plt.gca().invert_yaxis()
plt.xlabel('temperature [K]')
plt.ylabel('pressure [kPa]')
plt.savefig('./figs/THref_May.pdf')
plt.show()
plt.close()

# #plot transformed adiabats
# plt.title('TRANSFORMED MOIST ADIABATS $\\theta_{trans}$')
# plt.plot(THref,arrayTHnorm[0:int(Tmax),:].T,'b')
# plt.plot(THref,arrayTHnorm[int(Tmax):,:].T,'r')
# plt.xlabel('$\\theta_{ref}$ [K]')
# plt.ylabel('$\\theta_{trans}$')
# # plt.savefig('./figs/THtrans.pdf')
# plt.show()
# plt.close()

#subplot of fits for individual parameters
fig = plt.figure(figsize=(12, 10)) 
# plt.suptitle('FIT PARAMETERS')
import matplotlib.ticker as mtick
for iDeg in range(degree+1):
    plt.subplot(3,4,iDeg+1)
    plt.title(tags[iDeg])
    plt.xlabel('$\\theta_w$ [$^\circ$C]',fontsize=8)
    plt.plot(C0axis,store_args[iDeg,:],'g')
    plt.plot(C0axis,fitFCNs[iDeg](C0axis),'r')
    plt.gca().tick_params(labelsize=7)
# plt.tight_layout()
plt.subplots_adjust(top = .92, hspace=0.4, wspace=0.3, left=0.06, right=0.97, bottom=0.05)
plt.savefig('./figs/fit_params_May.pdf')
plt.show()
plt.close()




# #plot true and fitted normalized saturated adiabats
# plt.figure(figsize=(8,6))
# plt.title('TRUE AND MODELLED $\\theta_{norm}$')
# plt.plot(Prange,arrayTHnorm[1::10,:].T,color='g', label='directly computed $\\theta_{norm}$')
# plt.plot(PrangeFit,arrayTHfit[1::10,:].T,'r',label='modelled $\\theta_{norm}$')
# plt.gca().invert_xaxis()
# # plt.xlim([101,1])
# # plt.ylim([0,1.1])
# plt.grid()
# plt.ylabel("normalized moist adiabats")
# plt.xlabel("pressure [kPa]")
# # plt.savefig('./figs/THfit.pdf')
# plt.show()
# plt.close()

#plot error distribution contours
plt.figure(figsize=(8,6))
# plt.title('ERROR CONTOURS')
plt.imshow(arrayDiff.T,aspect='auto',origin='lower',cmap='RdBu_r',vmin=-0.1,vmax=0.1)
plt.xlabel("$\\theta_w$ [$^\circ$C]")
plt.ylabel("pressure [kPa]")
ax = plt.gca()
ax.set_xticks(np.arange(0,len(nThetaW),20))
ax.set_xticklabels(C0axis[::20].astype(int))
ax.set_yticks(PaxisIdx)
ax.set_yticklabels(PrangeFit[PaxisIdx].round())
cbar = plt.colorbar(format='%.2f')
cbar.set_label('temperature difference [K]')
plt.savefig('./figs/ErrorTHw_May.pdf')
plt.show()
plt.close()

# #plot error distribution contours in degrees
# plt.figure(figsize=(8,6))
# plt.title('ERROR (C)')
# plt.imshow(arrayDiffTemp.T,aspect='auto',origin='lower',cmap='RdBu_r',vmin=-0.05, vmax=0.05)
# plt.xlabel("temperature [C]")
# plt.ylabel("pressure [kPa]")
# ax = plt.gca()
# ax.set_xticks(np.arange(0,len(nThetaW),10))
# ax.set_yticks(np.arange(13,len(PrangeFit),200))
# ax.set_xticklabels(nThetaW[::10])
# ax.set_yticklabels(np.arange(100,1,-20))
# cbar = plt.colorbar()
# cbar.set_label('temperature difference [C]')
# # plt.savefig('./figs/ErrorTHw.pdf')
# plt.show()
# plt.close()
