#=======================================================
# This script is an updated version of the original tephigram work.
# Created by: nmoissee@eoas.ubc.ca Nov 2016
#=======================================================
#INPUT
Tmin=-40.
Tmax=40.
Pstep=0.001


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

#thermodynamic constants
Rd = 287.058    #[J K^-1 kg^-1] gas constant for dry air
Rv = 461.5  #[J K^-1 kg^-1] gas constant for water vapour
Cp = 1006.  #[J K^-1 kg^-1] specific heat of dry air at constant pressure
# Lv = 2.501e6  #latent heat of vapourization at standard temperature
T0 = 273.16     #standard temperature
e0 = 0.611657   #kPa: adjusted Clausius-Clayperon constant (Koutsoyiannis 2011)

#derived constants
Eps = Rd/Rv     #dimensionless
c1 = Rd/Cp      #dimensionless
# c2 = (Lv**2)/(Rv*Cp)  #[K^2]
# c3 = Lv/Cp        #[K]


Prange = np.arange(100,1, -Pstep)
nThetaW = np.arange(Tmin,Tmax)

arrayTHw = np.empty((len(Prange),len(nThetaW)))
arrayTHd = np.empty_like(arrayTHw)
arrayTHe = np.empty_like(arrayTHw)

def f_es(T):
    #saturation vapour pressure
    #REPLACING STANDARD EQUATION WITH Koutsoyiannis 2011
    return  e0*np.exp(24.921*(1.-(T0/T)))*((T0/T)**5.06)
def f_rs(P,es):
    #saturated mixing ratio of water at temperature
    return  Eps*es / (P - es)
def dTdP(P,T):
    return (c1*T + c3*rs)/(P*(1.+(c2*rs/T**2.)))
def f_thE(T,rs0):
    # return T * np.exp(c4*rs0/T)
    return T + Lv/Cp*rs0
    #second formula produces smaller error for positive ThetaW values


for nT, THw in enumerate(nThetaW):
    T = THw + T0
    Tz = np.copy(T) #save surface temperature which will be iterated up the pressure levels
    #------variable latent heat of vapourization constants-------
    # Lv = (2500.8 - 2.36*Temp + 0.0016*(Temp**2) - 0.00006*(Temp**3))*1000
    Lv = 3.139e6 - 2336 * T #(Koutsoyiannis 2011) - theoretically derived
    c2 = (Lv**2)/(Rv*Cp)    #[K^2]
    c3 = Lv/Cp      #[K]
    # print c3
    #------------------------------------------------------------
    print('Current adiabat: %s' %THw)
    
    for nP,P in enumerate(Prange):
        #get dry adiabat
        arrayTHd[nP,nT] = T*((P/100.)**c1) #Temp + T0 to avoid overwrite
        #get moist adiabat
        es = f_es(Tz)
        rs = f_rs(P,es)
        grad = dTdP(P,Tz)
        Tz = Tz - grad*Pstep
        arrayTHw[nP,nT] = Tz
        #qet equivalent adiabat
        rs0 = f_rs(100.,es)
        arrayTHe[nP,nT] = f_thE(Tz,rs0)

#plot stuve's diagram
plt.figure(figsize=(9,6))
plt.title('SANITY CHECK: "STUVE" PLOT')
plt.plot(arrayTHw[:,-1]-T0,Prange, 'b',label='moist adiabat $\\theta_w$')
plt.plot(arrayTHd[:,-1]-T0,Prange, 'r--',label='dry adiabat $\\theta_d$')
plt.plot(arrayTHe[:,-1]-T0,Prange, 'g:',label='equivalent potential temperature $\\theta_e}$')
plt.plot(arrayTHw[:,0::10]-T0,Prange, 'b')
plt.plot(arrayTHd[:,0::10]-T0,Prange, 'r--')
plt.plot(arrayTHe[:,0::10]-T0,Prange, 'g:')
plt.ylim([40,100])
plt.gca().invert_yaxis()
plt.xlim([-40,40])
plt.grid()
plt.xlabel("temperature [C]")
plt.ylabel("pressure [kPa]")
plt.legend(loc='upper right',fontsize=12)
plt.savefig('./figs/stuve.pdf')
plt.show()

#plot normalized adiabats
plt.title('NORMALIZED SATURATED ADIABATS $\\theta_{norm}$')
arrayTHnorm = (arrayTHw - arrayTHd)/(arrayTHe - arrayTHd)
plt.plot(arrayTHnorm[:,0::10],Prange, 'b')
plt.gca().invert_yaxis()
plt.xlim([0,1])
plt.xlabel('normalized temperature')
plt.ylabel('pressure [kPa]')
plt.savefig('./figs/THnorm.pdf')
plt.show()
plt.close()

#normailzing by one of the adiabats removes the non-linearity from the data
plt.title('$\\theta_{ref} = \\theta_{40C}$ POLYNOMIAL FIT')
THref = arrayTHnorm[:,-1]
THref_fit = np.poly1d(np.polyfit(Prange,THref,28))
plt.plot(THref,Prange,'g')
plt.plot(THref_fit(Prange),Prange,'r')
plt.gca().invert_yaxis()
plt.xlim([0,1])
plt.xlabel('normalized temperature')
plt.ylabel('pressure [kPa]')
plt.savefig('./figs/THref.pdf')
MAE_THref = np.mean(abs(THref-THref_fit(Prange)))
print('MAE for polynomial fit of +40C reference curve: %.2E' %MAE_THref)
plt.show()
plt.close()

plt.title('TRANSFORMED MOIST ADIABATS $\\theta_{trans}$')
plt.plot(THref,arrayTHnorm[:,0:Tmax],'b')
plt.plot(THref,arrayTHnorm[:,Tmax:],'r')
plt.xlim([0,1])
plt.xlabel('reference saturated adiabat')
plt.ylim([0,1])
plt.ylabel('transformated saturated adiabats')
plt.savefig('./figs/THtrans.pdf')
plt.show()
plt.close()

degree = 7 #degree of polinomial to model the curves
# Now model the fit parameters (for 8th degree polynomial)
store_args = np.zeros((degree+1,len(nThetaW)-1))
tags = ['k%s' %(i+1) for i in range(degree+1)]
for i in range(len(nThetaW)-1):
    main_pfit = np.poly1d(np.polyfit(THref,arrayTHnorm[:,i],degree))
    store_args[:,i] = main_pfit.coeffs
    plt.plot(THref,main_pfit(THref),'r')
    plt.plot(THref,arrayTHnorm[:,i],'b')
    plt.xlim([0,1])
    plt.ylim([0,1])
plt.show()

fig = plt.figure(figsize=(9, 9)) 
plt.suptitle('FIT PARAMETERS')
xvals = nThetaW[:-1]
#fits for individual parameters
fitFCNs = []
for iDeg in range(degree):
    pfit = np.poly1d(np.polyfit(xvals,store_args[iDeg,:],25))
    plt.subplot(3,3,iDeg+1)
    plt.title(tags[iDeg])
    plt.plot(xvals,store_args[iDeg,:],'g')
    plt.plot(xvals,pfit(xvals),'r')
    MAE = np.mean(abs(store_args[iDeg,:] - pfit(xvals)))
    print('%s MAE = %0.2E' %(tags[iDeg],MAE))
    fitFCNs.append(pfit)
plt.tight_layout()

plt.savefig('./figs/fit_params.pdf')
plt.show()

#+!!!!!!! STOPPED THE EDITS HERE

#TESTING THE METHOD======================================
fit_adiabats = np.empty_like(adiabats)
# for nT, Temp in enumerate(testvals):
for nT, Temp in enumerate(xvals):
    k1,k2,k3,k4,k5,k6,k7,k8,k9 = pfit1(Temp),pfit2(Temp),pfit3(Temp),pfit4(Temp),pfit5(Temp),pfit6(Temp),pfit7(Temp),pfit8(Temp),pfit9(Temp)
    for nP,Pres in enumerate(Prange):
        normP = pref_fit(Pres)
        normTH = k1*normP**8 + k2*normP**7 + k3*normP**6 + k4*normP**5 + k5*normP**4 + k6*normP**3 + k7*normP**2 + k8*normP + k9
        fit_adiabats[nP,nT] = normTH
plt.figure(figsize=(8,6))
plt.title('FINAL FIT RESULTS')
plt.plot(norm_adiabats[:,1::5],color='0.5')
plt.plot(fit_adiabats[:,1::5],'r')
plt.ylim([0,1.1])
plt.grid()
plt.ylabel("normalized moist adiabats")
plt.xlabel("pressure [kPa]")
plt.savefig('final_fit.pdf')
plt.show()

plt.figure(figsize=(8,6))
plt.title('ERROR DISTRIBUTION PLOT')
plt.plot(abs(norm_adiabats[:,:]- fit_adiabats[:,:]))

plt.grid()
plt.ylabel("normalized moist adiabats")
plt.xlabel("pressure [kPa]")
# plt.savefig('final_fit.pdf')
plt.show()


error = norm_adiabats - fit_adiabats
mean_error = np.mean(error,1)
plt.title('MEAN ERROR PROFILE')
plt.plot(mean_error, Prange)
plt.gca().invert_yaxis()
plt.grid()
plt.xlabel("normalized mean error")
plt.ylabel("pressure [kPa]")
plt.savefig('error_profile.pdf')
plt.show()
