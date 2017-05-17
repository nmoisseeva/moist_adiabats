#=======================================================
# This script is an updated version of the original tephigram work.
# Created by: nmoissee@eoas.ubc.ca May 2017
#=======================================================
#INPUT
Tmin=-40.
Tmax=-30.
# Pstep=0.0001
Ptop = 1    #kPa - upper atmosphere limit surface
degree =10      #degree of polinomial to model the curves

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
T0 = 273.16     #standard temperature
P0 = 101.3      #kPa
e0 = 0.611657   #kPa: adjusted Clausius-Clayperon constant (Koutsoyiannis 2011)

#derived constants
Eps = Rd/Rv     #dimensionless
c1 = Rd/Cp      #dimensionless


#create temp and pressure axes
Pbrange = np.arange(P0,10,-0.0001)
Pmrange = np.arange(10,2,-0.00001)
Ptrange = np.arange(2,Ptop,-0.000001)

PrangeList = np.concatenate((Pbrange,Pmrange,Ptrange))
Prange = PrangeList[:-1]


# nThetaW = np.arange(Tmin,Tmax, 0.1)
nThetaW = np.arange(Tmin,Tmax, 1)
nThetaWk = [i + T0 for i in nThetaW]

nThetaWlist = [round(i,2) for i in nThetaW]


# #from Thet fit
# popt = [  2.31265091e+06,   1.35750266e+03,   9.09425443e+05]
# #create a new equation for THe of THw at surface
# def func(T, L0,L1,K):
#     Tk = T + T0
#     es0 = f_es(Tk)
#     return Tk * np.exp((((L0 - L1*(T-T0)) + K*f_rs(P0,es0))*f_rs(P0,es0))/(Cp*Tk))

# func(nThetaW, *popt)





#create storage arrays
arrayTw = np.empty((len(nThetaW),len(Prange)))     #theta moist
arrayTd = np.empty_like(arrayTw)                  #theta dry
arrayTe = np.empty_like(arrayTw)                  #theta equivalent (from ref pressure level)

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
    # return T*np.exp((Lv*rs0)/(Cp*T))  #holton 
    return T * np.exp((rs0* (Lv+1.137e6*rs0) )/(Cp * T))   #David-Jones6.5

for nT, THw in enumerate(nThetaW):
    T = THw + T0
    Tz = np.copy(T)         #save surface temperature which will be iterated up the pressure levels
    #------variable latent heat of vapourization constants-------
    Lv = 3.139e6 - 2336 * T #(Koutsoyiannis 2011) - theoretically derived
    c2 = (Lv**2)/(Rv*Cp)    #[K^2]
    c3 = Lv/Cp      #[K]
    # es = f_es(T)
    grad0 = 0
    #------------------------------------------------------------
    print('Current temperature: %s' %THw)
    for nP,P in enumerate(Prange):
        #get dry adiabat
        arrayTd[nT,nP] = T*((P/P0)**c1) 

        # # #qet equivalent adiabat
        # arrayTe[nT,nP] = THe0*((P0/P)**c1)

        #get moist adiabat
        es = f_es(T+grad0)
        rs = f_rs(P,es)
        grad = dTdP(P,T+grad0)
        Pstep = P - PrangeList[nP+1]

       
        Tz = Tz + grad*Pstep
        arrayTw[nT,nP] = Tz

        # Tc = Tz-T0
        # if any((nThetaW-Tc)<0.0001):
        #     # Tidx = nThetaWlist.index(round(Tc,1))
        #     Tidx = nThetaWlist.index(int(Tc))
        #     WofTPreal[Tidx,nP] = T
        #     EofTPreal[Tidx,nP] = THe0*((P/P0)**c1)
        #     DofTPreal[Tidx,nP] = T*((P/P0)**c1)
        grad0 = grad*Pstep
    THe0 = Tz
    for nP,P in enumerate(Prange):
        #qet equivalent adiabat
        arrayTe[nT,nP] = THe0*((P/Ptop)**c1)

# arrayTnorm = (arrayTw - arrayTd)/(arrayTe - arrayTd)

# #monotonically select points every 0.1 kPa for fitting
# PrangeIdx = [np.argmin(abs(PrangeList - i)) for i in np.arange(P0,1,-0.1)]
# PrangeFit = Prange[PrangeIdx]

# #normailzing by one of the adiabats removes the non-linearity from the data
# THref = arrayTHnorm[-1,:]
# THref_fit = np.poly1d(np.polyfit(PrangeFit,THref[PrangeIdx],20))
# MAE_THref = np.mean(abs(THref[PrangeIdx]-THref_fit(PrangeFit)))
# print('MAE for polynomial fit of Tmax reference curve: %.2E' %MAE_THref)


#load array from other part
# arrayTHw = np.load('%s-%s.npy' %(Tmin,Tmax))
arrayTHw = np.load('-40.0-40.0.npy')
plt.contour(np.arange(-40,40),Prange[:-1799999],arrayTHw[:,:-1799999].T,nThetaWk[::5])
plt.plot(arrayTw[0::5,:].T-T0,Prange, 'r')
plt.gca().invert_yaxis()
plt.xlim([-40,100])
plt.show()





# # Now model,store coeffs and plot (for specified degree polynomial)
# print('Fitting polynomials to normalized curves')
# numterms = degree+1
# store_args = np.zeros((numterms,len(nThetaW)))
# tags = ['k%s' %(i+1) for i in range(numterms)]
# for i in range(len(nThetaW)):
#     main_pfit = np.poly1d(np.polyfit(PrangeFit,arrayTnorm[i,PrangeIdx],degree))
#     store_args[:,i] = main_pfit.coeffs
#     plt.plot(PrangeFit,main_pfit(PrangeFit),'r')
#     plt.plot(PrangeFit,arrayTnorm[i,PrangeIdx],'b')
#     # plt.xlim([0,1])
#     # plt.ylim([0,1])
# plt.show()
# plt.close()

# #now do fits for individual parameters
# # xvals = nThetaW[:-1]
# fitFCNs = []
# for iDeg in range(numterms):
#     # pfit = np.poly1d(np.polyfit(xvals,store_args[iDeg,:],25))
#     pfit = np.poly1d(np.polyfit(nThetaW,store_args[iDeg,:],25))
#     # MAE = np.mean(abs(store_args[iDeg,:] - pfit(xvals)))
#     MAE = np.mean(abs(store_args[iDeg,:] - pfit(nThetaW)))
#     print('%s MAE = %0.2E' %(tags[iDeg],MAE))
#     fitFCNs.append(pfit)

# #TESTING THE METHOD======================================
# print('Evaluating polynomial fit method....')
# arrayTHfit = np.empty_like(arrayTHw)
# # for nT, Temp in enumerate(xvals):

# for nT, T in enumerate(nThetaW):
#     print('.....current adiabat: %s' %T)
#     k = []
#     #calculate parameters 
#     for iDeg in range(numterms):
#         k.append(fitFCNs[iDeg](T))
#     #fit the moist adiabats
#     for nP,P in enumerate(Prange):
#         THrefm = THref_fit(P)        
#         THfit = 0.
#         #sum up the polynomial terms
#         for iDeg in range(numterms):
#             THfit = THfit + k[iDeg]*THrefm**(degree-iDeg)
#         arrayTHfit[nT,nP] = THfit

# arrayDiff = arrayTHnorm-arrayTHfit
# MAE = np.mean(abs(arrayDiff),0)



# #convert back to true adiabats
# arrayTHwm = np.empty_like(arrayTHw)
# arrayTHem = np.empty_like(arrayTHw)

# for nT, THw in enumerate(nThetaW):
#     T = THw + T0
#     #------variable latent heat of vapourization constants-------
#     Lv = 3.139e6 - 2336 * T #(Koutsoyiannis 2011) - theoretically derived
#     c2 = (Lv**2)/(Rv*Cp)    #[K^2]
#     c3 = Lv/Cp      #[K]
#     #------------------------------------------------------------
#     for nP,P in enumerate(Prange):
#         es = f_es(T)
#         rs0 = f_rs(P0,es)
#         # arrayTHem[nT,nP] = f_thE(T,rs0)
#         # arrayTHdm[nT,nP] = (arrayTHfit[nT,nP]*arrayTHem[nT,nP] - T)/(1+arrayTHfit[nT,nP])
#         # arrayTm[nT,nP] = arrayTHdm[nT,nP]*((P/P0)**c1)
#         THe = f_thE(T,rs0) * ((P/P0)**c1)
#         arrayTHem[nT,nP] = THe
#         THfit = arrayTHfit[nT,nP]
#         # THd = T + THfit*(THe-T)
#         THd = T*((P/P0)**c1)
#         # arrayTHwm[nT,nP] = THd*((P/P0)**c1)
#         arrayTHwm[nT,nP] = THfit*(THe-THd) + THd


# arrayDiffTemp = arrayTHw-arrayTHwm
# arrayDiffTHe = arrayTHe - arrayTHe_form

# #attempt to get THw from TP
# # for nT, THw in enumerate(nThetaW):
# #     T = THw + T0
# #     #------variable latent heat of vapourization constants-------
# #     Lv = 3.139e6 - 2336 * T #(Koutsoyiannis 2011) - theoretically derived
# #     c2 = (Lv**2)/(Rv*Cp)    #[K^2]
# #     c3 = Lv/Cp      #[K]
# #     #------------------------------------------------------------
# #     for nP,P in enumerate(Prange):
# #         es = f_es(T)
# #         rs0 = f_rs(P0,es)
# #         THe = f_thE(T,rs0)
# #         THfit = arrayTHfit[nT,nP]   ## THIS PART IS WRONG, GETS THnorm of THw not of T
# #         WofTPmodel[nT,nP] = ((P/P0)**c1)*(THfit*THe - T) / (THfit-1)
# # WofTPmodelC = WofTPmodel - T0
# # WofTPmodelC[np.isnan(WofTPreal)] = np.nan
# # WofTP_MAE = WofTPreal -  WofTPmodelC




# #=====================PLOTTING===========================
# #plot stuve's diagram
# plt.figure(figsize=(9,6))
# plt.title('SANITY CHECK: "STUVE" PLOT')
# plt.plot(arrayTw[-1,:]-T0,Prange, 'b',label='moist adiabat $\\theta_w$')
# plt.plot(arrayTd[-1,:]-T0,Prange, 'r--',label='dry adiabat $\\theta_d$')
# plt.plot(arrayTe[-1,:]-T0,Prange, 'g:',label='equivalent potential temperature $\\theta_e}$')
# plt.plot(arrayTw[0::10,:].T-T0,Prange, 'b')
# plt.plot(arrayTd[0::10,:].T-T0,Prange, 'r--')
# plt.plot(arrayTe[0::10,:].T-T0,Prange, 'g:')
# # plt.ylim([80,101])
# plt.gca().invert_yaxis()
# # plt.xlim([-100,100])
# plt.grid()
# plt.xlabel("temperature [C]")
# plt.ylabel("pressure [kPa]")
# plt.legend(loc='upper right',fontsize=12)
# # plt.savefig('./figs/stuve.pdf')
# plt.show()
# plt.close()

# #plot normalized adiabats
# plt.title('NORMALIZED SATURATED ADIABATS $\\theta_{norm}$')
# plt.plot(Prange,arrayTnorm[0::10,:].T, 'b')
# # plt.gca().invert_xaxis()
# # plt.ylim([0,1.1])
# plt.xlim([101,1])
# plt.ylabel('normalized temperature')
# plt.xlabel('pressure [kPa]')
# # plt.savefig('./figs/THnorm.pdf')
# plt.show()
# plt.close()

# #plot fit of single adiabat THref
# plt.title('$\\theta_{ref} = \\theta_{40C}$ POLYNOMIAL FIT')
# plt.plot(THref,Prange,'g')
# plt.plot(THref_fit(Prange),Prange,'r')
# plt.gca().invert_yaxis()
# plt.xlim([0,1.1])
# plt.ylim([101,1])
# plt.xlabel('normalized temperature')
# plt.ylabel('pressure [kPa]')
# plt.savefig('./figs/THref.pdf')
# plt.show()
# plt.close()

# #plot transformed adiabats
# plt.title('TRANSFORMED MOIST ADIABATS $\\theta_{trans}$')
# plt.plot(THref,arrayTHnorm[0:int(Tmax),:].T,'b')
# plt.plot(THref,arrayTHnorm[int(Tmax):,:].T,'r')
# plt.xlim([0,1])
# plt.xlabel('reference saturated adiabat')
# plt.ylim([0,1])
# plt.ylabel('transformated saturated adiabats')
# plt.savefig('./figs/THtrans.pdf')
# plt.show()
# plt.close()

# #subplot of fits for individual parameters
# fig = plt.figure(figsize=(10, 10)) 
# plt.suptitle('FIT PARAMETERS')
# for iDeg in range(degree):
#     plt.subplot(4,4,iDeg+1)
#     plt.title(tags[iDeg])
#     plt.xlabel('temperature [K]')
#     # plt.plot(xvals,store_args[iDeg,:],'g')
#     # plt.plot(xvals,fitFCNs[iDeg](xvals),'r')
#     plt.plot(nThetaW,store_args[iDeg,:],'g')
#     plt.plot(nThetaW,fitFCNs[iDeg](nThetaW),'r')
# # plt.tight_layout()
# plt.subplots_adjust(top = .92, hspace=0.3, wspace=0.2, left=0.05, right=0.97, bottom=0.05)
# plt.savefig('./figs/fit_params.pdf')
# plt.show()
# plt.close()

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