#import necessary packages 
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, minimize
from scipy import spatial

# http://www.researchgate.net/publication/229942047_A_new_formula_for_latent_heat_of_vaporization_of_water_as_a_function_of_temperature

#thermodynamic constants
Rd = 287. 	#[J K^-1 kg^-1] gas constant for dry air
Rv = 461. 	#[J K^-1 kg^-1] gas constant for water vapour
Cp = 1005. 	#[J K^-1 kg^-1] specific heat of dry air at constant pressure
# Lv = 2.5e6 	#latent heat of vapourization at standard temperature
T0 = 273.16		#standard temperature
e0 = 0.611657 	#kPa: adjusted Clausius-Clayperon constant (Koutsoyiannis 2011)


#derived constants
Eps = Rd/Rv 	#dimensionless
c1 = Rd/Cp 		#dimensionless
# c2 = (Lv**2)/(Rv*Cp) 	#[K^2]
# c3 = Lv/Cp 		#[K]
c4 = 2490. 		#[K kg_air kg_vapour^-1]


Prange = np.arange(100,1, -0.001)
ThetaW = np.arange(-45.,45.)
adiabats = np.empty((len(Prange),len(ThetaW)))
dry_adiabats = np.empty_like(adiabats)
eq_adiabats = np.empty_like(adiabats)

def f_es(T):
    #REPLACING STANDARD EQUATION WITH Koutsoyiannis 2011
    return  e0*np.exp(24.921*(1.-(T0/T)))*((T0/T)**5.06)
def f_rs(P,es):
    return  Eps*es / (P - es)
def dTdP(P,T):
    return (c1*T + c3*rs)/(P*(1+(c2*rs/T**2)))
def f_thE(T,rs0):
    # return T * np.exp(c4*rs0/T)
    return T + c4*rs0
    #second formula produces smaller error for positive ThetaW values


for nT, Temp in enumerate(ThetaW):
    T = Temp + T0
    #------variable latent heat of vapourization constants-------
    Lv = (2500.8 - 2.36*Temp + 0.0016*(Temp**2) - 0.00006*(Temp**3))*1000
    c2 = (Lv**2)/(Rv*Cp) 	#[K^2]
    c3 = Lv/Cp 		#[K]
    #------------------------------------------------------------
    print('Current adiabat: %s' %Temp)
    
    for nP,Pres in enumerate(Prange):
        #get dry adiabat
        dry_adiabats[nP,nT] = (Temp+T0)*((Pres/100.)**c1) 
        #get moist adiabat
        es = f_es(T)
        rs = f_rs(Pres,es)
        grad = dTdP(Pres,T)
        T = T - grad*0.001
        adiabats[nP,nT] = T
        #qet equivalent adiabat
        rs0 = f_rs(100.,es)
        eq_adiabats[nP,nT] = f_thE(T,rs0)

#plot stuve's diagram
plt.plot(adiabats[:,0::5]-T0,Prange, 'b')
plt.plot(dry_adiabats[:,0::5]-T0,Prange, 'r-')
plt.plot(eq_adiabats[:,0::5]-T0,Prange, 'g:')
plt.gca().invert_yaxis()
plt.xlim([-50,20])
plt.grid()
plt.xlabel("moist adiabats [C]")
plt.ylabel("pressure [kPa]")
plt.show()

#plot normalized adiabats
norm_adiabats = (adiabats - dry_adiabats)/((eq_adiabats - dry_adiabats))
plt.plot(norm_adiabats[:,0::10],Prange/100.)
plt.xlabel('normalized Theta')
plt.ylabel('pressure/100.')
plt.show()

#normailzing by one of the adiabats removes the non-linearity from the data
Pref = norm_adiabats[:,0]
def moist_fit(P, a,b,c,d,e,f,g):
    return a*P**6 + b*P**5 + c*P**4 + d*P**3 + e*P**2 + f*P + g
p0 = (5E-12,- 2E-09,4E-07,- 3E-05, 0.0005, 0.024,0.0053)

popt, pcov = curve_fit(moist_fit, Prange, Pref,p0=p0)
a,b,c,d,e,f,g = popt[:]
plt.plot(Prange,moist_fit(Prange, a,b,c,d,e,f,g),'r')
plt.plot(Prange,Pref,'b')
plt.xlabel('pressure range')
plt.ylabel('P_ref curve')
plt.show()
k0_params = popt

# plt.plot(Pref,norm_adiabats[:,0:50],'r')
# plt.plot(Pref,norm_adiabats[:,50:],'b')
plt.plot(norm_adiabats[:,0:50],Pref,'r')
plt.plot(norm_adiabats[:,50:],Pref ,'b')
plt.xlim([0,1])
plt.xlabel('normalized pressure')
plt.ylim([0,1])
plt.ylabel('ref adiabat')
plt.show()


#Fit biexponential to the family of curves
def biexp(X, a,b):
     return b*(1-np.exp(-a*X))
store_args = np.empty((2,len(ThetaW)-1))
for i in range(len(ThetaW)-1):
    xdata = norm_adiabats[:,i+1]
    popt, covp = curve_fit(biexp,xdata,Pref)
    store_args[:,i] = popt[:]
    plt.plot(xdata,biexp(xdata, popt[0],popt[1]),'r')
    plt.plot(xdata,Pref,'b')
    plt.xlim([0,1])
    plt.ylim([0,1])
    plt.show()
    # p0 = (popt)
    print popt
    k3 = popt

#-----------------------------

#FIRST VAR======================================

#this works a bit better (std = 0.004)
in_theta = ThetaW[1:]
pfit = np.poly1d(np.polyfit(in_theta,store_args[0,:],5))
plt.plot(in_theta,store_args[0,:],'g')
plt.plot(in_theta,pfit(in_theta),'r')
print np.std(store_args[0,:] - pfit(in_theta))
plt.show()
k1_params = popt


#SECOND VAR======================================
#this works ok (std = 0.06), test if biexponential works better
in_itheta = ThetaW[1:] + 45.
def fit3(X, a,b,c,d):
    return a/((X-d)**b) + c
p0 = (1.,1.,1.,-1)
popt, covp = curve_fit(fit3,in_itheta,store_args[1,:],p0=p0)
plt.plot(in_itheta,store_args[1,:],'g')
plt.plot(in_itheta,fit3(in_itheta, popt[0],popt[1],popt[2],popt[3]),'r')
print np.std(store_args[1,:] - fit3(in_itheta, popt[0],popt[1],popt[2],popt[3]))
k2_params = popt
plt.show()


#TESTING THE METHOD======================================

fit_adiabats = np.empty_like(adiabats)
for nT, Temp in enumerate(ThetaW[1:]):
    k1 = pfit(Temp)
    THi = Temp+45.
    k2 = fit3(THi,k2_params[0],k2_params[1],k2_params[2],k2_params[3])
    print Temp, k1, k2, store_args[:,nT]
    for nP,Pres in enumerate(Prange):
        normP = moist_fit(Pres, a,b,c,d,e,f,g)
        normTH = (1./-k1)*np.log(1-(normP/k2))
        fit_adiabats[nP,nT] = normTH

plt.plot(norm_adiabats[:,1::5],Pref, color='0.5')
plt.plot(fit_adiabats[:,0::5],Pref, 'r')
# plt.xlim([0,1])
plt.grid()
plt.xlabel("moist adiabats [C]")
plt.ylabel("pressure [kPa]")
plt.show()


