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
ThetaW = np.arange(-30.,40.)
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


def tanfit (P, a,b,c,d,e):
        # return a*np.arctan(b**P) * np.exp(-c*P)
        # return a * a*np.arctan(b**P)+ c*P + d
        return a* np.arctan(b**(P - d)) + c 
# p0 = (1.27133715e+00 ,  1.35999749e-03  , 6.99196638e+00 ,  1.00261648e+00)
p0 = (1.,1,1,1,1)
store_args = np.zeros((5,len(ThetaW)))
for i in range(len(ThetaW)):
# for i in range(6):
    popt, covp = curve_fit(tanfit,Prange,norm_adiabats[:,i],p0 = p0)
    store_args[:,i] = popt[:]
    plt.plot(Prange,tanfit(Prange,popt[0],popt[1],popt[2],popt[3],popt[4]),'r')
    plt.plot(Prange,norm_adiabats[:,i],'b')
    plt.show()
    p0 = popt[:]
    print popt


#TESTING THE METHOD======================================
fit_adiabats = np.empty_like(adiabats)
for nT, Temp in enumerate(xvals):
    k1,k2,k3,k4,k5,k6,k7,k8 = pfit1(Temp),pfit2(Temp),pfit3(Temp),pfit4(Temp),pfit5(Temp),pfit6(Temp),pfit7(Temp),pfit8(Temp)
    for nP,Pres in enumerate(Prange):
        normP = pref_fit(Pres)
        normTH = k1*normP**7 + k2*normP**6 + k3*normP**5 + k4*normP**4 + k5*normP**3 + k6*normP**2 + k7*normP + k8
        fit_adiabats[nP,nT] = normTH

plt.plot(norm_adiabats[:,1::5],color='0.5')
plt.plot(fit_adiabats[:,0::5],'r')
# plt.xlim([0,1])
plt.grid()
plt.xlabel("moist adiabats [C]")
plt.ylabel("pressure [kPa]")
plt.show()

