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
ThetaW = np.arange(-40.,40.)
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


store_args = np.empty((9,len(ThetaW)))
error = 0.
for i in range(len(ThetaW)):
    pfit_main = np.poly1d(np.polyfit(Prange,norm_adiabats[:,i],8))
    store_args[:,i] = pfit_main.coeffs
    plt.plot(Prange,norm_adiabats[:,i],'g')
    plt.plot(Prange,pfit_main(Prange),'r')
    plt.show()
    error = error + np.std(norm_adiabats[:,i] - pfit_main(Prange))
error = error/len(ThetaW)
print error



fig = plt.figure(figsize=(8, 6)) 
xvals = ThetaW
#fits for individual parameters
pfit1 = np.poly1d(np.polyfit(xvals,store_args[0,:],15))
plt.subplot(4,2,1)
plt.title('k1')
plt.plot(xvals,store_args[0,:],'g')
plt.plot(xvals,pfit1(xvals),'r')
print np.std(store_args[0,:] - pfit1(xvals))

plt.subplot(4,2,2)
plt.title('k2')
pfit2 = np.poly1d(np.polyfit(xvals,store_args[1,:],15))
plt.plot(xvals,store_args[1,:],'g')
plt.plot(xvals,pfit2(xvals),'r')
print np.std(store_args[1,:] - pfit2(xvals))

plt.subplot(4,2,3)
plt.title('k3')
pfit3 = np.poly1d(np.polyfit(xvals,store_args[2,:],15))
plt.plot(xvals,store_args[2,:],'g')
plt.plot(xvals,pfit3(xvals),'r')
print np.std(store_args[2,:] - pfit3(xvals))

plt.subplot(4,2,4)
plt.title('k4')
pfit4 = np.poly1d(np.polyfit(xvals,store_args[3,:],15))
plt.plot(xvals,store_args[3,:],'g')
plt.plot(xvals,pfit4(xvals),'r')
print np.std(store_args[3,:] - pfit4(xvals))

plt.subplot(4,2,5)
plt.title('k5')
pfit5 = np.poly1d(np.polyfit(xvals,store_args[4,:],15))
plt.plot(xvals,store_args[4,:],'g')
plt.plot(xvals,pfit5(xvals),'r')
print np.std(store_args[4,:] - pfit5(xvals))

plt.subplot(4,2,6)
plt.title('k6')
pfit6 = np.poly1d(np.polyfit(xvals,store_args[5,:],15))
plt.plot(xvals,store_args[5,:],'g')
plt.plot(xvals,pfit6(xvals),'r')
print np.std(store_args[5,:] - pfit6(xvals))

plt.subplot(4,2,7)
plt.title('k7')
pfit7 = np.poly1d(np.polyfit(xvals,store_args[6,:],15))
plt.plot(xvals,store_args[6,:],'g')
plt.plot(xvals,pfit7(xvals),'r')
print np.std(store_args[6,:] - pfit7(xvals))

plt.subplot(4,2,8)
plt.title('k8')
pfit8 = np.poly1d(np.polyfit(xvals,store_args[7,:],15))
plt.plot(xvals,store_args[7,:],'g')
plt.plot(xvals,pfit8(xvals),'r')
print np.std(store_args[7,:] - pfit8(xvals))
plt.show()

pfit9 = np.poly1d(np.polyfit(xvals,store_args[8,:],15))
plt.plot(xvals,store_args[8,:],'g')
plt.plot(xvals,pfit9(xvals),'r')
print np.std(store_args[8,:] - pfit9(xvals))

#-----------------------------


# #TESTING THE METHOD======================================
fit_adiabats = np.empty_like(adiabats)
for nT, Temp in enumerate(xvals):
    k1,k2,k3,k4,k5,k6,k7,k8,k9 = pfit1(Temp),pfit2(Temp),pfit3(Temp),pfit4(Temp),pfit5(Temp),pfit6(Temp),pfit7(Temp),pfit8(Temp),pfit9(Temp)
    for nP,Pres in enumerate(Prange):
        normP = Pres
        normTH = k1*normP**8 + k2*normP**7 + k3*normP**6 + k4*normP**5 + k5*normP**4 + k6*normP**3 + k7*normP**2 + k8*normP**1 + k9
        fit_adiabats[nP,nT] = normTH

plt.plot(norm_adiabats[:,1::5],color='0.5')
plt.plot(fit_adiabats[:,0::5],'r')
# plt.xlim([0,1])
plt.grid()
plt.xlabel("moist adiabats [C]")
plt.ylabel("pressure [kPa]")
plt.show()

