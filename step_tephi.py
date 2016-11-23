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
plt.title('SANITY CHECK: "STUVE" PLOT')
plt.plot(adiabats[:,-1]-T0,Prange, 'b',label='moist adiabat')
plt.plot(dry_adiabats[:,-1]-T0,Prange, 'r-',label='dry adiabat')
plt.plot(eq_adiabats[:,-1]-T0,Prange, 'g:',label='equivalent potential temperature')
plt.plot(adiabats[:,0::5]-T0,Prange, 'b')
plt.plot(dry_adiabats[:,0::5]-T0,Prange, 'r-')
plt.plot(eq_adiabats[:,0::5]-T0,Prange, 'g:')
# plt.ylim([60,100])
plt.gca().invert_yaxis()
plt.xlim([-20,20])
plt.grid()
plt.xlabel("moist adiabats [C]")
plt.ylabel("pressure [kPa]")
plt.legend(loc='upper right')
plt.savefig('stuve.pdf')
plt.show()

#plot normalized adiabats
plt.title('NORMALIZED MOIST ADIABATS')
norm_adiabats = (adiabats - dry_adiabats)/((eq_adiabats - dry_adiabats))
plt.plot(norm_adiabats[:,0::10],Prange, 'b')
plt.xlabel('normalized $\Theta_w$')
plt.ylabel('pressure [kPa]')
plt.savefig('norm_theta.pdf')
plt.show()


#normailzing by one of the adiabats removes the non-linearity from the data
plt.title('40C ADIABAT FIT')
Pref = norm_adiabats[:,-1]
pref_fit = np.poly1d(np.polyfit(Prange,Pref,28))
plt.plot(Prange,Pref,'g')
plt.plot(Prange,pref_fit(Prange),'r')
plt.savefig('pref_fit.pdf')
plt.show()
print(sum(abs(Pref-pref_fit(Prange))))


plt.title('TRANSFORMED MOIST ADIABATS')
plt.plot(Pref,norm_adiabats[:,0:40],'b')
plt.plot(Pref,norm_adiabats[:,40:],'r')
# plt.plot(norm_adiabats[:,0:50],Pref,'r')
# plt.plot(norm_adiabats[:,50:],Pref ,'b')
plt.xlim([0,1])
plt.xlabel('normalized pressure')
plt.ylim([0,1])
plt.ylabel('ref adiabat')
plt.savefig('trans_adiabats.pdf')
plt.show()


# def tanfit (P, a,b,c,d):
#         # return a*np.arctan(b**P) * np.exp(-c*P)
#         return a * a*np.arctan(b**P)+ c*P + d
# # p0 = (1.27133715e+00 ,  1.35999749e-03  , 6.99196638e+00 ,  1.00261648e+00)
# p0 = (1.,1,1,1)
# store_args = np.zeros((4,len(ThetaW)-1))
# for i in range(len(ThetaW)-1):
# # for i in range(6):
#     popt, covp = curve_fit(tanfit,Pref,norm_adiabats[:,i+1],p0 = p0)
#     store_args[:,i] = popt[:]
#     plt.plot(Pref,tanfit(Pref,popt[0],popt[1],popt[2],popt[3]),'r')
#     plt.plot(Pref,norm_adiabats[:,i+1],'b')
#     plt.xlim([0,1])
#     plt.ylim([0,1])
#     plt.show()
#     p0 = popt[:]
#     print popt


store_args = np.zeros((9,len(ThetaW)-1))
for i in range(len(ThetaW)-1):
# for i in range(6):
    main_pfit = np.poly1d(np.polyfit(Pref,norm_adiabats[:,i],8))
    store_args[:,i] = main_pfit.coeffs
    plt.plot(Pref,main_pfit(Pref),'r')
    plt.plot(Pref,norm_adiabats[:,i],'b')
    plt.xlim([0,1])
    plt.ylim([0,1])
    plt.show()


fig = plt.figure(figsize=(14, 14)) 
plt.suptitle('FIT PARAMETERS')
xvals = ThetaW[:-1]
#fits for individual parameters
pfit1 = np.poly1d(np.polyfit(xvals,store_args[0,:],23))
plt.subplot(3,3,1)
plt.title('k1')
plt.plot(xvals,store_args[0,:],'g')
plt.plot(xvals,pfit1(xvals),'r')
print sum(abs(store_args[0,:] - pfit1(xvals)))

plt.subplot(3,3,2)
plt.title('k2')
pfit2 = np.poly1d(np.polyfit(xvals,store_args[1,:],24))
plt.plot(xvals,store_args[1,:],'g')
plt.plot(xvals,pfit2(xvals),'r')
print sum(abs(store_args[1,:] - pfit2(xvals)))

plt.subplot(3,3,3)
plt.title('k3')
pfit3 = np.poly1d(np.polyfit(xvals,store_args[2,:],24))
plt.plot(xvals,store_args[2,:],'g')
plt.plot(xvals,pfit3(xvals),'r')
print sum(abs(store_args[2,:] - pfit3(xvals)))

plt.subplot(3,3,4)
plt.title('k4')
pfit4 = np.poly1d(np.polyfit(xvals,store_args[3,:],24))
plt.plot(xvals,store_args[3,:],'g')
plt.plot(xvals,pfit4(xvals),'r')
print sum(abs(store_args[3,:] - pfit4(xvals)))

plt.subplot(3,3,5)
plt.title('k5')
pfit5 = np.poly1d(np.polyfit(xvals,store_args[4,:],24))
plt.plot(xvals,store_args[4,:],'g')
plt.plot(xvals,pfit5(xvals),'r')
print sum(abs(store_args[4,:] - pfit5(xvals)))

plt.subplot(3,3,6)
plt.title('k6')
pfit6 = np.poly1d(np.polyfit(xvals,store_args[5,:],24))
plt.plot(xvals,store_args[5,:],'g')
plt.plot(xvals,pfit6(xvals),'r')
print sum(abs(store_args[5,:] - pfit6(xvals)))

plt.subplot(3,3,7)
plt.title('k7')
pfit7 = np.poly1d(np.polyfit(xvals,store_args[6,:],21))
plt.plot(xvals,store_args[6,:],'g')
plt.plot(xvals,pfit7(xvals),'r')
print sum(abs(store_args[6,:] - pfit7(xvals)))

plt.subplot(3,3,8)
plt.title('k8')
pfit8 = np.poly1d(np.polyfit(xvals,store_args[7,:],17))
plt.plot(xvals,store_args[7,:],'g')
plt.plot(xvals,pfit8(xvals),'r')
print sum(abs(store_args[7,:] - pfit8(xvals)))

plt.subplot(3,3,9)
plt.title('k9')
pfit9 = np.poly1d(np.polyfit(xvals,store_args[8,:],17))
plt.plot(xvals,store_args[8,:],'g')
plt.plot(xvals,pfit9(xvals),'r')
print sum(abs(store_args[8,:] - pfit9(xvals)))

plt.savefig('fit_params.pdf')
plt.show()



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
