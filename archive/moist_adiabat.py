import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.optimize import curve_fit


a = 0.285611 #R/cp
b = 1.35e7	#K2
c = 2488.4 	#K
LvRv = 5422. #K
T0 = 273.15
Eps = 0.622 #gv/gd
e0 = 0.611 	#kPa


Prange = np.arange(100,1, -0.01)		#defined until 1, otherwise potential temperature formula blows up
ThetaW = np.arange(-50,25)				#specific heat of air changes above 40

adiabats = np.empty((len(Prange),len(ThetaW)))
dry_adiabats = np.empty_like(adiabats)

def f_es(T):
	return  e0*np.exp(LvRv*(1./T0 - 1./T))

def f_rs(P,es):
	return  Eps*es / (P - es)

def dTdP(P,T):
	return (a*T + c*rs)/(P*(1+(b*rs/T**2)))



for nT, Temp in enumerate(ThetaW):
	T = Temp + 273.15
	for nP,Pres in enumerate(Prange):
		es = f_es(T)
		rs = f_rs(Pres,es)
		grad = dTdP(Pres,T)
		T = T - grad*0.01
		adiabats[nP,nT] = T
		#dry
		dry_adiabats[nP,nT] = T*((100./Pres)**a)	
	print Temp


plt.plot(adiabats[:,0::10]-273.15,Prange)
plt.gca().invert_yaxis()
plt.show()



#adjust the top part of the adiabat
a3 = 2600.
def f_thE(ThetaW,rs0):
	return ThetaW * np.exp(a3*rs0/ThetaW)

ThetaE = np.empty((len(ThetaW)))
for nT, Temp in enumerate(ThetaW):
	T = Temp + 273.15
	es0 = f_es(T)
	rs0 = f_rs(100,es0)
	ThetaE[nT] = f_thE(T,rs0)

# top_adiabats = adiabats/ThetaE
# plt.plot(top_adiabats[:,0::10]-273.15,Prange)
# plt.gca().invert_yaxis()
# plt.show()




norm_adiabats = (adiabats - dry_adiabats)/ThetaE
plt.plot(norm_adiabats[:,0::10])
plt.show()

spread = np.std(norm_adiabats, axis = 1)
plt.plot(spread)
plt.show()
print sum(spread)


# #find an equation for the bottom
# P0 = np.mean(top_adiabats[-1,:])


# c0 = [1.002,59.,-1.]


# testP = np.arange(0,9000,1000)
# for nP,Pval in enumerate(testP):

# 	def fix_bottom(xvals,base,xshift,ystretch):
# 		T = xvals
# 		return ystretch*(base**((T-xshift)/P0) - 1.)

# 	out = curve_fit(fix_bottom,ThetaW,top_adiabats[Pval,:],p0=c0)
# 	params = tuple(out[0])
# 	mean_fit = fix_bottom(ThetaW,params[0],params[1],params[2])
# 	plt.plot(mean_fit)
# 	plt.plot(top_adiabats[Pval,:])
# 	plt.show()
# 	print params


'''
error (sum of RMSE for all tested Tw values) for the similarity plot can be reduced by reducing the value of a3. 
The error shifts to top of the atmosphere values (minimizes in the center) - overall lower. Doesn't mean ideal - 
but could be useful in particular circumstances

All error is from higher temperature values - cause, constants no longer hold. 
'''
