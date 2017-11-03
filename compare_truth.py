#INPUT
Tmin=-100.
Tmax=100.
THmin =-68. #because we don't start at standard pressure -86 is ~70C moist adiabat
THmax= 42.
Pbot = 105
Ptop = 1    #kPa - upper atmosphere limit surface
Plim = 1
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
import os.path

P0 = 100.      #kPa

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

datafileB = '%s-%s_B.npy' %(Tmin,Tmax)
datafileT = '%s-%s_T.npy' %(Tmin,Tmax)

dataB = np.load(datafileO)
dataT = np.load(datafileT)

#monotonically select points every 0.1 kPa for fitting
PrangeIdx = [np.argmin(abs(PrangeList - i)) for i in np.arange(Pbot,Plim,-0.1)]
PrangeFit = Prange[PrangeIdx]
P0idx = np.argmin(abs(Prange - P0))
Tidx = [np.argmin(abs(nTw - i)) for i in nThetaW]

subarrayTHwB = dataB[Tidx,:]
subarrayTHwT = dataT[Tidx,:]

fitarrayB = subarrayTHwB[:,PrangeIdx]
fitarrayT = subarrayTHwT[:,PrangeIdx]

diff = fitarrayT - fitarrayB
truthMAE = MAE = np.mean(abs(diff.ravel()))
print('TRUTH MAE: %s' %truthMAE) 

plt.imshow(diff.T,aspect='auto',origin='lower',cmap='cubehelix_r',vmin=0,vmax=0.08)
plt.xlabel("$\\theta_w$ [$^\circ$C]")
plt.ylabel("pressure [kPa]")
ax = plt.gca()
ax.set_xticks(np.arange(0,len(nThetaW),20))
ax.set_xticklabels(C0axis[::20].astype(int))
ax.set_yticks(PaxisIdx)
ax.set_yticklabels(PrangeFit[PaxisIdx].round())
cbar = plt.colorbar(format='%.2f')
cbar.set_label('temperature difference [K]')
plt.title('SENSITIVITY TO CONSTANTS')
plt.savefig('./figs/TruthSensitivity.pdf')
plt.show()