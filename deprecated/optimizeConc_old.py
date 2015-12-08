import matplotlib
matplotlib.use('TKAgg');
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import stats

from lbl2xs import lbl_xs, get_lbl_data, NoneType, check_pT
from pairTypes import Interval
from cgsUnits import unitConversion
import numpy as np

datFile = '/media/win/Users/jaymz/Desktop/new_data/overnight-7-10_dataDump.csv'
vacFile = '/media/win/Users/jaymz/Desktop/new_data/vacuum-7-13_dataDump.csv'

fitScale = False

iFiles=['test/testH2O','test/testCO2']
oFile=None
initialGuess = [70000.0 ,5000.0]
commentChar='#'
molecule=None
xLimits=Interval(6240.0, 6255.0)
temperature=298.0
ptFile=None
lineShape='Voigt'
airWidth=0.1
sampling=5.0
wingExt=5.0
interpolate='3'
nGrids=1
gridRatio=8
nWidths=25.0,
format=None
verbose=False

###Load Experimental Data
dat = []

file = open(datFile,'r')

for line in file:
    l = str.split(line,',')
    ln = []
    for thing in l:
        ln.append(float(thing))
    dat.append(ln)
    
dat = np.matrix(dat)

data = {'x' : np.array((dat[:,0].flatten().tolist())[0]),
        'x(wn)' : np.array(((float(10**7)/dat[:,0]).flatten().tolist())[0]),
        'y' : np.array((dat[:,1].flatten().tolist())[0]), 
        'ringdownTime' : np.array(((-1.0/dat[:,0]).flatten().tolist())[0]) }
file.close()

dat = []

file = open(vacFile,'r')

for line in file:
    l = str.split(line,',')
    ln = []
    for thing in l:
        ln.append(float(thing))
    dat.append(ln)
    
dat = np.matrix(dat)

dataVac = {'x' : np.array((dat[:,0].flatten().tolist())[0]),
        'x(wn)' : np.array(((float(10**7)/dat[:,0]).flatten().tolist())[0]),
        'y' : np.array((dat[:,1].flatten().tolist())[0]), 
        'ringdownTime' : np.array(((-1.0/dat[:,0]).flatten().tolist())[0]) }
file.close()

bSlope, bIntercept, r_val, p_val, std_err = stats.linregress(dataVac['x(wn)'],dataVac['y'])
def baselineFit(wn):
    return bIntercept + bSlope*wn
c = 3.0*10**8
def alpha(decayConst,wavenumber):
    return 1.0/c*(-decayConst+baselineFit(wavenumber))

alph = []
for i in range(0,len(data['x'])):
    alph.append(alpha(data['y'][i],data['x(wn)'][i]))

data['alpha'] = np.array(alph)

###Load Stuff from files generated by extract.py
Line_Data, Mol_Data = [], []
for i in range(0,len(iFiles)):
    l, m = get_lbl_data (iFiles[i], xLimits, airWidth, wingExt, molecule, commentChar)
    if len(l['position'])<1:
        print 'no line data in ', xLimits; sys.exit();
    Line_Data.append(l)
    Mol_Data.append(m)

# interpolation method used in multigrid algorithms: only Lagrange!
if   interpolate in 'bBsScC':  lagrange=4
elif interpolate in 'qQ':      lagrange=3
elif interpolate in 'lL':      lagrange=2
else:                          lagrange=int(interpolate)

xMin, xMax = xLimits.limits()

#def fitScales():
#    A = np.column_stack((np.ones(len(Tnew)), Inew, Tnew))
#    c, resid,rank,sigma = np.linalg.lstsq(A,Wnew)
#    return c

def fitErr(pressures):
    crossSections = []
    sum = np.array([0.0]*len(data['x']))
    for i in range(0,len(Line_Data)):
        factor = pressures[i]/1013500.0/8.31/temperature*6.022*10**23
        crossSections.append(lbl_xs (Line_Data[i], Mol_Data[i], np.array([pressures[i]]), np.array([temperature]), xLimits, lineShape, sampling, wingExt, nGrids, gridRatio, nWidths, lagrange, verbose)[0]  )
        xarr = np.linspace(xMin,xMax,len(crossSections[-1]['y']))
        interp = interp1d(xarr,crossSections[-1]['y']*factor)
        sum = sum + interp(data['x(wn)'])
    
    diff = sum - data['alpha']
    return diff.dot(diff)

def fitErrScale(pressures):
    crossSections = []
    sum = np.array([0.0]*len(data['x']))
    interps = []
    for i in range(0,len(Line_Data)):
        factor = pressures[i]/1013500.0/8.31/temperature*6.022*10**23
        crossSections.append(lbl_xs (Line_Data[i], Mol_Data[i], np.array([pressures[i]]), np.array([temperature]), xLimits, lineShape, sampling, wingExt, nGrids, gridRatio, nWidths, lagrange, verbose)[0]  )
        xarr = np.linspace(xMin,xMax,len(crossSections[-1]['y']))
        interp = interp1d(xarr,crossSections[-1]['y']*factor)
        interps.append(interp(data['x(wn)']))
    
    A = np.column_stack(interps)
    c, resid,rank,sigma = np.linalg.lstsq(A,data['alpha'])
    
    for i in range(0,len(Line_Data)):
        sum = sum + interps[i]*abs(c[i])
    
    diff = sum - data['alpha']
    return diff.dot(diff)
    

from scipy.optimize import minimize
if(fitScale):
    minResult = minimize(fitErrScale,initialGuess,method='Nelder-Mead')
else:
    minResult = minimize(fitErr,initialGuess,method='Nelder-Mead')

#plot data    
plt.plot(data['x(wn)'],data['alpha'],'ro')
crossSections = []
xplt = np.linspace(xMin,xMax,2000)
interps = []
interpFuncs = []
sum = np.array([0.0]*len(xplt))
for i in range(0,len(Line_Data)):
    factor = minResult['x'][i]/1013500.0/8.31/temperature*6.022*10**23
    crossSections.append(lbl_xs (Line_Data[i], Mol_Data[i], np.array([minResult['x'][i]]), np.array([temperature]), xLimits, lineShape, sampling, wingExt, nGrids, gridRatio, nWidths, lagrange, verbose)[0]  )
    xarr = np.linspace(xMin,xMax,len(crossSections[-1]['y']))
    interp = interp1d(xarr,crossSections[-1]['y']*factor)
    interpFuncs.append(interp)
    interps.append(interp(data['x(wn)']))
if(fitScale):
    A = np.column_stack(interps)
    c, resid,rank,sigma = np.linalg.lstsq(A,data['alpha'])
else:
    c = [1.0]*len(Line_Data)
for i in range(0,len(Line_Data)):
    factor = minResult['x'][i]/1013500.0/8.31/temperature*6.022*10**23
    sum = sum + interpFuncs[i](xplt)*c[i]
    xarr = np.linspace(xMin,xMax,len(crossSections[i]['y']))
    plt.plot(xarr,crossSections[i]['y']*factor*c[i])
#plot total hitran absorbance
plt.plot(xplt,sum)
plt.axis([min(data['x(wn)']), max(data['x(wn)']), 0.0, max(data['alpha'])])

print('Pressures (mb): ' + str(minResult['x']/1000) + '\n')
print('Scale fits' + str(c) +'\n')

plt.show()