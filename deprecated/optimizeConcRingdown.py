import matplotlib
matplotlib.use('TKAgg');   #makes the plots work on some systems
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d
from scipy import stats
from scipy.optimize import minimize

from string import *
from types  import *
import sys
import os.path
import time 

from IO import open_outFile, commonExtension, paste, writeArray, grep_from_header
from pairTypes import Interval
from lbl2xs import lbl_xs
from lines import get_lbl_data
from xs2ac2od import sum_xsTimesDensity, integral_acTimesDistance
from lagrange_interpolation import *
from atmos import Atmos1D, get_profile_list, reGridAtmos1d, Profile
from cgsUnits import unitConversion
from molecules import molecules
molecNames   = molecules.keys()

##Set Parameters
datFile = '/media/win/Users/jaymz/Desktop/new_data/overnight-7-10_dataDump.csv'
vacFile = '/media/win/Users/jaymz/Desktop/new_data/vacuum-7-13_dataDump.csv'

calcPressure = False
calcTemperature = False
findOffset = True

lineFiles=['test/testH2O','test/testCO2']
gasSpecies = ['H2O','CO2']
concGuess = [.007, 330.0E-5]
temperature = 295.0
pressure = 1.0149E6*760/760
xLimits=Interval(6240.0, 6265.0)

##default stuff
outFile=None
commentChar='#'
zToA=None
lineShape='Voigt'
airWidth=0.1
sampling=5.0
wingExt=5.0
interpolate='3'
nGrids=2
gridRatio=8
nWidths=25.0
mode='d'
quad='T'
nanometer=False
yOnly=False
flipUpDown=False
plot=False
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

data['alpha'] = np.array(alph)/60

##
if(calcPressure & calcTemperature):
    concGuess = concGuess + [pressure] + [temperature]
elif(calcTemperature):
    concGuess = concGuess + [temperature]
elif(calcPressure):
    concGuess = concGuess + [pressure]

# parse file header and retrieve molecule names
species = [grep_from_header(file,'molecule') for file in lineFiles]			

print 'lbl2od for',  len(species), ' molecules:', join(species)


# interpolation method used in multigrid algorithms: only Lagrange!
if   interpolate in 'bBsScC':  lagrange=4
elif interpolate in 'qQ':      lagrange=3
elif interpolate in 'lL':      lagrange=2
else:                          lagrange=int(interpolate)

# loop over molecules:  read its line parameters (and some other data like mass) and compute cross sections for all p,T
lD, mD = [], []
for file in lineFiles:
    # read a couple of lines and convert to actual pressure and temperature
    Line_Data, Mol_Data = get_lbl_data (file, xLimits, airWidth, wingExt, commentChar=commentChar)
    lD.append(Line_Data)
    mD.append(Mol_Data)

def computeAbsorption(pressure, temperature, gasSpecies, concValues, Line_Data, Mol_Data):
    profs = [Profile([0.0],      what='altitude',    unit = 'km'),
                 Profile([pressure],  what='pressure',    unit = 'g/cm/s**2'),
                 Profile([temperature],    what='temperature', unit = 'K')]
    for i in range(0,len(gasSpecies)):
        profs.append(Profile([concValues[i]], what=gasSpecies[i], unit='ppV'))
    atmos = Atmos1D(profs)    
    xsMatrix = []    
    for i in range(0,len(gasSpecies)):
        # compute cross sections for selected line shape, returns a list with npT dictionaries
        molCrossSections = lbl_xs (Line_Data[i], Mol_Data[i], atmos.p, atmos.T, xLimits, lineShape, sampling, wingExt, nGrids, gridRatio, nWidths, lagrange, verbose)
        xsMatrix.append(molCrossSections)
        
    # compute absorption coefficients, also return uniform, equiidistant wavenumber grid corresponding to most dense cross section
    
    vGrid, absCoeff = sum_xsTimesDensity (xsMatrix, atmos.densities, interpolate=interpolate, verbose=verbose)
    return vGrid, absCoeff[:,0]
    
def computeFitError(pressure, temperature, gasSpecies, concValues, Line_Data, Mol_Data):
    vGrid, absCoeff = computeAbsorption(pressure, temperature, gasSpecies, concValues, Line_Data, Mol_Data)
    interp = interp1d(vGrid, absCoeff)
    diff = -interp(data['x(wn)']) + data['alpha']
    if(findOffset):
        diff = diff - sum(diff)/len(diff)
    return diff.dot(diff)
    
def funcToMinimize(x,Line_Data,Mol_Data):
    if(calcPressure & calcTemperature):
        out = computeFitError(x[-2], x[-1], gasSpecies, x[0:len(x)-2],Line_Data, Mol_Data)
    elif(calcTemperature):
        out = computeFitError(pressure, x[-1], gasSpecies, x[0:len(x)-1],Line_Data, Mol_Data)
    elif(calcPressure):
        out = computeFitError(x[-1], temperature, gasSpecies, x[0:len(x)-1],Line_Data, Mol_Data)
    else:
        out = computeFitError(pressure, temperature, gasSpecies, x,Line_Data, Mol_Data)
    if(any(n < 0 for n in x)):   #used to avoid negative values in the fit by giving them an artificially large error.
        return 100*out
    return out
    
res = minimize(funcToMinimize,concGuess,method='Nelder-Mead',args=(lD,mD))


if(calcTemperature & calcPressure):
    x, y = computeAbsorption(res['x'][-2], res['x'][-1], gasSpecies, res['x'][0:len(res['x'])-2],lD,mD)
    concVals = res['x'][0:len(res['x'])-2]
elif(calcTemperature):
    x, y = computeAbsorption(pressure, res['x'][-1], gasSpecies, res['x'][0:len(res['x'])-1],lD,mD)
    concVals = res['x'][0:len(res['x'])-1]
elif(calcPressure):
    x, y = computeAbsorption(res['x'][-1], temperature, gasSpecies, res['x'][0:len(res['x'])-1],lD,mD)
    concVals = res['x'][0:len(res['x'])-1]
else:
    x, y = computeAbsorption(pressure, temperature, gasSpecies, res['x'],lD,mD)
    concVals = res['x']

offset = 0
if(findOffset):
    interpFunc = interp1d(x,y)
    diff = data['alpha'] - interpFunc(data['x(wn)'])
    offset = sum(diff)/len(diff)
    data['alpha'] = data['alpha'] - offset

print('------------------\nConcentrations (ppm):\n' + repr(gasSpecies) + '\n' + repr((concVals*1e6).tolist()) + '\n')
if(calcTemperature & calcPressure):
    print('Computed Pressure:\t' + repr(res['x'][-2]/1000) + ' mb\n')
    print('Computed Temperature:\t' + repr(res['x'][-1]/1000) + ' K')
elif(calcPressure):
    print('Computed Pressure:\t' + repr(res['x'][-1]/1000) + ' mb')
elif(calcTemperature):
    print('Computed Temperature:\t' + repr(res['x'][-1]/1000) + ' K')
print('\n-------------------------\n')

##Plot using matplotlib
def plotAllSpecies():
    speciesGraphs = []
    for i in range(0,len(gasSpecies)):
        if(calcTemperature & calcPressure):
            xi, yi = computeAbsorption(res['x'][-2], res['x'][-1], [gasSpecies[i]], [res['x'][i]],[lD[i]],[mD[i]])
        elif(calcTemperature):
            xi, yi = computeAbsorption(pressure, res['x'][-1], [gasSpecies[i]], [res['x'][i]],[lD[i]],[mD[i]])
        elif(calcPressure):
            xi, yi = computeAbsorption(res['x'][-1], temperature, [gasSpecies[i]], [res['x'][i]],[lD[i]],[mD[i]])
        else:
            xi, yi = computeAbsorption(pressure, temperature, [gasSpecies[i]], [res['x'][i]],[lD[i]],[mD[i]])
        plt.plot(1.0E7/xi,yi,label=gasSpecies[i])
        speciesGraphs.append({'x' : xi, 'y' : yi, 'species' : gasSpecies[i], 'interpFunc' : interp1d(xi,yi)})
    return speciesGraphs

plt.plot(data['x'],data['alpha'],'bo',label = 'data')
speciesGraphs = plotAllSpecies()
plt.plot(1.0E7/x,y,label='total absorbance')
plt.axis([min(data['x']), max(data['x']), 0.0, max(data['alpha'])])
plt.xlabel("Wavelength(nm)")
plt.ylabel("Absorption Coef (cm^-1)")
plt.legend()

###File Outputs
fnOut = str.split(datFile,'.')
filenameOut = fnOut[0] + '_final'

filenameOutData = filenameOut + '_data.csv'
file = open(filenameOutData,'w')
file.write('#Data set absorption computed from the spectrum in ' + datFile + ' and the vacuum data set in' + vacFile + '\n')

file.write('\"wavelength(nm)\"\t\"data(optical depth)\"\n')
for i in range(0,len(data['alpha'])):
    file.write(str(data['x'][i]) + '\t' + str(data['alpha'][i]) + '\n')
file.close()

print("-----------------------------------\nMeasured absorption coefficients written to " + filenameOutData +"\n")

filenameOutHitran = filenameOut + '_hitranFit.csv'
file = open(filenameOutHitran,'w')
file.write('#Hitran data set fit for'+str(len(gasSpecies)) +' molecules\n')
if(calcPressure & calcTemperature):
    file.write('#Computed Pressure: \t'+str(res['x'][-2]) +' mb\n')
    file.write('#Computed Temperature:\t'+str(res['x'][-1]) +' K\n')
elif(calcPressure):
    file.write("#Temperature Set at:\t" + str(temperature) + '\n')
    file.write('#Computed Pressure: \t'+str(res['x'][-1]) +' mb\n')
elif(calcTemperature):
    file.write("#Pressure Set at:\t" + st(pressure) + '\n')
    file.write('#Computed Temperature: \t'+str(res['x'][-1]) +' K\n')
else:
    file.write("#Pressure Set at:\t" + str(pressure) + '\n')
    file.write('#Temperature Set at: \t'+str(temperature) +'\n')

file.write('\"Wavelength(nm)\"\t\"Total Absorbance\"')
for i in range(0,len(concVals)):
    file.write('\t\"' + gasSpecies[i] + ':' +  str(concVals[i]*1e6) + ' ppmV\"')
file.write('\n')
for i in range(0,len(x)):
    file.write(str(1.0E7/x[i]) + '\t' + str(y[i]))
    for j in range(0,len(speciesGraphs)):
        try:
            file.write('\t' + str(speciesGraphs[j]['interpFunc']( 1.0E7/x[i] ) ) )
        except ValueError:
            file.write('\t' + '0.0' )
    file.write('\n')
file.close()

print("Fitted HITRAN curves written to " + filenameOutHitran + "\n-----------------------------------\n")

##Show Final Plot (blocking statement, therfore it is last)
plt.show()