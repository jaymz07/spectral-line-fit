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

if (len(sys.argv) == 1):
    ##Set Parameters for running this script as standalone
    print("Standalone mode.\n")
    
    datFile = '/media/win/Users/jaymz/Desktop/new_data/reserviorSparklingWater_dataset2_-7-29-15_dataDump_final_data.csv'
    seperatorCharactor = '\t'
    nmWavelengthUnits = True
    
    calcPressure = True
    calcTemperature = False
    findOffset = True
    
    fitWavelengthTrue = True
    wavelengthMaxErr = 0.1
    
    lineFiles=['test/testH2O','test/testCO2']
    gasSpecies = ['H2O','CO2']
    concGuess = [.007, 330.0E-5]
    temperature = 295.0
    pressure = 1.0149E6*360/760
    
    wnCalculationExtension = 20
else:
    ##Scan command line parameters when executed from command line
    nmWavelengthUnits = True
    calcPressure = False
    calcTemperature = False
    fitWavelength = False
    findOffset = True
    seperatorCharactor = '\t'
    
    concGuess = None
    pressure, temperature = None, None
    
    for i in range(0,len(sys.argv)):
        if( ( str.lower(sys.argv[i]) == '-i' or str.lower(sys.argv[i]) == '--input' ) and i < len(sys.argv) - 1):
            datFile = sys.argv[i+1]
        if(str.lower(sys.argv[i]) == '--seperator' and i < len(sys.argv) - 1):
            speratorCharactor = sys.argv[i+1]
        if(str.lower(sys.argv[i]) == '--wavenumbers'):
            nmWavelengthUnits = False
        if(str.lower(sys.argv[i]) == '--find-pressure'):
            calcPressure = True
        if(str.lower(sys.argv[i]) == '--find-temperature'):
            calcTemperature = True
        if(str.lower(sys.argv[i]) == '--no-offset'):
            findOffset = False
        if(str.lower(sys.argv[i]) == '--line-files' and i < len(sys.argv) - 1):
            lineFiles = str.split(sys.argv[i+1],',')
        if(str.lower(sys.argv[i]) == '--species-labels' and i < len(sys.argv) - 1):
            gasSpecies = str.split(sys.argv[i+1],',')
        if(str.lower(sys.argv[i]) == '--conc-guess' and i < len(sys.argv) - 1):
            concGuessStr = str.split(sys.argv[i+1],',')
            concGuess = [float(val) for val in concGuessStr]
        if( ( str.lower(sys.argv[i]) == '-p' or str.lower(sys.argv[i]) == '--pressure' ) and i < len(sys.argv) - 1):
            pressure = float(sys.argv[i+1])
        if( ( str.lower(sys.argv[i]) == '-t' or str.lower(sys.argv[i]) == '--temperature' ) and i < len(sys.argv) - 1):
            temperature = float(sys.argv[i+1])
            
    if(not calcPressure and pressure == None or not calcTemperature and temperature == None):
        print("\nArguement error: temperature or pressure not defined. Set these or allow as open parameters\n")
        sys.exit()
    if(concGuess == None):
        concGuess = [0.001]*len(gasSpecies)

print(concGuess)
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
headers = []

file = open(datFile,'r')

for line in file:
    l = str.split(line,seperatorCharactor)
    ln = []
    append = False
    for thing in l:
        try:
            ln.append(float(thing))
            append = True
        except ValueError:
            headers.append(thing)
    if(append):
        dat.append(ln)
    
dat = np.matrix(dat)

if(nmWavelengthUnits):
    data = {'x' : np.array((dat[:,0].flatten().tolist())[0]),
            'x(wn)' : np.array(((float(10**7)/dat[:,0]).flatten().tolist())[0]),
            'y' : np.array((dat[:,1].flatten().tolist())[0]) }
else:
    data = {'x' : np.array(((float(10**7)/dat[:,0]).flatten().tolist())[0]),
            'x(wn)' : np.array((dat[:,0].flatten().tolist())[0]),
            'y' : np.array((dat[:,1].flatten().tolist())[0]) }
            
file.close()

xLimits=Interval(min(data['x(wn)'])-wnCalculationExtension, max(data['x(wn)'])+50)

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
#wavelengthOffset = -0.0011197771627370701 #0.0
wavelengthOffset = 0.0
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
    diff = -interp(data['x(wn)']+wavelengthOffset) + data['y']
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
    
##Main minimization
res = minimize(funcToMinimize,concGuess,method='Nelder-Mead',args=(lD,mD))

def computeFinalCurves(res):
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
        diff = data['y'] - interpFunc(data['x(wn)'])
        offset = sum(diff)/len(diff)
        y = y + offset
    
    
        print('------------------\nConcentrations (ppm):\n' + repr(gasSpecies) + '\n' + repr((concVals*1e6).tolist()) + '\n')
    if(calcTemperature & calcPressure):
        print('Computed Pressure:\t' + repr(res['x'][-2]/1000) + ' mb\n')
        print('Computed Temperature:\t' + repr(res['x'][-1]/1000) + ' K')
    elif(calcPressure):
        print('Computed Pressure:\t' + repr(res['x'][-1]/1000) + ' mb')
    elif(calcTemperature):
        print('Computed Temperature:\t' + repr(res['x'][-1]/1000) + ' K')
    print('\n-------------------------\n')
    
    return x,y
    
x,y = computeFinalCurves(res)
if(fitWavelength):
    plt.plot(1.0E7/x,y,label='first fit')

if(fitWavelength):
    print("---------------------------------------------\nAttempting wavelength correction\n---------------------------------------\n")
    edgeCutoff = 85
        
    def residue(dataX, dataY, refCurve):
        diff = refCurve(dataX) - np.array(dataY)
        return diff.dot(diff)
        
    def offsetFunctionResidue(offset):
        return residue(data['x'][::-1][edgeCutoff:-edgeCutoff] + offset,data['y'][::-1][edgeCutoff:-edgeCutoff], interp1d(1.0e7/x,y))
        
    offsetResult = minimize_scalar(offsetFunctionResidue,bounds = [ -0.1, 0.1], method = 'Bounded', tol = 1.0e-25)
    
    wavelengthOffset = wavelengthOffset + offsetResult['x']
    res2 = minimize(funcToMinimize,res['x'],method='Nelder-Mead',args=(lD,mD))
    
    xp, yp = computeFinalCurve(res2)

    

##Plot using matplotlib

def plotAllSpecies(res):
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
        yi = yi + offset
        plt.plot(1.0E7/xi,yi,label=gasSpecies[i])
        speciesGraphs.append({'x' : xi, 'y' : yi, 'species' : gasSpecies[i], 'interpFunc' : interp1d(xi,yi)})
    return speciesGraphs

plt.plot(data['x'],data['y'],'bo',label = 'data')
speciesGraphs = plotAllSpecies(res2)
if(not fitWavelength):
    plt.plot(1.0E7/x,y,label='total absorbance')
else:
    plt.plot(1.0E7/xp,yp,label='total absorbance (fit 2)')
plt.axis([min(data['x']), max(data['x']), 0.0, max(data['y'])])
plt.xlabel("Wavelength(nm)")
plt.ylabel("Absorption Coef (cm^-1)")
plt.legend()

###File Outputs
fnOut = str.split(datFile,'.')
filenameOut = fnOut[0] + '_final'

filenameOutHitran = filenameOut + '_hitranFit.txt'
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

print("-----------------------------------\nFitted HITRAN curves written to " + filenameOutHitran + "\n-----------------------------------\n")

##Show Final Plot (blocking statement, therfore it is last)
plt.show()