import matplotlib
matplotlib.use('TKAgg');   #makes the plots work on some systems
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d
from scipy import stats
from scipy.optimize import minimize, minimize_scalar

from string import *
from types  import *
import sys
import os.path
import time
import random 

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

useErrorValue = False

if (len(sys.argv) == 1):
    ##Set Parameters for running this script as standalone
    print("Standalone mode.\n")
    
    datFile = '/media/win/Users/jaymz/Desktop/new_data/reserviorSparklingWater-7-29-15_medians_absorptionValues.csv'
    seperatorCharactor = '\t'
    nmWavelengthUnits = True
    
    calcPressure = True
    calcTemperature = False
    findOffset = False
    
    resampleData = False
    numResamples = 50
    
    fitWavelength = False
    maxWavelengthOffset = 0.1
    
    lineFiles=['test/testH2O','test/CO2-C12','test/CO2-C13']
    gasSpecies = ['H2O', 'CO2', 'CO2-13']
    concGuess = [.007, 330.0E-5, 330.0E-7]
    isotopes = ['0', '626', '636']
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
    rasampleData = False
    seperatorCharactor = '\t'
    
    resampleData = False
    
    wnCalculationExtension = 30
    
    concGuess = None
    pressure, temperature = None, None
    isotopeInput = False
    
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
        if(str.lower(sys.argv[i]) == '--resample' and i < len(sys.argv) - 1):
            resampleData = True
            numResamples = int( sys.argv[i+1] )
        if(str.lower(sys.argv[i]) == '--isotopes'):
            splitStr = str.split(sys.argv[i+1],',')
            isotopes = [val for val in splitStr]
            isotopeInput = True
            
        if(str.lower( sys.argv[i] ) == '--help' or str.lower( sys.argv[i] ) == '-h'):
            print('\n----------------------------------\nScrript for the fitting of spectral data to the hitran line database\n\n')
            print('REQUIRED INPUTS:\n\n')
            print('-i , --input\t\t[input data file]\t\tInput data file')
            print('--line-files\t\t[input line files]\t\tInput line files generated from extract.py. Comma separated list. Ex: folder1/file1,folder2/file2')
            print('--species-labels\t[labels]\t\t\tComma separated list of gas names. Two cannot be identical. Ex: H2O,CO2')
            print('--conc-guess\t\t[Guess Values]\t\t\tComma separated list of gas concentration guesses. For a guess of 100 ppmV for two gasses: 1.0E-4,1.0E-4')
            print('-p , --pressure\t\t[pressure value]\t\tPressure value. If the --find-pressure option is enabled, this value is used as an initial guess in the fitting.')
            print('-t , --temperature\t[temperature value]\t\tTemperature value. If the --find-temperature option is enabled, this value is used as an initial guess in the fitting.')
            
            print('\n\nOptional Inputs:\n')
            print('--seperator\t\t[seperator string]\t\tCharactor which is used to seperate data collumns in input data file. (Default is set to tab charactor)')
            print('--wavenumbers\t\t\t\t\t\tUse wavenumber units (cm^-1) for input data file')
            print('--no-offset\t\t\t\t\t\tDo not fit a offset to aborption values (in addition to gas species)')
            print('--find-pressure\t\t\t\t\t\tAttempt finding total pressure of test gas by fitting spectrum to data.')
            print('--find-temperature\t\t\t\t\tAttempt finding the temperature of the test gas.')
            print('--resample\t\t[number of resamples]\t\tUse monte-carlo resampling of data to produce an error bar for all fitted values')
            print('--isotopes\t\t[isotope list]\t\t\tComma separated list of isotope identifiers. Used for mass and abundance parameters. Will still use all lines present in line file, however.')
            print('--help\t\t\t\t\t\t\tDisplay this help screen')
            
            print('\n--------------------------------------------------\n\n')
            sys.exit()
            
    if(not calcPressure and pressure == None or not calcTemperature and temperature == None):
        print("\nArguement error: temperature or pressure not defined. Set these or allow as open parameters\n")
        sys.exit()
    if(concGuess == None):
        concGuess = [0.001]*len(gasSpecies)
    if(not isotopeInput):
        abundances = ['0']*len(gasSpecies)

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

if(len(dat[0]) > 2):
    useErrorValue = True
dat = np.matrix(dat)

if(nmWavelengthUnits):
    data = {'x' : np.array((dat[:,0].flatten().tolist())[0]),
            'x(wn)' : np.array(((float(10**7)/dat[:,0]).flatten().tolist())[0]),
            'y' : np.array((dat[:,1].flatten().tolist())[0]) }
else:
    data = {'x' : np.array(((float(10**7)/dat[:,0]).flatten().tolist())[0]),
            'x(wn)' : np.array((dat[:,0].flatten().tolist())[0]),
            'y' : np.array((dat[:,1].flatten().tolist())[0]) }

if(useErrorValue):
    data['error'] = np.array( (dat[:,2].flatten().tolist())[0] )
            
file.close()

xLimits=Interval(min(data['x(wn)'])-wnCalculationExtension, max(data['x(wn)'])+wnCalculationExtension)


##load additional isotope information from text file
moleculeFile = open('molparam.txt')

moles = {}
mol = None
molData = {}

for line in moleculeFile:
    split = line.split()
    if( len(split) == 2):
        if(len(molData)>0):
            moles[mol] = molData
        mol = split[0]
        molData = {}
    if(len(split) == 5):
        molData[ split[0] ] =   { 'abundance' :  float(split[1]),
                                'Q(296K)' :    float(split[2]),
                                'gj' :         float(split[3]),
                                'mass' :       float(split[4])     }

moles[mol] = molData

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
index =0
for i in range(0,len(lineFiles)):
    # read a couple of lines and convert to actual pressure and temperature
    Line_Data, Mol_Data = get_lbl_data (lineFiles[i], xLimits, airWidth, wingExt, commentChar=commentChar)
    if(not isotopes[i] == '0'):
        if(moles[Line_Data['molecule']].has_key(isotopes[i])):
            Mol_Data['mass']  = moles[Line_Data['molecule']][isotopes[i]]['mass']
            Line_Data['strength']=Line_Data['strength']/moles[Line_Data['molecule']][isotopes[i]]['abundance']
        else:
            print("Invalid isotope identifier......\n")
            sys.exit()
    lD.append(Line_Data)
    mD.append(Mol_Data)
    index = index+1

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
    
def computeFitError(pressure, temperature, gasSpecies, concValues, Line_Data, Mol_Data, dataX, dataY, errorY):
    vGrid, absCoeff = computeAbsorption(pressure, temperature, gasSpecies, concValues, Line_Data, Mol_Data)
    dataX, dataY = np.array(dataX), np.array(dataY)
    interp = interp1d(vGrid, absCoeff)
    diff = -interp(dataX+wavelengthOffset) + dataY
    if(findOffset):
        diff = diff - sum(diff)/len(diff)
    if(errorY != None):
        diff = diff/np.array(errorY)
    return diff.dot(diff)
    
def funcToMinimize(x,Line_Data,Mol_Data, dataX, dataY, errorY):
    if(calcPressure & calcTemperature):
        out = computeFitError(x[-2], x[-1], gasSpecies, x[0:len(x)-2],Line_Data, Mol_Data, dataX, dataY, errorY)
    elif(calcTemperature):
        out = computeFitError(pressure, x[-1], gasSpecies, x[0:len(x)-1],Line_Data, Mol_Data, dataX, dataY, errorY)
    elif(calcPressure):
        out = computeFitError(x[-1], temperature, gasSpecies, x[0:len(x)-1],Line_Data, Mol_Data, dataX, dataY, errorY)
    else:
        out = computeFitError(pressure, temperature, gasSpecies, x,Line_Data, Mol_Data, dataX, dataY, errorY)
    if(any(n < 0 for n in x)):   #used to avoid negative values in the fit by giving them an artificially large error.
        return 100*out
    return out
    
##Main minimization


if(not resampleData):
    if(useErrorValue):
        res = minimize(funcToMinimize,concGuess,method='Nelder-Mead',args=(lD,mD, data['x(wn)'], data['y'], data['error']))
    else:
        res = minimize(funcToMinimize,concGuess,method='Nelder-Mead',args=(lD,mD, data['x(wn)'], data['y'], None))
else:
    avgConc = [0.0]*len(concGuess)
    avgConcSqr = [0.0]*len(concGuess)
    fitSuccessCounts =0
    for i in range(0,numResamples):
        syntheticX, syntheticY, syntheticE = [], [], []
        for j in range(0,len(data['y'])):
            randIndex = random.randrange(0,len(data['y']))
            syntheticX.append( data['x(wn)'][randIndex] )
            syntheticY.append( data['y'][randIndex] )
            if(useErrorValue):
                syntheticE.append( data['error'][randIndex] )
        if(useErrorValue):
            res = minimize(funcToMinimize,concGuess,method='Nelder-Mead',args=(lD,mD, syntheticX, syntheticY, syntheticE))
        else:
            res = minimize(funcToMinimize,concGuess,method='Nelder-Mead',args=(lD,mD, syntheticX, syntheticY, None))
        if(res['success']):
            print('-----------------------------------------------------\nSuccessful fit found at with values: ' + repr(res['x']))
            concGuess = res['x']
            avgConc = avgConc + res['x']
            avgConcSqr = avgConcSqr + res['x']**2
            fitSuccessCounts = fitSuccessCounts + 1
        else:
            print('Bad fit....\n')
    avgConc, avgConcSqr = avgConc/fitSuccessCounts, avgConcSqr/fitSuccessCounts
    concError = np.sqrt(avgConcSqr - avgConc**2)

concVals = []

def computeFinalCurves(res):
    global concVals
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
        
    offsetResult = minimize_scalar(offsetFunctionResidue,bounds = [ -maxWavelengthOffset, maxWavelengthOffset], method = 'Bounded', tol = 1.0e-25)
    
    wavelengthOffset = wavelengthOffset + offsetResult['x']
    if(useErrorValue):
        res2 = minimize(funcToMinimize,res['x'],method='Nelder-Mead',args=(lD,mD, data['x(wn)'], data['y'], data['error']) )
    else:
        res2 = minimize(funcToMinimize,res['x'],method='Nelder-Mead',args=(lD,mD, data['x(wn)'], data['y'], None) )
    
    xp, yp = computeFinalCurves(res2)

    

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
        if(findOffset):
            yi = yi + offset
        plt.plot(1.0E7/xi,yi,label=gasSpecies[i])
        speciesGraphs.append({'x' : xi, 'y' : yi, 'species' : gasSpecies[i], 'interpFunc' : interp1d(xi,yi)})
    return speciesGraphs
if(useErrorValue):
    plt.errorbar(data['x'],data['y'], yerr = data['error'], fmt = 'bo',label = 'data')
else:
    plt.plot(data['x'],data['y'],'bo',label = 'data')
if(fitWavelength):
    speciesGraphs = plotAllSpecies(res2)
else:
    speciesGraphs = plotAllSpecies(res)
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
    file.write("#Pressure Set at:\t" + str(pressure) + '\n')
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

if(resampleData):
    print("Resampling results:")
    if(calcTemperature and calcPressure):
        print("Average concentration result: " + str(avgConc[0:-2]) + "  +-  " + str(concError[0:-2]) )
        print("Average pressure result:" + repr(avgConc[-2]) + "  +-  " + repr(concError[-2]) )
        print("Average temperature result:" + repr(avgConc[-1]) + "  +-  " repr(concError[-1]) )
    elif(calcTemperature or calcPressure):
        print("Average concentration result: " + str(avgConc[0:-1]) + "  +-  " + str(concError[0:-1]) )
        if(calcPressure):
            print("Average pressure result: " + repr(avgConc[-1]) + " +- " + repr(concError[-1]) )
        elif(calcTemperature):
            print("Average temperature result: " + repr(avgConc[-1]) + " +- " + repr(concError[-1]) )
    else:
        print("Average concentration result: " + str(avgConc) + "  +-  " + str(concError) )
    print("-------------------------------------------------------------------------------

##Show Final Plot (blocking statement, therfore it is last)
plt.show()