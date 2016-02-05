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


#Test Arguments:
arguments = ['--line-files',
 '../tempFiles/H2O-161,../tempFiles/CO2-626',
 '--species-labels',
 'H2O-161,CO2-626',
 '--conc',
 '0.007,0.0003',
 '-p',
 '1013500.0',
 '-t',
 '296.0',
 '--x-range',
 '1598,1605',
 '--isotopes',
 '161,626']
 
 
sys.argv.extend(arguments)


useErrorValue = False

if (len(sys.argv) == 1):
    ##Set Parameters for running this script as standalone
    print("Standalone mode.\n")
    
    nmWavelengthUnits = True
    
    resampleData = False
    numResamples = 50
    
    maxWavelengthOffset = 0.1
    
    lineFiles=['../test/testH2O','../test/CO2-C12','../test/CO2-C13']
    gasSpecies = ['H2O', 'CO2', 'CO2-13']
    concGuess = [.007, 330.0E-5, 330.0E-7]
    isotopes = ['0', '626', '636']
    temperature = 295.0
    pressure = 1.0149E6*360/760
    
    wnCalculationExtension = 20
    xLimits=Interval(6000.0-wnCalculationExtension, 8000.0+wnCalculationExtension)
else:
    ##Scan command line parameters when executed from command line
    nmWavelengthUnits = True
    
    seperatorCharactor = '\t'
    
    resampleData = False
    
    wnCalculationExtension = 30
    
    concGuess = None
    pressure, temperature = None, None
    isotopeInput = False
    
    for i in range(0,len(sys.argv)):
        if(str.lower(sys.argv[i]) == '--wavenumbers'):
            nmWavelengthUnits = False
        if(str.lower(sys.argv[i]) == '--line-files' and i < len(sys.argv) - 1):
            lineFiles = str.split(sys.argv[i+1],',')
        if(str.lower(sys.argv[i]) == '--species-labels' and i < len(sys.argv) - 1):
            gasSpecies = str.split(sys.argv[i+1],',')
        if(str.lower(sys.argv[i]) == '--conc' and i < len(sys.argv) - 1):
            concGuessStr = str.split(sys.argv[i+1],',')
            concGuess = [float(val) for val in concGuessStr]
        if( ( str.lower(sys.argv[i]) == '-p' or str.lower(sys.argv[i]) == '--pressure' ) and i < len(sys.argv) - 1):
            pressure = float(sys.argv[i+1])
        if( ( str.lower(sys.argv[i]) == '-t' or str.lower(sys.argv[i]) == '--temperature' ) and i < len(sys.argv) - 1):
            temperature = float(sys.argv[i+1])
        if(str.lower(sys.argv[i]) == '--resample' and i < len(sys.argv) - 1):
            resampleData = True
            numResamples = int( sys.argv[i+1] )
        if( ( str.lower(sys.argv[i]) == '-x' or str.lower(sys.argv[i]) == '--x-range' ) and i < len(sys.argv) - 1):
            xrangeIn = sys.argv[i+1].split(',')
            print("X range:")
            print(xrangeIn)
            if(len(xrangeIn) ==2):
                xLimits = Interval(float(xrangeIn[0])-wnCalculationExtension, float(xrangeIn[1])+wnCalculationExtension)
            else:
                print("Invalid x range")
                sys.exit(1)
        if(str.lower(sys.argv[i]) == '--isotopes'):
            splitStr = str.split(sys.argv[i+1],',')
            isotopes = [val for val in splitStr]
            isotopeInput = True
            
        if(str.lower( sys.argv[i] ) == '--help' or str.lower( sys.argv[i] ) == '-h'):
            print('\n----------------------------------\nScrript for the fitting of spectral data to the hitran line database\n\n')
            print('REQUIRED INPUTS:\n\n')
            print('--line-files\t\t[input line files]\t\tInput line files generated from extract.py. Comma separated list. Ex: folder1/file1,folder2/file2')
            print('--species-labels\t[labels]\t\t\tComma separated list of gas names. Two cannot be identical. Ex: H2O,CO2')
            print('--conc\t\t[Concentration Values]\t\t\tComma separated list of gas concentration guesses. For a guess of 100 ppmV for two gasses: 1.0E-4,1.0E-4')
            print('-p , --pressure\t\t[pressure value]\t\tPressure value.')
            print('-t , --temperature\t[temperature value]\t\tTemperature value.')
            print('\n\nOptional Inputs:\n')
            print('--wavenumbers\t\t\t\t\t\tUse wavenumber units (cm^-1) for input data file')
            print('--isotopes\t\t[isotope list]\t\t\tComma separated list of isotope identifiers. Used for mass and abundance parameters. Will still use all lines present in line file, however.')
            print('-x,--x-range\t\t[X Range]\t\tLower bound, upper x bound')
            print('--help\t\t\t\t\t\t\tDisplay this help screen')
            
            print('\n--------------------------------------------------\n\n')
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
    

vGrid, absCoeff = computeAbsorption(pressure, temperature, gasSpecies, concGuess, lD, mD)

##Plot using matplotlib

plt.plot(vGrid,absCoeff,label = "HITRAN")

plt.xlabel("Wavelength(nm)")
plt.ylabel("Absorption Coef (cm^-1)")
plt.legend()

##Show Final Plot (blocking statement, therfore it is last)
plt.show()