import sys
from molecules import molecules
from subprocess import call
import numpy as np

databaseDirectory  = '../par/'
databaseNameString = '_hit12.par'
tempDirectory = '../tempFiles/'

dataSeperatorCharactor = '\t'
nmWavelengthUnits = True

temperature, pressure = None, None  #Specify default value here. Will prompt user if = None.

names = []
isos = []
isoNames = []
concGuess = []

lineExtension = 100.0

minX, maxX = None, None

###-----------------Scan command line arguments--------------------------------------
for i in range(0,len(sys.argv)):
    if(not isinstance(sys.argv[i],str)):
        continue
    if( ( str.lower(sys.argv[i]) == '-m' or str.lower(sys.argv[i]) == '--molecules' ) and i < len(sys.argv) - 1):
        inSplit = str.split(sys.argv[i+1], ',')
        for mol in inSplit:
            if(molecules.has_key(mol)):
                names.append(mol)
            else:
                print("Error! Invalid molecule name.\n")
                sys.exit()
                
    if( ( str.lower(sys.argv[i]) == '-is' or str.lower(sys.argv[i]) == '--isotopologues' ) and i < len(sys.argv) - 1):
        inSplit = str.split(sys.argv[i+1], ',')
        if(len(inSplit) == len(names)):
            for i in range(0,len(inSplit)):
                if(inSplit[i] == '0' or molecules[names[i]].has_key(inSplit[i])):
                    isos.append(inSplit[i])
                    isoNames.append(names[i] + '-' + inSplit[i])
                else:
                    print("Invalid isotopologue number!!")
                    sys.exit()
        else:
            print("Molecule names must be specified first in command line arguements!!")
            sys.exit()
    if( ( str.lower(sys.argv[i]) == '-c' or str.lower(sys.argv[i]) == '--concentrations' ) and i < len(sys.argv) - 1):
        inSplit = str.split(sys.argv[i+1],',')
        for conc in inSplit:
            concGuess.append(float(conc))
    if( (str.lower(sys.argv[i]) == '-t' or str.lower(sys.argv[i]) == '--temperature' ) and i < len(sys.argv) - 1):
        temperature = float(sys.argv[i+1])
    if( (str.lower(sys.argv[i]) == '-p' or str.lower(sys.argv[i]) == '--pressure' ) and i < len(sys.argv) - 1):
        pressure = float(sys.argv[i+1])
    if( ( str.lower(sys.argv[i]) == '-x' or str.lower(sys.argv[i]) == '--x-range' ) and i < len(sys.argv) - 1):
        xrangeIn = sys.argv[i+1].split(',')
        if(len(xrangeIn) ==2):
            minX = float(xrangeIn[0])
            maxX = float(xrangeIn[1])
        else:
            print("Invalid x range")
            sys.exit(1)
        
        
    if(str.lower(sys.argv[i]) == '--help'):
        print('\n-----------------------------------\n\nWrapper script for easy computation of HITRAN absorption simulation using generateHitranCurve.py script\n')
        print('Optional parameters (if not specified, they will be prompted):\n')
        print('-m,--molecules\t\t[Molecules]\t\tComma separated list of molecule names (H2O,CO2,...)')
        print('-is,--isotopologues\t[Isotopologues]\t\tComma separated list of Isotopologue codes (Ex. for CO2: 626,636). Zero for all species')
        print('-c,--concentrations\t[Concentrations]\tComma separated list of concentration values (1E-6,2E-5 for 1 ppmV and 10ppmV respectively)')
        print('-t,--temperature\t[Temperature]\t\tTemperature value or initial guess')
        print('-p,--pressure\t\t[Pressure]\t\tPressure value or initial guess')
        print('-x,--x-range\t\t[X Range]\t\tLower bound, upper x bound')
        print('\n------------------------------------------------------\n\n')
        sys.exit()


if(len(names) == 0 and len(isos) ==0 ):  #prompt for molecule name and isotopologue numbers
    doneInputting = False
    while not doneInputting:
        name = raw_input("Enter a molecule name, ex. H2O (or \"done\" when finished entering all molecules):")
        if(str.lower(name) != 'done'):
            if(name not in molecules):
                print("Invalid Molecule Name!\n")
            else:
                names.append(name)
                haveValidIsotope = False
                while not haveValidIsotope:
                    print("Isotopes of " + name + ': ' + str(molecules[name]['isotopes']) )
                    iso = raw_input('Enter an isotopologue number (not specifying one will fit all isotopologues): ')
                    if(iso != '' and iso != '0'):
                        if( iso not in molecules[name]['isotopes'] ):
                            print("Invalid Isotope number!\n")
                        else:
                            isos.append(iso)
                            isoNames.append(name + '-' + iso)
                            haveValidIsotope = True
                    else:
                        print("Using all isotopologues of " + name +". Using built-in hitran abundance scaling")
                        isoNames.append(name)
                        isos.append('0')
                        haveValidIsotope = True
                print('\n')
        else:
            doneInputting = True
if(minX == None or maxX == None):
    minX = raw_input("Enter minimum nm range: ")
    maxX = raw_input("Enter maximum nm range: ")

###--------------Make sure the correct HITRAN database files are read------------------
def get_hitran_name(moleculeName):
    molNumber = molecules[moleculeName]['hitran']
    molStr = str(molNumber)
    if(molNumber < 10):
        molStr = '0' + molStr
    return databaseDirectory + molStr + databaseNameString

inputHitranFiles = [get_hitran_name(mol) for mol in names]


##execute extract command for extracting info from HITRAN databases -- extract.py
print('\n\n---------------------------------------------\nExtracting lines from database files....\n')
for i in range(0,len(names)):
    extractCommand = [  'python',
                        'extract.py',
                        '-o',
                        tempDirectory+isoNames[i],
                        '-m',
                        str(names[i]),
                        '-x',
                        str(minX) + ',' + str(maxX),
                        '-X',
                        'nm',
                        inputHitranFiles[i]]
    call(extractCommand)
    
##Get values for concentrations, temperature, pressure
if(len(concGuess) != len(names)):
    concGuess = []
    defaultVal = 0.002
    print("\n\n--------------------------------------\n")
    for name in isoNames:
        try:
            val = float(raw_input("Enter concentration (1E-6 == ppmV) for " + name  + ": "))
        except ValueError:
            print("Using default of " + str(defaultVal))
            val = defaultVal
        concGuess.append(val)
pressDefault, tempDefault = 1013.2, 296.0
while pressure == None:
    
    try:
        pressure = raw_input("Enter pressure value (mb):")
    except ValueError:
        print("Must have pressure value!!!\n")
        continue

while temperature == None:
    try:
        temperature = raw_input("Enter temerature value (K):")
    except ValueError:
        print("Must have temperature value!!!\n")
        continue

pressure = float(pressure) * 1000
temperature = float(temperature)

##Run curve fitting script -- optimizeConc.py
print('\n\n---------------------------------------------\nGenerating plots of abosrption wavelength dependence....\n')

def argString(listIn):
    outString = ''
    for thing in listIn:
        outString = outString + str(thing) + ','
    return outString[0:-1]    
    

fitCommand = ['python',
              'generateHitranCurve.py',
              '--line-files',
              argString([(tempDirectory + nam) for nam in isoNames]),
              '--species-labels',
              argString(isoNames),
              '--conc',
              argString(concGuess),
              '-p',
              str(pressure),
              '-t',
              str(temperature),
              '--x-range',
              str(minX) + ',' + str(maxX),
              '--isotopes',
              argString(isos) ]    
call(fitCommand)
