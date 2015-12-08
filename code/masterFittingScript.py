import sys
from molecules import molecules
from subprocess import call
import numpy as np

databaseDirectory  = '../par/'
databaseNameString = '_hit12.par'
tempDirectory = 'tempFiles/'

dataSeperatorCharactor = '\t'
nmWavelengthUnits = True

calcPressure, calcTemperature = False, False
calcOffset = True

temperature, pressure = None, None  #Specify default value here. Will prompt user if = None.

names = []
isos = []
isoNames = []

###-----------------Scan command line arguments--------------------------------------
for i in range(0,len(sys.argv)):
    if(str.lower(sys.argv[i]) == '-t' or str.lower(sys.argv[i]) == '--find-temperature'):
        calcTemperature = True
        
    if(str.lower(sys.argv[i]) == '-p' or str.lower(sys.argv[i]) == '--find-pressure'):
        calcPressure = True
        
    if( ( str.lower(sys.argv[i]) == '-m' or str.lower(sys.argv[i]) == '--molecules' ) and i < len(sys.argv) - 1):
        inSplit = split(sys.argv[i+1], ',')
        for mol in inSplit:
            if(molecules.has_key(mol)):
                names.add(mol)
            else:
                print("Error! Invalid molecule name.\n")
                sys.exit()
                
    if( ( str.lower(sys.argv[i]) == '-is' or str.lower(sys.argv[i]) == '--isotopologues' ) and i < len(sys.argv) - 1):
        inSplit = split(sys.argv[i+1], ',')
        if(len(inSplit) == len(names)):
            for i in range(0,len(inSplit)):
                if(molecules[names[i]].has_key(insplit[i])):
                    isos.append(inSplit[i])
                    isoNames.append(name[i] + '-' + inSplit[i])
                else:
                    print("Invalid isotopologue number!!")
                    sys.exit()
        else:
            print("Molecule names must be specified first in command line arguements!!")
            sys.exit()
    
    if( ( str.lower(sys.argv[i]) == '--temperature' ) and i < len(sys.argv) - 1):
        temperature = float(sys.argv[i+1])
    if( ( str.lower(sys.argv[i]) == '--pressure' ) and i < len(sys.argv) - 1):
        temperature = float(sys.argv[i+1])
        
    if(str.lower(sys.argv[i]) == '--help'):
        print('\n-----------------------------------\n\nWrapper script for easy fitting of spectral data using the optimizeConc.py script\n')
        print('Optional parameters:\n')
        print('-t, --find-temperature\t\t\t\tFind the temperature of the test gas by fitting data.')
        print('-p, --find-pressure\t\t\t\tFind the pressure of the test gas by fitting data.')
        print('-i,--input\t\t[Data File name]\tFilename of tab separated list of wavelegnths, absorpbance coefficients (and optionally error)')
        print('-m,--molecules\t\t[Molecules]\t\tComma separated list of molecule names (H2O,CO2,...)')
        print('-is,--isotopologues\t[Isotopologues]\t\tComma separated list of Isotopologue codes (Ex. for CO2: 626,636). Zero for all species')
        print('--temperature\t\t[Temperature]\t\tTemperature value or initial guess')
        print('--pressure\t\t[Pressure]\t\tPressure value or initial guess (mb)')
        print('--no-offset\t\t\t\t\tDo not fit offset to absorption values. In other words, do not allow for wavelength independent absoprtion in the fit')
        print('\n------------------------------------------------------\n\n')
        sys.exit()

def tk_file_open_dialog():
    import Tkinter as tk
    import tkFileDialog
    
    root = tk.Tk()
    root.withdraw()
    fileName = tkFileDialog.askopenfilename()
    
    if(fileName == ()):
        sys.exit()
    
    return fileName

def get_path(wildcard):
    import wx
    app = wx.App(None)
    style = wx.FD_OPEN | wx.FD_FILE_MUST_EXIST
    dialog = wx.FileDialog(None, 'Open', wildcard = wildcard, style = style)
    if(dialog.ShowModal() == wx.ID_OK):
        path = dialog.GetPath()
    else:
        path = None
    dialog.Destroy()
    return path
    
print("Choose data file to open\n")

if(dataFileName == None):
    dataFileName = tk_file_open_dialog()        #File dialog that works on the windows Anaconda package by default
#dataFileName = get_path("*")               #File dialog that might work on other systems
#dataFileName = '/media/win/Users/jaymz/Desktop/new_data/reserviorSparklingWater-7-29-15_dataDump_final_data.csv'        #static path
if(dataFileName == None):           #exit if no file chosen
    sys.exit()

print("Using input data file: " + dataFileName)

print("--------------------------------------------\n")
if(len(names) == 0 and len(isos) ==0 ):  #prompt for molecule name and isotopologue numbers
    doneInputting = False
    while not doneInputting:
        name = raw_input("Enter a molecule name to fit, ex. H2O (or \"done\" when finished entering all molecules):")
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

###--------------Make sure the correct HITRAN database files are read------------------
def get_hitran_name(moleculeName):
    molNumber = molecules[moleculeName]['hitran']
    molStr = str(molNumber)
    if(molNumber < 10):
        molStr = '0' + molStr
    return databaseDirectory + molStr + databaseNameString

inputHitranFiles = [get_hitran_name(mol) for mol in names]

##-------------input data-----------
##Note this is redundant since it is done again in a subprocess.
##It is needed here to determine wavelength range of hitran computation.

dat = []
headers = []

file = open(dataFileName,'r')

for line in file:
    l = str.split(line,dataSeperatorCharactor)
    ln = []
    append = True
    for thing in l:
        try:
            ln.append(float(thing))
        except ValueError:
            headers.append(thing)
            append = False
    if(append):
        dat.append(ln)
    
dat = np.matrix(dat)

if(nmWavelengthUnits):
    data = {'x' :       np.array((dat[:,0].flatten().tolist())[0]),
            'x(wn)' :   np.array(((float(10**7)/dat[:,0]).flatten().tolist())[0]),
            'y' :       np.array((dat[:,1].flatten().tolist())[0]) }
else:
    data = {'x' :       np.array(((float(10**7)/dat[:,0]).flatten().tolist())[0]),
            'x(wn)' :   np.array((dat[:,0].flatten().tolist())[0]),
            'y' :       np.array((dat[:,1].flatten().tolist())[0]) }

minX, maxX = min(data['x']) - 75, max(data['x']) + 75
file.close()

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
    
##Get guesses for concentrations, temperature, pressure
concGuess = []
defaultVal = 0.002
print("\n\n--------------------------------------\n")
for name in isoNames:
    try:
        val = float(raw_input("Enter concentration guess (1E-6 == ppmV) for " + name  + ": "))
    except ValueError:
        print("Using default guess of " + str(defaultVal))
        val = defaultVal
    concGuess.append(val)
pressDefault, tempDefault = 1013.2, 296.0
while pressure == None:
    if(calcPressure):
        try:
            pressure = raw_input("Enter pressure guess (mb):")
        except ValueError:
            print("Using default pressure guess of:"+str(pressDefault) +"\n")
            pressure = pressDefault
            continue
    else:
        try:
            pressure = raw_input("Enter pressure value (mb):")
        except ValueError:
            print("Must have pressure value!!!\n")
            continue

while temperature == None:
    if(calcTemperature):
        try:
            temperature = raw_input("Enter temerature guess (K):")
        except ValueError:
            print("Using default temperature guess of"+tempDefault+"\n")
            temperature = tempDefault
            continue
    else:
        try:
            temperature = raw_input("Enter temerature value (K):")
        except ValueError:
            print("Must have temperature value!!!\n")
            continue

pressure = float(pressure) * 1000
temperature = float(temperature)

##Run curve fitting script -- optimizeConc.py
print('\n\n---------------------------------------------\nFitting data to line parameters....\n')

def argString(listIn):
    outString = ''
    for thing in listIn:
        outString = outString + str(thing) + ','
    return outString[0:-1]    
    

fitCommand = ['python',
              'optimizeConc.py',
              '-i',
              dataFileName,
              '--line-files',
              argString([(tempDirectory + nam) for nam in isoNames]),
              '--species-labels',
              argString(isoNames),
              '--conc-guess',
              argString(concGuess),
              '-p',
              str(pressure),
              '-t',
              str(temperature),
              '--isotopes',
              argString(isos) ]
if(calcPressure):
    fitCommand.append('--find-pressure')
if(calcTemperature):
    fitCommand.append('--find-temperature')
if(not calcOffset):
    fitCommand.append('--no-offset')
    
call(fitCommand)
