import matplotlib
matplotlib.use('TKAgg');
import matplotlib.pyplot as plt

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
from atmos import Atmos1D, get_profile_list, reGridAtmos1d
from cgsUnits import unitConversion
from molecules import molecules
molecNames   = molecules.keys()

atmFile = '../data/atmosDataFile'
lineFiles=['test/testH2O','test/testCO2']
outFile=None
commentChar='#'
zToA=None
lineShape='Voigt'
xLimits=Interval(6210.0, 6255.0)
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

# parse file header and retrieve molecule names
species = [grep_from_header(file,'molecule') for file in lineFiles]			

print 'lbl2od for',  len(species), ' molecules:', join(species)

# read atmospheric data from file,  returns an instance
atmos = Atmos1D(get_profile_list (atmFile, commentChar, ['z', 'p', 'T']+species, verbose))
# optinally remove high altitudes 
if zToA: atmos=reGridAtmos1d(atmos,zToA)
print '\n', atmFile, '   ==> atmos\n', atmos
# vertical column densities
vcd = atmos.vcd()
print str('   vcd [molec/cm**2]' + 20*' '+len(atmos.gases)*'%10.4g') % tuple(vcd)

# interpolation method used in multigrid algorithms: only Lagrange!
if   interpolate in 'bBsScC':  lagrange=4
elif interpolate in 'qQ':      lagrange=3
elif interpolate in 'lL':      lagrange=2
else:                          lagrange=int(interpolate)

# loop over molecules:  read its line parameters (and some other data like mass) and compute cross sections for all p,T
xsMatrix = []
for file in lineFiles:
    # read a couple of lines and convert to actual pressure and temperature
    Line_Data, Mol_Data = get_lbl_data (file, xLimits, airWidth, wingExt, commentChar=commentChar)
    # compute cross sections for selected line shape, returns a list with npT dictionaries
    molCrossSections = lbl_xs (Line_Data, Mol_Data, atmos.p, atmos.T, xLimits, lineShape, sampling, wingExt, nGrids, gridRatio, nWidths, lagrange, verbose)
    xsMatrix.append(molCrossSections)
# for xss in xsMatrix: print 'xs', type(xss), len(xss), [(xs.get('molecule'),len(xs.get('y'))) for xs in xss]

# compute absorption coefficients, also return uniform, equiidistant wavenumber grid corresponding to most dense cross section
vGrid, absCoeff = sum_xsTimesDensity (xsMatrix, atmos.densities, interpolate=interpolate, verbose=verbose)

for i in range(0,len(absCoeff[0])):
    plt.plot(vGrid,absCoeff[:,i])
plt.show()
