#!/usr/bin/env python

"""  lbl2od
  computation of line-by-line optical depth due to molecular absorption

  usage:
  lbl2od  [options]  atm_file  line_parameter_file(s)

  -h             help
  -c   char      comment character(s) used in input,output file (default '#')
  -o   string    output file for saving of cross sections (if not given: write to StdOut)
                 (in case of a list of input files this options specifies just the extension)

  -m   char      mode: 'c' ---> cumulative optical depth
                       'd' ---> difference (delta) optical depth (default)
                       'r' ---> reverse cumulative optical depth
		       anything else: total optical depth

  -i   char      interpolation method   [default: '3' three-point Lagrange,  choices are one of "234lqcsbk"]
  -L   char      Lineshape: V(oigt), L(orentz), G(auss)     [default: Voigt]
  -q   char      quadrature method: default 'T'=trapez,  'S'=simpson
  -s   float     sampling rate used for x-grid (default: 5.0 grid points per (mean) half width)
  -w   float     wing extension (cutoff wavenumber, default 10.0cm-1)

  -n   int       number of grids --- selects 'multigrid mode' for nGrids=2 or 3  (default: nGrids=3 three grid)
  -g   int       gridRatio = ratio of coarse to fine grid spacing  (only 2, 4, or 8, default 8)
  -W   float     transition from fine to coarse grid  (in units of half widths, default 25.0)

 --ToA float     top-of-atmosphere altitude (compute opt.depth only for levels below) 
  -x   Interval  lower,upper wavenumbers (comma separated pair of floats [no blanks!],
                                          default set according to range of lines in datafiles)
 --nm            on output write optical depth versus wavelength [nm] (default: wavenumber 1/cm)
  -r             on output reverse layer optical depth order:  top <--> bottom of atmosphere
  -y             on output write only optical depth without x-values (wavenumbers etc)


  line_parameter_file(s)
  The line parameter file(s) should contain a list of (preselected) lines
  that can be generated from HITRAN or GEISA database with extract.py 
  (See the documentation header of lbl2xs.py for more details)
  
  atm_file 
  The file containing an user's atmospheric profile data has to be in xy  format
  with columns for altitude, pressure, temperature and molecuar concentrations
  (see the documentation for details and the data directory for some examples)

  NOTES:
  interpolation required for multigrid approach and to regrid cross sections to the final dense grid (typically in the top levels)
  * Lagrange linear, quadratic, or cubic interpolation
  * B spline interpolation (uses scipy.interpolate modules splrep, splev)
  * Krogh (1970) interpolation (uses scipy.interpolate.krogh_interpolate)
  Lagrange two-point interpolation (selected by "2", "l", or "L") is clearly the least elaborate and most robust,
  even three-point interpolation can sometimes lead to oscillations associated with some negative xs values 
"""

####################################################################################################################################

# import some standard python modules
from string import *
from types  import *
import sys
import os.path
import time 

try:                      import numpy as np
except ImportError, msg:  raise SystemExit, str(msg) + '\nimport numeric python failed!'

from IO import open_outFile, commonExtension, paste, writeArray, grep_from_header #fgs120310 parse_file_header, readFileHeader
from pairTypes import Interval
from lbl2xs import lbl_xs
from lines import get_lbl_data
from xs2ac2od import sum_xsTimesDensity, integral_acTimesDistance
from lagrange_interpolation import *
from atmos import Atmos1D, get_profile_list, reGridAtmos1d
from cgsUnits import unitConversion
from molecules import molecules
molecNames   = molecules.keys()

####################################################################################################################################

def _lbl2od_ (atmFile, lineFiles, outFile=None, commentChar='#', zToA=None,
              lineShape='Voigt',
	      xLimits=None, airWidth=0.1, sampling=5.0, wingExt=5.0, interpolate='3', nGrids=2, gridRatio=0, nWidths=25.0,
	      mode='d', quad='T', nanometer=False, yOnly=False, flipUpDown=False, plot=False, verbose=False):

	# parse file header and retrieve molecule names
	species = [grep_from_header(file,'molecule') for file in lineFiles]			
	
	#print 'lbl2od for',  len(species), ' molecules:', join(species)

	# read atmospheric data from file,  returns an instance
	atmos = Atmos1D(get_profile_list (atmFile, commentChar, ['z', 'p', 'T']+species, verbose))
	# optinally remove high altitudes 
	if zToA: atmos=reGridAtmos1d(atmos,zToA)
	#print '\n', atmFile, '   ==> atmos\n', atmos
	# vertical column densities
	vcd = atmos.vcd()
	#print str('   vcd [molec/cm**2]' + 20*' '+len(atmos.gases)*'%10.4g') % tuple(vcd)

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
		molCrossSections = lbl_xs (Line_Data, Mol_Data, atmos.p, atmos.T, xLimits, lineShape,
					   sampling, wingExt, nGrids, gridRatio, nWidths, lagrange, verbose)
		xsMatrix.append(molCrossSections)
	# for xss in xsMatrix: print 'xs', type(xss), len(xss), [(xs.get('molecule'),len(xs.get('y'))) for xs in xss]

	# compute absorption coefficients, also return uniform, equiidistant wavenumber grid corresponding to most dense cross section
	vGrid, absCoeff = sum_xsTimesDensity (xsMatrix, atmos.densities, interpolate=interpolate, verbose=verbose)

	# integrate absorption coefficients along line of sight
	optDepth = integral_acTimesDistance (absCoeff, atmos.z, quad, mode, verbose)

	# save optical depth
	if outFile and os.path.splitext(outFile)[1].lower() in ('.nc', '.ncdf', '.netcdf'):
		save_optDepth_netcdf (vGrid, optDepth, outFile, atmos, nanometer)
	else:
		save_optDepth_xy (vGrid, optDepth, outFile, atmos, nanometer, yOnly, flipUpDown, commentChar)

	if plot:
		try:
			from pylab import plot, semilogy, show
			if mode in 'cdr':
				for j in range(optDepth.shape[1]):  semilogy (vGrid, optDepth[:,j])
				show()
			else:
				plot (vGrid, optDepth); show()
		except ImportError:
			raise SystemExit, 'ERROR --- lbl2od:  matplotlib not available, no quicklook!'

####################################################################################################################################

def save_optDepth_netcdf (vGrid, optDepth, outFile=None, atmos=None, nanometer=False):
	""" Save optical depth (vs wavenumber or wavelength) in netcdf format. """
	#print 'save_optDepth_netcdf', optDepth.shape

	# Open netcdf file for writing
	ncf = NetCDFFile(outFile,'w')

	# Create Dimensions
	ncf.createDimension('nlev',optDepth.shape[1]+1)  # fgs: changed, nLevels=len(zGrid)
	ncf.createDimension('nlyr',optDepth.shape[1])    #               nLayers=number of intervals of zGrid
	if nanometer:  ncf.createDimension('nwvl',optDepth.shape[0])
	else:          ncf.createDimension('nwvn',optDepth.shape[0])
	ncf.createDimension('none',1)

	# Create Variables:
	z      = ncf.createVariable('z', 'd', ('nlev',))
	if nanometer:
		tau    = ncf.createVariable('tau', 'd', ('nlyr', 'nwvl',))
		wvlMin = ncf.createVariable('wvlmin', 'd', ('none',))
		wvlMax = ncf.createVariable('wvlmax', 'd', ('none',))
		wvl    = ncf.createVariable('wvl', 'd', ('nwvl',))
	else:
		tau    = ncf.createVariable('tau', 'd', ('nlyr', 'nwvn',))
		wvnMin = ncf.createVariable('wvnmin', 'd', ('none',))
		wvnMax = ncf.createVariable('wvnmax', 'd', ('none',))
		wvn    = ncf.createVariable('wvn', 'd', ('nwvn',))

	# Assign data
	z.assignValue(atmos.z)
	if nanometer:
		depthOpt = np.flipud(optDepth)
		tau.assignValue(depthOpt.T)
		vGrid = 1e7/vGrid # convert wavenumber -> wavelength (note: simply overwrite variable!)
		wvl.assignValue(vGrid)
		wvlMin.assignValue(vGrid[0])
		wvlMax.assignValue(vGrid[-1])
	else:
		tau.assignValue(optDepth.T)
		wvn.assignValue(vGrid)
		wvnMin.assignValue(vGrid[0])
		wvnMax.assignValue(vGrid[-1])

	# Set units as attributes
	setattr(z, 'units', 'km')
	setattr(tau, 'units', '-')
	if nanometer:
		setattr(wvlMin, 'units', 'nm')
		setattr(wvlMax, 'units','nm')
		setattr(wvl, 'units', 'nm')
	else:
		setattr(wvnMin, 'units', '1/cm')
		setattr(wvnMax, 'units','1/cm')
		setattr(wvn, 'units', '1/cm')

	# Close netcdf file 
	ncf.close()



####################################################################################################################################

def save_optDepth_xy (vGrid, optDepth, outFile=None, atmos=None, nanometer=False, yOnly=False, flipUpDown=False, commentChar='#'):
	""" Save optical depth (vs wavenumber or wavelength) in ascii tabular format. """
	out = open_outFile (outFile, commentChar)
	
	nLevels = len(atmos.z)
	nGas    = len(atmos.gases)
	if flipUpDown:
		if len(optDepth.shape)==2:  optDepth = np.fliplr(optDepth) # swap columns
		comments = ['altitudes [km]:  ' + nLevels*' %10.1f' % tuple(unitConversion(atmos.z[::-1],'length',new='km')),
	                    'temperatures [K]:' + nLevels*' %10.2f' % tuple(atmos.T[::-1]),
	                    'pressures [mb]:  ' + nLevels*' %10.4g' % tuple(unitConversion(atmos.p[::-1],'p',new='mb'))]
	else:
		comments = ['altitudes [km]:  ' + nLevels*' %10.1f' % tuple(unitConversion(atmos.z,'length',new='km')),
	                    'temperatures [K]:' + nLevels*' %10.2f' % tuple(atmos.T),
	                    'pressures [mb]:  ' + nLevels*' %10.4g' % tuple(unitConversion(atmos.p,'p',new='mb'))]
	comments += ['gases:           ' +    nGas*' %10s'   % tuple(atmos.gases),
	             'vcd [molec/cm^2]:' +    nGas*' %10.2g' % tuple(atmos.vcd())]
	if nanometer:
		if yOnly: 
			writeArray (np.flipud(optDepth),                  out, comments=comments, format= '%10g')
		else:
			comments += ['%10s %10s' % ('wavelength','optical depth'),  '%10s' % 'nm']
			writeArray (np.flipud(paste(1e7/vGrid,optDepth)), out, comments=comments, format= '%10f %10g')
	else:
		if yOnly:
			writeArray (optDepth,              out, comments=comments, format= '%10g')
		else:	
			comments += ['%10s %10s' % ('wavenumber','optical depth'),  '%10s' % '1/cm']
			writeArray (paste(vGrid,optDepth), out, comments=comments, format= '%10f %10g')

	out.close()

####################################################################################################################################

if (__name__ == "__main__"):

	from command_parser import parse_command, standardOptions, multiple_outFiles 

        # parse the command, return (ideally) one file and some options
        opts = standardOptions + [  # h=help, c=commentChar, o=outFile
	       {'ID': 'm', 'name': 'mode', 'type': StringType, 'default': 'd', 'constraint': 'mode.lower() in "cdrt"'},
               {'ID': 'ToA', 'name': 'zToA', 'type': FloatType, 'constraint': 'zToA>0.0'},
               {'ID': 'a', 'name': 'airWidth', 'type': FloatType, 'constraint': 'airWidth>0.0'},
               {'ID': 'L', 'name': 'lineShape', 'type': StringType, 'constraint': 'lineShape in ["V","L","G"]', 'default': 'V'},
               {'ID': 'i', 'name': 'interpolate', 'type': StringType, 'default': '3', 'constraint': 'interpolate in "234lLqQcCbBsS"'},
               {'ID': 's', 'name': 'sampling', 'type': FloatType, 'constraint': 'sampling>0.0', 'default': 5.0},
               {'ID': 'w', 'name': 'wingExt', 'type': FloatType, 'constraint': 'wingExt>0.0', 'default': 5.0},
	       {'ID': 'x', 'name': 'xLimits', 'type': Interval, 'constraint': 'xLimits.lower>=0.0'},
               {'ID': 'n', 'name': 'nGrids', 'type': IntType, 'default': 3, 'constraint': 'nGrids>0'},
               {'ID': 'g', 'name': 'gridRatio', 'type': IntType, 'default': 8, 'constraint': 'gridRatio in [4,8]'},
               {'ID': 'q', 'name': 'quad', 'type': StringType, 'default': 'T', 'constraint': 'quad.upper() in ["T","S","H"]'},
               {'ID': 'W', 'name': 'nWidths', 'type': FloatType, 'default': 25.0,'constraint': 'nWidths>2.0'},
               {'ID': 'r', 'name': 'flipUpDown'},
               {'ID': 'y', 'name': 'yOnly'},
               {'ID': 'plot'},
               {'ID': 'nm', 'name': 'nanometer'},
               {'ID': 'v', 'name': 'verbose'}
               ]
 	
	inFiles, options, commentChar, outFile = parse_command (opts,(2,99))

	if options.has_key('h'):  print __doc__%globals();  raise SystemExit, " end of lbl2od help"

	# translate some options to boolean flags
	for key in ['yOnly', 'nanometer', 'plot', 'flipUpDown', 'verbose']:  options[key] = options.has_key(key)

	if commonExtension(inFiles):
		print '\nWARNING:  all input files have the same extension, first file probably NOT an atmospheric data file!!!\n'

	# unpack list of input files
	atmFile, lineFiles = inFiles[0], inFiles[1:]
	
	# optionally already import module used at very end (check availablity of netcdf before starting a lengthy job)
	if outFile and os.path.splitext(outFile)[1].lower() in ('.nc', '.ncdf', '.netcdf'):
		try:
			from Scientific.IO.NetCDF import *
		except ImportError:
			print 'INFO:  import Scientific.IO.NetCDF failed, trying pynetcdf'
			try:                 from pynetcdf import *
			except ImportError:  raise SystemExit, 'import "pynetcdf" failed, cannot find module!'
			else:                print 'INFO:  import pynetcdf'

	_lbl2od_ (atmFile, lineFiles, outFile, commentChar, **options)
