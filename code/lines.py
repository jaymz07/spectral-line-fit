#!/usr/bin/env python

"""  lines
  Read molecular spectroscopic line parameters (Hitran, Geisa, ... extract)  and convert to new pressure/temperature

  usage:
  lines [options] line_parameter_file(s)

  -h          help
  -c char     comment character(s) used in input, output file (default '#')
  -o string   output file for saving of line data (if not given: write to StdOut)

  -m string   molecule (no default, should be given in the line file, required otherwise)
  -p float    pressure (in mb,  default: p_ref of linefile, usually 1013.25mb=1atm)
  -T float    Temperature (in K, default: T_ref of linefile, usually 296K)
  -x Interval lower, upper wavenumbers (comma separated pair of floats [no blanks!],
                                        default set according to range of lines in datafile)

  NOTES:
  in:  The line parameter file(s) should contain a list of (preselected) lines
       that can be generated from HITRAN or GEISA database with extract.py 
       i.e., the original HITRAN or GEISA data bases cannot be used as input files
       (See the documentation header of lbl2xs.py for more details)

  out: The output file(s) are formatted automatically according to the extension given for the output file(s):
       if the extension is vSL or vSLD or vSLG, then three or four columns with position, strengths, Lorentz (and Doppler) width are written
       otherwise, three columns with position and strengths at T_ref and T are written
       (actually there is a fourth column with zeros to facilitate plotting of delta S as (one sided) error bar)
"""

####################################################################################################################################
import os
import sys
from string import *
from types import *
from math import pi, sqrt, log, ceil
import numpy as np
import time

from ir import *
from molecules import get_molec_data
from IO import readDataAndComments, parse_file_header, open_outFile
from cgsUnits import unitConversion
from pairTypes import Interval

# and some math constants
ln2     = log(2.0)
sqrtLn2 = sqrt(ln2)

####################################################################################################################################

def get_lbl_data (file, xLimits=None, airWidth=0.1, wingExt=5.0, molecule=None, commentChar='#'):
	""" Read spectroscopic line parameters and convert to appropriate pressure and temperature.  """
	# read line list and check consistency of molecule specs
	Line_Data = read_line_file (file, xLimits, airWidth, wingExt, molecule, commentChar) 

	# read molecular data
	Mol_Data = get_molec_data (Line_Data['molecule'])
	print ' %-s:  %.1famu  %i isotopes  %i modes\n' % \
		(Line_Data['molecule'], Mol_Data['mass'], len(Mol_Data['isotopes']), len(Mol_Data['VibFreq']))

	# line parameters at actual atmospheric conditions
	return Line_Data, Mol_Data

####################################################################################################################################

def read_line_file (lineFile, xLimits=None, airWidth=0.1, wingExt=0.0, molecule=None,  commentChar='#'):
	""" Read a simple line parameter list and return a dictionary:  lines and some attributes. """

	# wavenumber interval to be searched
	if xLimits:
		xLow, xHigh = xLimits.limits()
	else:
		xLow, xHigh = 0.0, 0.0

	# read entire line file and return a dictionary
	lines, comments = readDataAndComments (lineFile,commentChar)
	# if there is just a single line in the dataset, a 1dim array is returned
	lines = np.atleast_2d (lines)

	# parse comment header and extract some infos
	Line_Data = parse_file_header (comments, ['molecule','format','pressure','temperature'])

	# check if molecule is specified  (and consistent if specified in file and on command line)
	if molecule and Line_Data.has_key('molecule'):
		if not molecule==Line_Data['molecule']:
			raise SystemExit, 'ERROR:  inconsistent molecule specification!   ' + \
			                  repr(molecule) + ' <--> ' + repr(Line_Data['molecule'])
	elif Line_Data.has_key('molecule'):
		pass
	elif molecule:
		Line_Data['molecule'] = molecule
	else:
		print 'lineFile: ', lineFile
		raise SystemExit, 'ERROR:  molecule not specified (neither in line list header nor as command option)!'

	# also need reference pressure and temperature
	if Line_Data.has_key('temperature'):
		# remove unit spec 'Kelvin'
		Line_Data['temperature'] = float(split(Line_Data['temperature'])[0])
	else:
		raise SystemExit, 'ERROR:  reference temperature of line parameters not given!'
	if Line_Data.has_key('pressure'):
		try:    # remove unit spec 'millibar' and return pressure in cgs units!
			value_unit = split(Line_Data['pressure'])
			Line_Data['pressure'] = unitConversion(float(value_unit[0]), 'pressure', value_unit[1])
		except StandardError,errMsg: 
			raise SystemExit, str(errMsg) + '\nparsing pressure spec in line file failed ' + repr(Line_Data['pressure'])
	else:
		raise SystemExit, 'ERROR:  reference pressure of line parameters not given!'

	# copy columns into appropriate arrays and put into dictionary
	if Line_Data.has_key('format'):
		Line_Data.update (dataColumms_to_dict (lines, upper(Line_Data['format'])))
	elif data['columns']==2:
		Line_Data['position'] = lines[:,0]
		Line_Data['strength'] = lines[:,1]
		print '\n !!! WARNING: no column specification, ',\
		      'assuming line position and strength in first and second column!'
	else:
		raise SystemExit, 'ERROR:  no format specifier in line list!'

	# check if at least position and strengths are found
	if Line_Data.has_key('position') and Line_Data.has_key('strength'):
		firstLine = min(Line_Data['position'])
		lastLine  = max(Line_Data['position'])
		print '\n %-20s %8i lines in %10f ... %10f cm-1  with %8.2g < S < %8.2g   (T=%5.1fK)' % \
		       ( os.path.basename(lineFile), lines.shape[0], firstLine, lastLine, \
		         min(Line_Data['strength']), max(Line_Data['strength']), Line_Data['temperature'])
	else:
		raise SystemExit, 'ERROR:  Need at least line positions and strengths!'

	# set air broadening parameter to default if not found in line file
	if Line_Data.has_key('airWidth'):
		print '%s %8g < g < %8g   (p=%gmb)' % \
		       ( 75*' ', min(Line_Data['airWidth']), max(Line_Data['airWidth']), unitConversion(Line_Data['pressure'],'p', new='mb'))
	else:
		Line_Data['airWidth'] = np.zeros_like(Line_Data['position']) + airWidth
		print 'Air width initialized to ', airWidth

	if Line_Data.has_key('TempDep'):
		print '%s %8g < n < %8g' % ( 75*' ', min(Line_Data['TempDep']), max(Line_Data['TempDep']))

	if xLow and xHigh:
		# subset of lines requested, truncate lists
		print ' get_line_data:  xLow,xHigh specified: ', xLow, xHigh, '  (extension:  +/-', wingExt,')'
		Line_Data   = truncate_lineList (Line_Data,  xLow-wingExt, xHigh+wingExt)

	return Line_Data

####################################################################################################################################

def truncate_lineList (Line_Data, xLow, xHigh):
	""" Remove some lines at head and/or tail of line list. """
	# setup the mask
	freqMask = np.logical_and (Line_Data['position']>=xLow, Line_Data['position']<=xHigh)

	if sum(freqMask)<len(Line_Data['position']):
		print ' line data subset:  select', sum(freqMask), ' of', len(Line_Data['position']), ' lines in', xLow, xHigh
		# mandatory line parameters
		Line_Data['position'] = Line_Data['position'][freqMask]
		Line_Data['strength'] = Line_Data['strength'][freqMask]
		# optional line parameters
		if Line_Data.has_key('airWidth'):  Line_Data['airWidth']  = Line_Data['airWidth'][freqMask]
		if Line_Data.has_key('selfWidth'): Line_Data['selfWidth'] = Line_Data['selfWidth'][freqMask]
		if Line_Data.has_key('energy'):    Line_Data['energy']    = Line_Data['energy'][freqMask]
		if Line_Data.has_key('TempDep'):   Line_Data['TempDep']   = Line_Data['TempDep'][freqMask]

	return Line_Data

####################################################################################################################################

def dataColumms_to_dict (lines, format):
	""" Copy columns of line parameters into appropriate dictionary elements. """
	if len(format) == lines.shape[1]:
		Line_Dict = {}
		if find(format,'V')>-1: Line_Dict['position'] = lines[:,find(format,'V')]
		if find(format,'X')>-1: Line_Dict['position'] = lines[:,find(format,'X')]
		if find(format,'S')>-1: Line_Dict['strength'] = lines[:,find(format,'S')]
		if find(format,'E')>-1: Line_Dict['energy']   = lines[:,find(format,'E')] 
		if find(format,'A')>-1: Line_Dict['airWidth'] = lines[:,find(format,'A')]
		if find(format,'N')>-1: Line_Dict['TempDep']  = lines[:,find(format,'N')]
	else:
		raise SystemExit, 'ERROR: format spec inconsistent with number of columns ' + `lines.shape[1]`
	return Line_Dict

####################################################################################################################################

def line_widths (pressRef, TempRef, press, Temp, mass, positions, airWidths=0.1, TempExp=0.5):
	""" Convert pressure (air) broadening and set Doppler broadening half widths to actual p, T. """
	# Lorentzian half widths: only air broadening (self broadening is ignored!)
	if isinstance(airWidths,float):
		print 'Air widths initializing to ', airWidths
		airWidths = np.zeros_like(positions) + airWidths
	if isinstance(TempExp,float):
		print 'Air width temperature exponent initializing to ', TempExp
		TempExp = np.zeros_like(positions) + TempExp
	gammaL = airWidths * (press/pressRef) * (TempRef/Temp)**TempExp
	# Gaussian half widths (note: molecular mass in atomic mass units!)
	gammaD = positions * sqrt(2.*ln2*k*Temp/(mass*amu*c*c))
	return gammaL, gammaD

####################################################################################################################################

def strengths_and_widths (Line_Data, Mol_Data, pressure=1013.25e3, temperature=296.0):
	""" Convert line strength and Lorentzian width to pressure/temperature and set Doppler width. """
	Strengths = line_strengths (Line_Data['temperature'], temperature,
	                            Line_Data['position'], Line_Data['strength'], Line_Data.get('energy'),
	                            Mol_Data['NumDeg'], Mol_Data['VibFreq'], Mol_Data['TempExpQR'])
	gammaL, gammaD = line_widths  (Line_Data['pressure'], Line_Data['temperature'], pressure, temperature,
	                               Mol_Data['mass'], Line_Data['position'], # mandatory parameters
				       Line_Data.get('airWidth',0.1), Line_Data.get('TempDep',0.5))
	return Strengths, gammaL, gammaD

####################################################################################################################################

def write_strengths (positions, strengths, strRef, temperature, tempRef, molecule, outFile=None, commentChar='#'):
	""" Write line strengths (for two temperatures) vs line positions. """
	out = open_outFile (outFile, commentChar)
	out.write ('%s %s %s\n' % (commentChar, "molecule:", molecule))
	out.write ('%s %s %8.2f %s\n' % (commentChar, "temperature T_ref:" ,  tempRef, "K"))
	out.write ('%s %s %8.2f %s\n' % (commentChar, "temperature T:    ", temperature, "K"))
	out.write ('%s %10s %23s\n' % (commentChar, "position", "   S(T_ref)    |S(T)- S(T_ref)|"))
	format = '%12f  %11.3e %10.2e 0\n'
	for v,S0,S in zip(positions, strRef, strengths): 
		#out.write ( format % (v,S0,S) )
		if S>S0: out.write ( format % (v,S0,S-S0) )
		else:    out.write ( format % (v,S0,S0-S) )
	# close the output file (if its not stdout)
	if outFile: out.close()

def write_voigtLines (positions, strengths, gammaL, gammaD, pressure, temperature, molecule, outFile=None, commentChar='#'):
	""" Write Voigt line parameters (strengths, Lorentz and Gaussian widths vs line positionsp. """
	out = open_outFile (outFile, commentChar)
	out.write ('%s %s %s\n' % (commentChar, "molecule:", molecule))
	out.write ('%s %s %-12g\n' % (commentChar, "pressure  [mb]:     " ,  unitConversion(pressure,'p', new='mb')))
	out.write ('%s %s %8.2f\n' % (commentChar, "temperature  [K]: ", temperature))
	out.write ('%s %10s  %11s %11s %11s\n' % (commentChar, "position", "strength", "gammaL", "gammaG"))
	format = '%12f  %11.3e %11.3g %11.3g\n'
	for v,S, gL, gG in zip(positions, strengths, gammaL, gammaD):  out.write ( format % (v,S,gL,gG) )
	# close the output file (if its not stdout)
	if outFile: out.close()

def write_vSL (positions, strengths, gammaL, pressure, temperature, molecule, mass, outFile=None, commentChar='#'):
	""" Write Voigt line parameters (strengths, Lorentz and Gaussian widths vs line positionsp. """
	out = open (outFile, 'w')
	out.write ('%s %s %f %f %f\n' % (commentChar, molecule, pressure, temperature, mass))
	format = '%12f  %11.3e %11.3g\n'
	for v,S, gL in zip(positions, strengths, gammaL):  out.write ( format % (v,S,gL) )
	# close the output file (if its not stdout)
	if outFile: out.close()

####################################################################################################################################

def line_strengths (TempRef, Temp, positions, strengths, energies, numDeg, vibFreq, TempExpQR=1.0):
	""" Convert line strengths to actual p, T. """
	if abs(TempRef-Temp)>0.1:
		if energies.any():
			# ratio of rotational partition function
			ratioQR  =  (TempRef/Temp)**TempExpQR
			# ratio of vibrational partition function
			ratioQV  =  vibPartitionFunction (numDeg, vibFreq, TempRef) / vibPartitionFunction (numDeg, vibFreq, Temp)
			# Boltzmann factor
			deltaInvTemp = C2 * (Temp-TempRef) / (Temp*TempRef)
			sb = np.exp(deltaInvTemp*energies)
			# stimulated emission factor
			#se = (1.0 - np.exp(-C2*positions/Temp)) / (1.0 - np.exp(-C2*positions/TempRef))
			se = np.expm1(-C2*positions/Temp) / np.expm1(-C2*positions/TempRef)
			# multiply all conversion factors with original line strength
			return strengths * sb * se * ratioQR * ratioQV
		else:
			raise SystemExit, 'ERROR:  line strength temperature conversion impossible, no lower state energies!'
	else:
		return strengths

####################################################################################################################################

def  vibPartitionFunction (degeneracy, omega, temperature):
	""" Calculate vibrational partition function as a function of temperature.
	    Norton & Rinsland: ATMOS Data Processing and Science Analysis Methods; Appl. Opt. 30,389(1991) """
	c2T        = C2 / temperature
	factors    = 1.0 / (1.0 - np.exp(-c2T*omega))**degeneracy
	qVib      = np.product(factors)
	return qVib

####################################################################################################################################

def meanLineWidth (gammaL, gammaD):
	""" Evaluate mean pressure, Doppler, and combined broadening half width, return mean Voigt width. """
	# Voigt width (Whiting approximation, 1% accuracy)
	gammaV = 0.5 * (gammaL + np.sqrt(gammaL**2 + 4.*gammaD**2))
	# averages
	nLines = len(gammaL)
	meanGL = np.sum(gammaL)/nLines
	meanGD = np.sum(gammaD)/nLines
	meanGV = np.sum(gammaV)/nLines
	print  ' mean width (L, G, V): %10g %10g %10g   y %8g' % (meanGL,meanGD,meanGV, sqrtLn2*meanGL/meanGD)
	return  meanGL, meanGD , meanGV

####################################################################################################################################

if (__name__ == "__main__"):

	from command_parser import parse_command, standardOptions, multiple_outFiles 
        # parse the command, return (ideally) one file and some options
        opts = standardOptions + [  # h=help, c=commentChar, o=outFile
               {'ID': 'a', 'name': 'airWidth', 'type': FloatType, 'constraint': 'airWidth>0.0'},
               {'ID': 'm', 'name': 'molecule', 'type': StringType},
               {'ID': 'p', 'name': 'pressure', 'type': FloatType, 'constraint': 'pressure>0.0', 'default': 1013.25},
               {'ID': 'T', 'name': 'temperature', 'type': FloatType, 'constraint': 'temperature>0.0', 'default': 296.0},
	       {'ID': 'x', 'name': 'xLimits', 'type': Interval, 'constraint': 'xLimits.lower>=0.0'}
               ]
 	
	files, options, commentChar, outFile = parse_command (opts,(1,99))
	outFiles = multiple_outFiles (files, outFile)

	if options.has_key('h'): print __doc__%globals();  raise SystemExit, " end of lines help"

	# put named options into global variables
	for opt in opts:
		if opt.has_key('name') and opt.has_key('type'):
			exec opt['name'] + ' = ' + repr(options.get(opt['name']))

	wingExt = 0.0
	
	for inFile,outFile  in  zip(files, outFiles):
		# read a couple of lines from vSEan file
		Line_Data, Mol_Data = get_lbl_data (inFile, xLimits, airWidth, wingExt, molecule, commentChar)

		# extract line position
		positions = Line_Data['position']
		# adjust line parameters to p, T
		strengths, gammaL, gammaD = strengths_and_widths (Line_Data, Mol_Data, pressure, temperature)
		
		if outFile and outFile.lower().endswith('.vsl'):
			write_vSL (positions, strengths, gammaL, pressure, temperature, Line_Data['molecule'], Mol_Data['mass'], outFile, commentChar)
		elif outFile and (outFile.lower().endswith('.vsld') or outFile.lower().endswith('.vslg')):
			write_voigtLines (positions, strengths, gammaL, gammaD, pressure, temperature, Line_Data['molecule'], outFile, commentChar)
		else:
			write_strengths (positions, strengths, Line_Data['strength'], temperature, Line_Data['temperature'],
			                 Line_Data['molecule'], outFile, commentChar)
