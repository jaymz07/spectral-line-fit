#!/usr/bin/env python

"""
  cgsUnits
  convert (physical) data between (compatible) units

  usage:
  cgsUnits oldUnit newUnit value(s)

  NOTE:  the (physical) quantities and units supported are far from complete.
         (see the source file for quantities and units known)
"""

from ir import *

####################################################################################################################################

cgs_units = {'length': 'cm', 'pressure': 'g/(cm*s**2)', 'temperature': 'K', 'density': '1/cm**3', 'vmr': 'pp1', 'energy': 'erg', 'power': 'erg/s'} 

# conversion factors for cgs units
pressureUnits    = {'g/cm/s**2': 1.0, 'g/(cm.s**2)': 1.0, 'mb': 1.e3, 'mbar': 1.e3, 'hPa': 1.e3,  'atm': 1013250., 'Pa': 10., 'N/m**2': 10.}
frequencyUnits   = {'Hz': 1.0, 'kHz': 1.0e3, 'MHz': 1.0e6, 'GHz': 1.0e9, 'THz': 1.0e12, 'cm-1': c, '1/cm': c}
wavelengthUnits  = {'cm': 1.0, 'mm': 0.1, 'mue': 1.e-4, 'mum': 1.e-4, 'micrometer': 1.e-4,  'nm': 1.e-7, 'A': 1.e-8}
lengthUnits      = {'km': 1.e5, 'm': 1.e2, 'dm': 10., 'inch': 2.54};  lengthUnits.update(wavelengthUnits)
areaUnits        = dict([(nam+'**2', val**2) for nam,val in lengthUnits.items()])
volumeUnits      = dict([(nam+'**3', val**3) for nam,val in lengthUnits.items()]);  volumeUnits.update({'l': 1.e3, 'hl': 1.e5})
mixingRatioUnits = {'ppv': 1.0, 'ppV': 1.0, 'pp1': 1.e-2, '%': 1.e-6, 'ppm': 1.e-6, 'ppb': 1.e-9, 'ppt': 1.e-12, 'vmr': 1.0}
densityUnits     = dict( [('1/'+nam+'**3', val**-3) for nam,val in lengthUnits.items()] +
                         [(nam+'-3', val**-3) for nam,val in lengthUnits.items()])
temperatureUnits = {'Kelvin': 0.0, 'K': 0.0, 'k': 0.0, 'C': 273.15, 'Celsius': 273.15} # hmmm, lower case 'k' only to satisfy libradtran atmospheric data files
energyUnits      = {'erg': 1.0, 'g.cm**2/s**2': 1.0, 'kg.m**2/s**2': 1.e7, 'Nm': 1.e7, 'N.m': 1.e7, 'J': 1.e7, 'mJ': 1.e4}
powerUnits       = {'erg/s': 1.0, 'g.cm**2/s**3': 1.0, 'kg.m**2/s**3': 1.e7, 'W': 1.e7, 'mW': 1.e4, 'nW': 1.e-2}

allUnits = {'length':      lengthUnits,
            'altitude':    lengthUnits,
            'area':        areaUnits,
	    'volume':      volumeUnits,
	    'pressure':    pressureUnits, 'p':    pressureUnits,
	    'density':     densityUnits,
	    'mixingratio': mixingRatioUnits,
	    'frequency':   frequencyUnits,
	    'power':       powerUnits,
	    'energy':      energyUnits,
	    'wavelength':  wavelengthUnits}
 
####################################################################################################################################
def unitConversion (data, WHAT, old=None, new=None):
	""" Conversion of (scalar or array) physical values to different units.
	    If old (input) or new (output) unit is not given: assume cgs standard unit for this quantity. """
	# alternative approaches (see also discussion in H.P. Langtangen's book: 
	# UNUM     http://pypi.python.org/pypi/Unum/4.1.0
	# ScientificPython: http://dirac.cnrs-orleans.fr/ScientificPython/
	what = WHAT.lower()
	if what in allUnits.keys():
		xUnits = allUnits[what]
		if old==new:
			return data
		elif old in xUnits and new in xUnits:
			return data * xUnits[old] / xUnits[new]
		elif old in xUnits and not new:
			return data * xUnits[old] # convert to cgs 'base' unit
		elif new in xUnits and not old:
			return data / xUnits[new]  # convert from cgs 'base' unit
		else:
			print '\n', WHAT, xUnits.keys()
			raise SystemExit, '%s %s %s %s' % ('ERROR --- unitConversion:  unknown/unsupported units ', old, ' ---> ', new)
	elif WHAT=='T' or what[:4]=='temp':
		return cgsTemperature (data, old, new)
	else:
		raise SystemExit, 'ERROR --- unitConversion failed, unknown/unsupported quantity ' + WHAT

####################################################################################################################################
		
def cgsTemperature (data, old=None, new=None):
	if old==new:
		return data
	elif old in temperatureUnits and new in temperatureUnits:
		return data + temperatureUnits[old] - temperatureUnits[new]
	elif old in temperatureUnits and not new:
		return data - temperatureUnits[old]
	elif new in temperatureUnits and not old:
		return data + temperatureUnits[new]
	else:
		raise SystemExit, '%s %s %s %s' % ('ERROR --- cgsTemperature:  unit conversion failed, unknown/unsupported units ', old, ' ---> ', new)


####################################################################################################################################

# some 'old' conversion routines, still in use somewhere

def cgsLength (data, old=None, new=None):
	if old==new:
		return data
	elif old in lengthUnits and new in lengthUnits:
		return data * lengthUnits[old] / lengthUnits[new]
	elif old in lengthUnits and not new:
		return data * lengthUnits[old]
	elif new in lengthUnits and not old:
		return data / lengthUnits[new]
	else:
		raise SystemExit, '%s %s %s %s' % ('ERROR --- cgsLength:  unit conversion failed, unknown/unsupported units ', old, ' ---> ', new)

def cgsPressure (data, old=None, new=None):
	if old==new:
		return data
	elif old in pressureUnits and new in pressureUnits:
		return data * pressureUnits[old] / pressureUnits[new]
	elif old in pressureUnits and not new:
		return data * pressureUnits[old]
	elif new in pressureUnits and not old:
		return data / pressureUnits[new]
	else:
		raise SystemExit, '%s %s %s %s' % ('ERROR --- cgsPressure:  unit conversion failed, unknown/unsupported units ',  old,  ' ---> ', new)

def change_frequency_units (x, xUnitOld, xUnitNew):
	if xUnitOld=='cm-1':
		if   xUnitNew=='cm-1': pass
		elif xUnitNew=='THz':  x = x*c * 1e-12
		elif xUnitNew=='GHz':  x = x*c * 1e-9
		elif xUnitNew=='MHz':  x = x*c * 1e-6
		elif xUnitNew=='Hz':   x = x*c
		elif xUnitNew=='nm':   x = 10000000./x
		elif xUnitNew in ['mue','micro','mum']: x = 10000./x
		else: raise SystemExit, 'ERROR: unknown/invalid unit for wavenumber/frequency/wavelength!'
	elif xUnitOld.endswith('Hz'):
		if   xUnitOld=='Hz':   pass
		elif xUnitOld=='kHz':  x = x*1e3
		elif xUnitOld=='MHz':  x = x*1e6
		elif xUnitOld=='GHz':  x = x*1e9
		elif xUnitOld=='THz':  x = x*1e12
		if xUnitNew=='cm-1':  x = x/c
		else: raise SystemExit, 'ERROR: conversion ' + xUnitOld + ' --> ' + xUnitNew + ' not yet implemented!'
	elif xUnitOld in ['mue','micro','mum']:
		if   xUnitNew=='cm-1':  x = 10000./x
		else: raise SystemExit, 'ERROR: conversion ' + xUnitOld + ' --> ' + xUnitNew + ' not yet implemented!'
	elif xUnitOld=='nm':
		if   xUnitNew=='cm-1':  x = 1.e7/x
		else: raise SystemExit, 'ERROR: conversion ' + xUnitOld + ' --> ' + xUnitNew + ' not yet implemented!'
	else:
		raise SystemExit, 'ERROR: unknown/invalid unit for wavenumber/wavelength,frequency ' + xUnitOld
	return x

####################################################################################################################################

if (__name__ == "__main__"):
	import sys, os
	import numpy as np

	args = sys.argv[1:]

	if '-h' in sys.argv or '--help' in sys.argv:
		print __doc__%globals()
		if args[0]=="--help":  print 'default units in the cgs systems:\n', cgs_units
		raise SystemExit, "end of cgsUnits help"

	if len(args)>2:
		old, new = args[:2]
		values   = np.array(map(float,args[2:]))
		if old in temperatureUnits.keys() and new in temperatureUnits.keys():
			print 'temperature: ', cgsTemperature(values, old, new)
			sys.exit()
		what=''
		for name, units in allUnits.items():
			if old in units.keys() and new in units.keys():  what = name;  break
		if what:
			print '%s[%s]: ' % (what,new), unitConversion(values, what, old, new)
		else:
			raise SystemExit, '%s "%s" %s "%s" %s' % ('old', old, 'and/or new', new, 'unit unknown or incompatible')
	else:
		raise SystemExit, 'need at least three inputs:  oldUnit, newUnit, value(s)'
