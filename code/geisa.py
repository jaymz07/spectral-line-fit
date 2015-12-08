#!/usr/bin/env python

""" geisa

    usage:
    geisa [options] line_parameter_database

    command line options without arguments:
        -h  help

    optional command line options with string argument:
	-i  isotope number
	-m  molecule number
        -o  output file (default: standard output)

    optional command line options with float argument:
	-S  minimum line strength to accept
	-x  lower and upper end of wavenumber range (comma separated pair without blanks)

    Note: at least wavenumber range or molecule has to be selected!
"""
##############################################################################################################
#####  ToDo:                                                                                             #####
#####                                                                                                    #####
gIASI2003_fieldLength = map(int,'12 11 6 10 9 9 9 9 4 3 3 1 2 1 10 5 8 3 6 6 10 11 6 4 8 6 5 4 4 8 8 4 4'.split()) # sum 209
geisa2003_fieldLength = map(int,'12 11 6 10 9 9 9 9 4 3 3 3 2 1 10 5 8 3 6 6 10 11 6 4 8 6 5 4 4 8 8 4 4'.split()) # sum 211
geisa2008_fieldLength = map(int,'12 11 6 10 25 25 15 15 4 3 3 3 2 1 10 7 9 6 10 11 6 4 9 6 7 4 4 8 8 4 4'.split()) # sum 252
##############################################################################################################

import sys, os
from string import *
from types  import *
from time   import clock

from hitran import bisect_first_line
from pairTypes import Interval
 
##########################################################################################################################################

####################################################################################################################################

def extract_geisa (file, xLimits, molNr=0, isoNr=0, strMin=0.0):
	""" Read Geisa formatted database. """
	try:
	    geisa  = open (file)
	except IOError:
	    raise SystemExit, 'opening GEISA data file "' + file + '" failed!'
	# wavenumber interval to be searched
	xLow, xHigh = xLimits.limits()
	# initialize time and search first useful line
	time0  = clock()
 	if molNr<=0:
		if isoNr>0 or strMin>0:
			raise SystemExit, 'ERROR --- geisa:  no isotope or linestrength selection without molecule specification!'
	 	lines = extract_range (geisa, xLow, xHigh)
 	elif not isoNr and strMin<=0.0:
	 	lines = extract_Mol (geisa, xLow, xHigh, molNr)
 	elif not isoNr and strMin>0.0:
	 	lines = extract_MolStr (geisa, xLow, xHigh, molNr, strMin)
 	elif strMin<=0.0:
	 	lines = extract_MolIso (geisa, xLow, xHigh, molNr, isoNr)
 	else:
	 	lines = extract_MolIsoStr (geisa, xLow, xHigh, molNr, isoNr, strMin)
	geisa.close()
	return lines

##########################################################################################################################################

def set_geisa_fields (fileName):
	""" Use length of first accepted record to determine position (first and last indices) of wavenmber, moleculeID, isotopeID. """
	beginFreq = 0
	if '97' in fileName:
		endFreq = 10
		beginStr, endStr = 10, 20
		beginAir, endAir = 20, 25
		beginEne, endEne = 25, 35
		beginIso, endIso = 75, 79
		beginMol, endMol = 79, 82
		beginTEx, endTEx = 72, 75
		beginSlf, endSlf = 98,103
		#print 'assuming GEISA 97 format!'
	elif '2003' in fileName:
		endFreq=12
		beginStr, endStr = 12, 23
		beginAir, endAir = 23, 29
		beginEne, endEne = 29, 39
		beginTEx, endTEx = 75, 79
		beginIso, endIso = 79, 82
		beginMol, endMol = 82, 85
		beginSlf, endSlf =101,106
		#print 'assuming GEISA 03 format!'
	elif '08' in fileName or '09' in fileName:
		endFreq=12                 # A
		beginStr, endStr = 12, 23  # B
		beginAir, endAir = 23, 29  # C
		beginEne, endEne = 29, 39  # D
		beginTEx, endTEx =119,123  # F
		beginIso, endIso =123,126  # G
		beginMol, endMol =126,129  # I
		beginSlf, endSlf =145,152  # N
		#print 'assuming GEISA 08 format!'
	elif 'iasi' in fileName.lower():
		endFreq=12
		beginStr, endStr = 12, 23
		beginAir, endAir = 23, 29
		beginEne, endEne = 29, 39
		beginTEx, endTEx = 75, 79
		beginIso, endIso = 79, 82
		beginMol, endMol = 82, 85
		beginSlf, endSlf = 99,105
		#print 'assuming GEISA IASI format!'
	else:
		raise SystemExit, repr(fileName)+'\nGEISA database --- unknown version!?!\n(Cannot find release year and/or "IASI" in filename to identify proper format)'
	return {'iw': beginFreq, 'lw': endFreq, 'iM': beginMol, 'lM': endMol, 'iI': beginIso, 'lI': endIso, 
	        'iS': beginStr,  'lS': endStr,  'iE': beginEne, 'lE': endEne,
		'iA': beginAir, 'lA': endAir, 'isw': beginSlf, 'lsw': endSlf, 'iT': beginTEx, 'lT': endTEx}

##########################################################################################################################################

def extract_range (geisa, xLow, xHigh):
	""" Read all lines up to a upper wavenumber limit from Geisa formatted database. """
	# determine position (first and last indices) of wavenmber, moleculeID, isotopeID
	fields = set_geisa_fields (geisa.name)
	iw, lw, lM = fields['iw'], fields['lw'], fields['lM']
	# proceed to first requested line
	record = bisect_first_line (geisa, xLow, xHigh, iw, lw)
	# initialize list if lines
	lines = []
	# collect lines
	while record:
		wvn = float(record[:lw])
		if wvn<=xHigh:  lines.append(record)
		else:           break
		# read next record
		record = geisa.readline()

	if len(lines)>0:  print '# last  line     accepted \n', lines[-1][:lM]
	if record:        print '# first line not accepted \n', record[:lM]     # empty string returned at end-of-file

	return lines

##########################################################################################################################################

def extract_Mol (geisa, xLow, xHigh, getMol):
	""" Read lines of a given molecule up to a upper wavenumber limit from Geisa formatted database. """
	# determine position (first and last indices) of wavenmber, moleculeID, isotopeID
	fields = set_geisa_fields (geisa.name)
	iw, lw, iM, lM = fields['iw'], fields['lw'], fields['iM'], fields['lM']
	# proceed to first requested line
	record = bisect_first_line (geisa, xLow, xHigh, iw, lw)
	# initialize list if lines
	lines = []
	# collect lines
	while record:
		wvn = float(record[:lw])
		if wvn>xHigh: break
		mol = int(record[iM:lM])
		if mol==getMol: lines.append(record)
		# read next record
		record = geisa.readline()

	if len(lines)>0:  print '# last  line     accepted \n', lines[-1][:lM]
	if record:        print '# first line not accepted \n', record[:lM]     # empty string returned at end-of-file

	return lines

##########################################################################################################################################

def extract_MolIso (geisa, xLow, xHigh, getMol, getIso):
	""" Read lines of a given molecule/isotope up to a upper wavenumber limit from Geisa formatted database. """
	# determine position (first and last indices) of wavenmber, moleculeID, isotopeID
	fields = set_geisa_fields (geisa.name)
	iw, lw, iM, lM, iI, lI = fields['iw'], fields['lw'], fields['iM'], fields['lM'], fields['iI'], fields['lI']
	# proceed to first requested line
	record = bisect_first_line (geisa, xLow, xHigh, iw, lw)
	# initialize list if lines
	lines = []
	# collect lines
	while record:
		wvn = float(record[:lw])
		if wvn>xHigh: break
		mol = int(record[iM:lM])
		iso = int(record[iI:lI])
		if mol==getMol and iso==getIso: lines.append(record)
		# read next record
		record = geisa.readline()

	if len(lines)>0:  print '# last  line     accepted \n', lines[-1][:lM]
	if record:        print '# first line not accepted \n', record[:lM]     # empty string returned at end-of-file

	return lines

##########################################################################################################################################

def extract_MolStr (geisa, xLow, xHigh, getMol, strMin):
	""" Read strong lines of a given molecule up to a upper wavenumber limit from Geisa formatted database. """
	# determine position (first and last indices) of wavenmber, moleculeID, isotopeID
	fields = set_geisa_fields (geisa.name)
	iw, lw, iM, lM, iS, lS = fields['iw'], fields['lw'], fields['iM'], fields['lM'], fields['iS'], fields['lS']
	# proceed to first requested line
	record = bisect_first_line (geisa, xLow, xHigh, iw, lw)
	# initialize list if lines
	lines = []
	# collect lines
	while record:
		wvn = float(record[:lw])
		if wvn>xHigh: break
		mol = int(record[iM:lM])
		str = float(replace(record[iS:lS],'D','e'))
		if mol==getMol and str>=strMin: lines.append(record)
		# read next record
		record = geisa.readline()

	if len(lines)>0:  print '# last  line     accepted \n', lines[-1][:lM]
	if record:        print '# first line not accepted \n', record[:lM]     # empty string returned at end-of-file

	return lines

##########################################################################################################################################

def extract_MolIsoStr (geisa, xLow, xHigh, getMol, getIso, strMin):
	""" Read strong lines of a given molecule/isotope up to a upper wavenumber limit from Geisa formatted database. """
	# determine position (first and last indices) of wavenmber, moleculeID, isotopeID
	fields = set_geisa_fields (geisa.name)
	iw, lw, iM, lM, iI, lI, iS, lS = fields['iw'], fields['lw'], fields['iM'], fields['lM'], fields['iI'], fields['lI'], fields['iS'], fields['lS']
	# proceed to first requested line
	record = bisect_first_line (geisa, xLow, xHigh, iw,lw)
	# initialize list if lines
	lines = []
	# collect lines
	while record:
		wvn = float(record[:lw])
		if wvn>xHigh: break
		mol = int(record[iM:lM])
		iso = int(record[iI:lI])
		str = float(replace(record[iS:lS],'D','e'))
		if mol==getMol and iso==getIso and str>=strMin: lines.append(record)
		# read next record
		record = geisa.readline()

	if len(lines)>0:  print '# last  line     accepted \n', lines[-1][:lM]
	if record:        print '# first line not accepted \n', record[:lM]     # empty string returned at end-of-file

	return lines

##########################################################################################################################################

if (__name__ == "__main__"):
	from IO import open_outFile

	from command_parser import parse_command, standardOptions
        opts = standardOptions + [  # h=help, c=commentChar, o=outFile
	       {'ID': 'i', 'name': 'isoNr',  'type': IntType, 'default': 0},
               {'ID': 'm', 'name': 'molNr', 'type': IntType, 'default': 0},
               {'ID': 'S', 'name': 'strMin', 'type': FloatType,   'default': 0.0, 'constraint': 'strMin>0'},
	       {'ID': 'x', 'name': 'xLimits', 'type': Interval, 'default': Interval(0.0,99999.9), 'constraint': 'xLimits.lower>=0.0'}
               ]

	files, options, commentChar, outFile = parse_command (opts, 1)

	for opt in opts:
		if opt.has_key('name') and opt.has_key('type'): exec opt['name'] + ' = ' + repr(options.get(opt['name']))

    	if options.has_key('h'):
    		print __doc__%globals();  raise SystemExit, " end of geisa help"
	elif molNr or xLimits.size()<50000.0:
		# Read lines from GEISA line parameter file
		lines = extract_geisa (files[0], xLimits, molNr, isoNr, strMin)
	else:
		# at least molecule or wavenumber range needed
		raise SystemExit, ' ERROR: neither molecule nor wavenumber range specified!'

	# open an output file if explicitely specified (otherwise use standard output)
	out = open_outFile (outFile, commentChar='#')
	# print selected lines
	for line in lines: out.write (line)
	# close the output file (if its not stdout)
	if outFile: out.close()
