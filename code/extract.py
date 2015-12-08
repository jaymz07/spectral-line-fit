#!/usr/bin/env python

"""
  extract

  extract (grep/select) line parameters from spectroscopic data base file

  usage:
  extract  [options]  line_parameter_database

  command line options:
    -h            help 
   --help         help extended
    -c   char(s)  comment character used in output file (default #)
    -o   file     output file (default: standard output, see last note of extended help)
    
    -d   string   output directory (default: save file(s) in current directory)
    -f   string   output format:  original or simple lists (position vs strengths etc)
                                  default "vSEan" (see notes of extended help)
    -i   integer  isotope name (e.g. 4 or 162 for heavy water HDO, the fourth most abundant isotope)
    -m   string   molecule name  (one molecule only!)
    -S   float    minimum line strength to accept
    -x   interval lower and upper end of spectral (default: wavenumber) range (comma separated pair without blanks)
    -X   string   unit used for setting spectral range with -x option (does not affect output listing!)
                  (default: "cm-1",  other choices: "Hz", "MHz", "GHz", "mue", "nm")

  For more information use
  extract --help
"""

more_help = \
"""

  OUTPUT FORMAT:
    *  "simple lists" --- for each molecule a (tabular ascii) file with columns for position etc is generated
                          (where the filename is set automatically by molecule name with the format indicated by the extension)

       The actual format is defined by a combination of single letters
       "v" --- wavenumber/frequency position   (Hint: the letter "v" looks like the greek nu)
       "S" --- Strengths
       "E" --- Energy (lower state)
       "a" --- air broadening half width (pressure, collision broadening)
       "s" --- self broadening half width
       "n" --- temperature exponent n of pressure broadening
       "i" --- isotope number
       "b" --- broadening parameters, equivalent to "asni"
       
       use "vSEan" or "vSEasni"="vSEb" to produce a line list acceptable by lbl2xs and lbl2od
       use "o" or "h" or "g" to save the extract in the original format

       Valid formats are:   "o", "g", "h",  "vS","vSE","vSEa","vSEb","vSEan","vSEasn","vSEasni"


  NOTES:
    *  at least wavenumber range OR molecule has to be selected!
    *  molecule names are case sensitive!
    *  extracting lines of some selected molecules simultaneously is (currently) not supported:
       either specify one molecule or none
    *  the database filename must include either the string "hit" or "geisa" or "sao" (case insensitive)
       in order to give extract a chance to read with the proper format!
    *  currently only HITRAN, GEISA or SAO are supported
    *  format conversion hitran <---> geisa not yet implemented
    *  if lines for all gases in a spectral range are to be selected and if the output is to be written
       to a simple line list (format 'vS' etc), separate files are produced for each molecule individually
       (an output file specified with the -o option then lists ALL selected lines in the original format)

"""

####################################################################################################################################
#####  ToDo:                                                                                                                   #####
#####                                                                                                                          #####
#####  allow isotope aliases like HDO                                                                                          #####
#####                                                                                                                          #####
####################################################################################################################################

from string import *
from types import *
import os
import sys
import os.path

from ir import c, h, k
from cgsUnits import change_frequency_units 
from molecules import molecules, get_mol_id_nr, isotope_id
from pairTypes import Interval

####################################################################################################################################

def print_lines (outFile, Lines, dataFile, molecule, format='', commentChar='#'):
    	if   'hit'   in lower(dataFile):    xUnit, TRef, pRef = 'cm-1', 296.0, 1013.25
    	elif 'sao'   in lower(dataFile):    xUnit, TRef, pRef = 'cm-1', 296.0, 1013.25
    	elif 'geisa' in lower(dataFile):    xUnit, TRef, pRef = 'cm-1', 296.0, 1013.25
    	elif 'jpl'   in lower(dataFile):    xUnit, TRef, pRef  = 'MHz', 300.0, 1013.0
    	else:                               raise SystemExit, 'unknown database type!'
	
	# check if output file extension indicates the format
	if outFile:
		outFileExt = os.path.splitext(outFile)[1]
		if not format and lower(outFileExt).startswith('.vs'):
			format = outFileExt[1:]
			print 'format automatically determined from output file extension ', format

	if lower(format).startswith('vs'):
		if molecule: 
			write_lines_xy (outFile, Lines, dataFile, format, TRef, pRef, molecule, commentChar)
		else:
			if outFile:  print '\nWARNING:  specified output filename will be ignored!\n          (lines will be saved to files with names defined by molecule and format)\n'
			lineLists = split_molecules (Lines, dataFile) # actually this returns a dictionary!
			for molec,lines in lineLists.items():
				outFile = molec + os.path.extsep + format
				write_lines_xy (outFile, lines, dataFile, format, TRef, pRef, molec, commentChar)
	else:
		# save line extract in original format
		commentChar = '00'
		if   'hit' in dataFile:
			if lower(format).startswith('g'):  raise SystemExit, '\nERROR --- extract:  no format conversion hitran -> geisa'
		elif 'geisa' in dataFile:
			commentChar= 5*commentChar
			if lower(format).startswith('h'):  raise SystemExit, '\nERROR --- extract:  no format conversion geisa -> hitran'
		else:
			if lower(format).startswith('g') or lower(format).startswith('h'):
				raise SystemExit, '\nERROR --- extract:  no format conversion ??? -> hitran'

		out = open_outFile (outFile, commentChar)
		if 'hit' in dataFile.lower():
			if len(Lines[0].rstrip())==100:  # old hitran versions <=2000
				out.write ('%3s%12s%10s%10s%5s%5s%10s%4s%8s%3s%3s%9s%9s%3s%6s\n' % \
		    		          ('000','wavenumber','S','A','air','self','Energy','n','pShift','uV','lV','uL','lL','er','ref'))
			else:  # new hitran versions >=2004
				out.write ('%3s%12s%10s%10s%5s%5s%10s%4s%8s%15s%15s%15s%15s%6s%12s %s7%7s\n' % \
		    		          ('000','wavenumber','S','A','air','self','Energy','n','pShift','upVib','loVib','upLocal','loLocal','err','ref','usw','lsw'))
		for line in Lines: out.write (line)
		if outFile: out.close()

    	return

####################################################################################################################################

def write_lines_xy (outFile, lines, dataFile, job='vS', TRef=0, pRef=0, molecule='', commentChar='#'):
	""" Print 'core' line parameters, i.e., positions vs strengths, and optionally energies, airWidths, tempExponents. """
	# open an output file if explicitely specified (otherwise use standard output)
	out = open_outFile (outFile, commentChar)
	if molecule: out.write ('%s %s %s\n' % (commentChar, "molecule:", molecule))
	# reference pressure and temperature of database
	if TRef: out.write ('%s %s %8.2f %s\n' % (commentChar, "temperature:", TRef, "K"))
	if pRef: out.write ('%s %s %8.2f %s\n' % (commentChar, "pressure:   ", pRef, "mb"))
	# extract most important numeric line parameters
	positions, strengths, energies, airWidths, selfWidths, tempDep, isoNr = core_parameters (lines, dataFile)

	# some statistical information
	nLines = len(positions)
	out.write ('%s %-30s %12i\n' % (commentChar, "number of lines:     ", nLines))
	out.write ('%s %-30s %12.3g%13.3g\n' % (commentChar, "min, max line strength:     ", min(strengths),  max(strengths)))
	out.write ('%s %s %s\n' % (commentChar, "format:", job))
	if 'a' in job:
		out.write ('%s %-30s %12.3f%13.3f\n' % (commentChar, "min, max airbroad.  widths: ", min(airWidths),  max(airWidths)))
		out.write ('%s %-30s %12.3f%13.3f\n' % (commentChar, "min, max selfbroad. widths: ", min(selfWidths), max(selfWidths)))
		out.write ('%s %-30s %12.3f%13.3f\n' % (commentChar, "min, max temp. exponent:    ", min(tempDep),    max(tempDep)))
	if outFile:
		print '%-8s %8i %s %11.3g %s %-11.3g %s %s' % (molecule, nLines, ' lines with', min(strengths),  '< S <', max(strengths), ' writing to file', outFile)
	#
	if job=="vSEani":
		format = '%12f %11.3e %11.5f %8.5f %7.4f %3i\n'
		out.write ('%1s %10s %11s %11s %8s %7s %3s\n' % (commentChar,'position','strength', 'energy', 'airWidth', 'Tdep', 'iso'))
		out.write ('%1s %10s %11s %11s %8s %7s\n' % (commentChar,'cm-1','cm-1/cm-2', 'cm-1', 'cm-1', ''))
		for l in xrange(nLines):  out.write ( format % (positions[l],strengths[l],energies[l],airWidths[l],tempDep[l],isoNumber[l]))
	elif job=="vSEan":
		format = '%12f %11.3e %11.5f %8.5f %7.4f\n'
		out.write ('%1s %10s %11s %11s %8s %7s\n' % (commentChar,'position','strength', 'energy', 'airWidth', 'Tdep'))
		out.write ('%1s %10s %11s %11s %8s %7s\n' % (commentChar,'cm-1','cm-1/cm-2', 'cm-1', 'cm-1', ''))
		for l in xrange(nLines):  out.write ( format % (positions[l],strengths[l],energies[l],airWidths[l],tempDep[l]))
	elif job=="vSEasni":
		format = '%12f %11.3e %11.5f %8.5f %8.5f %7.4f %4i\n'
		out.write ('%1s %10s %11s %11s %8s %8s %7s %4s\n' % (commentChar,'position','strength', 'energy', 'airWidth', 'selfWidth', 'Tdep', 'iso'))
		out.write ('%1s %10s %11s %11s %8s %8s %7s\n' % (commentChar,'cm-1','cm-1/cm-2', 'cm-1', 'cm-1',  'cm-1', ''))
		for l in xrange(nLines):
			out.write ( format % (positions[l],strengths[l],energies[l],airWidths[l],selfWidths[l],tempDep[l], isoNr[l]))
	elif job=="vSEasn":
		format = '%12f %11.3e %11.5f %8.5f %8.5f %7.4f\n'
		out.write ('%1s %10s %11s %11s %8s %8s %7s\n' % (commentChar,'position','strength', 'energy', 'airWidth', 'selfWidth', 'Tdep'))
		out.write ('%1s %10s %11s %11s %8s %8s %7s\n' % (commentChar,'cm-1','cm-1/cm-2', 'cm-1', 'cm-1',  'cm-1', ''))
		for l in xrange(nLines):  out.write ( format % (positions[l],strengths[l],energies[l],airWidths[l],selfWidths[l],tempDep[l]))
	elif job=="vSEa":
		format = '%12f %11.3e %11.5f %8.5f\n'
		for l in xrange(nLines):  out.write ( format % (positions[l],strengths[l],energies[l],airWidths[l]) )
	elif job=="vSE":
		format = '%12f %11.3e %11.5f\n'
		for l in xrange(nLines):  out.write ( format % (positions[l],strengths[l],energies[l]) )
	else:
		format = '%12f %11.3e\n'
		for l in xrange(nLines):  out.write ( format % (positions[l],strengths[l]) )
	# close the output file (if its not stdout)
	if outFile: out.close()
    	return

####################################################################################################################################

def check_database_file (files):
	""" Check command line arguments supplied to extract. """
	if len(files)==0:
		return 'ERROR: no line parameter file given!'
	elif len(files)==1:
		if os.path.isfile(files[0]):
			lineDataFile = files[0]
		else:
			return 'ERROR: ' + `files[0]` + ' --- line parameter file invalid, nonexisting, ... ?'
	else:
		return 'ERROR: too many arguments (need just one line parameter database file)!'

	# parse line database filename to identify type
	count_hit   = int('HIT'   in upper(lineDataFile))
	count_geisa = int('GEISA' in upper(lineDataFile))
	count_sao   = int('SAO'   in upper(lineDataFile))
	count_jpl   = int('JPL'   in upper(lineDataFile))
	countAll    = count_hit+count_geisa+count_sao+count_jpl
	#
	if   countAll>1:  return 'ERROR: line parameter database filename is ambiguous!\n(filename should include either "HIT" or "GEISA" or "JPL" (case insensitive))'
	elif countAll<1:  return 'ERROR: type of line parameter database invalid or unknown!\n(either "HITRAN" or "GEISA")'
	else:             return


####################################################################################################################################

def extract_lines (files, xLimits=None, molecule=None, isotope=None, strMin=0.0, xUnit='cm-1'):

	errMsg = check_database_file (files)
	if errMsg: raise SystemExit, errMsg
	file = files[0]

	if not (molecule or xLimits):
		raise SystemExit, ' ERROR: neither molecule nor wavenumber range specified!'
	else:
		if isotope and not molecule:
 		    if isotope<10: raise SystemExit, '%s %i%s' % ('selecting', isotope,'. abundant isotope of all molecules not implemented!')
		    else:          raise SystemExit, 'ERROR: searching isotopes requires specification of molecule!'

	if xLimits: 
		if xUnit!='cm-1': 
			print xLimits, xUnit, '-->',
			xLimits = change_frequency_units (xLimits, xUnit, 'cm-1')
			print xLimits
	else:
		xLimits = Interval(0.0,999999.9)

	# read all lines in interval and return dictionary
	if   count(upper(file),'HIT'):
		if molecule:
			try:             MolNr = molecules[molecule]['hitran']
			except KeyError: raise SystemExit, 'ERROR: invalid/unknown molecule ' + repr(Mol)
			IsoNr = isotope_id (molecule, isotope, 'hitran')
		else:   MolNr = 0; IsoNr=0
		from hitran import extract_hitran
		lines = extract_hitran (file, xLimits, MolNr, IsoNr, strMin)
	elif count(upper(file),'GEISA'):
		if molecule:
			try:             MolNr = molecules[molecule]['geisa']
			except KeyError: raise SystemExit, 'ERROR: invalid/unknown molecule ' + repr(molecule)
			if isotope:
				if isotope<10:
					print str(isotope) + '. most abundant isotope',
					try:               isotope=int(molecules[molecule]['isotopes'][isotope-1])
					except IndexError: raise SystemExit, 'ERROR: invalid/unknown isotope ' + repr(isotope)
					else:              print '-->', isotope
		else:
			MolNr = 0; IsoNr=0
		from geisa import extract_geisa
		lines = extract_geisa (file, xLimits, MolNr, isotope, strMin)
	elif   count(upper(file),'SAO'):
		if molecule:
			try:             MolNr = molecules[molecule]['sao']
			except KeyError: raise SystemExit, 'ERROR: invalid/unknown molecule ' + repr(molecule)
			IsoNr = isotope_id (molecule, isotope, 'hitran')
		else:   MolNr = 0; IsoNr=0
		from sao import extract_sao
		lines = extract_sao (file, xLimits, MolNr, IsoNr, strMin)
	else:
		raise SystemExit, 'sorry, currently only hitran and geisa atlas!'
	print len(lines), ' lines extracted from ', file
	return lines

####################################################################################################################################

def core_parameters (lines, dataFile):
	""" Given a list of data base records (one entry per transition) return python lists of the most important spectrocopic line parameters. """
	# column start/stop for all types of parameters
	if 'hit' in dataFile.lower() or 'sao' in dataFile.lower():
	 	iw,lw, iS,lS, iE, lE, iA,lA, isw, lsw, iT,lT, iI,lI = 3,15, 15,25, 45,55, 35,40, 40,45, 55,59, 2,3
	else:
		from geisa import set_geisa_fields
		fields = set_geisa_fields (dataFile)
		for key,val in fields.items(): exec key + '=' + repr(val)
	# now extract columns, convert to appropriate type (int/float)
	positions = [float(line[iw:lw]) for line in lines]
	strengths = [float(replace(line[iS:lS],'D','e')) for line in lines]
	energies  = [float(line[iE:lE]) for line in lines]
	airWidths = [float(line[iA:lA]) for line in lines]
	tempDeps  = [float(line[iT:lT]) for line in lines]
	isoNr     = [int(line[iI:lI]) for line in lines]
	if 'hit' in dataFile.lower() or 'geisa' in dataFile.lower():
		selfWidths = [float(line[isw:lsw]) for line in lines]
	else:
		selfWidths = len(lines)*[0] # a list with nLines zeros
	return positions, strengths, energies, airWidths, selfWidths, tempDeps, isoNr

####################################################################################################################################

def split_molecules (lines, dataFile):
	""" Given the list of database records extracted from Hitran/Geisa, distribute the entries in separate lists for each molecule. """
	# initialize dictionary
	lineLists = {}
	if 'hit' in dataFile.lower():
		mol_id = get_mol_id_nr (molecules, 'hitran') # translation dictionary: hitran molecular ID numbers --> names
		im, lm =  0, 2                               # set indices for molecular id
	elif 'geisa' in dataFile.lower():
		from geisa import set_geisa_fields
		fields = set_geisa_fields (dataFile)
		im, lm = fields['iM'], fields['lM']
		mol_id = get_mol_id_nr (molecules, 'geisa') # translation dictionary: geisa molecular ID numbers --> names
	# split list of lines
	for line in lines:
		molNr = int(line[im:lm])
		molec = mol_id.get(molNr)
		if molec in lineLists: lineLists[molec].append(line)
		elif molec>0:          lineLists[molec]   = [line]
		else:                  raise SystemExit, 'ERROR --- extract:  unknown molecule with ID number ' + `molNr`
	return lineLists

####################################################################################################################################

if (__name__ == "__main__"):
 	
	from IO import open_outFile
	import re
	from command_parser import parse_command, standardOptions
        opts = standardOptions + [  # h=help, c=commentChar, o=outFile
	       {'ID': 'help'},
	       {'ID': 'f', 'name': 'format', 'type': StringType, 'default': 'vSEan',
	                   'constraint': 'lower(format) in ["o", "g", "h", "vs","vse","vsea","vseb","vsean","vseasn","vseasni"]'},
               {'ID': 'd', 'name': 'outDir', 'type': StringType},
               {'ID': 'i', 'name': 'isotope', 'type': IntType, 'constraint': 'isotope>=0', 'default': 0},
               {'ID': 'm', 'name': 'molecule', 'type': StringType, 'constraint': 'len(re.split("[,;\s]",molecule))==1'},
               {'ID': 'S', 'name': 'strMin', 'type': FloatType, 'constraint': 'strMin>=0.0', 'default': 0.0},
	       {'ID': 'x', 'name': 'xLimits', 'type': Interval, 'constraint': 'xLimits.lower>=0.0'},
	       {'ID': 'X', 'name': 'xUnit', 'type': StringType, 'default': 'cm-1', 'constraint': "xUnit in ['cm-1', 'mue', 'nm', 'Hz', 'kHz', 'MHz', 'GHz', 'THz']"}
               ]

	files, options, commentChar, outFile = parse_command (opts, 1)

	if options.has_key('h'):  print __doc__;  raise SystemExit
	if options.has_key('help'):  print __doc__[:-42] + more_help;  raise SystemExit, " end of extract help"

	outDir = options.pop('outDir',None)
	format = options.pop('format',None)
	# replace format shortcut with 'full' name
	if format=="vSEb":  format='vSEasni'

	lines = extract_lines (files, **options)

	if len(lines)>0:
		if outDir:
			# if not outFile: raise SystemExit, 'directory specified, but no output file name!'
			if not os.path.isdir(outDir): os.mkdir(outDir); print 'mkdir ', outDir
			os.chdir(outDir); print '\nchanging directory -->', outDir, '\n'
		# print extracted lines
		print_lines (outFile, lines, files[0], options.get('molecule'), format, commentChar)
	else:
		raise SystemExit, 'no lines'
