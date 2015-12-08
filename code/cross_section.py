#!/usr/bin/env python

""" cross_section

  Read and write / plot molecular cross sections (e.g., to reformat or interpolate).

  usage:
  cross_section [options] xs_file(s)

  command line options:
   -c   char(s)  comment character in input file(s) (default #)
   -C   ints     sequence of integers: cross section (levels/layers) to extract
                 (in the xy ascii file 0 corresponds to the first column=wavenumber)
   -f   string   format for output file [default: 'xy', 'h' for hitran format]
   -h            help
   -i   string   interpolation method for spectral domain (cross section vs wavenumber)
                 "2", "3", "4" for Lagrange interpolation, "s" for spline
                 default: '3' three-point Lagrange
		 '0' in combination with 'xy' tabular output format generates individual files for each p, T, molecule
  --plot         matplotlib for quicklook of cross sections
   -o   file     output file (default: standard output)

  Cross Section Files:
  *   xy formatted ascii file with wavenumbers in column 1 and cross section(s) (for some p,T) in the following column(s).
  *   hitran formatted cross section file
  *   pickled cross section file (default output of lbl2xs)

  NOTE:  
  *  If no output file is specified, only a summary 'statistics' is given !!!
  *  xy tabular output format:
     For each molecule cross sections of all p,T pairs will be interpolated to a common wavenumber grid and saved as 'matrix';
     To write an individual file for each p/T, suppress interpolation with '-i0' option.
  *  When reading files,  cross_section tries to determine the format from the first line (record)
  *  When reading Hitran xs files, the header record is sometimes incorrectly formatted (missing blanks), and cross_section might fail.
     (Reading files with xs data distributed over several records is not yet implemented!).
"""

import os
import sys
from string import *
from types import *
import time

import numpy as np

from cgsUnits import unitConversion, cgsTemperature
from pairTypes import Interval
from IO import parse_file_header, readDataAndComments, open_outFile
from lagrange_interpolation import *

####################################################################################################################################

def write_crossSections (crossSections, outFile=None, format='', linePositions=None, interpolate='3', commentChar='#'):
	""" Save cross sections of one molecule and several p,T levels on file. """
	if format=='xy':
		write_crossSections_xyy (crossSections, outFile, linePositions, interpolate, commentChar)
	elif isinstance(format,str) and format.startswith('h'):
		write_crossSections_hitran (crossSections, outFile)
	else:
		write_crossSections_pickled (crossSections, outFile)

####################################################################################################################################

def write_crossSections_pickled (crossSections, outFile=None):
	""" Write cross sections (typically a list of xs instances for some p,T) to output (file) using Python's pickle module. """
	import pickle
	if outFile: out = open (outFile, 'w')
	else:       raise SystemExit, '\n ERROR --- write_crossSections_pickled:  no cross section pickling for standard out!' + \
	                              '\n                                         (choose different format and/or specify output file(s))'
	# NOTE:  the sequence of loads has to corrrespond to the sequence of dumps!
	pickle.dump(join([os.path.basename(sys.argv[0])] + sys.argv[1:]), out)
	pickle.dump([len(xs['y']) for xs in crossSections], out)
	for xs in crossSections:
		pickle.dump(xs, out)
    	out.close()

####################################################################################################################################

def write_crossSections_hitran (crossSections, outFile=None, lineLength=-1, format=' %12g'):
	""" Write cross section to output (file) in hitran format. """
	if outFile: out = open (outFile, 'w')
	else:       out = sys.stdout
	for xs in crossSections:
		nxs = len(xs['y'])
		xLow, xHigh = crossSections[0]['x'].limits()
		maxXS = max(xs['y'])
		out.write ('%s %f %f %i %g %8.3f %g\n' % (xs['molecule'], xLow, xHigh, nxs, unitConversion(xs['p'],'p',new='mb'), xs['T'], maxXS))
		if lineLength<1:
			# all data in a single long line
    			out.write ( str(nxs*format+'\n') % tuple(xs['y']) )
		elif lineLength>1:
			for i in xrange(0,nxs,lineLength):
    				lo, hi = i, min(i+lineLength,nxs)
    				out.write ( (hi-lo)*format % tuple(xs['y'][i:min(i+lineLength,nxs)]) + '\n' )
		else:
			# just one value per line
			frmt = format+'\n'
			for i in xrange(0,nxs):  out.write ( frmt % xs['y'][i])
    	if outFile: out.close()

####################################################################################################################################

def write_crossSections_xy (crossSections, outFile=None, linePositions=None, commentChar='#'):
	""" Write cross section(s) to output (file) in tabular xy format.
	    (Cross sections of a given molecule (i.e., data for all p,T levels) are written to a individual files) """

	# extract attributes like pressure and temperatures
	nXS = len(crossSections)
	molecules    = tuple([xs.get('molecule') for xs in crossSections])
	pressures    = tuple([unitConversion(xs.get('p',0.0),'p',new='mb') for xs in crossSections])
	temperatures = tuple([xs.get('T',0.0) for xs in crossSections])
	lineShapes   = tuple([xs.get('lineShape') for xs in crossSections])
	if isinstance(linePositions,np.ndarray):
		lineInfo =  '%i %s %f ... %f' % (len(linePositions), 'lines in ', min(linePositions),max(linePositions))
	else:
		lineInfo = ''

	for l, xs in enumerate(crossSections):
		if isinstance(outFile,str) and '.' in outFile:
			outRoot, outExt = os.path.splitext(outFile)
			oFile = '%s_%2.2i_%gmb_%.1fK%s' % (outRoot, l+1, pressures[l], temperatures[l], outExt)
		else:   
			oFile = '%s_%gmb_%.1fK.%s' % (molecules[l], pressures[l], temperatures[l], 'xs')
		out = open_outFile (oFile, commentChar)
		# print a summary
		out.write ( '%1s  %-16s %s\n' % (commentChar, 'molecule:', molecules[l]))
		out.write ( '%1s  %-16s %s\n' % (commentChar, 'lineshape:', lineShapes[l]))
		if lineInfo:  out.write ( '%1s  %s\n' % (commentChar, lineInfo))
		out.write ( commentChar+'  pressure [mb]:  ' + ' %10g' % pressures[l] + '\n')
		out.write ( commentChar+'  temperature [K]:' + ' %10.2f' % temperatures[l] + '\n')
		out.write ( commentChar + '\n' )

		# extract cross section spectra and save in list
		yyy  = xs['y']
		nxy  = len(yyy)-1
		xLow, xHigh = crossSections[0]['x'].limits()
		xLo,  xHi   = float(xLow/nxy), float(xHigh/nxy)
		# and write cross section vs wavenumber
		out.write ( '%1s%12s %15s\n' % (commentChar, 'wavenumber', 'cross section') )
		format = ' %12f %12.6g\n'
		for i,y in enumerate(yyy): 
			x = i*xHi + (nxy-i)*xLo
			out.write (format % (x,y))
		
		if outFile: out.close()

####################################################################################################################################

def write_crossSections_xyy (crossSections, outFile=None, linePositions=None, interpolate='3', commentChar='#'):
	""" Write cross section(s) to output (file) in tabular xy format.
	    (All cross sections of a given molecule (i.e., data for all p,T levels) 
	     are interpolated to a common wavenumber grid and written to a single file) """

	out = open_outFile (outFile, commentChar)
	if outFile:
		# print a summary
		molecules = tuple([xs.get('molecule') for xs in crossSections])
		sameMolec = np.alltrue([molecules[0] == other for other in molecules[1:]])
		if sameMolec:      out.write ( '%1s  %-16s %s\n' % (commentChar, 'molecule:', molecules[0]))

		lineShapes = tuple([xs.get('lineShape') for xs in crossSections])
		sameShapes = np.alltrue([lineShapes[0] == other for other in lineShapes[1:]])
		if sameShapes:     out.write ( '%1s  %-16s %s\n' % (commentChar, 'lineshape:', lineShapes[0]))

		if isinstance(linePositions,np.ndarray):
			out.write ( '%1s  %i %s %f ... %f\n' % (commentChar, len(linePositions), 'lines in ', min(linePositions),max(linePositions)) )
		# extract pressure and temperatures
		nXS = len(crossSections)
		pressures    = tuple([unitConversion(xs.get('p',0.0),'p',new='mb') for xs in crossSections])
		temperatures = tuple([xs.get('T',0.0) for xs in crossSections])
		out.write ( commentChar+'  pressure [mb]:  ' + nXS*' %10g' % pressures + '\n')
		out.write ( commentChar+'  temperature [K]:' + nXS*' %10.2f' % temperatures + '\n')
		out.write ( commentChar + '\n' )

	# select an interpolation method
	if   interpolate.lower() in 'sb':
                try:                      from scipy.interpolate import splrep, splev
	        except ImportError, msg:  raise SystemExit, msg
                else:                     intMethod =  'Spline interpolation (scipy splrep/splev)'
	elif   interpolate.lower()=='k':
                try:                      from scipy.interpolate import krogh_interpolate as kInt
	        except ImportError, msg:  raise SystemExit, msg
                else:                     intMethod = 'Krogh interpolation (scipy)';  print 'WARNING --- use with care!?! ', intMethod
	elif interpolate.lower() in '234lqc':
                interpolate = interpolate.translate(maketrans('lqcLQC','234234')) # first replace letters by th numeric order
		intMethod = '%s %s %s' % ('Lagrange', interpolate, 'interpolation')
		if   interpolate=='2': interpolation = lagrange2_regularGrid 
		elif interpolate=='3': interpolation = lagrange3_regularGrid 
		elif interpolate=='4': interpolation = lagrange4_regularGrid 
	else:                  raise SystemExit, 'invalid interpolation scheme' + repr(interpolate)

	# extract cross section spectra and save in list
	yyy  = [xs['y'] for xs in crossSections]
	lenY = [len(y) for y in yyy]
	# select densest grid as final for output
	xLow, xHigh = crossSections[0]['x'].limits()
	dx  = (xHigh-xLow)/(max(lenY)-1)
	x   = np.arange(xLow,xHigh+0.1*dx,dx)

	# interpolate 
	tStart = time.clock()
	if interpolate in '234':
		for j,y in enumerate(yyy):  yyy[j] = interpolation (y, len(x))
	elif interpolate in 'kK': 
		for j,y in enumerate(yyy):  yyy[j] = kInt (np.linspace(xLow, xHigh, len(y)), y, x)
	else:
		for j,y in enumerate(yyy):
			xOld = np.linspace (xLow, xHigh, len(y))
			tck  = splrep (xOld, y)
			yyy[j] = splev(x,tck)
	tStop  = time.clock();  print  ' ', tStop-tStart, 'sec ', intMethod

	# transform into numerical array
	yyy = np.transpose(np.array(yyy))
	# and write cross sections vs wavenumber
	out.write ( '%1s%12s %15s\n' % (commentChar, 'wavenumber', 'cross section') )
	xFormat = ' %12f '
	yFormat = yyy.shape[1]*'%12.5g' + '\n'
	for i in xrange(len(x)):
		out.write (xFormat % x[i])
		out.write (yFormat % tuple(yyy[i,:]))
	
    	if outFile: out.close()

####################################################################################################################################

def read_cross_sections_pickled (file):
	""" Read one cross section file (python pickle format format; one molecule, several p,T pairs). """
	try:    f = open(file)
	except: raise SystemExit, 'ERROR:  opening pickled cross section file ' + repr(file) + ' failed (check existance!?!)'
	# initialize list of cross sections to be returned
	xsList = []
	# NOTE:  the sequence of loads has to corrrespond to the sequence of dumps!
	import pickle
	info  = pickle.load(f); print file, info
	lenXS = pickle.load(f); print file, lenXS,
	while 1:
		try: xsList.append(pickle.load(f))
		except EOFError:
			f.close(); print len(xsList); break
	return xsList

####################################################################################################################################

def read_cross_sections_hitran (file):
	""" Read one cross section file (hitran format; one molecule, several p,T pairs). """
	# initialize list of cross sections to be returned
	xsList = []
	# read entire cross section file (incl. commented header lines)
	try:    data = open(file).readlines()
	except: raise SystemExit, 'ERROR:  opening hitran cross section file ' + repr(file) + ' failed (check existance!?!)'
	else:   print len(data), ' records in file, probably for', len(data)/2, ' header and data pairs'
		
	# parse comment header and extract some infos
	npT = len(data)/2
	for i in range(0,npT):
		header      = split(data[2*i])
		try:
			molecule    = header[0]
			xLow, xHigh = map(float,header[1:3])
			nxy         = int(header[3])
			temp, press = map(float,header[4:6])
			restOfLine  = header[6]
		except ValueError, msg:
			raise SystemExit, '%s\n%s\n%s' % (msg, 'reading hitran formatted xs file failed\nproblems with header line', repr(data[2*i]))
		else:
			print molecule, xLow, xHigh, nxy, temp, press, restOfLine, len(split(data[2*i+1]))
		y           = np.array(map(float,split(data[2*i+1])))
		print y.shape, min(y), max(y), '\n'
		xs = {'molecule': molecule, 'p': press, 'T': temp, 'nxy': nxy, 'x': Interval(xLow, xHigh), 'y': y}
		xsList.append(xs)
	return xsList

####################################################################################################################################

def read_cross_sections_xy (file, commentChar):
	""" Read one cross section file (xy ascii format; typically one molecule, several p,T pairs). """
	# initialize list of cross sections to be returned
	xsList = []
	# read entire cross section file (incl. commented header lines)
	try:
		xyData, comments = readDataAndComments(file,commentChar)
	except:
		raise SystemExit, 'ERROR:  reading cross section file failed (check format etc!?!)\n' + repr(file)

	# parse comment header and extract some infos
	need = ('pressure', 'temperature', 'molecule')
	headerDict = parse_file_header (comments, need)

	npT = xyData.shape[1]-1
	xLimits = Interval(xyData[0,0], xyData[-1,0])
	for l in range(1,npT+1):
		xsList.append({'molecule': headerDict['molecule'],  'x': xLimits,  'y': xyData[:,l]})

	pValues, pUnit = headerDict['pressure'][0], headerDict['pressure'][1]
	pValues = unitConversion (pValues, 'pressure', pUnit)
	if len(pValues)==npT:
		for l in range(npT): xsList[l]['p'] = pValues[l]

	TValues, TUnit = headerDict['temperature'][0], headerDict['temperature'][1]
	TValues = unitConversion (TValues, 'temperature', TUnit)
	if len(TValues)==npT:
		for l in range(npT): xsList[l]['T'] = TValues[l]

	return xsList

####################################################################################################################################

def read_cross_sections (file, commentChar='#', format=''):
	""" Read cross section data from a file;  return a list of dictionaries, each with entries for molecule, pressure, temperature, (wavenumber) interval, and data. """
	# enforce to read data with a certain file format  (well, currently format is not used in the call!!!)
	if   format.startswith('h'):  firstLine = 'ABC'
	elif format.startswith('xy'): firstLine = commentChar
	else:                         firstLine = strip(open(file).readline())

	# try to determine filetype automatically from first nonblank character in file
	if firstLine[0] in uppercase and firstLine[1] in letters+digits: # probably name of a molecule
		crossSections = read_cross_sections_hitran (file)
	elif firstLine.startswith(commentChar):
		crossSections = read_cross_sections_xy (file, commentChar)
	else:
		crossSections = read_cross_sections_pickled (file)

	return crossSections

####################################################################################################################################

def pyplot_xs (crossSections):
	
	if isinstance(crossSections,dict):
		figure()
		xs  = crossSections
		xGrid=np.linspace(xs['x'].lower,xs['x'].upper,len(xs['y']))
		semilogy (xGrid, xs['y'], label='%10.3gmb %7.2fK' % (unitConversion(xs['p'],'p',new='mb'), xs['T']))
		#show()
	elif isinstance(crossSections,list):
		figure()
		for xs in crossSections:
			xGrid=np.linspace(xs['x'].lower,xs['x'].upper,len(xs['y']))
			semilogy (xGrid, xs['y'], label='%10.3gmb %7.2fK' % (unitConversion(xs['p'],'p',new='mb'), xs['T']))
		title('%s  %10.3gmb %7.2fK' % (xs['molecule'], unitConversion(xs['p'],'p','mb'), xs['T']) )
		legend()
		#show()

####################################################################################################################################

def _cross_section_ (iFile, oFile=None, commentChar='#', format=None, interpolate='3', columns=[], plot=False):
	# read cross sections from file
	crossSections = read_cross_sections (iFile, commentChar)
	
	# remove unwanted cross sections
	if columns:
		columns = [int(C)-1 for C in columns.split(",")]
		newCrossSections = [xs for l,xs in enumerate(crossSections) if l in columns]
		crossSections = newCrossSections

	# print summary
	for xs in crossSections:
		nv     = len(xs['y']) - 1 # number of intervals !!!
		deltaX = xs['x'].size() / nv
		iMin   = np.argmin(xs['y'])
		iMax   = np.argmax(xs['y'])
		print '%-10s %10i wavenumber %10.3g %s %10gmb %8.2fK %12.3g < xs < %10.3g  @ %10i %10f   avg %10.3g' % \
		      (xs['molecule'], nv+1, deltaX, xs['x'], unitConversion(xs['p'], 'p', new='mb'), xs['T'],
		       min(xs['y']), max(xs['y']), iMax, xs['x'].lower+iMax*deltaX, np.mean(xs['y']))

	if plot:  pyplot_xs(crossSections)
	
	# save cross sections
	if oFile:
		if format=='xy':
			if interpolate=='0':
				write_crossSections_xy (crossSections, oFile, commentChar=commentChar)
			else:
				write_crossSections_xyy (crossSections, oFile, interpolate=interpolate, commentChar=commentChar)
		elif format.startswith('h'):
			write_crossSections_hitran (crossSections, oFile)
		else:
			write_crossSections_pickled (crossSections, oFile)

####################################################################################################################################

if (__name__ == "__main__"):

	from command_parser import parse_command, standardOptions, multiple_outFiles 
        opts = standardOptions + [  # h=help, c=commentChar, o=outFile
               {'ID': 'f', 'name': 'format', 'default': 'xy'},
	       {'ID': 'C', 'name': 'columns', 'type': StringType, 'constraint': 'all([digit.isdigit() for digit in split(columns,",")])'},
               {'ID': 'i', 'name': 'interpolate', 'type': StringType, 'default': '3', 'constraint': 'len(interpolate)==1 and interpolate.lower() in "0234lqcbhks"'},
               {'ID': 'plot'} ]

	Files, options, commentChar, outFile = parse_command (opts,(1,99))
 	
	if 'h' in options: print __doc__%globals(); raise SystemExit, "End of cross_section help"

	outFiles    = multiple_outFiles (Files, outFile)
	options['plot'] = options.has_key('plot')

	if options['plot']:
		try: from matplotlib.pyplot import figure, semilogy, plot, legend, title, show
		except ImportError: raise SystemExit, 'ERROR --- cross_section:  matplotlib not available, no quicklook!'

	for iFile,oFile in zip(Files,outFiles):
		_cross_section_ (iFile, oFile, commentChar, **options)
	
	if options['plot']: show()
