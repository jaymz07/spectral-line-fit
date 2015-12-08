import os
import sys
from string import lower, strip, split, count, whitespace, replace, find, rfind, upper, join
import re
import numpy as np

####################################################################################################################################

def parse_file_header (header, keywords, sep=':', commentChar='#'):
	""" Scan thru list of records (typically read as file header) and search for certain keys;
	    returns a dictionary with (ideally) len(keywords) key,value pairs
	    where each value is either a string or a pair (tuple) of value(s) and unit. """
	# if there is just one keyword (given as a string) put it into a list nevertheless
	if isinstance(keywords,str): keywords=[keywords]
	# get rid of leading or trailing blanks in keywords
	for keyword in keywords: keyword = strip(keyword)
	# initialize
	dict = {}
	# loop over all file header lines
	for record in header:
		# get rid of comment character(s) and leading blanks
		record = re.sub ('^'+commentChar+'* *','',record)
		if count(record,sep)==0: continue
		key,val = split(record,sep,1)
		#print '>> ', key, ' <<', val, '>>'
		for keyword in keywords:
			if find(upper(key),upper(keyword))>-1:
				# check if there is an unit specification [in square brackets]
				mo = re.search(r'\[.*\]',key)
				if mo:
					# only use first part as key without units
					name   = strip(key[:mo.start()])
					values = split(strip(val))
					unit   = key[mo.start()+1:mo.end()-1]
					dict[name] = (np.array(map(float,values)), unit)
				else:
					dict[strip(key)] = strip(val)
	return dict

####################################################################################################################################

def readFileHeader (file, commentChar='#'):
	""" Read a tabular (xy) formatted ascii file and return list of comments (without commentChar) found in header. """
	# open file, initialize list of comment and read first line
	f = open(file)
	comments = []
	record   = strip(f.readline())
	# loop over file header and move all comment records to a separate list
	while record.startswith(commentChar): 
	        # get rid of comment character(s) and leading blanks
		record = re.sub ('^'+commentChar+'* *','',record)
		if len(record)>0:  comments.append(record)
		record   = strip(f.readline())
	f.close()
	return comments

####################################################################################################################################

def grep_from_header (file, keyword, sep=':', commentChar='#'):
	""" Scan thru list of records (typically read as file header), search for ONE keyword, and return its 'value'.
	    (Equivalent to parse_file_header(readFileHeader(file),keyword)[keyword], but returns only the entry). """
	keyword = strip(keyword)
	try: f = open(file)
	except IOError, errMsg:  raise SystemExit, str(errMsg) + '\ncheck your input file list!'
	while True:
		record   = strip(f.readline())
		if len(record)==0: return
		# get rid of comment character(s) and leading blanks
		record = re.sub ('^'+commentChar+'* *','',record)
		if count(record,sep)==0: continue
		key,val = split(record,sep,1)
		# comparison of given and wanted keyword case insensitive!
		if find(upper(key),upper(keyword))>-1:
			# check if there is an unit specification [in square brackets]
			mo = re.search(r'\[.*\]',key)
			if mo:
				# only use first part as key without units
				name   = strip(key[:mo.start()])
				values = split(strip(val))
				unit   = key[mo.start()+1:mo.end()-1]
				return (np.array(map(float,values)), unit)
			else:
				return strip(val)

####################################################################################################################################

def readDataAndComments (file, commentChar='#', delimiter=None, converters=None, skiprows=0, usecols=None, unpack=False):
	""" Read tabular (xy) formatted ascii file, return data as numpy array and list of comments in file header.
	    (Most options are simply passed thru to numpy.loadtxt) """
	if not os.path.isfile(file):
		raise SystemExit, 'ERROR --- readDataAndComments: input file not existing!?!\n' + repr(file)
	# note different naming convention for comment character
	comments = readFileHeader (file, commentChar)
	try:
		data     = np.loadtxt (file, comments=commentChar, converters=converters, skiprows=skiprows, usecols=usecols, unpack=unpack)
	except ValueError, msg:
		raise SystemExit, 'reading (numeric) data failed\n' + str(msg)
	return data, comments

####################################################################################################################################

def cstack (*arrays):
	""" Shorthand robust version of numpy.column_stack: 'paste' arrays side-by-side. """
	try:
		return np.column_stack(arrays)
	except ValueError, msg:
		print 'array dimensions mismatch\n', [a.shape for a in arrays]
		print str(msg)

def paste (*arrays):
	""" Merge a set of arrays with the same number of rows to a new array with nColumns = sum(arrays.shape[1]) . """
	# for an alternative approach using concatenate, see function get_cira in module atmos.py
	# paste is !almost! equivalent to numpy.column_stack
	shapes  = [np.shape(a) for a in arrays]
	nRows   = [s[0] for s in shapes]
	if min(nRows)==max(nRows):
		maxCols = max(map(len,shapes))
		if maxCols==1:
			new = np.array(arrays)
		else:
			arrays = list(arrays)
			for i,a in enumerate(arrays):
				if len(a.shape)==1: arrays[i] = np.expand_dims(a,1)
			new = [a[:,j] for a in arrays for j in range(a.shape[1])]
		return np.transpose(new)
	else:
		print 'ERROR:  cannot paste arrays of different number of rows!'

####################################################################################################################################

def commonExtension (files):
	""" Return the common extension of all files, if identical; otherwise return None. """
	extensions = [os.path.splitext(file)[1] for file in files]
	ext0 = extensions[0]
	for ext in extensions:
		if not ext==ext0: return
	return ext0

####################################################################################################################################

def writeArray (xy, file=None, format='%g ', comments=[], commentChar='#'):
	""" Write a numeric array xy to file (or stdout if unspecified).
	    (format must have one, two, or xy.shape[1] specifiers, line feed is appended if necessary.) """
	if file:
		if hasattr(file,'closed'):
			if file.mode=='w':  out=file
			else:               print 'ERROR --- IO.writeArray:  file "' + file + '" opened in readmode!'
		else:
			out=open(file,'w')
			# print command line as first line to output file
			if not (sys.argv[0]=='' or 'ipython' in sys.argv[0]):
				sys.argv[0] = os.path.basename(sys.argv[0])
				out.write (commentChar + ' ' + join(sys.argv) + '\n' + commentChar + '\n')
	else:   out = sys.stdout

	for com in comments: out.write ( '%s %s\n' % (commentChar, strip(com)) )

	if len(xy.shape)==1:
		if count(format,'\n')==0: format = format.rstrip()+'\n'
		for i in range(xy.shape[0]): out.write (format % xy[i] )
	elif len(xy.shape)==2:
		npc = count(format,'%')
		if npc==1:
			format = xy.shape[1] * format
		elif npc==2 and xy.shape[1]>2:
			f2 = rfind(format,'%')
			format = format[:f2] + (xy.shape[1]-1) * (' '+format[f2:])
		elif npc!=xy.shape[1]:
			print "ERROR --- IO.writeArray:  check format (number of format specs does'nt match number of columns in data)"
			return
		if count(format,'\n')==0: format = format.rstrip()+'\n'
		for i in range(xy.shape[0]): out.write (format % tuple(xy[i,:]) )
	else:
		print 'ERROR --- IO.writeArray:  writing arrays with more than 2 dimensions not supported!'
	if not out.closed: out.close()

####################################################################################################################################

def open_outFile (outFile, commentChar='#'):
	""" Open output file and write job  specification (command line). """
	if outFile:
		try:
			out = open(outFile,'w')
		except IOError, errMsg:
			raise SystemExit, str(errMsg) + '\nERROR --- opening output file failed!'
		# print command line as very first record to file header
		sys.argv[0] = os.path.basename(sys.argv[0])
		if len(sys.argv)<20: sysArgv = join(sys.argv)
		else:                sysArgv = join(sys.argv[:10]) + ' ....'
		out.write (commentChar + ' ' + sysArgv + '\n' + commentChar + '\n')
	else:
		out = sys.stdout
	return out
