import getopt
import sys
import os
import re
from string import replace, strip, lower, upper, capitalize, swapcase, split, splitfields
from types import *
import keyword
#import Numeric as Num
import numpy   as np
from IO import commonExtension
from pairTypes import *

standardOptions = [{'ID': 'h'},
	           {'ID': 'c', 'type': StringType, 'default': '#'},
                   {'ID': 'o', 'type': StringType}]

####################################################################################################################################

def prepare4getopt (knownOptions):
	""" Extract short and long option id's as an input for getopt. """
	ShortOptions = ''
	LongOptions  = []
	for option in knownOptions:
		# retrieve type from default value
		if  not option.has_key('type') and option.has_key('default'): option['type'] = type(option['default'])

		if option.has_key('ID') or option.has_key('id'):
		        ID = option.setdefault('ID',option.get('id'))
			if len(ID)==1:
				if option.has_key('type'): ShortOptions = ShortOptions + ID + ':' 
				else:                      ShortOptions = ShortOptions + ID
			else:
				if option.has_key('type'): LongOptions.append(ID+'=')
				else:                      LongOptions.append(ID)
		elif option.has_key('name'):
			option['ID'] = option['name']
			if option.has_key('type'): LongOptions.append(ID+'=')
			else:                      LongOptions.append(ID)
		else:
			raise SystemExit, 'option specification requires "ID" and/or "name"'
	return ShortOptions, LongOptions

####################################################################################################################################

def getopt_parser (ShortOptions, LongOptions):
	""" Parse command line string using getopt.  Return list of files and a dictionary of options! """
    	args = sys.argv[1:]
    	if len(args)==1 and 'help=' in LongOptions and '--help' in args:
		options = {'h': None}  # print standard help message
		files   = []
    	elif len(args)>0:
        	try:
            		if LongOptions:
	    			OptionsList, files = getopt.getopt(args, ShortOptions, LongOptions)
	    		else:
               			OptionsList, files = getopt.getopt(args, ShortOptions)
        	except getopt.error, errMsg:
            		print "\ncheck your options list!!!"
			print errMsg
            		print "valid options are ",
            		for i in range(len(ShortOptions)):
                		if not ShortOptions[i]==':':  print '-' + ShortOptions[i],
            		if LongOptions:
            			for option in LongOptions: print '--' + replace(option,'=',''),
            		raise SystemExit, "parsing input arguments failed!"
        	# return options as a dictionary (getopt returns a double list!)
        	options = {}
        	for i in range(len(OptionsList)):
            		key          = replace(OptionsList[i][0],'-','') # remove leading dash(es)
            		options[key] = OptionsList[i][1] 
    	else:
        	files   = []
       		options = {}
	return files, options	

####################################################################################################################################

def check_constraint (value, name, constraint):
	if name in constraint:
		if keyword.iskeyword(name):
			raise SystemExit, 'name conflict: ' + name + ' is a reserved word in PYTHON!'
		if type(value) in (IntType,FloatType):
			statement = name + '=' + `value`
			exec  statement
			if  not eval(constraint):
				raise SystemExit, statement + '\nconstraint ' + repr(constraint) + ' violated!'
		elif type(value) is np.ndarray:
			statement = name + ' = np.' + `value`
			exec  statement
			if not np.alltrue(eval(constraint)):
				raise SystemExit, statement + '\nconstraint ' + repr(constraint) + ' violated (comparison elementwise)!'
		#elif type(value) is Num.ArrayType:
		#	statement = name + ' = Num.' + `value`
		#	exec  statement
		#	if not Num.alltrue(eval(constraint)):
		#		raise SystemExit, statement + '\nconstraint ' + repr(constraint) + ' violated (comparison elementwise)!'
		elif type(value) is StringType:
			statement = name + '=' + `value`
			exec  statement
			if  not eval(constraint):
				raise SystemExit,  statement + '\nconstraint ' + repr(constraint) + ' violated!'
		elif isinstance(value,(Interval,PairOfInts,PairOfFloats)):
			statement = name + '=' + `value`
			exec  statement
			if  not eval(constraint):
				raise SystemExit,  statement + '\nconstraint ' + repr(constraint) + ' violated!'
		else:
			raise SystemExit, 'unknown/unsupported type ' + repr(type(value)) + ' for check_constraint'
	else:
		raise SystemExit, 'Variable name ' + repr(name) + ' not used in constraint expression ' + repr(constraint)

####################################################################################################################################

def check_type (id, name, given, oType):
	try:
		if   oType==StringType:             typeChecked = strip(given)
		elif oType==FloatType:              typeChecked = float(given)
		elif oType==IntType:                typeChecked = int(given)
		#elif oType in (ListType,TupleType,np.ndarray,Num.ArrayType, Interval, PairOfInts, PairOfFloats):
		elif oType in (ListType,TupleType,np.ndarray, Interval, PairOfInts, PairOfFloats):
			typeChecked = re.split('[,;\s]',given)
			if   oType==np.ndarray:
				typeChecked=np.array(map(float,typeChecked))
			#elif   oType==Num.ArrayType:
			#	typeChecked=Num.array(map(float,typeChecked))
			elif oType==Interval: 
				low, hi = map(float,typeChecked)
				typeChecked=Interval(low,hi)
			elif oType==PairOfInts: 
				left, right = map(int,typeChecked)
				typeChecked=PairOfInts(left,right)
			elif oType==PairOfFloats: 
				left, right = map(float,typeChecked)
				typeChecked=PairOfFloats(left,right)
			#else: print 'unrecognized type ', id, name, given, typeChecked, oType, oType==np.ndarray
		else:
			raise SystemExit, 'type ' + repr(option['type']) + ' not yet supported, sorry'
	except ValueError, errMsg:
		raise SystemExit, 'ERROR:  option ' + id + ' = ' + name + '   ' + str(errMsg)
	except:	
		print id, name, given, oType
		raise SystemExit, 'type checking of options failed!'
	return typeChecked

####################################################################################################################################

def check_options (optionsGiven, knownOptions, verbose=0):
	""" Check the options specified on the command line wrt type and constraits, add unspecified options with defaults if available."""
	for option in knownOptions:
		id = option.get('ID',option.get('id'))
		name = option.get('name',id)
		if optionsGiven.has_key(id):
			given = optionsGiven.pop(id)
			if option.has_key('type'):
				typeChecked = check_type (id, name, given, option['type'])
			else:
				typeChecked = given
			if option.has_key('constraint'):
				check_constraint (typeChecked, option['name'], option['constraint'])
			optionsGiven[name]=typeChecked
		else:	
			if option.has_key('default'):
				optionsGiven[name]=option['default']
				if verbose and not id=='c': print id, name, 'set to default:', optionsGiven[name]
	return optionsGiven

####################################################################################################################################

def parse_command (knownOptions, numFiles=None, env4defaults='', verbose=0):
    	""" Parse command line arguments or interactively ask directly for files,  return options as a dictionary. """

	if env4defaults:
        	userDefaults = os.environ.get(env4defaults)
		if userDefaults:  knownOptions = change_defaults (knownOptions, userDefaults)
	
	# translate options specs into string (for short options) and optionally list (for long options) appropriate for standard getopt
	ShortOptions, LongOptions = prepare4getopt (knownOptions)

	# parse the command using getopt, but return a dictionary instead of getopt's list of two-element tuples (only with those options specified in the command line!)
	files, optionsGiven = getopt_parser (ShortOptions, LongOptions)

	# check options for type etc, return dictionary, now with defaults added!
	optionsGiven = check_options (optionsGiven, knownOptions, verbose)

	# check number of (input) files
	if not (optionsGiven.has_key('h') or optionsGiven.has_key('help')):
		commandName = os.path.basename(sys.argv[0])
		if isinstance(numFiles, int) and numFiles>0: 
			if not len(files)==numFiles:
				raise SystemExit, 'incorrect number of files/arguments, need exactly ' + str(numFiles) + \
				                  '\n(see "' + commandName + ' -h" for a detailed usage help)'
		elif isinstance(numFiles, (list,tuple)) and len(numFiles)==2: 
			if not numFiles[0]<=len(files)<=numFiles[1]:
				raise SystemExit, 'incorrect number of files/arguments, need ' + `numFiles[0]` + ' ... ' + `numFiles[1]` + \
				                  '\n(see "' + commandName + ' -h" for a detailed usage help)'

	# extract common options
	commentChar = optionsGiven.pop('c','#')
	outFile     = optionsGiven.pop('o',None)

    	return files, optionsGiven, commentChar, outFile

####################################################################################################################################

def multiple_outFiles (inFiles, outFile):
	""" Given a list of input file names and a 'template' (e.g. extension) for the output files, return a list of output file names. """
	if outFile:
		if len(inFiles)>1:
			if not outFile.startswith('.'): outFile='.'+outFile
			# -o option specifies extension of output files
			commonExt = commonExtension(inFiles)
			print 'commonExt:', commonExt
			# replace extension only when all input files have the same extension
			if commonExt:  outFiles = [os.path.splitext(file)[0]+outFile for file in inFiles]
			else:          outFiles = [file+outFile for file in inFiles] 
		else:
			if outFile.startswith('.'): return [os.path.splitext(inFiles[0])[0]+outFile] # do not write to a hidden file!
			else:                       return [outFile]
	else:
		outFiles = [None for file in inFiles]
	return outFiles

####################################################################################################################################

def change_defaults (opts, userDefaults):
	""" Replace standard options by user specified options. """
	newDefaults = dict([splitfields(strip(kv),'=') for kv in splitfields(userDefaults,';')])
	for opt in opts:
		#if opt.has_key('default') and newDefaults.has_key(opt['ID']):
		if opt['ID'] in newDefaults:
		    if 'default' in opt:
			if isinstance(opt['default'],int):      opt['default'] = int(newDefaults[opt['ID']])
			elif isinstance(opt['default'],float):  opt['default'] = float(newDefaults[opt['ID']])
			else:                                   opt['default'] = newDefaults[opt['ID']]
		    elif 'type' in opt:
			if   opt['type']==StringType:           opt['default']= newDefaults[opt['ID']]
			elif opt['type']==FloatType:            opt['default']= float(newDefaults[opt['ID']])
			elif opt['type']==IntType:              opt['default']= int(newDefaults[opt['ID']])
	print '\n---> new defaults <---'
	for opt in opts:
		if newDefaults.has_key(opt['ID']):  print opt['ID'], opt.get('default')
	print
	return opts

