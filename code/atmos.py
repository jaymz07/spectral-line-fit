#!/usr/bin/env python
"""
atmos

Read atmospheric data file(s) in xy tabular ascii format:
extract profiles, interpolate to new altitude grid, save in xy or namelist format.

usage:

  atmos [options] file(s)

  -h          help
  -c char     comment character(s) used in input,output file (default '#')
  -o string   output file for saving of atmospheric data
  -f format   to be used in output file:
              'xy'        simple tabular form
	      'namelist'  for use in MIRART, GARLIC, BIRRA, ...
	      otherwise:  a prettified representation of the dictionary 
  -p float    pressure ratio for regridding (constraint: pRatio = p(lower)/p(upper) > 1)
  -x string   eXtract profiles from file(s) (a string with comma seperated entries, e.g., 'p,H2O,T')
  -z string   change altitude levels to (evtl. equidistant) grid:
              either the grid step,   i.e. simply specify just one number
              or     the grid limits, i.e. lower and upper bounds seperated with a comma
              or     an uniform grid, i.e. a triple of floats zLow,zHigh,deltaZ
              or (in apostrophes) 'start[step]stop', e.g. '20[2.5]50[5]80' to specify a grid with variable steps
 --pT         do NOT save altitude as first column in xy file
	      (use pressure as first column, e.g., to save a pressure/temp file as input for lbl2xs)
 --saveVMR    save mixing ratios [ppm] instead of number densities
 --ToA float  top-of-atmosphere altitude in km

"""

import os, sys
import re
from string import lower, upper, strip, index, split, join, count, ascii_uppercase
from types import *
import numpy as np

from ir import k as kBoltzmann
from cgsUnits import unitConversion, cgsTemperature

from molecules import molecules
molecNames   = molecules.keys()
molecNamesLo = [lower(mol) for mol in molecNames]

####################################################################################################################################
atmos_short_names = ['user', 'TRO', 'MLS', 'MLW', 'SAS', 'SAW', 'US']
atmos_long_names  = ['tropical', 'midlatitude summer', 'midlatitude winter', 'subarctic summer', 'subarctic winter', 'US standard']

def atmos_short2long_name (short):
	""" replace common acronyms used for AFGL atmospheric model profiles by long names used in namelist file. """ 
	SHORT = short.upper()
	if SHORT in atmos_short_names: return  atmos_long_names[atmos_short_names.index(SHORT)-1]
	else:                          return  short

####################################################################################################################################

class Atmos1D:
	def __init__ (self, atmData):
		""" Given a list of profile data (dictionaries with entries for 'what', 'data' etc) generate an atmos instance. """
		if not isinstance(atmData,list):  raise SystemExit, 'Atmos1D initialization failed, need a list of profiles!'
		self.gases = []

		for prof in atmData:
			what = prof.what
			if lower(what).startswith('alt') or lower(what)=='z': 
				self.z = unitConversion (prof.data, 'length', prof.unit)
			elif lower(what).startswith('press') or what=='p': 
				if hasattr(self,'p'):  raise SystemExit, 'ERROR Atmos1D init:  pressure already present!'
				self.p = unitConversion (prof.data, 'pressure', prof.unit)
			elif lower(what).startswith('temp') or what=='T': 
				if hasattr(self,'T'):  raise SystemExit, 'ERROR Atmos1D init:  temperature already present!'
				self.T = cgsTemperature (prof.data, prof.unit)
			elif lower(what).startswith('air') or lower(what).startswith('rho'): 
				if hasattr(self,'air'):  raise SystemExit, 'ERROR Atmos1D init:  "air" already present!'
				self.air = unitConversion (prof.data, 'density', prof.unit)
			elif what[0] in ascii_uppercase: 
			    	if what not in self.gases:  self.gases.append(what)
			    	else:                       raise SystemExit, '%s %s %s' % ('ERROR Atmos1D init:  gas ', what, 'already defined')
				
			    	if prof.unit.startswith('pp') or prof.unit=='%' or lower(prof.unit)=='vmr':
					densities = unitConversion (prof.data, 'MixingRatio', prof.unit)
			    	elif prof.unit.endswith('m-3') or prof.unit.endswith('m**3') or prof.unit.endswith('m^3'): 
					densities = unitConversion (prof.data, 'density', prof.unit)
				else:
					raise SystemExit, '%s %s %s %s' % ('unknown profile unit ', repr(prof.unit), ' for molecule',  prof.what)
				if not hasattr(self,'densities'): self.densities = densities
				else:                             self.densities = np.column_stack((self.densities,densities))
		# pressure <---> temperature <---> air number density
		if   hasattr(self,'p')+hasattr(self,'T')+hasattr(self,'air')==3:  print 'Atmos1d: got p, T, air'
		elif hasattr(self,'p') and hasattr(self,'T'):       self.air = self.p / (kBoltzmann*self.T)
		elif hasattr(self,'T') and hasattr(self,'air'):     self.p   = self.air * (kBoltzmann*self.T)
		elif hasattr(self,'p') and hasattr(self,'air'):     self.T   = self.p / (kBoltzmann*self.air)
		else:  raise SystemExit, 'ERROR --- Atmos1D init: need at least two of three data: p, T, air'
		# finally convert all mixing ratios to number densities
		if len(self.gases)>0:
			if len(self.densities.shape)==1:  self.densities = np.expand_dims(self.densities,1)
			for m,gas in enumerate(self.gases):
				if max(self.densities[:,m])<1.0:
					self.densities[:,m] = self.air * self.densities[:,m]
					print 'gas #', m+1, gas, ' vmr converted to number density'
	def __str__(self):
		format = '%10g %10.4g %7.2f %10.4g' + len(self.gases)*'%10.4g' +'\n'
		atmString = str('%10s %10s %7s %10s' + len(self.gases)*'%10s' + '\n') % (('z', 'p', 'T', 'air') + tuple(self.gases))
		for l,z in enumerate(self.z): atmString += format % ((z, self.p[l], self.T[l], self.air[l]) + tuple(self.densities[l,:]))
		return atmString
	def __len__(self): 
		return len(self.z) * (len(self.gases) + hasattr(self,'p')+hasattr(self,'T')+hasattr(self,'air'))
	def shape(self): 
		return len(self.z), len(self.gases) + hasattr(self,'p')+hasattr(self,'T')+hasattr(self,'air')
	def get(self,what):
		if hasattr(self,what):
			return getattr(self,what)
		else:
			if   what=='air':                                  return self.p / (kBoltzmann*self.T)
			elif what=='p' or lower(what).startswith('pres'):  return self.air * (kBoltzmann*self.T)
			elif what=='T' or lower(what).startswith('temp'):  return self.p * (kBoltzmann*self.air)
			elif what in self.gases:                           return self.densities[:,self.gases.index(what)]
	def vcd(self, gas=None, zBottom=0.0, zTop=None):
		""" Vertical column densities. """
		columns = np.zeros(len(self.gases))
		deltaZ = self.z[1:] - self.z[:-1] 
		# for m,gas in enumerate(self.gases):
		#	columns[m] = sum(0.5*deltaZ*(self.densities[1:,m] + self.densities[:-1,m]))
		return np.array([sum(0.5*deltaZ*(self.densities[1:,m] + self.densities[:-1,m])) for m in range(len(self.gases))])
	def perturb(self, what, amount=1.0, how='%', where=None):
				raise SystemExit, 'no perturbation yet'

####################################################################################################################################

class Profile:
	what=''; unit=''; where=''; when=''
	def __init__ (self, data, what=None, unit  = '', when  = '', where = ''):
		if   isinstance(data,np.ndarray):
			self.data = data
		elif   isinstance(data,(list,tuple)):
			self.data = np.array(data,np.float)
		elif isinstance(data,(float,int)):
			self.data = np.array([data],np.float)
		elif isinstance(data,str):
			if ',' in data: self.data = np.array(split(data,','),np.float)
			else:           self.data = np.array(split(data),np.float)
		if what:  self.what = what
		else:     raise SystemExit, 'ERROR:  Profile needs at least an identifier name "what"'
		if unit:  self.unit  = strip(unit)
		if where: self.where = strip(where)
		if when:  self.when  = strip(when)
	def __repr__ (self): 
		text  = "Profile(" + repr(len(self.data)*' %g' % tuple(self.data) )
		if self.what:   text = text + ', what='+repr(self.what)
		if self.unit:   text = text + ', unit='+repr(self.unit)
		if self.where:  text = text + ', where='+repr(self.where)
		if self.when:   text = text + ', when='+repr(self.when)
		return text + ")"                    
	def __str__ (self): 
		if len(self.data)>10: dataString = '%g %g %g ..... %g %g %g  (%i)' \
		                                   % tuple(np.concatenate((self.data[:3],self.data[-3:],[len(self.data)])))
		else:                 dataString = len(self.data)*' %g' % tuple(self.data)
		if self.unit:  return '%s [%s]: %s' % (self.what, self.unit, dataString)
		else:          return '%s: %s'      % (self.what,            dataString)
	def __len__(self): 
		return len(self.data)
	def __add__ (self, other):
		self.data = self.data + other
		return self
	def __mul__ (self, other):
		self.data = self.data * other
		return self
	def __radd__ (self, other):
		self.data = self.data + other
		return self
	def __rmul__ (self, other):
		self.data = self.data * other
		return self
	def __rdiv__ (self, other):
		self.data = self.data / other
		return self



##############################################################################################################
def get_namelist_atmos (file, extract='*', where=None):	
	""" Read a namelist formatted atmospheric data file (as used, e.g., by Mirart) and return some profiles. """
	import nameList
	# open first file, read everything as a string, and split into list of strings
    	file = open(file)
	data = file.read()
	file.close()
	# Split dataset into a list of namelist records, returned as strings
	records = nameList.findall_namelists(data)
	# Parse string containing a namelist record into its identifier (group name) and key/value(s) pair(s)
	parsedRecords = [nameList.parse_namelist (record) for record in records]
	# first record should be altitude grid
	if parsedRecords[0]['namelist']=='altitude':
		profileList = [Profile(parsedRecords[0]['z'], what='altitude', unit=parsedRecords[0]['unit'])]
	else:
		raise SystemExit, 'ERROR --- get_namelist_atmos: first record in file is not altitude!'
	# add further records to list
	if extract in ('*', 'all'):
		for record in parsedRecords[1:]:
			print  record['what'], record.keys(), (not where), ('where' not in record), record.get('where')==where
			if (not where) or ('where' not in record or strip(record.get('where'))=='' or record.get('where')==where):
				profileList.append( Profile(record['prof'], what=record['what'], unit=record['unit']) )
	else:
		for record in parsedRecords[1:]:
			print 'where:', record.get('where'), '-->', where
			if record['what'] in extract:
				if (not where) or record.get('where')==where:
					profileList.append( Profile(record['prof'], what=record['what'], unit=record['unit']) )
	return profileList

####################################################################################################################################

def read_atmos_xy (file, extract='*', commentChar='#'):
	""" Read tabular ascii atmospheric data files and return header lines AND ALL data saved as instance of class Atmos1D. """
	try:             input  = open(file)
	except IOError:  raise SystemExit, 'ERROR --- read_atmos_xy:  Opening file "' +  file +  '" failed, check file specification'
	# get comments from file
	lcc = len(commentChar)
	where = os.path.split(file)[1] # possibly overwritten by file header content
	whatsGiven = ''
	while 1:
		record = input.readline()
		if record.startswith(commentChar):
			headerLine = lower(strip(record[lcc:]))
			if   headerLine.startswith("where:"):    where      = record[index(record,':')+1:]
			elif headerLine.startswith("when:"):     when       = record[index(record,':')+1:]
			elif headerLine.startswith("what:"):     whatsGiven = split(record[index(record,':')+1:])
			elif headerLine.startswith("height "):   whatsGiven = split(record[lcc:]) # just a short version of "what: ..."
			elif headerLine.startswith("altitude "): whatsGiven = split(record[lcc:]) # just a short version of "what: ..."
			elif headerLine.startswith("z "):        whatsGiven = split(record[lcc:]) # just a short version of "what: ..."
			elif headerLine.startswith("units:"):    unitsGiven = split(record[index(record,':')+1:])
			elif headerLine.startswith("km "):       unitsGiven = split(record[lcc:]) # just a short version of "units: ..."
			elif re.findall('\w.*?[[(]\w.*?[])]',headerLine):
				whatsGivenAndUnits = [item.groups() for item in re.finditer('(\w.*?)[[(](\w.*?)[])]', headerLine)]
				whatsGiven = [what[0] for what in whatsGivenAndUnits]
				unitsGiven = [what[1] for what in whatsGivenAndUnits]
		else:
			break
	input.close()

	if whatsGiven:  whatsGiven = correct_names (whatsGiven)
	else:           raise SystemExit, "%s %s\n%s" % ("ERROR --- read_atmos_xy:  reading file ", file,
	                                                 "                          no information on data entries (column id's) found in file header!")

		
	# get numeric data from file
	allProfiles = np.loadtxt(file, comments=commentChar)

	# finally check consistency
	if not len(whatsGiven)==len(unitsGiven):
		raise SystemExit, '%s\n%i %s %s\n%i %s %s' % ('ERROR --- read_atmos_xy: check profile identifiers and units',
		                  len(whatsGiven), ' profiles:', join(whatsGiven), len(unitsGiven), ' units:   ', join(unitsGiven))

	if len(whatsGiven)==allProfiles.shape[1]:
		pass
	elif len(whatsGiven)==allProfiles.shape[1]-1:
		# probably very first column is an integer numbering of the levels
		if max(abs(allProfiles[:,0]-np.arange(allProfiles.shape[0])-1))<0.001:
			print 'WARNING: ignoring the very first column,  seems to be a counter of the levels!'
			allProfiles = allProfiles[:,1:]
	else:
		raise SystemExit, '%s %i %i' % ('number of profiles (columns) and number of identifiers in file header ("#what:" record) inconsistent', len(whatsGiven), allProfiles.shape[1])
	print len(whatsGiven), ' atmospheric data profiles for ', allProfiles.shape[0], ' levels:', whatsGiven

	if extract in ('*', 'all', 'All', 'ALL'):
		profileList = [Profile(allProfiles[:,m], what=what, unit=unitsGiven[m]) for m,what in enumerate(whatsGiven)]
	else:
		extract = correct_names(extract)
		if not 'altitude' in extract: extract.insert(0,'altitude')
		profileList = []
		for m,need in enumerate(extract):
			try:
				iz = whatsGiven.index(need)
				profileList.append( Profile(allProfiles[:,iz], what=need, unit=unitsGiven[iz]) )
			except ValueError:
				raise SystemExit, 'ERROR --- read_atmos_xy:  cannot find '+ repr(need) + \
				                  ' in list of profiles given: ' + join(whatsGiven)

	return profileList

####################################################################################################################################

def correct_names (names):
    	""" Loop thru a list of names and correct some names. """
    	for j,name in enumerate(names):
    		given = strip(upper(name))
		if   given=='T'   or given[:4]=="TEMP":                      names[j] = 'temperature'
		elif given=='P'   or given[:4]=="PRES":                      names[j] = 'pressure'
		elif given=='AIR' or given=='RHO' or given[:4]=="DENS":      names[j] = 'air'
		elif given=='Z'   or given[:3]=="ALT" or given=="HEIGHT":    names[j] = 'altitude'
		else:
			if name[0] in ascii_uppercase:
				names[j] = strip(name)
			elif name.lower() in molecNamesLo:
				names[j] = molecNames[molecNamesLo.index(name)]
				print 'WARNING:  replaced molecular name: ', name, ' --->', names[j]
			else:
				raise SystemExit, 'ERROR --- atmos:  unknown molecule!?! ' + repr(name)
			
	return names

####################################################################################################################################

def write_atmProfiles_namelist (out, profileList):
	""" Write atmospheric profiles in namelist format suitable for MIRART. """
	nz = len(profileList[0].data)
	# print altitude information
	write_namelist_altitude (out, profileList[0].data, profileList[0].unit, format='%10.2f')
	for profile in profileList[1:]:
		if nz==len(profile.data):  write_namelist_profile (out, profile.data, profile.what, profile.unit, format='%10.4g')
		else:                      raise SystemExit, 'ERROR --- write_atmProfiles_namelist:  profiles have different length!'
	
def write_atmos_namelist (out, data, saveVMR=False):
	""" Write atmospheric data in namelist format suitable for MIRART. """
	nz = len(data.z)
	# print altitude information
	write_namelist_altitude (out, unitConversion(data.z,'length',new='km'), unit='km', format='%10.2f')
	# print profiles
	#write_namelist_profile (out, data.p, 'pressure', unit='g/cm/s**2', format='%10.4g')
	write_namelist_profile (out, unitConversion(data.p, 'p', new='mb'), 'pressure', unit='mb', format='%10.4g')
	write_namelist_profile (out, data.T, 'temperature', unit='K', format='%10.2f')
	write_namelist_profile (out, data.air, 'air', unit='1/cm**3', format='%10.4g')
	for m,gas in enumerate(data.gases):
		vmr = max(data.densities[:,m])<1.0 # very primitive test !!!
		if vmr:
			write_namelist_profile (out, data.densities[:,m], data.gases[m], unit='ppv', format='%10.4g')
		else:
			if saveVMR and hasattr(data,'densities'):
				VMR = 1.e6*data.densities[:,m]/data.air
				write_namelist_profile (out, VMR, data.gases[m], unit='ppm', format='%10.4g')
			else:
				write_namelist_profile (out, data.densities[:,m], data.gases[m], unit='1/cm**3', format='%10.4g')
		
def write_namelist_altitude (out, altitude, unit, format='%10.2f', nCols=8):
	nz = len(altitude)
	out.write (" &ALTITUDE  NZ  = %i\n"  % nz )
	out.write ("            Unit='%s'\n" % unit )
	out.write ("            Z = ")
	indent = ''
	for l in xrange(0,nz,nCols):
    		lo, hi = l, min(l+nCols,nz)
    		out.write ( indent + (hi-lo)*" %6.2f" % tuple(altitude[lo:hi]) + '\n' )
    		indent = 16*' '
	out.write ( " &END\n\n" )

def write_namelist_profile (out, profile, what, unit=None, where=None, when=None, format='%10g', nCols=8):
    	out.write ( " &PROFILE  What  = '" + what + "'\n" )
    	if unit:  out.write ( "           Unit  = '" + unit + "'\n" )
    	if where: out.write ( "           Where  = '" + where + "'\n" )
    	if when:  out.write ( "           When  = '" + when + "'\n" )
    	out.write ( "           Prof  = " )
       	indent = ''; format = ' '+strip(format)
	nLevels = len(profile)
	for l in xrange(0,nLevels,nCols):
		lo, hi = l, min(l+nCols,nLevels)
		out.write ( indent + (hi-lo)*format % tuple(profile[lo:hi]) + '\n' )
		indent = 19*' '
	out.write ( " &END\n\n" )

####################################################################################################################################
def write_atmos_xy (out, data, pT=False, saveVMR=False):
	""" print profiles in xy format suitable for pfui or xmgr. """
	nz, ny = data.shape()

	if saveVMR:
		if not hasattr(data,'air'):  raise SystemExit, 'ERROR --- write_atmos_xy:  no air number density ---> no conversion to VMR'

	saveAltitudes = not pT # do not print altitudes, use pressure as first column (e.g. to save a pressure/temp file as input for lbl2xs)

	# print column identifiers
	if saveAltitudes:
		what  = '%7s%3s' % ('#what: ', 'z')
		units = '%7s%3s' % ('#units:', 'km')
		zGrid = unitConversion (data.z, 'length', new='km')
	else:
		what=''; units=''

	# print file header (column id's)
	if hasattr(data,'p'):
		#what  = what  + ' %11s' % 'pressure';     units = units + ' %11s' % 'g/cm/s**2'
		what  = what  + ' %11s' % 'pressure';     units = units + ' %11s' % 'mb'
	if hasattr(data,'T'):
		what  = what  + ' %11s' % 'temperature';  units = units + ' %11s' % 'K'
	if hasattr(data,'air'):
		what  = what  + ' %11s' % 'air';          units = units + ' %11s' % '1/cm**3'
	if hasattr(data,'densities'):
		what  = what  + ' ' + join(['%11s' % gas for gas in data.gases])
		if saveVMR:  units = units + ' ' + join(['%11s'%'ppm' for gas in data.gases])
		else:        units = units + ' ' + join(['%11s'%'1/cm**3' for gas in data.gases])
	out.write (what + '\n' + units + '\n')

	# print data in nxy format
	for l in range(nz):
		if saveAltitudes:  out.write ( '%10.1f' % zGrid[l] )
		#if hasattr(data,'p'):   out.write ( ' %11g' % data.p[l])
		if hasattr(data,'p'):   out.write ( ' %11g' % unitConversion(data.p[l],'p',new='mb'))
		if hasattr(data,'T'):   out.write ( ' %11.3f' % data.T[l])
		if hasattr(data,'air'): out.write ( ' %11.4g' % data.air[l])
		if saveVMR:
			if hasattr(data,'densities'): out.write ( data.densities.shape[1]*' %11.3g' % tuple(1.e6*data.densities[l,:]/data.air[l]) )
		else:
			if hasattr(data,'densities'): out.write ( data.densities.shape[1]*' %11.3g' % tuple(data.densities[l,:]) )
		out.write('\n')

####################################################################################################################################

def change_pressure_grid (zGrid, pressures, pRatio=2.0):
	""" Setup a new 'optimized' altitude grid according to a given pressure ratio.
	    (assumes a nearly exponential pressure decrease with altitude) """
	# check pressure array
	if min(pressures)<=0.0:
		raise SystemExit, 'pressure values not all positive!'
	elif max(pressures[1:]-pressures[:-1])>=0.0:
		raise SystemExit, 'pressure array not monotone decreasing!'
	else:
		# everything fine, build up new p array stepwise
		pMin = pressures[-1]
		pNew = [pressures[0]]
		while pNew[-1]>=pMin: pNew.append(pNew[-1] / pRatio)

		# so far a standard list, convert to a numeric array
		pNew = np.array(pNew)
		# interpolate: needs monotone increasing abscissa, hence revert x and y before and afterwards
		zNew = np.flipud(np.interp (np.flipud(pNew), np.flipud(pressures),np.flipud(zGrid)))
		# slightly shrink z to avoid extrapolation
		zFac = (zGrid[0]-zGrid[-1]) / (zNew[0]-zNew[-1])
		zNew = ((zNew-zNew[0]) * zFac) + zNew[0]
    	return zNew, pNew

####################################################################################################################################

def change_altitude_grid (zOld, gridSpec):
	""" Setup a new altitude grid and interpolate profiles to new grid. """
    	zFirst, zLast =  zOld[0], zOld[-1]
    	#specs = re.split ('[-\s:,]+',gridSpec)
    	if count(gridSpec,'[')+count(gridSpec,']')==0:
		if count(gridSpec,',')==0:
			try:                deltaZ = float(gridSpec)
			except ValueError:  raise SystemExit, 'z grid spacing not a number!'
			# set up new altitude grid
			zNew = np.arange(zFirst, zLast+deltaZ, deltaZ)
		elif count(gridSpec,',')==1:
			try:                zLow,zHigh = map(float,split(gridSpec,','))
			except ValueError:  raise SystemExit, 'z grid spacing not a pair of floats!'
			# for new grid simply extract old grid points within given bounds (also include altitudes slightly outside)
			eps  = min( zOld[1:]-zOld[:-1] ) / 10.
			zNew = np.compress(np.logical_and(np.greater_equal(zOld,zLow-eps), np.less_equal(zOld,zHigh+eps)), zOld)
		elif count(gridSpec,',')==2:
			try:                zLow,zHigh,deltaZ = map(float,split(gridSpec,','))
			except ValueError:  raise SystemExit, 'z grid spacing not a triple of floats (zLow.zHigh,deltaZ)!'
			# set up new altitude grid
			zNew = np.arange(max(zLow,zFirst), min(zHigh,zLast)+deltaZ, deltaZ)
		elif count(gridSpec,',')>2:
			try:                zNew = np.array(map(float, split(gridSpec,',')))
			except ValueError:  raise SystemExit, 'z grid not a set of floats separated by commas!'
    	elif count(gridSpec,'[')==count(gridSpec,']') > 0:
      		zNew = parseGridSpec (gridSpec)
		if not zFirst <= zNew[0] < zNew[-1] <= zLast:
			raise SystemExit, '%s  %f %f  %s  %f %f' % ('ERROR: new zGrid', zNew[0],zNew[-1], ' outside old grid', zFirst, zLast)
    	else:
       		raise SystemExit, 'New altitude not specified correctly\n' + \
       		      'either simply give altitude step size, a pair of lower,upper limits,  or "start(step)stop"!'
    	return zNew

####################################################################################################################################

def parseGridSpec (gridSpec):
	""" Set up (altitude) grid specified in format 'start[step1]stop1[step2]stop' or similar. """
	# get indices of left and right brackets
	lp = [];  rp = []
	for i in xrange(len(gridSpec)):
		if   (gridSpec[i]=='['):  lp.append(i)
		elif (gridSpec[i]==']'):  rp.append(i)
		else:                     pass
	if len(lp)==len(rp):
		gridStart = [];  gridStop = [];  gridStep = []
		for i in xrange(len(lp)):
			if i>0:  start=rp[i-1]+1
			else:    start=0
			if i<len(lp)-1: stop=lp[i+1]
			else:           stop=len(gridSpec)

			try:
				gridStart.append(float(gridSpec[start:lp[i]]))
			except ValueError:
				print 'cannot parse grid start specification\nstring not a number!'
				raise SystemExit
			try:
				gridStep.append(float(gridSpec[lp[i]+1:rp[i]]))
			except ValueError:
				print 'cannot parse grid step specification\nstring not a number!'
				raise SystemExit
			try:
				gridStop.append(float(gridSpec[rp[i]+1:stop]))
			except ValueError:
				print 'cannot parse grid stop specification\nstring not a number!'
				raise SystemExit
			if i==0:
				if gridStop[0]<=gridStart[0]: print 'incorrect grid specification:  Stop < Start'; raise SystemExit
				newGrid = np.arange(gridStart[0], gridStop[0]+gridStep[0], gridStep[0])
			else:
				if gridStop[i]<=gridStart[i]: print 'incorrect grid specification:  Stop < Start'; raise SystemExit
				newGrid = np.concatenate((newGrid[:-1],np.arange(gridStart[i], gridStop[i]+gridStep[i], gridStep[i])))
	else:
		print 'cannot parse grid specification\nnumber of opening and closing braces differs!\nUse format start[step]stop'
		raise SystemExit
	# set up new altitude grid
	return newGrid

####################################################################################################################################
def interpolate_atmos1d (data, zNew, interpolate=None):
	""" Interpolate atmospheric profiles to a new altitude grid. """
	if interpolate=='h':
		try: from pchi import hermite
		except ImportError: interpolate=''; print 'import pchi failed, switching to default interpolation'

    	# initialize new profile matrix
    	densityNew  = np.zeros([len(zNew),len(data.gases)], np.float)
    	# copy original altitude grid
    	zOld = data.z
    	# loop over all original profiles
	if interpolate=='h':
		data.p= np.exp(hermite(zOld, np.log(data.p), zNew)); print 'interpolating log(p)'
		data.T=        hermite(zOld,        data.T,  zNew)
		data.air=      hermite(zOld,      data.air,  zNew)
		for m in range(data.densities.shape[1]): densityNew[:,m] = hermite (zOld, data.densities[:,m], zNew)
	else:
		data.p= np.exp(np.interp(zNew, zOld, np.log(data.p))); print 'interpolating log(p)'
		data.T=        np.interp(zNew, zOld,        data.T)
		data.air=      np.interp(zNew, zOld,      data.air)
		for m in range(data.densities.shape[1]): densityNew[:,m] = np.interp (zNew, zOld, data.densities[:,m])
    	# replace altitudes and profiles in dictionary
    	data.z = zNew
    	data.densities = densityNew
    	return data

####################################################################################################################################

def reGridAtmos1d (atmos1d, zToA=None, pRatio=None, zGrid=None):
	if bool(zGrid)+bool(pRatio)+bool(zToA)>1:
		raise SystemExit, 'either specify ToA altitude OR new zGrid OR pressure ratio for atmospheric profile interpolation!'
	elif pRatio>0.0 or zGrid:
		zOld = atmos1d.z
		if pRatio>0.0:
			zNew, pNew = change_pressure_grid (atmos1d.z, atmos1d.p, pRatio)
		else:           
			zNew = change_altitude_grid (atmos1d.z, zGrid)
			pNew = np.exp(np.interp(zNew, zOld, np.log(atmos1d.p))); atmos1d.p = pNew
		atmos1d.z   = zNew
		atmos1d.T   = np.interp(zNew, zOld, atmos1d.T)
		atmos1d.air = np.exp(np.interp(zNew, zOld, np.log(atmos1d.air)))
		for m,mol in enumerate(atmos1d.gases):
    			atmos1d.densities[:,m] = np.interp(zNew, zOld, atmos1d.densities[:,m])
    	elif zToA:
		if zToA<200.:
			zToA = unitConversion(zToA, 'length', old='km')
			print 'WARNING --- reGridAtmos1d:  zToA very small, assuming kilometer units'
		nz = sum([1 for z in atmos1d.z if z<=zToA])
		atmos1d.z = atmos1d.z[:nz]
		atmos1d.p = atmos1d.p[:nz]
		atmos1d.T = atmos1d.T[:nz]
		atmos1d.air = atmos1d.air[:nz]
		atmos1d.densities = atmos1d.densities[:nz,:]
	return atmos1d

####################################################################################################################################

def reGridProfiles (profileList, zToA=None, zGrid=None, pRatio=None):
	""" Interpolate atmospheric profiles to new grid either using altitude or pressure specifications. """
        # find altitude grid array
	if profileList[0].what=='altitude':
		zOld = profileList[0].data; zUnit = profileList[0].unit
	else: raise SystemExit, 'very first profile in profileList is not "altitude"!'

    	# set up a new grid
	if zGrid and pRatio:
		raise SystemExit, 'either specify new zGrid OR pressure ratio for atmospheric profile interpolation!'
	elif pRatio>0.0:
		# check if pressure array is unique and extract (unique) pressure data
		if not sum([prof.what=='pressure' for prof in profileList])==1:
			print [prof.what for prof in profileList]
			raise SystemExit, 'ERROR --- reGridProfiles: need exactly one pressure profile'
		for prof in profileList:
			if prof.what=='pressure':  pressures = prof.data; break
		zNew, pNew = change_pressure_grid (zOld, pressures, pRatio)
		print 'change_pressure_grid', len(zNew), len(pNew), zNew, pNew
		if zToA:
			zNew = zNew[:sum(zNew<=zToA)]
			pNew = pNew[:sum(zNew<=zToA)]
	elif zGrid and zToA:
		raise SystemExit, 'either specify new zGrid for atmospheric profile interpolation OR top-of-atmosphere altitude!'
    	elif zGrid:
		zNew = change_altitude_grid (zOld, zGrid)
    	elif zToA:
		nz = sum([1 for z in zOld if z<=zToA])
		for prof in profileList:  prof.data = prof.data[:nz]
		return profileList
    	else:
		print 'z Grid o.k.!!!'

	# interpolate profiles
	for prof in profileList[1:]:
		if   prof.what=='altitude': 
			if not prof.unit==zUnit:  raise SystemExit, 'ERROR --- reGridProfiles:  altitude grids have different units!'
			zOld=prof.data; profilelist.remove(prof); print 'removed further altitude grid '+ str(prof.data)
		elif prof.what in ('pressure','air'):
			prof.data = np.exp(np.interp(zNew, zOld, np.log(prof.data)))
		else:
			prof.data = np.interp(zNew, zOld, prof.data)
	# finally replace altitude grid
	profileList[0].data = zNew

	return profileList

####################################################################################################################################

def get_profile_list (atmFiles, commentChar='#', extract='*', where=None, verbose=False):
	""" Read a list of atmospheric data files and return a list of profiles (not necessarily with the same z grid). """
	print 'get_profile_list:', extract
	if isinstance(atmFiles,str): atmFiles = [atmFiles]
	if isinstance(extract,str) and not extract.lower() in ('*', 'all'):
		print 'extracting subset of profiles: ', extract, 
		extract = correct_names(split(extract,',')); print extract

	profileList = []
	for file in atmFiles:
		if os.path.splitext(file)[1]=='.nml':  profileList += get_namelist_atmos (file, extract, where)
		else:                                  profileList += read_atmos_xy (file, extract, commentChar)

	if verbose:
		print len(profileList), ' profiles accepted'
		for prof in profileList:  print prof

	return profileList

####################################################################################################################################

def _atmos_ (atmFiles, outFile=None, commentChar='#', extract='*', where=None, zToA=None, pRatio=None, zGrid=None, format=None, saveVMR=False, pT=False, verbose=False):

	if isinstance(where, str):
		# replace acronyms with names used in namelist file
		where = atmos_short2long_name(where)

 	profileList = get_profile_list (atmFiles, commentChar, extract, where, verbose)

	if zToA or pRatio or zGrid:
		profileList = reGridProfiles (profileList, zToA, zGrid, pRatio)
		if verbose:
			for prof in profileList:  print prof

	atmos1d = Atmos1D(profileList)

	out = open_outFile (outFile, commentChar)
	if format in ('namelist', 'nml'):
		write_atmos_namelist (out, atmos1d, saveVMR)
	elif format=='xy':
		write_atmos_xy (out, atmos1d, pT, saveVMR)
	else:
		from pprint import pprint
		pprint(atmos1d.__dict__, out)
	# close the output file (if its not stdout)
	if outFile: out.close()

####################################################################################################################################

if (__name__ == "__main__"):
 	
	from IO import open_outFile
	from command_parser import parse_command, standardOptions
	from pairTypes import Interval

        opts = standardOptions + [  # h=help, c=commentChar, o=outFile
               {'ID': 'ToA', 'name': 'zToA', 'type': FloatType, 'constraint': 'zToA>0.0'},
               {'ID': 'f', 'name': 'format', 'type': StringType},
               {'ID': 'p', 'name': 'pRatio', 'type': FloatType, 'constraint': 'pRatio>1.0'},
               {'ID': 'x', 'name': 'extract', 'type': StringType, 'default': '*'},
               {'ID': 'z', 'name': 'zGrid', 'type': StringType},
               {'ID': 'where', 'type': StringType},
               {'ID': 'v', 'name': 'verbose'},
               {'ID': 'saveVMR'},
               {'ID': 'pT'}
               ]

	atmFiles, options, commentChar, outFile = parse_command (opts,[1,99])

	if options.has_key('h'):
		print __doc__%globals();  raise SystemExit, "end of atmos help"
	else:
		options['verbose'] = options.has_key('verbose')
		options['saveVMR'] = options.has_key('saveVMR')
		options['pT']      = options.has_key('pT')
		_atmos_ (atmFiles, outFile, commentChar, **options)
