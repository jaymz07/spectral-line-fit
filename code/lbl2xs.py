#!/usr/bin/env python

""" lbl2xs
  computation of line-by-line molecular absorption cross sections

  usage:
  lbl2xs [options] line_parameter_file(s)

  -h            help
  -c  char      comment character(s) used in input,output file (default '#')
  -o  string    output file for saving of cross sections (if not given: write to StdOut)
                (in case of a list of input files this options specifies just the extension)

  -p  float(s)  pressure(s) (in mb,  default: p_ref of linefile, usually 1013.25mb=1atm)
  -T  float(s)  Temperature(s) (in K, default: T_ref of linefile, usually 296K)
 --pT file      with list of pressures in mb (first column) and temperatures in K (second column)

  -n  int       number of grids --- selects 'multigrid mode' for nGrids=2 or 3 (default) (nGrids=1 brute force)
  -g  int       gridRatio = ratio of coarse to fine grid spacing (only 2, 4, or 8, default 8)
  -W  float     transition from fine to coarse grid (in units of half widths, default 25.0)

  -f  string    format for output file: 'a'='t'='xy' tabular ascii OR 'h' hitran OR pickled (default) 
  -i  string    interpolation method   [default: '3' three-point Lagrange,  choices are one of "234lqcsbk"]
                (only required for multigrid approach or when cross sections for several p,T pairs are saved in xy format)
  -L  char      Lineshape: V(oigt), L(orentz), G(auss)     [default: Voigt]
  -m  string    molecule (no default, should be given in the line file, required otherwise)
  -s  float     sampling rate used for x-grid (default: 5.0 grid points per (mean) half width)
  -w  float     wing extension (cutoff wavenumber, default 10.0cm-1)
  -x  Interval  lower,upper wavenumbers (comma separated pair of floats [no blanks!],
                                        default set according to range of lines in datafile)


  If several line parameter files are given (usually for several molecules)
  AND if an output file (extension) has been specified (-o option)
  a cross section file will be generated for each line file:
  * if all line files have the same extension, the cross section files will have the old basename with the extension as specified by the -o option
  * otherwise the input file name will be augmented by the string specified as -o option
  
  The line parameter file should contain a header section indicating the molecule, pressure, temperature,
  and type of columns followed by a list of (preselected) lines in the format
  # molecule: XYZ
  # temperature: 296 K
  # pressure:    1013.25 mb
  # format: vSEan
  position1  strength1  energy1  airWidth1  tDep1
  position2  strength2  energy2  airWidth2  tDep2
  position3  strength3  energy3  airWidth3  tDep3
  .........  .........  .......  .........  .....

  This file can be generated from HITRAN, SAO, or GEISA database with extract.py 
  The line list should contain at least 2 columns with line positions (in cm-1) and line strengths.
  If lower state energy is missing, cross sections can be calculated only at the reference temperature of the line list.
  If air broadened half width is missing, it can be set to a (constant) default value with the -a option.
  If temperature dependence of the air broadened half width is missing, the value in the molecular data file is used.
"""

####################################################################################################################################

# import some standard python modules
from string import *
import sys
import os.path
import time 

try:                      import numpy as np
except ImportError, msg:  raise SystemExit, str(msg) + '\nimport numeric python failed!'

from pairTypes import Interval
from lagrange_interpolation import *
from lineshapes import *
from lines import *
from cgsUnits import unitConversion

####################################################################################################################################

interpolationDict = {(2,2): lagrange2_interpolate2,
                     (3,2): lagrange3_interpolate2,
                     (4,2): lagrange4_interpolate2,
                     (2,4): lagrange2_interpolate4,
		     (3,4): lagrange3_interpolate4,
		     (4,4): lagrange4_interpolate4,
                     (2,8): lagrange2_interpolate8,
		     (3,8): lagrange3_interpolate8,
		     (4,8): lagrange4_interpolate8}

####################################################################################################################################

def lbl_3grids (positions, strengths, gammaL, gammaD,  uLow, uHigh,  lineShape='Voigt',
                sampling=5.0, gridRatio=4, nWidths=25.0, lagrange=3, verbose=False):
	""" Compute lbl cross sections using three grids: fine grid near line center, medium grid, and coarse grid everywhere.
	    (xs to be computed on monochromatic uniform grid in interval (uLow,uHigh) with spacing defined by half widths. """

	# check optional arguments, reset if not given
	if not lineShape.startswith('V'):
		raise SystemExit, 'ERROR --- lbl_3grids:  unknown/unsupported lineShape ' + lineShape

	# check optional spectral range, reset if not given (here just for safety, usually done by calling routine lbl_xs)
	if uLow and uHigh: pass
	elif uLow:         uHigh = positions[-1]
	elif uHigh:        uLow  = positions[0]
	else:              uLow  = positions[0];  uHigh = positions[-1]; print 'INFO --- lbl_3grids:  uLow, uHigh:', uLow, uHigh

	# Lagrange interpolation method
	interpolation = interpolationDict [(lagrange, gridRatio)]; #print ' lbl_3grids:  ', interpolation.__doc__

	# average line width for Lorentz, Gaussian, and Voigt line shapes
	meanGL, meanGD, meanGV  =  meanLineWidth (gammaL, gammaD)

	nGrids=3
	c2fhalf = gridRatio/2
	ipl = lagrange

	# fine monochromatic wavenumber grid 'u' with u[0]=uLow and u[-1]=uHigh
	uGrid = wavenumber_grid4xs (uLow, uHigh, meanGV, sampling, gridRatio, nGrids, verbose)
	nu    = len(uGrid)-1     # number of grid intervals
	du    = (uHigh-uLow)/nu  # grid spacing

	# medium resolution grid 'v' with v[1]=u[0] and v[-2]=u[-1]
	dv    = gridRatio*du
	recdv = 1.0/dv
	vGrid = np.arange(uLow-c2fhalf*dv,uHigh+(c2fhalf+0.9)*dv,dv)
	nv  = len(vGrid)-1  # number of grid intervals = index of last element

	# coarse resolution grid 'w' with w[1]=v[0] and w[-2]=v[-1]
	dw    = gridRatio*dv
	recdw = 1.0/dw
	wGrid = np.arange(vGrid[0]-c2fhalf*dw,vGrid[-1]+(c2fhalf+0.9)*dw,dw)
	nw  = len(wGrid)-1  # number of grid intervals = index of last element
	if verbose:
		print '\n lbl_3grids: ', '  uGrid', len(uGrid), du, '  vGrid', len(vGrid), dv, '  w', len(wGrid), dw
		print ' u:', 5*'%12.6f' % tuple(uGrid[:5]), ' ...', 5*'%12.6f' % tuple(uGrid[-5:])
		print ' v:', 5*'%12.6f' % tuple(vGrid[:5]), ' ...', 5*'%12.6f' % tuple(vGrid[-5:])
		print ' w:', 5*'%12.6f' % tuple(wGrid[:5]), ' ...', 5*'%12.6f' % tuple(wGrid[-5:])

	if len(wGrid)<ipl:
		raise SystemExit, 'ERROR --- lbl_3grids:  coarse grid W has less grid points than required for interpolation!'

	# and some frequently used wavenumber grid points
	u0      = uGrid[0]
	v0      = vGrid[0]
	vEnd    = vGrid[-1]
	w0      = wGrid[0]
	wEnd    = wGrid[-1]
	vLftMin = vGrid[-ipl]
	vRgtMax = vGrid[ipl-1]
	wLftMin = wGrid[-ipl]
	wRgtMax = wGrid[ipl-1]

	# allocate coarse, medium, and fine cross section arrays
	xsCrude   = np.zeros_like(wGrid)
	xsMedium  = np.zeros_like(vGrid)
	xsFine    = np.zeros_like(uGrid)
	vgtMedium = np.zeros_like(vGrid)
	vgtFine   = np.zeros_like(uGrid)

	# sum over lines
	for l in xrange(len(positions)):
		# compute cross section on entire coarse grid
		vgtCoarse = Voigt_Kuntz_Humlicek1 (wGrid, positions[l], strengths[l], gammaL[l], gammaD[l])
		xsCrude = xsCrude + vgtCoarse

		# set limits for medium resolution region
		centerExtension = nWidths*max(gammaL[l],gammaD[l])
		vLeft  = positions[l]-gridRatio*centerExtension
		vRight = positions[l]+gridRatio*centerExtension
		## need at least 2, 3, or 4 coarse grid points (according to degree of interpolation)
		#dbg if verbose: print ('# %12.6f  med: %f %f' % (lines[l,0], vLeft,vRight)), 
		if vRight<wRgtMax or vLeft>wLftMin:
			#dbg if verbose: print ' far wing only'
			continue   # line outside spectral interval: next line
		
		kLeft  = max(int((vLeft-w0)*recdw),c2fhalf)             # left and right coarse grid indices
		kRight = min(nw-int((wEnd-vRight)*recdw),nw-c2fhalf)   
		if kLeft>kRight: continue

		jLeft  = gridRatio*(kLeft-c2fhalf)                             # left and right medium grid indices
		jRight = gridRatio*(kRight-c2fhalf)+1                        # add 1 because of pythons range end
		#dbg if verbose:
		#	dbg print ('   %6i %6i  %f %f' % (kLeft,kRight, wGrid[kLeft], wGrid[kRight])), 
		#	dbg print ('   %6i %6i' % (jLeft,jRight)),
		#	dbg print (' %f %f' % (vGrid[jLeft], vGrid[jRight-1]))
		#	dbg print ('   %6i %6i  %f %f' % (jLeft,jRight, vGrid[jLeft], vGrid[jRight-1]))
		# compute cross section on medium grid in extended center region
		vgtMedium[jLeft:jRight] = Voigt (vGrid[jLeft:jRight], positions[l], strengths[l], gammaL[l], gammaD[l])
		# interpolate coarse grid cross section near extended line center
		li = interpolation (vgtCoarse[kLeft-1:kRight+2])
		# accumulate medium (near extended line center only) cross section
		xsMedium[jLeft:jRight] = xsMedium[jLeft:jRight] + vgtMedium[jLeft:jRight] - li[gridRatio:-gridRatio]

		# set limits for fine resolution region (for Voigt better use voigt width!!!)
		uLeft  = positions[l]-centerExtension
		uRight = positions[l]+centerExtension
		#dbg if verbose: print ('#              fine: %f %f' % (uLeft,uRight)), 
		if uRight<vRgtMax or uLeft>vLftMin:
			#dbg if verbose: print 'near wing only'
			continue   # line outside spectral interval: next line
		
		iLeft  = max(int((uLeft -v0)*recdv),c2fhalf)           # left and right medium grid indices
		iRight = min(nv-int((vEnd-uRight)*recdv),nv-c2fhalf)   
		if iLeft>iRight: continue                                                                        # bug fix 16jun2004
		#Right = min(int((uRight-v0)*recdv),nv-1) + 1  # add 1 because of python's range end
		jLeft  = gridRatio*(iLeft-c2fhalf)                           # left and right fine grid indices
		jRight = gridRatio*(iRight-c2fhalf) +1
		#dbg if verbose:
		#dbg 	print ('   %6i %6i  %f %f' % (iLeft,iRight, vGrid[iLeft], vGrid[iRight])), 
		#dbg 	print ('   %6i %6i  %f %f' % (jLeft,jRight, u[jLeft], u[jRight-1]))
		# compute cross section on fine grid in center region
		vgtFine[jLeft:jRight] = Voigt (uGrid[jLeft:jRight], positions[l], strengths[l], gammaL[l], gammaD[l])
		# interpolate medium grid cross section near line center
		li = interpolation (vgtMedium[iLeft-1:iRight+2])
		# accumulate fine (near line center only) and medium (everywhere) cross section
		xsFine[jLeft:jRight] = xsFine[jLeft:jRight] + vgtFine[jLeft:jRight] - li[gridRatio:-gridRatio]

	# finally interpolate coarse cross section to medium grid and add to medium
	IXSc  = interpolation(xsCrude)
	XScm  = xsMedium + IXSc[c2fhalf*gridRatio:-c2fhalf*gridRatio]

	# and interpolate coarse cross section to fine grid and add to fine
	IXSm   = interpolation(XScm)
	XS    = xsFine + IXSm[c2fhalf*gridRatio:-c2fhalf*gridRatio]
	
	return uGrid, XS

####################################################################################################################################

def lbl_2grids (positions, strengths, gammaL, gammaD,  vLow, vHigh, lineShape='Voigt',
                sampling=5.0, gridRatio=4, nWidths=25.0, lagrange=3, verbose=False):
	""" Compute lbl cross sections using two grids: fine grid near line center, coarse grid everywhere. """

	# check optional arguments, reset if not given
	if not lineShape.startswith('V'):
		raise SystemExit, 'ERROR --- lbl_2grids:  unknown/unsupported lineShape ' + lineShape

	# check optional spectral range, reset if not given (here just for safety, usually done by calling routine lbl_xs)
	if vLow and vHigh: pass
	elif vLow:         vHigh = positions[-1]
	elif vHigh:        vLow  = positions[0]
	else:              vLow  = positions[0];  vHigh = positions[-1]; #print 'INFO --- lbl_2grids:  vLow, vHigh:', vLow, vHigh

	# Lagrange interpolation method
	interpolation = interpolationDict [(lagrange, gridRatio)]; #print ' lbl_2grids:  ', interpolation.__doc__

	# average line width for Lorentz, Gaussian, and Voigt line shapes
	meanGL, meanGD, meanGV  =  meanLineWidth (gammaL, gammaD)

	# setup fine and coarse wavenumber grids
	nGrids = 2
	vGrid  = wavenumber_grid4xs (vLow, vHigh, meanGV, sampling, gridRatio, nGrids, verbose)
	nv     = len(vGrid)    # number of points of fine grid

	dw    = gridRatio*(vHigh-vLow)/(nv-1)    # coarse resolution grid spacing
	recdw = 1.0/dw
	wGrid = np.arange(vLow-dw,vHigh+1.9*dw,dw) # coarse grid
	nW    = len(wGrid)    # number of points of coarse grid
	iMax  = nW-1          # number of intervals of coarse grid
	if verbose: print '%s %8i%s %8f %8f %8f %8f %s %8f %8f %s %g%s' % ('coarse w  grid: ', 
	              len(wGrid)-1, '+1 points: ', wGrid[0], wGrid[1], wGrid[2], wGrid[3], ' ... ', wGrid[-2], wGrid[-1], '  (delta ', dw,')')

	# allocate coarse and fine cross section arrays and initialize
	xsFine    = np.zeros_like(vGrid)  # cross section on fine grid: sum over all lines
	vgtFine   = np.zeros_like(vGrid)  # cross section on fine grid for a single line
	xsCoarse  = np.zeros_like(wGrid)  # cross section on coarse grid: sum over all lines
	vgtCoarse = np.zeros_like(wGrid)  # cross section on coarse grid for a single line

	# some frequently used variables
	wBegin  = wGrid[0]
	wEnd    = wGrid[-1]
	wLftMin = wGrid[-lagrange]
	wRgtMax = wGrid[lagrange-1]

	# sum over all lines
	for l in xrange(len(positions)):
		#gtCoarse = Voigt_Kuntz_Humlicek1 (wGrid, positions[l], strengths[l], gammaL[l], gammaD[l])
		vgtCoarse = Voigt                 (wGrid, positions[l], strengths[l], gammaL[l], gammaD[l])

		xsCoarse  = xsCoarse + vgtCoarse

		# set limits for line center region
		centerExtension = nWidths*max(gammaL[l],gammaD[l])
		vLeft  = positions[l]-centerExtension
		vRight = positions[l]+centerExtension
		if vRight<wRgtMax or vLeft>wLftMin:  continue   #  next line

		iLeft  = max(int((vLeft-wBegin)*recdw),1)          # left and right coarse grid indices
		iRight = min(iMax-int((wEnd-vRight)*recdw),iMax-1)   
		jLeft  = gridRatio*(iLeft-1)                             # left and right fine grid indices
		jRight = gridRatio*(iRight-1) + 1                        # add 1 because of pythons range end
		# compute cross section on fine grid in center region
		vgtFine[jLeft:jRight] = Voigt (vGrid[jLeft:jRight], positions[l], strengths[l], gammaL[l], gammaD[l])
		# interpolate coarse grid cross section near line center
		vgtCoarseInt = interpolation (vgtCoarse[iLeft-1:iRight+2])  # add 1 because of pythons range end
		# accumulate fine (near line center only) and coarse (over entire line domain) cross section
		xsFine[jLeft:jRight] = xsFine[jLeft:jRight] + vgtFine[jLeft:jRight] - vgtCoarseInt[gridRatio:-gridRatio]

	# finally interpolate coarse cross section to fine grid and add to fine
	xsCoarseInt = interpolation (xsCoarse)
	XS    = xsFine + xsCoarseInt[gridRatio:-gridRatio]

	return vGrid, XS

####################################################################################################################################

def lbl_brute (positions, strengths, gammaL, gammaD, vLow, vHigh, lineShape='Voigt', sampling=5.0, verbose=False):
	""" Compute lbl cross sections the usual way (brute force, no cutoff in wings). """

	# check optional spectral range, reset if not given (here just for safety, usually done by calling routine lbl_xs)
	if vLow and vHigh: pass
	elif vLow:         vHigh = positions[-1]
	elif vHigh:        vLow  = positions[0]
	else:              vLow  = positions[0]; vHigh = positions[-1]; print 'INFO --- lbl_brute:  vLow, vHigh:', vLow, vHigh

	# average line width for Lorentz, Gaussian, and Voigt line shapes
	meanGL, meanGD, meanGV  =  meanLineWidth (gammaL, gammaD)

	if lineShape.startswith('V'):	
		vGrid = wavenumber_grid4xs (vLow, vHigh, meanGV, sampling)
		# allocate cross section array and initialize
		XS = np.zeros(len(vGrid),np.float)
		# sum over all lines
		for l in xrange(len(positions)):
			XS = XS + Voigt(vGrid, positions[l], strengths[l], gammaL[l], gammaD[l])
	elif lineShape.startswith('L'):	
		vGrid = wavenumber_grid4xs (vLow, vHigh, meanGL, sampling)
		# allocate cross section array and initialize
		XS = np.zeros(len(vGrid),np.float)
		# sum over all lines
		for l in xrange(len(positions)):
			XS = XS + Lorentz(vGrid, positions[l], strengths[l], gammaL[l])
	elif lineShape.startswith('G'):	
		vGrid = wavenumber_grid4xs (vLow, vHigh, meanGD, sampling)
		# allocate cross section array and initialize
		XS = np.zeros(len(vGrid),np.float)
		# sum over all lines
		for l in xrange(len(positions)):
			XS = XS + Gauss(vGrid, positions[l], strengths[l], gammaD[l])
	else:
		raise SystemExit, 'ERROR --- lbl_brute:  unknown lineShape ' + lineShape
	return vGrid, XS

####################################################################################################################################

def lbl_xs (Line_Data, Mol_Data,  pressures, temperatures, xLimits=None, lineShape="Voigt",
            sampling=5.0, wingExt=5.0, nGrids=1, gridRatio=1, nWidths=25.0, lagrange=3, verbose=False):
	""" Compute cross sections for a given molecule and some p,T pairs by summation of lines.
	    Returns a list of dictionaries with cross sections given on an equidistant wavenumber grid. """

	# extract line position
	positions = Line_Data['position']

	# if no limits specified, add some points to the left and right of the first/last line
	if xLimits:
		vLow, vHigh = xLimits.limits()
	else:
		vLow    = max(positions[ 0]-0.1*wingExt,0.0)
		vHigh   =     positions[-1]+0.1*wingExt
		xLimits = Interval(vLow,vHigh)
	     	print ' no wavenumber range set, added some points to left and right of the first, last line', xLimits

	# initialize list of cross sections
	nLevels = len(pressures); lvl=0
	crossSections = []
	# loop over all levels or layers
	for p, T in zip(pressures, temperatures):
		lvl = lvl+1
		print '%-10s%-8s lvl %2i/%2i  %20s %10g %s %g\n%50s %10.2f %s %6.2f' % (' lbl_xs:',  Line_Data['molecule'], lvl, nLevels,
		                                                           'pressure [g/cm/s**2]', Line_Data['pressure'], '--->', p,
		                                                           'temperature      [K]', Line_Data['temperature'], ' --->',  T)
		# adjust line parameters to p, T
		strengths, gammaL, gammaD = strengths_and_widths (Line_Data, Mol_Data, p, T)
 
		# compute cross sections by summing profiles of all contributing lines
		tStart = time.clock()
		if nGrids==3:
			v, xs  =  lbl_3grids (positions, strengths, gammaL, gammaD, vLow, vHigh, lineShape,
			                      sampling, gridRatio, nWidths, lagrange, verbose=verbose)
		elif nGrids==2:
			v, xs  =  lbl_2grids (positions, strengths, gammaL, gammaD, vLow, vHigh, lineShape,
			                      sampling, gridRatio, nWidths, lagrange, verbose)
		else:
			v, xs  =  lbl_brute (positions, strengths, gammaL, gammaD, vLow, vHigh, lineShape, sampling, verbose)
		tStop = time.clock(); deltaTime = tStop-tStart

		print '%s %-8s (%i %s, %8gmb, %5.1fK, %10i %.2fsec %.2fns):  %8g %s %8g\n' % \
			(' cross section ', Line_Data['molecule'], len(positions),'lines', unitConversion(p,'p',new='mb'), T, \
			 len(xs), deltaTime, 1.e9*deltaTime/(len(xs)*len(positions)), min(xs),'< xs <',max(xs))

		# brute force remedy to get rid of negative values most likely due to higer order interpolation
		if lagrange>2 and min(xs)<0.0:
			print 'truncate negative xs:', sum(np.where(xs<0,1,0))
			np.putmask(xs, xs<0.0, 0.0)

		# save new cross sections along with attributes in dictionary and append to list of cross sections
		crossSections.append({'molecule': Line_Data['molecule'], 'lineShape': lineShape,  'x': xLimits, 'p': p,  'T': T,  'y': xs})

	return crossSections

####################################################################################################################################

def wavenumber_grid4xs (vLow, vHigh, gammaMean, sampling=5.0,  gridRatio=1, nGrids=1, verbose=0):
	""" Set up equidistant wavenumber grid appropriate to typical line width. """
	# grid point spacing
	dv   = gammaMean/sampling
	# number of grid point intervals
	nv   = int(round( (vHigh-vLow)/dv ))
	
	if nGrids<1:  raise SystemExit, 'ERROR --- wavenumber_grid4xs: nonpositive number of grids!'
	if not gridRatio in [1,2,4,8]:  raise SystemExit, 'ERROR --- wavenumber_grid4xs: invalid grid ratio!'

	# adjust number of grid point intervals to an integer multiple of gridRatio^(nGrids-1)
	mm = gridRatio**(nGrids-1)
        #fgs120310 nv = mm*int(nv/mm)
	nv = mm*max(int(nv/mm),1) # make sure that there are at least mm grid points	

	# adjust spacing to provide an integer number of intervals
	dv   = (vHigh-vLow)/nv
	# set up array of grid points
	v  = np.arange(vLow,vHigh+0.9*dv,dv)
	if verbose: print '%s %8i%s %8f %8f %8f %s %8f %8f %s %g%s %i' % \
	      ('waveumber grid: ', nv, '+1 points: ', v[0], v[1], v[2], ' ... ', v[-2], v[-1], '  (delta ', dv,')',mm)
	return v

####################################################################################################################################

def check_pT (pressures, temperatures):
	""" check presence and consistency of pressure and temperatures. """

	#print 'p', len(pressures), pressures, '\nT', len(temperatures), temperatures
	if   len(pressures)==len(temperatures):
		pass
	elif len(pressures)>len(temperatures)==1:
		temperatures = np.zeros_like(pressures) + temperatures[0]
	elif len(temperatures)>len(pressures)==1:
		pressures = np.zeros_like(temperatures) + pressures[0]
	else:
		raise SystemExit, 'ERROR: invalid/wrong pressure/temperature!!!\n(number of temperatures != number of pressures)'
	return pressures, temperatures

####################################################################################################################################

def _lbl2xs_ (iFile, oFile=None, commentChar='#', molecule=None, pressures=None, temperatures=None, ptFile=None,
              lineShape='Voigt',
              xLimits=None, airWidth=0.1, sampling=5.0, wingExt=5.0, interpolate='3', nGrids=1, gridRatio=8, nWidths=25.0,
	      format=None, verbose=False):
	from cross_section import write_crossSections

	# read a couple of lines and convert to actual pressure and temperature
	Line_Data, Mol_Data = get_lbl_data (iFile, xLimits, airWidth, wingExt, molecule, commentChar)

	if len(Line_Data['position'])<1:
		print 'no line data in ', xLimits; return

	# pressure(s) and temperature(s)
	if ptFile and (pressures or temperatures):
		raise SystemExit, 'either give p,T file or pressure(s) and temperature(s) as command line option'
	elif ptFile:
		try:
                        pT = np.loadtxt(ptFile,comments=commentChar)
		except:
                        raise SystemExit, ' ERROR: reading ' + ptFile + ' failed!'
                pressures      = unitConversion(pT[:,0], 'p', 'mb')
                temperatures   = pT[:,1]
		# some simple tests
		if pT.shape[1]!=2 or min(temperatures)<100 or max(temperatures)>1000:
			print 'WARNING --- check your pT file:  more than 2 columns OR strange (very low, high) temperatures!'
	else:
		# if unset on command line, use line parameter reference values (redo for every line file, ref. values might change!)
		if isinstance(pressures,NoneType):     pressures    = np.array([Line_Data['pressure']])
		else:                                  pressures    = unitConversion(pressures, 'p', 'mb')
		if isinstance(temperatures, NoneType): temperatures = np.array([Line_Data['temperature']])
		pressures, temperatures = check_pT (pressures, temperatures)

	# interpolation method used in multigrid algorithms: only Lagrange!
	if   interpolate in 'bBsScC':  lagrange=4
	elif interpolate in 'qQ':      lagrange=3
	elif interpolate in 'lL':      lagrange=2
	else:                          lagrange=int(interpolate)

	# compute cross sections for selected line shape
	crossSections = lbl_xs (Line_Data, Mol_Data, pressures, temperatures, xLimits, lineShape,
	                        sampling, wingExt, nGrids, gridRatio, nWidths, lagrange, verbose)

	# save cross sections (and lines)
	print ' %-8s %5i %s %s' % (Line_Data['molecule'], len(crossSections), ' cross sections:  save to ', oFile)
	write_crossSections (crossSections, oFile, format, Line_Data['position'], interpolate, commentChar)

####################################################################################################################################

if (__name__ == "__main__"):

	from IO import open_outFile, parse_file_header
	from command_parser import parse_command, standardOptions, multiple_outFiles 

        # parse the command, return (ideally) one file and some options
        opts = standardOptions + [  # h=help, c=commentChar, o=outFile
               {'ID': 'a', 'name': 'airWidth', 'type': FloatType, 'constraint': 'airWidth>0.0'},
               {'ID': 'L', 'name': 'lineShape', 'type': StringType, 'constraint': 'lineShape in ["V","L","G"]', 'default': 'V'},
               {'ID': 'f', 'name': 'format', 'type': StringType,    'constraint': 'format in [" ","a","h","t","xy"]', 'default': ''},
               {'ID': 'i', 'name': 'interpolate', 'type': StringType, 'default': '3', 'constraint': 'interpolate in "234lLqQcCbBsS"'},
               {'ID': 'm', 'name': 'molecule', 'type': StringType},
               {'ID': 'p', 'name': 'pressures', 'type': np.ndarray, 'constraint': 'pressures>0.0'},
               {'ID': 'T', 'name': 'temperatures', 'type': np.ndarray, 'constraint': 'temperatures>0.0'},
               {'ID': 's', 'name': 'sampling', 'type': FloatType, 'constraint': 'sampling>0.0', 'default': 5.0},
               {'ID': 'w', 'name': 'wingExt', 'type': FloatType, 'constraint': 'wingExt>0.0', 'default': 5.0},
	       {'ID': 'x', 'name': 'xLimits', 'type': Interval, 'constraint': 'xLimits.lower>=0.0'},
               {'ID': 'n', 'name': 'nGrids', 'type': IntType, 'default': 3, 'constraint': 'nGrids>0'},
               {'ID': 'g', 'name': 'gridRatio', 'type': IntType, 'default': 8, 'constraint': 'gridRatio in [1,2,4,8]'},
               {'ID': 'W', 'name': 'nWidths', 'type': FloatType, 'default': 25.0,'constraint': 'nWidths>2.0'},
               {'ID': 'pT', 'name': 'ptFile', 'type': StringType},
               {'ID': 'A',  'name': 'approxAsymptotic'},
               {'ID': 'v',  'name': 'verbose'}
               ]
 	
	files, options, commentChar, outFile = parse_command (opts,(1,99))
	outFiles = multiple_outFiles (files, outFile)

	if options.has_key('h'): print __doc__%globals();  raise SystemExit, " end of lbl2xs help"
	options['verbose'] = options.has_key('verbose')
	
	# loop over line parameter files each with data for one molecule
	for iFile,oFile in zip(files,outFiles):
		 print 'lbl2xs', iFile, '-->', oFile
		 _lbl2xs_ (iFile, oFile, commentChar, **options)
