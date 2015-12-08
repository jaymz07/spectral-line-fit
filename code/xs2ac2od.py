import numpy as np

####################################################################################################################################
def integral_acTimesDistance  (absCoeff, sGrid=None, quad='T', mode=0, verbose=0):
	""" Integrate absorption coefficient along path thru atmosphere (starting at 'observer') and return optical depth. """
	if isinstance(sGrid,np.ndarray):  deltaS = sGrid[1:] - sGrid[:-1]
	else:                             deltaS = np.ones(absCoeff.shape[1]-1)

	# first try to import 'nonstandard' quadrature modules and switch to trapez in case of problems
	if quad.lower()=='s':
		try:
			from scipy.integrate import simps
		except ImportError, msg:
			quad='t'; 
			print 9*'WARNING   ' + '\nWARNING   import scipy.integrate failed, using trapez quadrature instead        WARNING\n'+9*'WARNING   '
	elif quad.lower()=='h':
		try:
			from pchi import quad, quads # Piecewise Cubic Hermite Interpolation
		except ImportError, msg:
			quad='t'; 
			print 8*'WARNING   ' + '\nWARNING   import pchip failed, using trapez quadrature instead        WARNING\n'+8*'WARNING   '
	print '\n integral_acTimesDistance:  ', absCoeff.shape, '   mode=', mode, '   quad=', quad

	nv = absCoeff.shape[0]
	if quad.lower()=='s': # Simpson quadrature
		print 8*'WARNING   ' + '\nWARNING   simpson quadrature --- use with care --- can be very slow!  WARNING\n'+8*'WARNING   '
		if mode in 'cdr':
			# cumulative difference optical depth:  integral for all path segments
			optDepth = []
			for i in range(nv):
				ac = absCoeff[i,:]
				od = np.array([simps(ac[:l+1], sGrid[:l+1]) for l in range(1,len(sGrid))])
				optDepth.append(od)
			optDepth = np.column_stack(optDepth).T
			if   mode=='d':
				optDepth = optDepth[:,1:]-optDepth[:,:-1]
			elif mode=='r':
				deltaOptDepth = optDepth[:,1:]-optDepth[:,:-1]
				totOptDepth   = optDepth[:,-1]
				optDepth      = np.column_stack((np.add.accumulate(deltaOptDepth[:,::-1],1),totOptDepth))
		else:
			# total optical depth
			optDepth = np.array([simps(absCoeff[i,:], sGrid) for i in range(nv)])
	elif quad.lower()=='h':
		if mode in 'cdr':
			intervals = np.array([(sGrid[l-1],sGrid[l]) for l in range(1,len(sGrid))])
			# difference optical depth:  integral for all path segments
			optDepth = np.array([quads(sGrid, absCoeff[i,:], intervals)  for i in range(nv)])
			# stepwise accumulated optical depth
			if   mode=='c':  optDepth = np.add.accumulate(optDepth,1)
			elif mode=='r':  optDepth = np.add.accumulate(optDepth[:,::-1],1)
		else:
			# total optical depth
			optDepth = np.array([quad(sGrid, absCoeff[i,:]) for i in range(nv)])
	else: # trapez quadrature
		if mode=='c':
			optDepth = 0.5*deltaS * (absCoeff[:,1:]+absCoeff[:,:-1])
			optDepth = np.add.accumulate(optDepth,1)
		elif mode=='r':
			optDepth = 0.5*deltaS * (absCoeff[:,1:]+absCoeff[:,:-1])
			optDepth = np.add.accumulate(optDepth[:,::-1],1)
		elif mode=='d':
			optDepth = 0.5*deltaS * (absCoeff[:,1:]+absCoeff[:,:-1])
		else:
			# total optical depth
			optDepth = 0.5*np.sum(deltaS*(absCoeff[:,1:] + absCoeff[:,:-1]),1)
	if verbose and mode in 'cdr':
		for l in range(optDepth.shape[1]):  print '%3i %10f   %10g < od < %10g' % (l, 1e-5*sGrid[l], min(optDepth[:,l]),max(optDepth[:,l]))

	return optDepth

####################################################################################################################################

def sum_xsTimesDensity (xsMatrix, densities, commonFreq=True, interpolate='3', verbose=False):
	""" Compute absorption coefficients as product cross section times molecular density, summed over all molecules.
	    xsMatrix is a "matrix" of cross sections with npT=nLevels rows and nGas columns,
	    each element is a dictionary with two entries: wavenumber interval and cross section array. """

	# print information about input data: both matrices should have the same shape
	nMol, npT = np.array(xsMatrix).shape
	nLevels, nGas = densities.shape
	print '\n nLevels, nMol:', nMol, npT, ' xsMatrix  ', nLevels, nGas, ' densities'

	if not (npT,nMol)==(nLevels,nGas):
		raise SystemExit, 'ERROR --- sum_xsTimesDensity: cross section "matrix" and density matrix have different shape!' 
	if verbose:  print ' densities:', densities.shape, densities	

	# four (npT*nGas) matrices with number of cross section values, first and last wavenumber grid point, and grid point spacing
	nvMatrix    = np.array([[len(xs['y']) for xs in xss] for xss in xsMatrix] )
	vLowMatrix  = np.array([[xs['x'].lower for xs in xss] for xss in xsMatrix])
	vHighMatrix = np.array([[xs['x'].upper for xs in xss] for xss in xsMatrix])
	deltaV      = np.array([[xs['x'].size()/(len(xs['y'])-1) for xs in xss] for xss in xsMatrix])
	# interval (vLowMax,vHighMin) defines the wavenumber range common to all data
	vLowMin   = np.min(vLowMatrix);  vLowMax  = np.max(vLowMatrix)
	vHighMin  = np.min(vHighMatrix); vHighMax = np.max(vHighMatrix)
	minDeltaV = np.min(deltaV)
	if verbose:
		print ' nv=len(xs_array)\n', nvMatrix
		print ' deltaV min', minDeltaV, '\n', deltaV
		print ' vLow   min', vLowMin, ' max', vLowMax, '\n', vLowMatrix
		print ' vHigh  min', vHighMin,' max', vHighMax,'\n', vHighMatrix

	if approx(vLowMin,vLowMax) and approx(vHighMin,vHighMax):
		print ' sum_xsTimesDensity:  x limits identical for all xs'
		# select an interpolation method: old and new xGrid have the same limits, new grid is only denser
		from lagrange_interpolation import lagrange2_regularGrid, lagrange3_regularGrid, lagrange4_regularGrid 
		if   interpolate=='3':       lagrange = lagrange3_regularGrid 
		elif interpolate=='4':       lagrange = lagrange4_regularGrid 
		elif interpolate in '2lL':   lagrange = lagrange2_regularGrid # identical results as np.interp = piecewise linear interpolation
		elif interpolate in 'bBsS':  from scipy.interpolate import splrep, splev
		else:                        raise SystemExit, 'invalid interpolation scheme' + repr(interpolate)

		nx    = np.max(nvMatrix)
		vGrid = np.linspace(vLowMax, vHighMin, nx)
		absCo = np.zeros((nx, npT));  print ' initialized absCo ', absCo.shape
		for m,xss in enumerate(xsMatrix): # loop over molecules (xs files)
			for l,xs in enumerate(xss): # loop over levels/layers (press, temp pairs)
				vGridLimits = xs['x']
				yValues     = xs['y']
				if verbose:  print m,l, densities[l,m], vGridLimits, len(yValues), min(yValues), max(yValues),
				if interpolate in '234':
					yyy = lagrange (yValues, nx)
				else:
					vOld  = np.linspace(vGridLimits.lower, vGridLimits.upper, len(yValues))
					vGrid = np.linspace(vLowMax, vHighMin, nx)
					tck   = splrep (vOld, yValues)
					yyy   = splev(vGrid,tck)
				absCo[:,l] = absCo[:,l] + densities[l,m]*yyy
				if verbose:  print ' --->', min(absCo[:,l]), max(absCo[:,l])
	elif vLowMax<vHighMin:
	    if commonFreq: # sum cross sections to abs.coefficients only in wavenumber interval common to all data
	     	print ' sum_xsTimesDensity:  x limits xLowMax<vHighMin'
		nx    = int((vHighMin-vLowMax)/minDeltaV)
		vGrid = np.linspace(vLowMax, vHighMin, nx+1)
		absCo = np.zeros((nx+1, npT));  print ' initialized absCo ', absCo.shape, ' for v', vLowMax, vHighMin
		for m,xss in enumerate(xsMatrix): # loop over molecules (xs files)
			for l,xs in enumerate(xss): # loop over levels/layers (press, temp pairs)
				vOldLimits = xs['x']; yValues=xs['y']
				vOld = np.linspace(vOldLimits.lower, vOldLimits.upper, len(yValues))
				if verbose:  print m,l, densities[l,m], vOldLimits, len(yValues), min(yValues), max(yValues),
				yyy = np.interp (vGrid, vOld, yValues) # piecewise linear interpolation
				# accumulate absorption coefficient
				absCo[:,l] = absCo[:,l] + densities[l,m]*yyy
				if verbose:  print ' --->', min(absCo[:,l]), max(absCo[:,l])
	    else:
	     	print 'sum_xsTimesDensity:  x limits xLowMax<vHighMin, extrapolate to 0 outside'
		nx    = int((vHighMax-vLowMin)/minDeltaV) # number of grid intervals!!!
		vGrid = np.linspace(vLowMin, vHighMax, nx+1)
		absCo = np.zeros((nx+1, npT));  print ' initialized absCo ', absCo.shape, ' for v', vLowMin, vHighMax, len(vGrid)
		for m,xss in enumerate(xsMatrix): # loop over molecules (xs files)
			for l,xs in enumerate(xss): # loop over levels/layers (press, temp pairs)
				vOldLimits = xs['x']; yValues=xs['y']
				vOld = np.linspace(vOldLimits.lower, vOldLimits.upper, len(yValues))
				if verbose:  print m,l, densities[l,m], vOldLimits, len(yValues), min(yValues), max(yValues),
				# piecewise linear interpolation (with y set to zero outside given domain) and accumulate
				#?# yyy = np.interp (vGrid, vOld, yValues, 0.0, 0.0)
				#?# absCo[:,l] = absCo[:,l] + densities[l,m]*yyy
				# indices of first and last old wavenumber grid point w.r.t. new grid
				iLow  = int((vOldLimits.lower-vLowMin)/minDeltaV) + 1 # add one to avoid extrapolation
				iHigh = int((vOldLimits.upper-vLowMin)/minDeltaV)
				yyy = np.interp (vGrid[iLow:iHigh+1], vOld, yValues, nx) # piecewise linear interpolation
				## # accumulate absorption coefficient
				absCo[iLow:iHigh+1,l] = absCo[iLow:iHigh+1,l] + densities[l,m]*yyy
				if verbose:  print ' --->', min(absCo[:,l]), max(absCo[:,l])
	else:
		raise SystemExit, 'ERROR --- sum_xsTimesDensity:  wavenumber grids do not overlap (no common subinterval)'

	if not verbose:
		print '\n Absorption coefficient for ', nMol, ' molecules and', npT, ' levels (or layers):'
		frmt = '   <ac> = %9.3g    %10.3g < ac < %9.3g  @ %10f'
		for l in range(npT): print '%3i' % (l+1), nMol*'%10.3g ' % tuple(densities[l,:]), \
		                            frmt % (np.mean(absCo[:,l]), min(absCo[:,l]), max(absCo[:,l]), vGrid[np.argmax(absCo[:,l])])
	return vGrid, absCo

####################################################################################################################################

def approx (self,other,eps=0.001):
	return abs(self-other) < eps*self
