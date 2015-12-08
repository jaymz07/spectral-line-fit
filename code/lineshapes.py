"""
Line profile functions normalized to 1.0 with wavenumber (array) and position, strengths, and width(s) as input variables.
"""

from math import pi, sqrt, log
# some math constants
recSqrtPi  = 1./sqrt(pi)
ln2     = log(2.0)
sqrtLn2 = sqrt(ln2)
sqrtPi  = sqrt(pi)
sqrtLn2overPi  = sqrtLn2/sqrtPi

import numpy as np

####################################################################################################################################

def voigtWidth (gammaL=1.0, gammaD=1.0):
	""" Half width half maximum (HWHM) of Voigt profile (Whiting's approximation). """
	return 0.5*(gammaL+sqrt(gammaL**2+4.0*gammaD**2)) 

####################################################################################################################################

def Lorentz (v, v0=0.0, S=1.0, gammaL=1.0):
	""" Pressure broadening: Lorentz profile normalized to one, multiplied with line strength. """
	return (S*gammaL/pi)/((v-v0)**2+gammaL**2)

####################################################################################################################################

def Gauss (v, v0=0., S=1.0, gammaD=1.0):
	""" Doppler broadening:  Gauss profile normalized to one, multiplied with line strength. """
	t = ln2*((v-v0)/gammaD)**2
	# don't compute exponential for very large numbers > log(1e308)=709.
	# (in contrast to python.math.exp the numpy exp overflows!!!)
	t = np.where(t>200.,200.0,t)
	return S * sqrtLn2/(sqrtPi*gammaD) * np.exp(-t)

####################################################################################################################################

def Voigt (v, v0=0.0, S=1.0, gammaL=1.0, gammaD=1.0):
	""" Voigt profile normalized to one, multiplied with line strength. """
	x = sqrtLn2 * (v-v0)/ gammaD
	y = sqrtLn2 * gammaL / gammaD
       	#K, L = hui6r (x,y)
	K = hum1wei24(x,y).real 
	vgt = (sqrtLn2overPi/gammaD)*S*K
	return vgt

####################################################################################################################################

def Voigt_Kuntz_Humlicek1 (v, v0=0.0, S=1.0, gammaL=1.0, gammaD=1.0):
	""" Voigt profile normalized to one, multiplied with line strength.
	    Rational approximation for asymptotic region for large |x|+y.
	    Real part only using Kuntz (JQSRT, 1997) implementation with Ruyten's (JQSRT, 2003) correction. """
	xx = (sqrtLn2 * (v-v0)/ gammaD)**2
	y  = sqrtLn2 * gammaL / gammaD
	yy = y*y
	kha1 = 0.5+yy
	kha2 = kha1 + yy*yy; khb2 = 2.*yy-1.0
	return  (sqrtLn2overPi/gammaD)*S * y*recSqrtPi * (kha1+xx) / (kha2 + xx*(khb2+xx))

####################################################################################################################################

# constants for Hui's rational approximations
a6 = np.array([ 122.607931777104326, 214.382388694706425, 181.928533092181549,
                  93.155580458138441,  30.180142196210589,   5.912626209773153, 0.564189583562615])
b6 = np.array([ 122.607931773875350, 352.730625110963558, 457.334478783897737,
                 348.703917719495792, 170.354001821091472,  53.992906912940207, 10.479857114260399])

def hui6r (x,y):
	"""
	Hui rational function approximation for complex error function w(z)=Kx,y)+iL(x,y) with z=x+i*y.
   
	K(x,y)  =  y/pi * integral exp(-t*t)/((x-t)*(x-t)+y*y) dt
	L(x,y)  =  1/pi * integral (x-t)*exp(-t*t)/((x-t)*(x-t)+y*y) dt
	(integration interval (-infinity,+infinity)

	rational approximation with 6,7 coefficients:   w = K + iL = P / Q
	where P and Q are complex polynomials of order 6 and 7, respectively

	Hui, Armstrong, Wray:
	Rapid computation of the Voigt and complex error function
	JQSRT 19, pp.509-516 (1978)

	Optimization as discussed in JQSRT 48, pp. 743--762, 1992
	"""

	# coefficients Re(P)
	cr0 = (((((a6[6]*y+a6[5])*y+a6[4])*y+a6[3])*y+a6[2])*y+a6[1])*y+a6[0]
	cr2 = -((((15.*a6[6]*y +10.*a6[5])*y +6.*a6[4])*y +3.*a6[3])*y +a6[2])
	cr4 = (15.*a6[6]*y +5.*a6[5])*y +a6[4]
	cr6 = -a6[6]
	# coefficients Im(P)
	ci1 = -(((((6.*a6[6]*y +5.*a6[5])*y +4.*a6[4])*y +3.*a6[3])*y +2.*a6[2])*y+a6[1])
	ci3 = ((20.*a6[6]*y +10.*a6[5])*y +4.*a6[4])*y +a6[3]
	ci5 = -(6.*y*a6[6] +a6[5])
	# coefficients Re(Q)
	dr0 = ((((((y+b6[6])*y +b6[5])*y +b6[4])*y +b6[3])*y +b6[2])*y +b6[1])*y +b6[0]
	dr2 = -(((((21.*y +15.*b6[6])*y +10.*b6[5])*y +6.*b6[4])*y +3.*b6[3])*y+b6[2])
	dr4 = ((35.*y +15.*b6[6])*y +5.*b6[5])*y +b6[4]
	dr6 = -(7.*y +b6[6])
	# coefficients Im(Q)
	di1 = -((((((7.*y+6.*b6[6])*y+5.*b6[5])*y+4.*b6[4])*y+3.*b6[3])*y+2.*b6[2])*y+b6[1])
	di3 = (((35.*y +20.*b6[6])*y +10.*b6[5])*y +4.*b6[4])*y +b6[3]
	di5 = -((21.*y +6.*b6[6])*y +b6[5])
	#
        xx = x*x
	# nominator    P(z) = a6*z**6 + ... + a1*z + a0  with z=y-i*x
        pRe = ( (cr6*xx+cr4) *xx+cr2) *xx+cr0
        pIm = ( (ci5*xx+ci3) *xx+ci1) *x
	# denominator  Q(z) = z**7 + b6*z**6 + ... + b1*z + b0
        qRe =  ((dr6*xx+dr4)*xx+dr2)*xx+dr0
        qIm = (((xx+di5)*xx+di3)*xx+di1)*x
	
	# inverse of squared denominator |Q|**2
        q2inv = 1./ (qRe*qRe + qIm*qIm)
	# real and imaginary part of P/Q
        vgtK =  q2inv * (pRe*qRe + pIm*qIm)
        vgtL =  q2inv * (pIm*qRe - pRe*qIm)
	#
	return vgtK, vgtL

####################################################################################################################################

wB = np.array([  0.000000000000e+00,  -1.513746165453e-10,   4.904821733949e-09,   1.331046162178e-09,  -3.008282275137e-08,
                -1.912225850812e-08,   1.873834348705e-07,   2.568264134556e-07,  -1.085647578974e-06,  -3.038893183936e-06,
                 4.139461724840e-06,   3.047106608323e-05,   2.433141546260e-05,  -2.074843151142e-04,  -7.816642995614e-04,
                -4.936426901280e-04,   6.215006362950e-03,   3.372336685532e-02,   1.083872348457e-01,   2.654963959881e-01,
                 5.361139535729e-01,   9.257087138589e-01,   1.394819673379e+00,   1.856286499206e+00,   2.197858936532e+00])

wL24 = sqrt(24.)/2**0.25 # sqrt(N)/2^(1/4)

def hum1wei24 (x,y):
 	""" Complex error function using Humlicek's and Weideman's rational approximation:

	    |x|+y>15:  Humlicek (JQSRT, 1982) rational approximation for region I;
	    else:      J.A.C. Weideman: SIAM-NA 1994); equation (38.I) and table I.
	    
	    F. Schreier, JQSRT 112, pp, 1010-1025, 2011:  doi: 10.1016/j.jqsrt.2010.12.010 """
	t   = y - 1j*x
	w = t * recSqrtPi / (0.5 + t*t)  # Humlicek (1982) approx 1 for s>15
	if y<15.0:
		mask = abs(x)+y<15.       # returns true for interior points
		iz  = -t[np.where(mask)]  # returns small complex array covering only the interior region
		# the following five lines are only evaluated for the interior grid points
		lpiz = wL24 + iz  # wL - y + x*complex(0.,1.)
		lmiz = wL24 - iz  # wL + y - x*complex(0.,1.)
		recLmiZ  = 1.0 / lmiz
		r       = lpiz * recLmiZ
		w24 = recLmiZ * (recSqrtPi + 2.0*recLmiZ*(wB[24]+(wB[23]+(wB[22]+(wB[21]+(wB[20]+(wB[19]+(wB[18]+(wB[17]+(wB[16]+(wB[15]+(wB[14]+(wB[13]+(wB[12]+(wB[11]+(wB[10]+(wB[9]+(wB[8]+(wB[7]+(wB[6]+(wB[5]+(wB[4]+(wB[3]+(wB[2]+wB[1]*r)*r)*r)*r)*r)*r)*r)*r)*r)*r)*r)*r)*r)*r)*r)*r)*r)*r)*r)*r)*r)*r)*r))
		# replace asympotitic Humlicek approximation in interior center region
		np.place(w, mask, w24)
	return w

