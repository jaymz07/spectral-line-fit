####################################################################################################################################

class PairOfInts:
	def __init__ (self, *pair):
		try:
			left, right = pair
			if isinstance(left,(float,int)) and isinstance(right,(float,int)):
				self.left  = left
				self.right = right
			else:
				print "ERROR:  PairOfInts must be floats or integers!"
		except ValueError:
			raise SystemExit, 'ERROR:  PairOfInts initialization requires a pair of floats (or integers)!'
		except:
			raise SystemExit, "ERROR:  PairOfInts initialization failed!"
	def __str__ (self):
		return  'PairOfInts [%g,%g]' % (self.left,self.right) 
	def __repr__ (self):
		return 'PairOfInts(%s,%s)' % (self.left,self.right) 
	def __add__ (self,other):
		if isinstance(other,PairOfInts):
			return PairOfInts(self.left+other.left, self.right+other.right)
		elif isinstance(other,(float,int)):
			return PairOfInts(self.left+other,self.right+other)
		else:
			print 'ERROR: other is not a PairOfInts (instance) or an integer or float!'
	def __radd__ (self, other):
		return PairOfInts(self.left+other, self.right+other)
	def __sub__ (self,other):
		if isinstance(other,PairOfInts):
			return PairOfInts(self.left-other.left, self.right-other.right)
		elif isinstance(other,(float,int)):
			return PairOfInts(self.left-other,self.right-other)
		else:
			print 'ERROR: other is not a PairOfInts (instance) or an integer or float!'
	def __rsub__ (self, other):
		return PairOfInts(other-self.left, other-self.right)
	def __mul__ (self, other):
		if isinstance(other,PairOfInts):
			return PairOfInts(self.left*other.left, self.right*other.right)
		#fgs150310 elif isinstance(self,PairOfInts) and type(other) in [FloatType,IntType]:
		elif isinstance(self,PairOfInts) and isinstance(other,(float,int)):
			return PairOfInts(other*self.left, other*self.right)
		else:
			print 'ERROR: other is not a PairOfInts (instance) or an integer or float!'
	def __rmul__ (self, other):
		return PairOfInts(other*self.left, other*self.right)
	def __div__ (self, other):
		if isinstance(other,PairOfInts):
			return PairOfInts(self.left/other.left, self.right/other.right)
		# fgs150310 elif isinstance(self,PairOfInts) and type(other) in [FloatType,IntType]:
		elif isinstance(self,PairOfInts) and isinstance(other,(float,int)):
			return PairOfInts(self.left/other, self.right/other)
		else:
			print 'ERROR: other is not a PairOfInts (instance) or an integer or float!'
	def __rdiv__ (self, other):
		return PairOfInts(float(other)/self.left, float(other)/self.right)
	def __eq__ (self,other):
		if isinstance(other,PairOfInts):
			return self.left==other.left and self.right==other.right
		else:
			print 'ERROR: other is not a PairOfInts (instance)!'
	def list (self):
		return [self.left,self.right] 

####################################################################################################################################

class PairOfFloats:
	def __init__ (self, *pair):
		try:
			left, right = pair
			if isinstance(left,(float,int)) and isinstance(right,(float,int)):
				self.left  = left
				self.right = right
			else:
				print "ERROR:  PairOfFloats must be floats or integers!"
		except ValueError:
			raise SystemExit, 'ERROR:  PairOfFloats initialization requires a pair of floats (or integers)!'
		except:
			raise SystemExit, "ERROR:  PairOfFloats initialization failed!"
	def __str__ (self):
		return  'PairOfFloats [%g,%g]' % (self.left,self.right) 
	def __repr__ (self):
		return 'PairOfFloats(%s,%s)' % (self.left,self.right) 
	def __add__ (self,other):
		if isinstance(other,PairOfFloats):
			return PairOfFloats(self.left+other.left, self.right+other.right)
		elif isinstance(other,(float,int)):
			return PairOfFloats(self.left+other,self.right+other)
		else:
			print 'ERROR: other is not a PairOfFloats (instance) or an integer or float!'
	def __radd__ (self, other):
		return PairOfFloats(self.left+other, self.right+other)
	def __sub__ (self,other):
		if isinstance(other,PairOfFloats):
			return PairOfFloats(self.left-other.left, self.right-other.right)
		elif isinstance(other,(float,int)):
			return PairOfFloats(self.left-other,self.right-other)
		else:
			print 'ERROR: other is not a PairOfFloats (instance) or an integer or float!'
	def __rsub__ (self, other):
		return PairOfFloats(other-self.left, other-self.right)
	def __mul__ (self, other):
		if isinstance(other,PairOfFloats):
			return PairOfFloats(self.left*other.left, self.right*other.right)
		#fgs150310 elif isinstance(self,PairOfFloats) and type(other) in [FloatType,IntType]:
		elif isinstance(self,PairOfFloats) and isinstance(other,(float,int)):
			return PairOfFloats(other*self.left, other*self.right)
		else:
			print 'ERROR: other is not a PairOfFloats (instance) or an integer or float!'
	def __rmul__ (self, other):
		return PairOfFloats(other*self.left, other*self.right)
	def __div__ (self, other):
		if isinstance(other,PairOfFloats):
			return PairOfFloats(self.left/other.left, self.right/other.right)
		#fgs150310 elif isinstance(self,PairOfFloats) and type(other) in [FloatType,IntType]:
		elif isinstance(self,PairOfFloats) and isinstance(other,(float,int)):
			return PairOfFloats(self.left/other, self.right/other)
		else:
			print 'ERROR: other is not a PairOfFloats (instance) or an integer or float!'
	def __rdiv__ (self, other):
		return PairOfFloats(float(other)/self.left, float(other)/self.right)
	def __eq__ (self,other):
		if isinstance(other,PairOfFloats):
			return self.left==other.left and self.right==other.right
		else:
			print 'ERROR: other is not a PairOfFloats (instance)!'
	def approx (self,other,eps=0.001):
		if isinstance(other,PairOfFloats):
			return abs(self.left-other.left)<eps*self.left and abs(self.right-other.right)<eps*self.right
		else:
			print 'ERROR: other is not a PairOfFloats (instance)!'
	def list (self):
		return [self.left,self.right] 

####################################################################################################################################

class Interval:
	def __init__ (self, *limits):
		try:
			lower, upper = limits
			if isinstance(lower,(float,int)) and isinstance(upper,(float,int)):
				self.lower = min(lower,upper)
				self.upper = max(lower,upper)
			else:
				print "ERROR:  Interval bounds must be floats or integers!"
		except ValueError:
			raise SystemExit, 'ERROR:  Interval initialization requires a pair of floats (or integers)!'
		except:
			raise SystemExit, "ERROR:  Interval initialization failed!"
	def __str__ (self):
		return  'Interval [%g,%g]' % (self.lower,self.upper) 
	def __repr__ (self):
		return 'Interval(%s,%s)' % (self.lower,self.upper) 
	def limits (self):
		return self.lower, self.upper
	def contains (self, other):
		if isinstance(Interval):  return self.lower <= other.lower <= other.upper <= self.upper
	def member (self, value):
		if isinstance(value,(float,int)):
			return self.lower <= value <= self.upper
		print "ERROR:  value must be an integer or float!"
	def inside (self, value):
			return self.lower < value < self.upper
			#raise SystemExit, "ERROR:  value must be an integer or float!"
	def intersect (self, other):
		return Interval(max(self.lower,other.lower),
		                min(self.upper,other.upper))
	def __add__ (self,other):
		if isinstance(other,Interval):
			if min(self.upper,other.upper)>max(self.lower,other.lower):
				return Interval(min(self.lower,other.lower),
				                max(self.upper,other.upper))
			else:
				return
		elif isinstance(other,(float,int)):
			if other>=0.0:
				return Interval(self.lower-other,self.upper+other)
			else:
				print 'ERROR:  positive float or integer required!'
		else:
			print 'ERROR: other is not an Interval (instance)!'
	def __radd__ (self, other):
		return Interval(self.lower-other, self.upper+other)
	def shift (self, other):
		if isinstance(other,(float,int)):
			return Interval(self.lower+other, self.upper+other)
		else:
			print 'ERROR: argument of Interval shift must be float or integrer!'
	def __sub__ (self,other):
		if isinstance(other,Interval):
			if min(self.upper,other.upper)>max(self.lower,other.lower):
				return Interval(max(self.lower,other.lower),
				                min(self.upper,other.upper))
			else:
				return
		elif isinstance(other,(float,int)):
			if self.upper-self.lower>2.*other: return Interval(self.lower+other,self.upper-other)
		else:
			print 'ERROR: other is not an Interval (instance)!'
	def __eq__ (self,other):
		if isinstance(other,Interval):
			return self.lower==other.lower and self.upper==other.upper
		else:
			print 'ERROR: other is not an Interval (instance)!'
	def approx (self,other,eps=0.001):
		if isinstance(other,Interval):
			return abs(self.lower-other.lower)<eps*self.lower and abs(self.upper-other.upper)<eps*self.upper
		else:
			print 'ERROR: other is not an Interval (instance)!'
	def __len__ (self):
		return int(self.upper-self.lower)
	def __nonzero__ (self):
		return bool(self.upper-self.lower)
	def size (self):
		return float(self.upper-self.lower)
	def __cmp__(self, other):
		return cmp(self.upper-self.lower,other.upper-other.lower)
	def __mul__ (self, other):
		#fgs150310 if isinstance(self,Interval) and type(other) in [FloatType,IntType]:
		if isinstance(self,Interval) and isinstance(other,(float,int)):
			return Interval(other*self.lower, other*self.upper)
		else:
			print 'ERROR: no multiplication of Intervals!'
	def __rmul__ (self, other):
		return Interval(other*self.lower, other*self.upper)
	def __div__ (self, other):
		#fgs150310 if isinstance(self,Interval) and type(other) in [FloatType,IntType]:
		if isinstance(self,Interval) and isinstance(other,(float,int)):
			return Interval(self.lower/other, self.upper/other)
		else:
			print 'ERROR: no division of Intervals!'
	def __rdiv__ (self, other):
		return Interval(other/self.upper, other/self.lower)

