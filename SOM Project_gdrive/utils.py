#from brian import *
#from numpy import *
from pylab import *

def cos_weight_vector(N2, cycles):
	"""
	Just for 1 cosine distribution, 1 cell to N2 cells
	"""
	x = zeros(N2)
	for i in arange(N2):
		x[i] = sin((cycles*2*pi/N2)*i)
	return x

def cos_weight_matrix(N1, N2, cycles):
	"""
	cos_weight_matrix
	@param: N1 number of neurons in layer 1
	@param: N2 number of neurons in layer 2
	@param: cycles width of cos func.
	@return: N1xN2 weight matrix
	"""

	x = zeros((N1,N2))
	j = arange(N2)
	
	for i in range(N1):
		a = cycles * (i * 2 * pi / N1 - j * 2 * pi / N2)
		a[ logical_and(a > pi/2, a < (2*pi*cycles - pi/2) )] = pi/2
		a[ logical_and(a < -pi/2, a > -(2*pi*cycles - pi/2) )] = -pi/2
		x[i, :] = abs( cos(a) )
	return x

def inh_cos_weight_matrix(N1, N2, cycles):
	"""
	inh_cos_weight_matrix
	@param: N1 number of neurons in layer 1
	@param: N2 number of neurons in layer 2
	@param: cycles width of cos func.
	@return: N1xN2 weight matrix shifted 1/2 N1
	"""

	x = zeros((N1,N2))
	j = arange(N2)
	
	for i in range(N1):
		a = cycles * ((i-N1/2) *2 * pi / N1 - j * 2 * pi / N2)
		a[ logical_and(a > pi/2, a < (2*pi*cycles - pi/2) )] = pi/2
		a[ logical_and(a < -pi/2, a > -(2*pi*cycles - pi/2) )] = -pi/2
		x[i, :] = cos(a)
	return x

def cos_weight(deg, rf, cycles=1):
    """
    converts difference of angles to cosine weights. This is the full cosine. 
    
    @param: deg the angle of inputs.
    @param: rf the preferred angles.
    @param: cycles the width of the cosine distribution.
    @return: weights the cosine weights of the relative angles.
    """
    
    weights = zeros((len(rf), len(deg)))
    
    for i in range(len(rf)):
        diff = 2 * pi * abs((rf[i] - deg)/360.0)
        # Wrap the angles to [0 2 pi] by doing this. Can't find a better function...
        a = arccos(cos(diff))
        a = cycles * a
        a[a > pi] = pi
        a[a < -pi] = -pi
        weights[i, :] = 0.5 * (cos(a) + 1)
        
    return weights


def touch_to_spikes(deg, N, cycles=1):
    """
    converts an angle of touch to a rate for cosine receptive fields.
    @param: deg the angle of stimulation in degrees
    @param: N the number of receptive fields to partition.
    @param: cycles the width of the cos function.
    """
    j = ones(1)
    cells = zeros((N, len(deg)))
	
    for i in range(N):
        a = cycles * (i *2 * pi / N - j * deg * 2 * pi / 360 )
        a[ logical_and(a > pi/2, a < (2*pi*cycles + pi/2) )] = pi/2
        a[ logical_and(a < -pi/2, a > -(2*pi*cycles + pi/2) )] = -pi/2
        cells[i] = cos(a)[0]
    
    return cells
    
def angle2color(deg, zangle=0):
    """
    angle2color(deg): returns an hsv color for a given angle.
    @param: deg the angle in degrees.
    @param: zangle relative angle that indicates 0.
    @return: col the element vector of the colors.
    """
    
    H = mod(deg - zangle, 360.) / 360.
    
    S = ones_like(H)
    V = ones_like(H)
    
    col = matplotlib.colors.hsv_to_rgb(dstack((H,S,V)))
    
    return col.reshape((-1,3))