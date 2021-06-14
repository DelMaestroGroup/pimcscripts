"""
pyutils.py

Useful python routines that we have come up with over the years.

JBLux, AGDelma 2009
"""

from __future__ import print_function 
import re
import gzip
import math
import sys

#TODO - add a dump file method that is complementary

# ----------------------------------------------------------------------------
def now():
    ''' Print the current date and time.'''
    from datetime import datetime
    print(datetime.now().strftime("Last Updated %Y/%m/%d at %H:%M:%S"))

# ----------------------------------------------------------------------------
def linear(r1,r2):
    ''' Return a function that performs a linear interpolation between two
    points. '''
    m = (r1[1]-r2[1])/(r1[0]-r2[0])
    b = r1[1] - m*r1[0]

    return lambda x: m*x + b

# ----------------------------------------------------------------------------
def loadFile(filename,format='float',comment='#'):
    """ Loads filename into an array of data with same dimensions 
    as rows and columns of file and ignores comments.  
    Reads in all numbers as the specified format (string) 
    (int : integer, float : float, scientific : scientific notation).
    """
    
    data_array = []
    
    if format == 'int':
        type_converter = lambda x: int(x)
    elif (format == 'float' or format == 'scientific'):
        type_converter = lambda x: float(x) 
    
    if re.match('.*\.gz',filename):
        f = gzip.GzipFile(filename,'r')
    else:
        f = open(filename,'r')
        
    for line in f.readlines():
        if line[0] != comment:
            numbers = line.split()
            if len(numbers) == 1:
                data_array.append([type_converter(numbers[0])])
            else:
                data_array.append([type_converter(num) for num in numbers])
    f.close()
    return data_array

# ----------------------------------------------------------------------------
def getDimensions(data):
    """ Get the dimensions of a matrix. """
    N = []
    N.append(len(data))
    try:
        N.append(len(data[0]))
    except:
        N.append(0)

    return N

# ----------------------------------------------------------------------------
def vector2Matrix(vec):
    """ Convert a numpy array to a matrix. """
    import numpy as n
    mat = n.zeros([1,len(vec)],float)
    mat[0,:] = vec
    return mat

# ----------------------------------------------------------------------------
def smooth(x, window_len=10, window='hanning'):
    """smooth the data using a window with requested size.
    
    from: http://scipy.org/Cookbook/SignalSmooth

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    import numpy as np    
    t = np.linspace(-2,2,0.1)
    x = np.sin(t)+np.random.randn(len(t))*0.1
    y = smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    import numpy as np
    s=np.r_[2*x[0]-x[window_len:1:-1], x, 2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))
    
    if window == 'flat': #moving average
        w = np.ones(window_len,'d')
    else:
        w = getattr(np, window)(window_len)
    y = np.convolve(w/w.sum(), s, mode='same')
    return y[window_len-1:-window_len+1]

# ----------------------------------------------------------------------------
def extrema(x, max = True, min = True, strict = False, withend = False):
    """
    This function will index the extrema of a given array x.
    
    Options:
        max         If true, will index maxima
        min         If true, will index minima
        strict      If true, will not index changes to zero gradient
        withend     If true, always include x[0] and x[-1]
    
    This function will return a tuple of extrema indexies and values
    """
    
    # This is the gradient
    from numpy import zeros
    dx = zeros(len(x))
    from numpy import diff
    dx[1:] = diff(x)
    dx[0] = dx[1]
    
    # Clean up the gradient in order to pick out any change of sign
    from numpy import sign
    dx = sign(dx)
    
    # define the threshold for whether to pick out changes to zero gradient
    threshold = 0
    if strict:
        threshold = 1
        
    # Second order diff to pick out the spikes
    d2x = diff(dx)
    
    if max and min:
        d2x = abs(d2x)
    elif max:
        d2x = -d2x
    
    # Take care of the two ends
    if withend:
        d2x[0] = 2
        d2x[-1] = 2
    
    # Sift out the list of extremas
    from numpy import nonzero
    ind = nonzero(d2x > threshold)[0]
    
    return ind, x[ind]


# ----------------------------------------------------------------------------
def average(data,dim=1):

    """ Average a data matrix over the dimension given in dim. 
    
    Only works with 2d lists. dim=0 means average over rows, while dim=1 means
    average over columns (default)."""

    # We first determine whether or not we are dealing with a 1d or 2d python list
    N = getDimensions(data)

    if N[dim] == 0:
        raise NameError('Dimension to average over does not match with dimension of array!')

    ave = []
    if dim == 0:
        for i in range(N[0]):
            temp = 0.0
            for j in range(N[1]):
                temp += data[i][j]
            ave.append(temp/(1.0*N[1]))
    elif dim == 1:
        for j in range(N[1]):
            temp = 0.0
            for i in range(N[0]):
                temp += data[i][j]
            ave.append(temp/(1.0*N[0]))

    return ave

# ----------------------------------------------------------------------------
def error(data,dim=1):
    """ Computes the error in a data matrix over the dimension given in dim. 
    
    Only works with 2d lists. dim=0 means average over rows, while dim=1 means
    average over columns (default)."""

    # we get the dimensions of our data matrix
    N = getDimensions(data)

    if N[dim] == 0:
        raise NameError('Dimension to average over does not match with dimension of array!')

    # Create the data^2 matrix
    data2 = []
    for i in range(N[0]):
        data2.append([d*d for d in data[i]])

    # now we get the averages of data and data2
    ave  = average(data,dim)
    ave2 = average(data2,dim)

    # try and get the error matrix
    norm = N[not dim] - 1
    if norm == 0:
        err = [0.0 for n in range(len(ave))]
    else:
        err = [math.sqrt(abs(ave2[n]-ave[n]*ave[n])/(1.0*norm)) for n in range(len(ave))]

    return err

# ----------------------------------------------------------------------------
def bootstrap(data,dim=1):
    """ Computes the bootstrap error in a data matrix over the dimension given in dim. 
    
    Only works with 2d lists."""

    import random

    # the number of bootstrap samples
    numBoot = 500

    # we get the dimensions of our data matrix
    N = getDimensions(data)

    # make sure we have at least 2 elements of data
    if N[not dim] > 1:
        sampleAve = []

        # we need to compute things differently based on which dimension
        # the average is being performed over

        # Average over rows (non-standard)
        if dim == 0:
            for i in range(N[0]):
                ave = []
                for b in range(numBoot):
                    temp = 0.0
                    for j in range(N[1]):
                        k = random.randint(0,N[1]-1)
                        temp += data[i][k]
                    ave.append(temp/(1.0*N[1]))
                sampleAve.append(ave)

        # average over columns (standard)
        elif dim == 1:
            for j in range(N[1]):
                ave = []
                for b in range(numBoot):
                    temp = 0.0
                    for i in range(N[0]):
                        k = random.randint(0,N[0]-1)
                        temp += data[k][j]
                    ave.append(temp/(1.0*N[0]))
                sampleAve.append(ave)

        sampleAve2 = []
        for i in range(len(sampleAve)):
            sampleAve2.append([d*d for d in sampleAve[i]])

        bootAve = average(sampleAve,0)
        bootAve2 = average(sampleAve2,0)

        norm = 1.0*numBoot/(1.0*(numBoot-1))
        bootErr = [math.sqrt(norm*abs(bootAve2[n]-bootAve[n]*bootAve[n])) for n in range(N[dim])]

    else:
        bootErr= [0.0 for n in range(N[dim])]

    return bootErr

# ----------------------------------------------------------------------------
def isList(x):
    '''Determine if x is a list.'''
    try:
        len(x)
        return True
    except:
        return False
#    return hasattr(x,'__iter__')

# ----------------------------------------------------------------------------
def DedekindEta(tau):
    ''' Compute the Dedekind Eta function for imaginary argument tau. 
        See: http://mathworld.wolfram.com/DedekindEtaFunction.html. '''
    try:
        import mpmath as mp
        mpmath_loaded = True
    except ImportError:
        mpmath_loaded = False 
    qbar = mp.exp(-2.0*mp.pi*tau)

    return mp.nsum(lambda n: ((-1)**n)*(qbar**((6*n-1)*(6.0*n-1)/24)),[-mp.inf,mp.inf])

# ----------------------------------------------------------------------------
def dDedekindEta(tau):
    ''' Compute the derivative of the Dedekind Eta function for imaginary argument tau. 
        See: http://mathworld.wolfram.com/DedekindEtaFunction.html. '''
    try:
        import mpmath as mp
        mpmath_loaded = True
    except ImportError:
        mpmath_loaded = False 
    qbar = mp.exp(-2.0*mp.pi*tau)

    res = mp.nsum(lambda n: ((-1)**n)*(6.0*n-1.0)*(6.0*n-1.0)*(qbar**((6.0*n-1.0)*(6.0*n-1.0)/24.0)),[-mp.inf,mp.inf])
    return mp.pi*mp.j*res / 12.0

# ----------------------------------------------------------------------------
def dThetaRatio(z,q,n):
    ''' Compute the ratio of the nth derivative of the third Jacobi theta function
        with itself.'''
    try:
        import mpmath as mp
        mpmath_loaded = True
    except ImportError:
        mpmath_loaded = False 
        
    #return (mp.djtheta(3,z,q,n)/mp.jtheta(3,z,q))
    return (mp.jtheta(3,z,q,n)/mp.jtheta(3,z,q))

# ----------------------------------------------------------------------------
class PyFit:
    ''' A class which performs non-linear regression fits to 1D data. '''

    # Test that we can load the desired modules
    try:
        import numpy as n
        numpy_loaded = True
    except ImportError:
        numpy_loaded = False 

    try:
        import scipy.optimize.minpack as min
        scipy_loaded = True
    except ImportError:
        scipy = False 

    # ------------------------------------------------------------------------
    def __init__(self,initialFitPar,fitFunc,inX,inY,err=None,funcPar=None,dFunc=None):
        ''' We setup all local class variables, which may be overwritten
            usint initFit. '''

        # Setup all class variables
        self.params    = funcPar
        self.function  = fitFunc
        self.dfunction = dFunc
        self.x         = inX
        self.y         = inY
        self.err       = err 
        self.aFit      = initialFitPar

        self.dataLength = len(inX)
        self.fitLength  = len(self.aFit)

        # Initialize and update the weights
        self.weight = self.n.ones([self.dataLength],float)
        self.getWeights()

    # ------------------------------------------------------------------------
    def objective(self,a): 
        '''Calculate the residual deviation with the fitting function. '''

        obj = self.n.zeros([self.dataLength],float)

        for i in range(self.dataLength):
            obj[i] = self.weight[i]*(self.y[i] - self.function(self.x[i],a,self.params))

        return obj

    # ------------------------------------------------------------------------
    def dobjective(self,a): 
        '''Calculate the residual deviation with the fitting function. '''

        dobj = self.n.zeros([self.dataLength,len(a)],float)

        for i in range(self.dataLength):
            dobj[i,:] = -self.weight[i]*self.dfunction(self.x[i],a,self.params)

        return dobj

    # ------------------------------------------------------------------------
    def getWeights(self):
        ''' Get the appropriate weights. ''' 
        if self.err is None:
            self.weight = self.n.ones([self.dataLength],float)
        else:
            for i in range(self.dataLength):
                if (self.err[i] == 0.0):
                    self.weight[i] = 1.0
                else:
                    self.weight[i] = 1.0/self.err[i]

    # ------------------------------------------------------------------------
    def updateFit(self,initx,inity,inita,err=None,funcPar=None): 
        ''' Update the data that will be fit. '''

        self.params = funcPar
        self.x = initx
        self.y = inity
        self.err = err
        self.aFit = inita

        # Update the weights
        self.getWeights()

    # ------------------------------------------------------------------------
    def getFit(self): 
        ''' Perform the actual least-squares fit. '''
        
        # Have we supplied a Jacobian for our fit function?
        if self.dfunction is not None:
            df = self.dobjective
        else:
            df = None

        # Perform the fit
        self.aFit,covMat,infodict,mesg,ierr = self.min.leastsq(self.objective, self.aFit,\
                Dfun=df, epsfcn=1.0E-6,full_output=True) 

        # Get chi^2
        self.yFit = self.n.zeros([self.dataLength],float)
        self.chi2 = 0.0
        dof = 0
        for i in range(self.dataLength):
            self.yFit[i] = self.function(self.x[i],self.aFit,self.params)
            if self.err[i] > 1.0E-12:
                dof += 1
                self.chi2 += (self.y[i]-self.yFit[i])**2/(self.err[i]**2)
        self.chi2 /= 1.0*(dof - self.fitLength - 1)

        # this is the adjusted covariance matrix
        s_sq = (self.objective(self.aFit)**2).sum()/(self.dataLength-self.fitLength) 
        try:
            covMat = covMat * s_sq
        except:
            covMat = self.n.zeros([self.fitLength,self.fitLength])

        # Find the error in the fit parameters
        if not isList(self.aFit):
            self.aErr = self.n.sqrt(covMat[0,0])
        else:
            self.aErr = []
            for i in range(self.fitLength):
                self.aErr.append(self.n.sqrt(covMat[i,i]))
