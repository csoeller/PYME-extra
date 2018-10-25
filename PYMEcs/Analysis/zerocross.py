# this code is borrowed from PyAstronomy
# I include their license below under which this code was made available

# Copyright (c) 2011, PyA group

# Permission is hereby granted, free of charge, to any person obtaining a copy of this
# software and associated documentation files (the "Software"), to deal in the Software
# without restriction, including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons
# to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies
# or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
# PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
# FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

from __future__ import division
import numpy as np
# from PyAstronomy.pyaC import pyaErrors as PE

def zerocross1d(x, y, getIndices=False):
  """
    Find the zero crossing points in 1d data.
    
    Find the zero crossing events in a discrete data set.
    Linear interpolation is used to determine the actual
    locations of the zero crossing between two data points
    showing a change in sign. Data point which are zero
    are counted in as zero crossings if a sign change occurs
    across them. Note that the first and last data point will
    not be considered whether or not they are zero. 
    
    Parameters
    ----------
    x, y : arrays
        Ordinate and abscissa data values.
    getIndices : boolean, optional
        If True, also the indicies of the points preceding
        the zero crossing event will be returned. Default is
        False.
    
    Returns
    -------
    xvals : array
        The locations of the zero crossing events determined
        by linear interpolation on the data.
    indices : array, optional
        The indices of the points preceding the zero crossing
        events. Only returned if `getIndices` is set True.
  """
  
  # Check sorting of x-values
  if np.any((x[1:] - x[0:-1]) <= 0.0):
    raise( RuntimeError("zerocross1d: The x-values must be sorted in ascending order! " +
                        "Sort the data prior to calling zerocross1d."))
  
  # Indices of points *before* zero-crossing
  indi = np.where(y[1:]*y[0:-1] < 0.0)[0]
  
  # Find the zero crossing by linear interpolation
  dx = x[indi+1] - x[indi]
  dy = y[indi+1] - y[indi]
  zc = -y[indi] * (dx/dy) + x[indi]
  
  # What about the points, which are actually zero
  zi = np.where(y == 0.0)[0]
  # Do nothing about the first and last point should they
  # be zero
  zi = zi[np.where((zi > 0) & (zi < x.size-1))]
  # Select those point, where zero is crossed (sign change
  # across the point)
  zi = zi[np.where(y[zi-1]*y[zi+1] < 0.0)]
  
  # Concatenate indices
  zzindi = np.concatenate((indi, zi)) 
  # Concatenate zc and locations corresponding to zi
  zz = np.concatenate((zc, x[zi]))
  
  # Sort by x-value
  sind = np.argsort(zz)
  zz, zzindi = zz[sind], zzindi[sind]
  
  if not getIndices:
    return zz
  else:
    return zz, zzindi
