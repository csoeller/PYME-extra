from PYME.recipes.base import ModuleBase, register_module, Filter
from PYME.recipes.traits import Input, Output, Float, Enum, CStr, Bool, Int,  File

import numpy as np
import skimage.filters as skf

@register_module('FlexiThreshold') 
class FlexiThreshold(Filter):
    """Chose a threshold using a range of available thresholding methods.
       Currently we can chose from: simple, fractional, otsu, isodata
    """
    method = Enum('simple','fractional','otsu','isodata')
    parameter = Float(0.5)

    def fractionalThreshold(self, data):
        N, bins = np.histogram(data, bins=5000)
        #calculate bin centres
        bin_mids = (bins[:-1] )
        cN = np.cumsum(N*bin_mids)
        i = np.argmin(abs(cN - cN[-1]*(1-self.parameter)))
        threshold = bins[i]
        return threshold

    def applyFilter(self, data, chanNum, frNum, im):

        if self.method == 'fractional':
            threshold = self.fractionalThreshold(data)
        if self.method == 'simple':
            threshold = self.parameter
        if self.method == 'otsu':
            threshold = skf.threshold_otsu(data)
        if self.method == 'isodata':
            threshold = skf.threshold_isodata(data)

        mask = data > threshold
        return mask

    def completeMetadata(self, im):
        im.mdh['Processing.ThresholdParameter'] = self.parameter
        im.mdh['Processing.ThresholdMethod'] = self.method
