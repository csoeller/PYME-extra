from PYME.recipes.base import ModuleBase, register_module, Filter
from PYME.recipes.traits import Input, Output, Float, Enum, CStr, Bool, Int,  File

import numpy as np
import skimage.filters as skf
from scipy import ndimage

from PYME.IO.image import ImageStack
from PYME.IO import MetaDataHandler


@register_module('FlexiThreshold') 
class FlexiThreshold(Filter):
    """Chose a threshold using a range of available thresholding methods.
       Currently we can chose from: simple, fractional, otsu, isodata
    """
    method = Enum('simple','fractional','otsu','isodata',
                  'li','yen') # newer skimage has minimum, mean and triangle as well
    parameter = Float(0.5)
    clipAt = Float(2e6) # used to be 10 - increase to large value for newer PYME renderings

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
            threshold = self.fractionalThreshold(np.clip(data,None,self.clipAt))
        elif self.method == 'simple':
            threshold = self.parameter
        else:
            method = getattr(skf,'threshold_%s' % self.method)
            threshold = method(np.clip(data,None,self.clipAt))

        mask = data > threshold
        return mask

    def completeMetadata(self, im):
        im.mdh['Processing.ThresholdParameter'] = self.parameter
        im.mdh['Processing.ThresholdMethod'] = self.method

@register_module('LabelRange')        
class LabelRange(Filter):
    """Asigns a unique integer label to each contiguous region in the input mask.
    Throws away all regions which are outside of given number of pixel range.
    Also uses the number of sites from a second input channel to decide if region is retained,
    retaining only those with the number sites in a given range.
    """
    inputSitesLabeled = Input("sites") # sites and the main input must have the same shape!
    minRegionPixels = Int(10)
    maxRegionPixels = Int(100)
    minSites = Int(4)
    maxSites = Int(6)
    sitesAsMaxima = Bool(False)

    def filter(self, image, imagesites):
        #from PYME.util.shmarray import shmarray
        #import multiprocessing
        
        if self.processFramesIndividually:
            filt_ims = []
            for chanNum in range(image.data.shape[3]):
                filt_ims.append(np.concatenate([np.atleast_3d(self.applyFilter(image.data[:,:,i,chanNum].squeeze().astype('f'), imagesites.data[:,:,i,chanNum].squeeze().astype('f'), chanNum, i, image)) for i in range(image.data.shape[2])], 2))
        else:
            filt_ims = [np.atleast_3d(self.applyFilter(image.data[:,:,:,chanNum].squeeze().astype('f'), imagesites.data[:,:,:,chanNum].squeeze().astype('f'), chanNum, 0, image)) for chanNum in range(image.data.shape[3])]
            
        im = ImageStack(filt_ims, titleStub = self.outputName)
        im.mdh.copyEntriesFrom(image.mdh)
        im.mdh['Parent'] = image.filename
        
        self.completeMetadata(im)
        
        return im
        
    def execute(self, namespace):
        namespace[self.outputName] = self.filter(namespace[self.inputName],namespace[self.inputSitesLabeled])

    def applyFilter(self, data, sites, chanNum, frNum, im):

        # siteLabels = self.recipe.namespace[self.sitesLabeled]
        
        mask = data > 0.5
        labs, nlabs = ndimage.label(mask)
        
        rSize = self.minRegionPixels
        rMax = self.maxRegionPixels

        minSites = self.minSites
        maxSites = self.maxSites
        
        m2 = 0*mask
        objs = ndimage.find_objects(labs)
        for i, o in enumerate(objs):
            r = labs[o] == i+1
            #print r.shape
            area = r.sum()
            if (area >= rSize) and (area <= rMax):
                if self.sitesAsMaxima:
                    nsites = sites[o][r].sum()
                else:
                    nsites = (np.unique(sites[o][r]) > 0).sum() # count the unique labels (excluding label 0 which is background)
                if (nsites >= minSites) and (nsites <= maxSites):
                    m2[o] += r

        labs, nlabs = ndimage.label(m2 > 0)

        return labs

    def completeMetadata(self, im):
        im.mdh['Labelling.MinSize'] = self.minRegionPixels
        im.mdh['Labelling.MaxSize'] = self.maxRegionPixels
        im.mdh['Labelling.MinSites'] = self.minSites
        im.mdh['Labelling.MaxSites'] = self.maxSites


@register_module('LabelByArea')        
class LabelByArea(Filter):
    """Asigns a unique integer label to each contiguous region in the input mask.
    Optionally throws away all regions which are smaller than a cutoff size.
    """
    
    def applyFilter(self, data, chanNum, frNum, im):
        mask = data > 0.5
        labs, nlabs = ndimage.label(mask)
        
        m2 = 0*mask
        objs = ndimage.find_objects(labs)
        for i, o in enumerate(objs):
            r = labs[o] == i+1
            #print r.shape
            area = r.sum()
            m2[o] += r*area

        return m2

    def completeMetadata(self, im):
        im.mdh['Labelling.Property'] = 'area'


import skimage.measure
import math

@register_module('LabelByRegionProperty')        
class LabelByRegionProperty(Filter):
    """Asigns a region property to each contiguous region in the input mask.
    Optionally throws away all regions for which property is outside a given range.
    """
    regionProperty = Enum(['area','circularity','aspectratio'])
    filterByProperty = Bool(False)
    propertyMin = Float(0)
    propertyMax = Float(1e6)
    
    def applyFilter(self, data, chanNum, frNum, im):
        mask = data > 0.5
        labs, nlabs = ndimage.label(mask)
        rp = skimage.measure.regionprops(labs,None,cache=True)
        
        m2 = np.zeros_like(mask,dtype='float')
        objs = ndimage.find_objects(labs)
        for region in rp:
            oslices = objs[region.label-1]
            r = labs[oslices] == region.label
            #print r.shape
            if self.regionProperty == 'area':
                propValue = region.area
            elif self.regionProperty == 'aspectratio':
                propValue = region.major_axis_length / region.minor_axis_length
            elif self.regionProperty == 'circularity':
                propValue = 4 * math.pi * region.area / (region.perimeter*region.perimeter)
            if self.filterByProperty:
                if (propValue >= self.propertyMin) and (propValue <= self.propertyMax):
                    m2[oslices] += r*propValue
            else:
                m2[oslices] += r*propValue

        return m2

    def completeMetadata(self, im):
        im.mdh['Labelling.Property'] = self.regionProperty


from scipy.signal import fftconvolve
def circle_kernel(diam=112.0, sigma=10.0, extent_nm=150.0,voxelsize_nm=5.0,eps=15.0):
    x = np.arange(-extent_nm/2.0, extent_nm/2.0+1.0, voxelsize_nm, dtype='f')
    y = x.copy()
    x2d,y2d = np.meshgrid(x,y)
    circ2d = (x2d**2 + y2d**2 -0.25*diam**2 <= eps**2) & (x2d**2 + y2d**2 -0.25*diam**2 >= -eps**2)
    g2d = np.exp(-(x2d**2+y2d**2)/2.0/(sigma**2))
    ckernel = np.clip(fftconvolve(circ2d,g2d,mode='same'),0,None)
    ckernel /= ckernel.sum()

    return ckernel

@register_module('CircleConvolution')
class CircleConvolution(ModuleBase):
    inputImage = Input('input')
    outputImage = Output('circle_convol')
    outputKernel = Output('circle_kernel')

    circleDiameter = Float(112.0)
    circleBlurRadius = Float(15.0)

    def run(self, inputImage):

        img = inputImage.data_xyztc[:,:,0,0,0].squeeze()
        ckern = circle_kernel(diam=self.circleDiameter,extent_nm=250.0,eps=30.0,
                              sigma=self.circleBlurRadius,voxelsize_nm=inputImage.voxelsize_nm.x)

        conv = np.clip(fftconvolve(img,ckern,mode='same'),0,None)
        
        imconvol = ImageStack(conv,mdh=MetaDataHandler.NestedClassMDHandler(inputImage.mdh), titleStub = self.outputImage)
        imconvol.mdh['CircleConvolution.Diameter'] = self.circleDiameter
        imconvol.mdh['CircleConvolution.BlurRadius'] = self.circleBlurRadius
        circleKernel = ImageStack(ckern, titleStub = self.outputKernel)
        circleKernel.mdh['Diameter'] = self.circleDiameter
        circleKernel.mdh['BlurRadius'] = self.circleBlurRadius
        
        return {'outputImage' : imconvol, 'outputKernel' : circleKernel}
    

# split single channel image to two channels for FRC
# following paper from Riegher et al:
#      Rieger, B., Droste, I., Gerritsma, F., Ten Brink, T. & Stallinga, S.
#      Single image Fourier ring correlation. Opt. Express 32, 21767 (2024).
@register_module('Split1FRC')
class CircleConvolution(ModuleBase):
    inputImage = Input('input')
    outputImage = Output('split_1frc')

    splitProb = Float(0.5)

    def run(self, inputImage):
        im0 = inputImage.data_xyztc[:,:,0,0,0]
        rng = np.random.default_rng()
        im0_1 = rng.binomial(im0,self.splitProb)
        im0_2 = im0 - im0_1
        im_1frc = np.stack([im0_1,im0_2],axis=2).astype(inputImage.data_xyztc.dtype)
        imstack1frc = ImageStack(im_1frc[:,:,None,None,:],mdh=MetaDataHandler.NestedClassMDHandler(inputImage.mdh),titleStub = self.outputImage)

        return imstack1frc

