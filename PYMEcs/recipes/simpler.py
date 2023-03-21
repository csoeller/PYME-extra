from PYME.recipes.base import register_module, ModuleBase, Filter
from PYME.recipes.traits import Input, Output, Float, Enum, CStr, Bool, Int, List, DictStrStr, DictStrList, ListFloat, ListStr, FileOrURI

import numpy as np
from PYME.IO import tabular

import logging
logger = logging.getLogger(__file__)

from scipy.stats import binned_statistic
def mystd(vec):
    return np.std(vec,ddof=0)

# note that we are using our manual implementation of stddev (mystd) rather than the builtin 'std'
# our tests showed that using 'std' can result in NaNs or other non-finite values
# hopfully this will be fixed at some stage
def get_stddev_property(ids, prop, use_builtin_std=False):
    maxid = int(ids.max())
    edges = -0.5+np.arange(maxid+2)
    idrange = (0,maxid)

    if use_builtin_std:
        statistic = 'std'
    else:
        statistic = mystd
        
    propstd, bin_edge, binno = binned_statistic(ids, prop, statistic=statistic, bins=edges, range=idrange)
    propstd[np.isnan(propstd)] = 1000.0 # (mark as huge error_z)
    std_events = propstd[ids]
    binno_events = binno[ids]
    return std_events, binno_events

def get_values_from_image(values_image, points, normalise=False, n0max=1.0):
    """
    Function to extract values from a segmented image (2D or 3D) at given locations.
    
    Parameters
    ----------
    values_image: PYME.IO.image.ImageStack instance
        an image containing object labels
    points: tabular-like (PYME.IO.tabular, np.recarray, pandas DataFrame) containing 'x', 'y' & 'z' columns
        locations at which to extract labels

    Returns
    -------
    ids: Label number from image, mapped to each localization within that label

    """
    from PYME.Analysis.points.coordinate_tools import pixel_index_of_points_in_image

    pixX, pixY, pixZ = pixel_index_of_points_in_image(values_image, points)

    values_data = values_image.data_xyztc

    if values_data.shape[2] == 1:
        # disregard z for 2D images
        pixZ = np.zeros_like(pixX)

    ind = (pixX < values_data.shape[0]) * (pixY < values_data.shape[1]) * (pixX >= 0) * (pixY >= 0) * (pixZ >= 0) * (
        pixZ < values_data.shape[2])

    vals = np.zeros_like(pixX)
    imgdata = np.clip(values_data[:,:,:,0,0].squeeze(),0,None)  # we assume no time sequence, only 1 colour

    if normalise:
        maxval = np.percentile(imgdata,97.5)
        imgdata *= n0max/maxval
    
    # assume there is only one channel
    vals[ind] = np.atleast_3d(imgdata)[pixX[ind], pixY[ind], pixZ[ind]].astype('i')

    return vals


@register_module('ClusterModes')
class ClusterModes(ModuleBase):
    
    inputName = Input('dbscanClustered')
    IDkey = CStr('dbscanClumpID')
    outputName = Output('with_clusterModes')
    PropertyKey = CStr('nPhotons')

    def execute(self, namespace):
        from PYMEcs.Analysis.Simpler import clusterModes
        
        inp = namespace[self.inputName]
        cmodes = tabular.MappingFilter(inp)

        ids = inp[self.IDkey] # I imagine this needs to be an int type key
        props = inp[self.PropertyKey]

        cm, ce, ccx, ccy = clusterModes(inp['x'],inp['y'],ids,props)
        cmodes.addColumn('clusterMode',cm)
        cmodes.addColumn('clusterModeError',ce)
        cmodes.addColumn('clusterCentroid_x',ccx)
        cmodes.addColumn('clusterCentroid_y',ccy)
        
        # propogate metadata, if present
        try:
            cmodes.mdh = inp.mdh
        except AttributeError:
            pass
        
        namespace[self.outputName] = cmodes

    @property
    def _key_choices(self):
        #try and find the available column names
        try:
            return sorted(self._parent.namespace[self.inputName].keys())
        except:
            return []

    @property
    def default_view(self):
        from traitsui.api import View, Group, Item
        from PYME.ui.custom_traits_editors import CBEditor

        return View(Item('inputName', editor=CBEditor(choices=self._namespace_keys)),
                    Item('_'),
                    Item('IDkey', editor=CBEditor(choices=self._key_choices)),
                    Item('PropertyKey', editor=CBEditor(choices=self._key_choices)),
                    Item('_'),
                    Item('outputName'), buttons=['OK'])


@register_module('N0FromImage')
class N0FromImage(ModuleBase):
    """
    Maps each point in the input table to a pixel in a labelled image, and extracts the pixel value at that location to
    use as a label for the point data. 

    Inputs
    ------
    inputName: Input
        name of tabular input containing positions ('x', 'y', and optionally 'z' columns should be present)
    inputImage: Input
        name of image input containing N0 data

    Outputs
    -------
    outputName: Output
        name of tabular output. A mapped version of the tabular input with one extra column
    value_key_name : CStr
        name of new column which will contain the label number from image, mapped to each localization within that label

    """
    inputName = Input('input')
    inputImage = Input('n0data')
    normaliseN0 = Bool(False)
    maxN0 = Float(1.0)

    value_key_name = CStr('N0')

    outputName = Output('with_N0')

    def execute(self, namespace):
        from PYME.IO import tabular

        inp = namespace[self.inputName]
        img = namespace[self.inputImage]

        n0 = get_values_from_image(img, inp, normalise=self.normaliseN0, n0max=self.maxN0)

        withN0 = tabular.MappingFilter(inp)
        withN0.addColumn(self.value_key_name, n0)

        # propagate metadata, if present
        try:
            withN0.mdh = namespace[self.inputName].mdh
        except AttributeError:
            pass

        namespace[self.outputName] = withN0



@register_module('N0FromInterpolationMap')
class N0FromInterpolationMap(ModuleBase):

    inputName = Input('filtered')
    outputName = Output('with_N0')
    # keywords inherited from FILE, see https://docs.enthought.com/traits/traits_user_manual/defining.html
    # note that docs do not emphasize that filter keyword value must be an array of wildcard strings!
    N0_map_file = FileOrURI(filter=['*.n0m'], exists=True)    
    normaliseN0 = Bool(False)
    maxN0 = Float(1.0)

    def execute(self, namespace):
        inp = namespace[self.inputName]
        mapped = tabular.MappingFilter(inp)

        # we may want to cache the read! traits has a way to do this, see
        #    traits.has_traits.cached_property in https://docs.enthought.com/traits/traits_api_reference/has_traits.html
        # but this may not be compatible with the recipes use of traits
        # in that case will have to use one of the ways in which existing recipe modules achieve this
        from six.moves import cPickle
        with open(self.N0_map_file, 'rb') as fid:
            n0m,bb,origin = cPickle.load(fid)        
        try:
            N0 = n0m(inp['x'],inp['y'],grid=False) # this should ensure N0 is floating point type
        except TypeError:
            N0 = n0m(inp['x'],inp['y'])

        if self.normaliseN0:
            maxval = np.percentile(N0,97.5)
            N0 *= self.maxN0/maxval           

        mapped.addColumn('N0', N0)

        # propagate metadata, if present
        try:
            mapped.mdh = namespace[self.inputName].mdh
        except AttributeError:
            pass

        namespace[self.outputName] = mapped


# here we process the data so that we keep one fixed N0max for each event and
# rescale all nPhotons values to reflect the local change in N0
# this is then more similar to what is described in the SIMPLER paper
# and allows us to pass the filtered events to our students for further processing!
@register_module('ScaleNPhotonsFromN0')
class ScaleNPhotonsFromN0(ModuleBase):
    inputName = Input('with_N0')
    outputName = Output('nPhotonsScaled')
    
    def execute(self, namespace):
        inp = namespace[self.inputName]
        mapped = tabular.MappingFilter(inp)

        N0 = inp['N0']
        N0maxval = np.percentile(N0,97.5)
        N0max = N0maxval * np.ones_like(N0, dtype='f')
        nPhotonsScaled = inp['nPhotons'] * N0maxval / N0
    
        mapped.addColumn('N0', N0max)
        mapped.addColumn('nPhotons', nPhotonsScaled)

        # propagate metadata, if present
        try:
            mapped.mdh = namespace[self.inputName].mdh
        except AttributeError:
            pass

        namespace[self.outputName] = mapped

        
@register_module('SIMPLERzgenerator')
class SIMPLERzgenerator(ModuleBase):

    inputName = Input('filtered')
    outputName = Output('with_simpler_z')
    df_in_nm = Float(88.0)
    alphaf = Float(0.9)
    N0_scale_factor = Float(1.0)
    N0_is_uniform = Bool(False)
    with_error_z = Bool(False)
    use_builtin_std = Bool(True)
    
    def execute(self, namespace):
        inp = namespace[self.inputName]
        mapped = tabular.MappingFilter(inp)
        if self.N0_is_uniform:
            N0 = np.ones_like(inp['x'])
        else:
            N0 = inp['N0']
        N = inp['nPhotons']
        NoverN0 = N/(N0*self.N0_scale_factor)
        simpler_z = self.df_in_nm*np.log(self.alphaf/(NoverN0 - (1 - self.alphaf)))
        simpler_z[np.isnan(simpler_z)] = -100.0
        simpler_z[np.isinf(simpler_z)] = -100.0

        
        mapped.addColumn('NoverN0', NoverN0)
        mapped.addColumn('z', simpler_z)
        if self.with_error_z:
            error_z, ezn = get_stddev_property(inp['clumpIndex'], simpler_z,
                                               use_builtin_std=self.use_builtin_std)
            mapped.addColumn('error_z', error_z)
            mapped.addColumn('error_zN', ezn)

        # propogate metadata, if present
        try:
            mapped.mdh = inp.mdh
        except AttributeError:
            pass
        
        namespace[self.outputName] = mapped
