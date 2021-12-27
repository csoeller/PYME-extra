from PYME.recipes.base import register_module, ModuleBase, Filter
from PYME.recipes.traits import Input, Output, Float, Enum, CStr, Bool, Int, List, DictStrStr, DictStrList, ListFloat, ListStr, FileOrURI

import numpy as np
from PYME.IO import tabular

import logging
logger = logging.getLogger(__file__)

def get_values_from_image(label_image, points, normalise=False, n0max=1.0):
    """
    Function to extract values from a segmented image (2D or 3D) at given locations.
    
    Parameters
    ----------
    label_image: PYME.IO.image.ImageStack instance
        an image containing object labels
    points: tabular-like (PYME.IO.tabular, np.recarray, pandas DataFrame) containing 'x', 'y' & 'z' columns
        locations at which to extract labels

    Returns
    -------
    ids: Label number from image, mapped to each localization within that label

    """
    from PYME.Analysis.points.coordinate_tools import pixel_index_of_points_in_image

    pixX, pixY, pixZ = pixel_index_of_points_in_image(label_image, points)

    label_data = label_image.data_xyztc

    if label_data.shape[2] == 1:
        # disregard z for 2D images
        pixZ = np.zeros_like(pixX)

    ind = (pixX < label_data.shape[0]) * (pixY < label_data.shape[1]) * (pixX >= 0) * (pixY >= 0) * (pixZ >= 0) * (
        pixZ < label_data.shape[2])

    ids = np.zeros_like(pixX)
    imgdata = np.clip(label_data[:,:,:,0,0].squeeze(),0,None)  # we assume no time sequence, only 1 colour

    if normalise:
        maxval = np.percentile(imgdata,97.5)
        imgdata *= n0max/maxval
    
    # assume there is only one channel
    ids[ind] = np.atleast_3d(imgdata)[pixX[ind], pixY[ind], pixZ[ind]].astype('i')

    return ids


@register_module('ClusterModes')
class ClusterModes(ModuleBase):
    
    inputName = Input('dbscanClustered')
    IDkey = CStr('dbscanClumpID')
    outputName = Output('with_clusterModes')
    PropertyKey = CStr('nPhotons')

    def execute(self, namespace):
        from PYMEcs.Analysis.Simpler import clusterModes
        
        inp = namespace[self.inputName]
        cmodes = tabular.mappingFilter(inp)

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
        name of tabular output. A mapped version of the tabular input with 2 extra columns
    label_key_name : CStr
        name of new column which will contain the label number from image, mapped to each localization within that label
    label_count_key_name : CStr
        name of new column which will contain the number of localizations within the label that a given localization
        belongs to
    minimum_localizations: Int
        threshold for the number of localizations required to propagate a label through to localizations

    """
    inputName = Input('input')
    inputImage = Input('n0data')
    normaliseN0 = Bool(False)
    maxN0 = Float(1.0)

    label_key_name = CStr('N0')

    outputName = Output('with_N0')

    def execute(self, namespace):
        from PYME.IO import tabular

        inp = namespace[self.inputName]
        img = namespace[self.inputImage]

        n0 = get_values_from_image(img, inp, normalise=self.normaliseN0, n0max=self.maxN0)

        withN0 = tabular.MappingFilter(inp)
        withN0.addColumn(self.label_key_name, n0)

        # propagate metadata, if present
        try:
            withN0.mdh = namespace[self.inputName].mdh
        except AttributeError:
            pass

        namespace[self.outputName] = withN0

@register_module('SIMPLERzgenerator')
class SIMPLERzgenerator(ModuleBase):

    inputName = Input('filtered')
    outputName = Output('with_simpler_z')
    df_in_nm = Float(88.0)
    alphaf = Float(0.9)
    N0_scale_factor = Float(1.0)
    N0_is_uniform = Bool(False)
    N0_map_file = FileOrURI(filter=['*.n0m'], exists=True) # keywords inherited from FILE, see https://docs.enthought.com/traits/traits_user_manual/defining.html
    # note that docs do not emphasize that filter keyword value must be an array of wildcard strings!
    
    def execute(self, namespace):
        inp = namespace[self.inputName]
        mapped = tabular.mappingFilter(inp)

        if not self.N0_is_uniform:
            # we may want to cache the read! traits has a way to do this, see
            #    traits.has_traits.cached_property in https://docs.enthought.com/traits/traits_api_reference/has_traits.html
            # but this may not be compatible with the recipes use of traits
            # in that case will have to use one of the ways in which existing recipe modules achieve this
            from six.moves import cPickle
            with open(self.N0_map_file, 'rb') as fid:
                n0m,bb = cPickle.load(fid)        
            try:
                N0 = n0m(inp['x'],inp['y'],grid=False) # this should ensure N) is floating point type
            except TypeError:
                N0 = n0m(inp['x'],inp['y'])
        else:
            N0 = np.ones_like(inp['x'])
        N0 *= self.N0_scale_factor
        N = inp['nPhotons']
        NoverN0 = N/N0
        simpler_z = self.df_in_nm*np.log(self.alphaf/(NoverN0 - (1 - self.alphaf)))
        simpler_z[np.isnan(simpler_z)] = -100.0

        mapped.addColumn('N0', N0)
        mapped.addColumn('NoverN0', NoverN0)
        mapped.addColumn('z', simpler_z)

        # propogate metadata, if present
        try:
            mapped.mdh = inp.mdh
        except AttributeError:
            pass
        
        namespace[self.outputName] = mapped
