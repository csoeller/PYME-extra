from PYME.recipes.base import register_module, ModuleBase, Filter
from PYME.recipes.traits import HasTraits, Float, List, Bool, Int, CStr, Enum, on_trait_change, Input, Output

from PYME.IO.image import ImageStack
import numpy as np

import logging
logger = logging.getLogger(__name__)

@register_module('ExtractChannelByName')    
class ExtractChannelByName(ModuleBase):
    """Extract one channel from an image using regular expression matching to image channel names - by default this is case insensitive"""
    inputName = Input('input')
    outputName = Output('filtered_image')

    channelNamePattern = CStr('channel0')
    caseInsensitive = Bool(True)
    
    def _pickChannel(self, image):
        import re
        flags = 0
        if self.caseInsensitive:
            flags |= re.I
        channelNames = image.mdh['ChannelNames']
        idx = -1
        for i, c in enumerate(channelNames):
            if re.search(self.channelNamePattern,c,flags):
                idx = i
        if idx < 0:
            raise RuntimeError("Expression %s did not maych any channel names" % self.channelNamePattern)
        
        chan = image.data[:,:,:,idx]
        
        im = ImageStack(chan, titleStub = 'Filtered Image')
        im.mdh.copyEntriesFrom(image.mdh)
        im.mdh['ChannelNames'] = [channelNames[idx]]
        im.mdh['Parent'] = image.filename
        
        return im
    
    def execute(self, namespace):
        namespace[self.outputName] = self._pickChannel(namespace[self.inputName])


