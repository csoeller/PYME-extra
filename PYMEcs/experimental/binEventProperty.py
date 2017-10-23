import numpy as np
import matplotlib.pyplot as plt
import wx

import logging
logger = logging.getLogger(__file__)

from traits.api import HasTraits, Str, Int, CStr, List, Enum, Float
from traitsui.api import View, Item, Group
from traitsui.menu import OKButton, CancelButton, OKCancelButtons

from PYME.DSView.dsviewer import ViewIm3D, ImageStack

class propertyChoice(HasTraits):
    clist = List([])
    EventProperty = Enum(values='clist')
    BinByProperty = Enum(values='clist')
    BinWidth = Float(100.0)

    traits_view = View(Group(Item(name = 'EventProperty'),
                             Item(name = 'BinByProperty'),
                             Item('_'),
                             Item(name = 'BinWidth'),
                             label = 'Select Properties and Binning',
                             show_border = True),
                       buttons = OKCancelButtons)

    def add_channels(self,chans):
        for chan in chans:
            if chan not in self.clist:
                self.clist.append(chan)


class PropertyBinning:
    """

    """
    def __init__(self, visFr):
        self.visFr = visFr
        self.pipeline = visFr.pipeline

        visFr.AddMenuItem('Experimental', 'Bin Event Property by coordinate', self.OnBinProperty,
                          helpText='Bin average an event property, e.g. background,  as a function of another property, e.g. x')

    def OnBinProperty(self, event=None):
        # these are the defaults
        avProperty = 'nPhotons'
        byProperty = 'x'
        binwidth = 100

        # traitsui dialog to set the event properties etc
        pChoice = propertyChoice()
        pChoice.add_channels(sorted(self.pipeline.keys()))
        # select defaults if these event properties exist
        if avProperty in pChoice.clist:
            pChoice.EventProperty = avProperty
        if byProperty in pChoice.clist:
            pChoice.BinByProperty = byProperty

        if pChoice.configure_traits(kind='modal'):
            avProperty = pChoice.EventProperty
            byProperty = pChoice.BinByProperty
            binwidth = pChoice.BinWidth
        else:
            return

        avp = self.pipeline[avProperty]
        byp = self.pipeline[byProperty]

        bypmin = byp.min()
        bypmax = byp.max()
        nbins = int((bypmax-bypmin)/binwidth)
        bins = bypmin + np.arange(nbins+1)*binwidth

        binctrs = 0.5*(bins[0:-1]+bins[1:])
        
        avbin, abins = np.histogram(byp, bins=bins, weights=avp)
        bybin, bybins = np.histogram(byp, bins=bins)

        good = bybin > 0
        bad = bybin <= 0

        avmean = np.zeros_like(avbin,dtype='float64')
        avmean[good] = avbin[good]/bybin[good]
        avmean[bad] = np.nan

        # direct matplotlib plotting scrapped for ImageStack approach
        # which allows saving of data
        #plt.figure()
        #plt.plot(binctrs, avmean)
        #plt.xlabel(byProperty)
        #plt.ylabel('mean of %s' % avProperty)

        plots = []
        plots.append(avmean.reshape(-1, 1,1))

        im = ImageStack(plots, titleStub='mean of %s' % avProperty)
        im.xvals = binctrs
        im.xlabel = byProperty

        im.ylabel = 'mean of %s' % avProperty
        im.defaultExt = '*.txt'

        im.mdh['voxelsize.x'] = binwidth
        im.mdh['ChannelNames'] = [avProperty]
        im.mdh['Profile.XValues'] = im.xvals
        im.mdh['Profile.XLabel'] = im.xlabel
        im.mdh['Profile.YLabel'] = im.ylabel

        im.mdh['OriginalImage'] = self.pipeline.filename

        ViewIm3D(im, mode='graph', parent=wx.GetTopLevelParent(self.visFr))

        
def Plug(visFr):
    """Plugs this module into the gui"""
    visFr.propBin = PropertyBinning(visFr)

