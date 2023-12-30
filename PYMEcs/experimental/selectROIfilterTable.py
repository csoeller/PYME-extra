import numpy as np
from PYME.recipes.tablefilters import FilterTable
import wx

def selectWithDialog(choices, message='select image from list', caption='Selection'):
    dlg = wx.SingleChoiceDialog(None, message, caption, choices, wx.CHOICEDLG_STYLE)
    if dlg.ShowModal() == wx.ID_OK:
        item = dlg.GetStringSelection()
    else:
        item = None
    dlg.Destroy()
    return item


class SelectROIFT:
    def __init__(self, visFr):
        self.visFr = visFr

        visFr.AddMenuItem('Experimental>View', "Add ROI FilterTable module from selection", self.OnSelROIFT)
        visFr.AddMenuItem('Experimental>View', "Add ROI FilterTable module from image", self.OnImROIFT)

    def OnSelROIFT(self, event):
        try:
            #old glcanvas
            x0, y0 = self.visFr.glCanvas.selectionStart[0:2]
            x1, y1 = self.visFr.glCanvas.selectionFinish[0:2]
        except AttributeError:
            #new glcanvas
            x0, y0 = self.visFr.glCanvas.selectionSettings.start[0:2]
            x1, y1 = self.visFr.glCanvas.selectionSettings.finish[0:2]

        filters = {}
        filters['x'] = [float(min(x0, x1)), float(max(x0, x1))] # must ensure all values are eventually scalars to avoid issue with recipe yaml output
        filters['y'] = [float(min(y0, y1)), float(max(y0, y1))] # ditto

        recipe = self.visFr.pipeline.recipe
        ftable = FilterTable(recipe, inputName=self.visFr.pipeline.selectedDataSourceKey,
                                  outputName='selectedROI', filters=filters)
        if not ftable.configure_traits(kind='modal'):
            return

        recipe.add_module(ftable)
        recipe.execute()

    def OnImROIFT(self, event):
        from PYME.DSView import dsviewer
        selection = selectWithDialog(list(dsviewer.openViewers.keys()))
        if selection is not None:
            img = dsviewer.openViewers[selection].image
        else:
            return
        if img.mdh.getOrDefault('Filter.Keys',None) is None:
            logger.debug('no Filter.Keys in image metadata')
            return

        filters = {}
        filters['x'] = list(img.mdh['Filter.Keys']['x'])
        filters['y'] = list(img.mdh['Filter.Keys']['y'])

        recipe = self.visFr.pipeline.recipe
        ftable = FilterTable(recipe, inputName=self.visFr.pipeline.selectedDataSourceKey,
                                  outputName='selectedROI', filters=filters)
        if not ftable.configure_traits(kind='modal'):
            return

        recipe.add_module(ftable)
        recipe.execute()
        

def Plug(visFr):
    """Plugs this module into the gui"""
    SelectROIFT(visFr)
