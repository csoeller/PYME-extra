import logging
logger = logging.getLogger(__file__)
import wx

class dsviewerUtils:
    def __init__(self, dsviewer):
        self.dsviewer = dsviewer
        dsviewer.AddMenuItem('Experimental>Utils',
                             'Save recipe from metadata',
                             self.OnSaveRecipeFromMDH,
                             helpText='try to save a recipe as yaml from image metadata')

    def OnSaveRecipeFromMDH(self, event=None):
        with wx.FileDialog(self.dsviewer, "Choose recipe file name", wildcard='YAML (*.yaml)|*.yaml',
                           style=wx.FD_SAVE) as dialog:

            if dialog.ShowModal() == wx.ID_CANCEL:
                return

            filename = dialog.GetPath()

        from PYMEcs.misc.utils import recipe_from_mdh
        recipe = recipe_from_mdh(self.dsviewer.image.mdh)
        if recipe is not None:
            with open(filename,"w") as fi:
                print(recipe,file=fi)


def Plug(dsviewer):
    """Plugs this module into the gui"""
    dsviewerUtils(dsviewer)
