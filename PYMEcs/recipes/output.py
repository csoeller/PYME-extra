from PYME.recipes.base import register_module, ModuleBase, OutputModule
from PYME.recipes.traits import Input, Output, Float, Enum, CStr, Bool, Int, DictStrStr, File, Button
#from traitsui.file_dialog import save_file # this call seems to shut down the gui of VisGUI

import numpy as np
import pandas as pd

@register_module('CSVOutputFileBrowse2')
class CSVOutputFileBrowse2(OutputModule):
    """
    Save tabular data as csv. This module uses a File Browser to set the fileName

    Parameters
    ----------

    inputName : basestring
        the name (in the recipe namespace) of the table to save.

    fileName : File
        a full path to the file

    Notes
    -----

    We convert the data to a pandas `DataFrame` and uses the `to_csv`
    method to save.

    This version of the output module uses the pyface.FileDialog.
    """

    inputName = Input('output')
    fileName = File('out.csv')
    saveAs = Button('Save as...')
    
    def save(self, namespace, context={}):
        """
        Save recipes output(s) to CSV

        Parameters
        ----------
        namespace : dict
            The recipe namespace
        context : dict
            Information about the source file to allow pattern substitution to generate the output name. At least
            'basedir' (which is the fully resolved directory name in which the input file resides) and
            'filestub' (which is the filename without any extension) should be resolved.

        Returns
        -------

        """

        out_filename = self.fileName
        v = namespace[self.inputName]

        if not isinstance(v, pd.DataFrame):
            v = v.toDataFrame()
                
        v.to_csv(out_filename)

    def _saveAs_changed(self):
        """ Handles the user clicking the 'Save as...' button.
        """
        from pyface.api import FileDialog, OK
        import os
        dirname = os.path.dirname(self.fileName)
        filename = os.path.basename(self.fileName)
        if not dirname:
            dirname = os.getcwd()
              
        dlg = FileDialog(action='save as', default_directory=dirname, default_filename=filename)
        if dlg.open() == OK:
            self.fileName = dlg.path
        
    @property
    def default_view(self):
        from traitsui.api import View, Group, Item, HGroup
        from PYME.ui.custom_traits_editors import CBEditor

        editable = self.class_editable_traits()
        inputs = [tn for tn in editable if tn.startswith('input')]
        return View(
            Group(Item('inputName', editor=CBEditor(choices=self._namespace_keys)),
                  HGroup(
                      Item('saveAs', show_label=False),
                      '_',
                      Item('fileName', style='readonly', springy=True)
                  )
            ), buttons=['OK'])

@register_module('CSVOutputFileBrowse')
class CSVOutputFileBrowse(OutputModule):
    """
    Save tabular data as csv. This module uses a File Browser to set the fileName

    Parameters
    ----------

    inputName : basestring
        the name (in the recipe namespace) of the table to save.

    fileName : File
        a full path to the file

    Notes
    -----

    We convert the data to a pandas `DataFrame` and uses the `to_csv`
    method to save.

    This version of the output module uses the wx.FileDialog.

    """

    inputName = Input('output')
    fileName = File('output.csv')
    saveAs = Button('Save as...')
    
    def save(self, namespace, context={}):
        """
        Save recipes output(s) to CSV

        Parameters
        ----------
        namespace : dict
            The recipe namespace
        context : dict
            Information about the source file to allow pattern substitution to generate the output name. At least
            'basedir' (which is the fully resolved directory name in which the input file resides) and
            'filestub' (which is the filename without any extension) should be resolved.

        Returns
        -------

        """

        out_filename = self.fileName
        v = namespace[self.inputName]

        if not isinstance(v, pd.DataFrame):
            v = v.toDataFrame()
                
        v.to_csv(out_filename)

    def _saveAs_changed(self):
        """ Handles the user clicking the 'Save as...' button.
        """
        import wx
        import os
        dirname = os.path.dirname(self.fileName)
        filename = os.path.basename(self.fileName)
        if not dirname:
            dirname = os.getcwd()
        dlg = wx.FileDialog(None, "Save as...", dirname, filename, "*.csv",
                            wx.SAVE|wx.OVERWRITE_PROMPT)
        result = dlg.ShowModal()
        inFile = dlg.GetPath()
        dlg.Destroy()

        if result == wx.ID_OK:          #Save button was pressed
            self.fileName = inFile
        
    @property
    def default_view(self):
        from traitsui.api import View, Group, Item, HGroup
        from PYME.ui.custom_traits_editors import CBEditor

        editable = self.class_editable_traits()
        inputs = [tn for tn in editable if tn.startswith('input')]
        return View(
            Group(Item('inputName', editor=CBEditor(choices=self._namespace_keys)),
                  HGroup(
                      Item('saveAs', show_label=False),
                      '_',
                      Item('fileName', style='readonly', springy=True)
                  )
            ), buttons=['OK'])
