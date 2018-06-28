from PYME.recipes.base import register_module, ModuleBase, OutputModule
from PYME.recipes.traits import Input, Output, Float, Enum, CStr, Bool, Int, DictStrStr, File, Button
from traitsui.file_dialog import save_file

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
               
        dlg = FileDialog(action='save as')
        if dlg.open() == OK:
            self.fileName = dlg.path
        
    @property
    def default_view(self):
        from traitsui.api import View, Group, Item, HGroup
        return View(
            Group(Item('inputName'),
                  HGroup(
                      Item('saveAs', show_label=False),
                      '_',
                      Item('fileName', style='readonly', springy=True)
                  )
            ), buttons=['OK'])
