## Copy and paste this code in PYME/PYME-extra/PYMEcs/experimental/selectROIfilterTable.py
## Don't forget that if this is the first time installing then you will also need to pip install openpyxl
## else none of this will work. No error messages will appear in the console to warn this is an import problem. 

import numpy as np
import pandas as pd
import openpyxl
import wx # new for save events
from pandas import DataFrame, read_csv
from PYME.recipes.tablefilters import FilterTable


def pipelineSaveCSV(pipeline,filename,keys):  #Required to sit outside of class?
    pipeline.save_txt(filename,keys)


class SaveROI:


    def __init__(self, visFr):
        self.visFr = visFr
        self.pipeline = visFr.pipeline

        visFr.AddMenuItem('Alex', "2. Save coordinates\tF1", self.FindMe)


    def FindMe(self, event):
        # New idea to get this out. Open a .csv or xlsx
        # Append new entry. #Save
        try:
            #old glcanvas
            x0, y0 = self.visFr.glCanvas.selectionStart[0:2]
            x1, y1 = self.visFr.glCanvas.selectionFinish[0:2]
        except AttributeError:
            #new glcanvas
            x0, y0 = self.visFr.glCanvas.selectionSettings.start[0:2]
            x1, y1 = self.visFr.glCanvas.selectionSettings.finish[0:2]

            
        try: 
            openfile1 = 'coords.csv'
            df = pd.read_csv(openfile1)
        except:
            df = pd.DataFrame(columns=["x0", "x1", "y0", "y1"])
            df.to_csv('coords.csv', index=True)
            openfile1 = 'coords.csv'
            df = pd.read_csv(openfile1)

        x0list, x1list, y0list, y1list = [], [], [], [] 
        for i in range (df.x0.size):#
            x0list.append(df.x0[i])
            x1list.append(df.x1[i])
            y0list.append(df.y0[i])
            y1list.append(df.y1[i])

        x0list.append(x0)
        x1list.append(x1)
        y0list.append(y0)
        y1list.append(y1)

        exportdata = pd.DataFrame({
            'x0' : x0list,
            'x1' : x1list,
            'y0' : y0list,
            'y1' : y1list,
            })

        exportdata.to_csv('coords.csv', index=True)


class SelectandfixFiducial:
    

    def __init__(self, visFr):
        self.visFr = visFr
        self.pipeline = visFr.pipeline

        visFr.AddMenuItem('Alex', "1. Select and correct the fiducial \tF4", self.OnSelROIFT)
    

    def OnSelROIFT(self, event):
        print("Howdy")
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
                                  outputName='Fiducial', filters=filters)
        if not ftable.configure_traits(kind='modal'):
            return

        recipe.add_module(ftable)
        recipe.execute()

    
        from PYMEcs.recipes import localisations
        recipe = self.pipeline.recipe
        # change defaults back to median filter
        ftrack = localisations.FiducialTrack(recipe, inputName='Fiducial',
                                             filterMethod='Median',
                                             timeWindow=350,
                                             outputName='fiducialAdded')
        if not ftrack.configure_traits(kind='modal'):
            return
        recipe.add_module(ftrack)
        recipe.add_module(localisations.FiducialApplyFromFiducials(recipe, inputData=self.pipeline.selectedDataSourceKey,
                                                                   inputFiducials='fiducialAdded',
                                                                   outputName='fiducialApplied',
                                                                   outputFiducials='corrected_fiducials'))
        recipe.execute()
        self.pipeline.selectDataSource('fiducialApplied')

class SelectROIFT:
    

    def __init__(self, visFr):
        self.visFr = visFr

        visFr.AddMenuItem('Experimental>View', "Add ROI FilterTable module from selection", self.OnSelROIFT)


class SelectROIFT:
    def __init__(self, visFr):
        self.visFr = visFr

        visFr.AddMenuItem('Experimental>View', "Add ROI FilterTable module from selection", self.OnSelROIFT)

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


class EventsSavetoCSV:
    """

    """
    def __init__(self, visFr):
        self.visFr = visFr
        self.pipeline = visFr.pipeline

        visFr.AddMenuItem('Alex', "5. Save MINFLUX events",self.OnSaveCaml)

    def OnSaveCaml(self,event):
        import os
        pipeline = self.pipeline
        keys = ['x','y']
        if 't' in pipeline.keys():
            keys.append('t')
            keys.append('cfr')
            keys.append('clumpIndex')
            keys.append('clumpSize')
            keys.append('efo')
            keys.append('error_x')
            keys.append('error_y')
            keys.append('fbg')
            keys.append('focus')
            keys.append('x_raw')
            keys.append('y_raw')

        if 'z' in pipeline.keys():
            keys.append('z')
            keys.append('error_z')
            keys.append('z_raw')
        if 'NNclass' in pipeline.keys():
            keys.append('NNclass')
            keys.append('NNdist')
            keys.append('NNdist2')
        if 'sig' in pipeline.keys():
            keys.append('sig')
        if 'dbscanClumpID' in pipeline.keys():
            keys.append('dbscanClumpID')
            keys.append('dbscanClumpSize')
        if 'nPhotons' in pipeline.keys():
            keys.append('nPhotons')
        if 'nchi2' in pipeline.keys():
            keys.append('nchi2')
        if 'dyeID' in pipeline.keys():
            keys.append('dyeID')
        if 'qid' in pipeline.keys():
            keys.append('qid')
            keys.append('ryrid')
            keys.append('su_num')
            keys.append('suid')
            keys.append('muid')

        filename = wx.FileSelector('Save MINFLUX events (select basename)...',
                                   wildcard="CSV files (*.csv)|*.csv", 
                                   flags = wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
                                   
        if not filename == '':
            base, ext = os.path.splitext(filename)
            pipelineSaveCSV(pipeline,base+".csv",keys)






def Plug(visFr):
    """Plugs this module into the gui"""
    SelectROIFT(visFr)
    SelectandfixFiducial(visFr)
    SaveROI(visFr)
    EventsSavetoCSV(visFr)
