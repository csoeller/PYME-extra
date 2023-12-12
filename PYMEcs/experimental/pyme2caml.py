import wx
import numpy as np

def pipelineSaveCSV(pipeline,filename,keys):
    pipeline.save_txt(filename,keys)


def p_fileinfo(pipeline,keys):
    xlim = [float(pipeline['x'].min()),float(pipeline['x'].max())]
    ylim = [float(pipeline['y'].min()),float(pipeline['y'].max())]
    fileinfo = dict(xCol = keys.index('x'),
                    yCol= keys.index('y'),
                    ChanIDCol= None,
                    ClusMembershipIDCol = None,
                    LabelIDCol = None,
                    UIDCol = None)
    fileinfo.update(dict(InputFileDelimiter=',',InputFileExt='.csv'))
    fileinfo.update(dict(AutoAxes=False,
                         AutoAxesNearest=1000,
                         ImageSize=[xlim[1]-xlim[0],ylim[1]-ylim[0],0],
                         DataScale=1,
                         xMin=xlim[0],
                         xMax=xlim[1],
                         yMin=ylim[0],
                         yMax=ylim[1]))
    fileinfo.update(dict(ClosestFriend=1,
                         FurthestFriend=100))
    return fileinfo


def pipelineSaveJSON(pipeline,filename,keys):
    import json
    fileinfo = p_fileinfo(pipeline,keys)
    
    with open(filename,"w+") as file:
        json.dump(fileinfo, file,skipkeys=True,indent=3)


class IOcaml:
    """

    """
    def __init__(self, visFr):
        self.visFr = visFr
        self.pipeline = visFr.pipeline

        visFr.AddMenuItem('Experimental>Deprecated>IO', "Save event data to CSV and JSON for CAML",self.OnSaveCaml)

    def OnSaveCaml(self,event):
        import os
        pipeline = self.pipeline
        keys = ['x','y','t','A','error_x','error_y']
        if 'sig' in pipeline.keys():
            keys.append('sig')
        if 'dbscanClumpID' in pipeline.keys():
            keys.append('dbscanClumpID')
        if 'nPhotons' in pipeline.keys():
            keys.append('nPhotons')
        if 'nchi2' in pipeline.keys():
            keys.append('nchi2')
        
        filename = wx.FileSelector('Save events for CAML (select basename)...',
                                   wildcard="CSV files (*.csv)|*.csv", 
                                   flags = wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
                                   
        if not filename == '':
            base, ext = os.path.splitext(filename)
            pipelineSaveCSV(pipeline,base+".csv",keys)
            pipelineSaveJSON(pipeline,base+".json",keys)

def Plug(visFr):
    """Plugs this module into the gui"""
    visFr.iocaml = IOcaml(visFr)
