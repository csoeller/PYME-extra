import wx
import os.path


def on_map(image, parentWindow=None, glCanvas=None):
    from PYME.Analysis import gen_sCMOS_maps
    from PYME.DSView import ViewIm3D
    from PYMEcs.Analysis.MapUtils import combine_maps

    # combine maps with dialogue here
    # also show valid map

    with wx.FileDialog(parentWindow, "Choose maps", wildcard='TIFF (*.tif)|*.tif',
                       style=wx.FD_OPEN | wx.FD_MULTIPLE) as dialog:

        if dialog.ShowModal() == wx.ID_CANCEL:
            return

        filelist = dialog.GetPaths()

        combinedMap, vMap = combine_maps(filelist,return_validMap=True)

        if combinedMap.mdh['CameraMap.Type'] == 'mean':
            mapType = 'dark'
        elif combinedMap.mdh['CameraMap.Type'] == 'variance':
            mapType = 'variance'
            
    if mapType == 'dark':
        ViewIm3D(combinedMap, title='Dark Map', parent=parentWindow, glCanvas=glCanvas)
    else:
        ViewIm3D(combinedMap, title='Variance Map', parent=parentWindow, glCanvas=glCanvas)

    ViewIm3D(vMap, title='Valid Regions', parent=parentWindow, glCanvas=glCanvas)
    
    if mapType == 'dark':
        mapname = gen_sCMOS_maps.mkDefaultPath('dark', image.mdh)
    else:
        mapname = gen_sCMOS_maps.mkDefaultPath('variance', image.mdh)

    map_dlg = wx.FileDialog(parentWindow, message="Save dark map as...",
                             defaultDir=os.path.dirname(mapname),
                             defaultFile=os.path.basename(mapname), 
                             wildcard='TIFF (*.tif)|*.tif', 
                             style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
    
    if map_dlg.ShowModal() == wx.ID_OK:
        mapfn = map_dlg.GetPath()
        combinedMap.Save(filename=mapfn)


def Plug(dsviewer):
    dsviewer.AddMenuItem(menuName='Experimental>Map Tools', itemName='Combine tiled ROI Maps',
                         itemCallback = lambda e : on_map(dsviewer.image, dsviewer, dsviewer.glCanvas))
