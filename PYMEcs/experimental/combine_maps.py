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
        mapname = gen_sCMOS_maps.mkDefaultPath('dark', combinedMap.mdh)
    else:
        mapname = gen_sCMOS_maps.mkDefaultPath('variance', combinedMap.mdh)

    # on windows we may need to pass the full path to defaultFile to force selecting the
    # directory we want; otherwise the last used directory may be used on some
    # windows 7+ installs, see also https://forums.wxwidgets.org/viewtopic.php?t=44404
    import platform
    if platform.system() == 'Windows':
        fname = mapname
    else:
        fname = os.path.basename(mapname)
    map_dlg = wx.FileDialog(parentWindow, message="Save dark map as...",
                            defaultDir=os.path.dirname(mapname),
                            defaultFile=fname, 
                            wildcard='TIFF (*.tif)|*.tif', 
                            style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
    
    if map_dlg.ShowModal() == wx.ID_OK:
        mapfn = map_dlg.GetPath()
        combinedMap.Save(filename=mapfn)

from PYME.recipes.batchProcess import bake
from PYME.recipes import modules
import os

def on_bake(image, parentWindow=None, glCanvas=None):
    with wx.FileDialog(parentWindow, "Choose all series", wildcard='H5 (*.h5)|*.h5',
                       style=wx.FD_OPEN | wx.FD_MULTIPLE) as dialog:

        if dialog.ShowModal() == wx.ID_CANCEL:
            return

        filelist = dialog.GetPaths()
        
    inputGlobs = {'input' : filelist}
    map_dir = os.path.dirname(filelist[0])
    output_dir = os.path.join(map_dir,'analysis')
     
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    recipe_str = """
- processing.DarkAndVarianceMap:
    input: input
    output_dark: dark
    output_variance: variance
- output.ImageOutput:
    filePattern: '{output_dir}/{file_stub}_dark.tif'
    inputName: dark
- output.ImageOutput:
    filePattern: '{output_dir}/{file_stub}_variance.tif'
    inputName: variance
"""

    try:
        recipe = modules.ModuleCollection.fromYAML(recipe_str)
    except AttributeError:
        from PYME.recipes.recipe import Recipe
        recipe = Recipe.fromYAML(recipe_str)
        
    bake(recipe, inputGlobs, output_dir)


def Plug(dsviewer):
    dsviewer.AddMenuItem(menuName='Experimental>Map Tools', itemName='Analyse tiled ROI Map Series',
                         itemCallback = lambda e : on_bake(dsviewer.image, dsviewer, dsviewer.glCanvas))
    dsviewer.AddMenuItem(menuName='Experimental>Map Tools', itemName='Combine tiled ROI Maps',
                         itemCallback = lambda e : on_map(dsviewer.image, dsviewer, dsviewer.glCanvas))
