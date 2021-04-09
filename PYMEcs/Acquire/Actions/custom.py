from traits.api import HasTraits, Str, Int, CStr, List, Enum, Float, Bool
from traitsui.api import View, Item, Group
from traitsui.menu import OKButton, CancelButton, OKCancelButtons

from PYMEcs.misc.guiMsgBoxes import YesNo
from PYME.IO.FileUtils.nameUtils import numToAlpha


# could be called in init file as:

# @init_gui('ROI Calibration')
# def roi_calibration(MainFrame, scope):
#
#     def roi_action_callback(event=None):
#         from PYMEcs.Acquire.Actions.custom import queue_calibration_series
#         queue_calibration_series(scope)
    
#     MainFrame.AddMenuItem('Calibration', 'Camera Maps>Sub ROIs', roi_action_callback)

class ChipROI(HasTraits):
    roiSize = Int(256)
    overlap = Int(20)
    numberOfFrames = Int(500)
    checkBeforeQueuingActions = Bool(True)


# commenting out the below which was just a quick trial how to do this
# keeping in for reference right now

# def queue_roi_series(scope):
#     cam = scope.cam

#     args = {'state': {'Camera.ROI' : [50, 50, 200, 200]}}
#     scope.actions.QueueAction('state.update', args)
#     args = {'maxFrames': 500, 'stack': False}
#     scope.actions.QueueAction('spoolController.StartSpooling', args)
#     args = {'state': {'Camera.ROI' : [100, 100, 250, 250]}}
#     scope.actions.QueueAction('state.update', args)
#     args = {'maxFrames': 500, 'stack': False}
#     scope.actions.QueueAction('spoolController.StartSpooling', args)    
#     args = {'state': {'Camera.ROI' : [0, 0, 256, 256]}}
#     scope.actions.QueueAction('state.update', args)
    
#     # in future we might code this as:
#     #
#     # calib = [actions.SpoolSeries(maxFrames=500, stack=False, 
#     #                              state={'Camera.ROI' : [50, 50, 200, 200]}),
#     #          actions.SpoolSeries(maxFrames=500, stack=False,
#     #                              state={'Camera.ROI' : [100, 100, 250, 250]}),
#     # ]

#     # scope.actions.queue_actions(calib)


def check_roi(x0,y0,x1,y1,width=None, height=None):
    if x0<0:
        x0=0
    if y0<0:
        y0=0
    if x1>(width):
        x1 = width
    if y1>(height):
        y1 = height
    return [x0,y0,x1,y1]


def spoolSeries(scope, maxFrames=500, stack=False, state=None,
                method=None, subdirectory=None):
    if state is not None:
        args = {'state': state}
        scope.actions.QueueAction('state.update', args)

    args = {'settings': {'max_frames': maxFrames,
                         'stack': stack,
                         'subdirectory': subdirectory,
                         'method': method}}
    scope.actions.QueueAction('spoolController.start_spooling', args)


def setState(scope,state):
    args = {'state': state}
    scope.actions.QueueAction('state.update', args)



def get_mapdir(scope):
    import os
    seriesCounter = 0
    subdirectory =  'map_' + numToAlpha(seriesCounter)
    while os.path.exists(scope.spoolController.get_dirname(subdirectory)):
        seriesCounter += 1
        subdirectory = 'map_' + numToAlpha(seriesCounter)
    return subdirectory
        
        
# ToDo - refine implementation as action interface improves
#      - add online help, using traits help infrastructure (i.e. in class ChipROI)
def camera_chip_calibration_series(scope):
    import os

    chipROI = ChipROI()
    if not chipROI.configure_traits(kind='modal'):
            return

    cam = scope.cam
    chipWidth  = cam.GetCCDWidth()
    chipHeight = cam.GetCCDHeight()
    curROI = cam.GetROI()
    
    stepsize = chipROI.roiSize - chipROI.overlap
    rsz = chipROI.roiSize
    
    xsteps = int(chipWidth / stepsize)
    ysteps = int(chipHeight / stepsize)

    rois = []

    # note PYME cam ROI conventions: 1) chip origin starts at 0,0
    #                                2) ROI with coordinates [x0,y0,x1,y1] starts and includes point x0,y0
    #                                   but extends to and *excludes* point x1, y1
    #                                This implies that the ROI width is x1-x0 and height is y1-y0 (used further below)
    x0 = 0
    y0 = 0
    x1 = rsz # note *not* rsz-1, since it excludes x1, y1, see above
    y1 = rsz
    
    for iy in range(0,ysteps):
        for ix in range(0,xsteps):
            rois.append(check_roi(x0+ix*stepsize,y0+iy*stepsize,
                                  x1+ix*stepsize,y1+iy*stepsize,
                                  width = chipWidth, height=chipHeight))

    # show tiling on chip
            
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    import numpy as np

    cols = ['r','g']
 
    fig,ax = plt.subplots(1)
    ax.imshow(np.ones((chipWidth,chipHeight)))
    for i, roi in enumerate(rois):
        rect = patches.Rectangle((roi[0],roi[1]),
                                 roi[2]-roi[0], # remember we get ROI width by difference of corner points
                                 roi[3]-roi[1], # same for height
                                 linewidth=1,edgecolor=cols[i %2],facecolor='none')
        ax.add_patch(rect)
        cx = 0.5*(roi[0]+roi[2])
        cy = 0.5*(roi[1]+roi[3])
        plt.text(cx,cy,'%d' % i,c='w') # plot in white to ensure we see on top of darkish chip plot


    if chipROI.checkBeforeQueuingActions:
        if not YesNo(None, "Will use %d ROIS.\nProceed with running ROI actions?" % len(rois), caption='Proceed'):
            return

    mapdir = get_mapdir(scope)
    # actually queue series
    for roi in rois:
        spoolSeries(scope, maxFrames=chipROI.numberOfFrames, stack=False,
                    state={'Camera.ROI' : roi}, subdirectory=mapdir,
                    method='FILE')

    # set back to original ROI
    setState(scope,state={'Camera.ROI' : curROI})
    
