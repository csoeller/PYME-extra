from traits.api import HasTraits, Str, Int, CStr, List, Enum, Float
from traitsui.api import View, Item, Group
from traitsui.menu import OKButton, CancelButton, OKCancelButtons


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


def queue_roi_series(scope):
    cam = scope.cam

    args = {'state': {'Camera.ROI' : [50, 50, 200, 200]}}
    scope.actions.QueueAction('state.update', args)
    args = {'maxFrames': 500, 'stack': False}
    scope.actions.QueueAction('spoolController.StartSpooling', args)
    args = {'state': {'Camera.ROI' : [100, 100, 250, 250]}}
    scope.actions.QueueAction('state.update', args)
    args = {'maxFrames': 500, 'stack': False}
    scope.actions.QueueAction('spoolController.StartSpooling', args)    
    args = {'state': {'Camera.ROI' : [0, 0, 256, 256]}}
    scope.actions.QueueAction('state.update', args)
    
    # in future we might code this as:
    #
    # calib = [actions.SpoolSeries(maxFrames=500, stack=False, 
    #                              state={'Camera.ROI' : [50, 50, 200, 200]}),
    #          actions.SpoolSeries(maxFrames=500, stack=False,
    #                              state={'Camera.ROI' : [100, 100, 250, 250]}),
    # ]

    # scope.actions.queue_actions(calib)


def check_roi(x0,y0,x1,y1,width=None, height=None):
    if x0<0:
        x0=0
    if y0<0:
        y0=0
    if x1>(width-1):
        x1 = width-1
    if y1>(height-1):
        y1 = height-1
    return [x0,y0,x1,y1]


def spoolSeries(scope, maxFrames=500, stack=False, state=None):
    if state is not None:
        args = {'state': state}
        scope.actions.QueueAction('state.update', args)
        
    args = {'maxFrames': maxFrames, 'stack': stack}
    scope.actions.QueueAction('spoolController.StartSpooling', args)


def setState(scope,state):
    args = {'state': state}
    scope.actions.QueueAction('state.update', args)

    
# ToDo - refine implementation as action interface improves
#      - add help explaining what is going on, using traits help infrastructure (i.e. in class ChipROI)
def camera_chip_calibration_series(scope):

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

    x0 = 0
    y0 = 0
    x1 = rsz-1
    y1 = rsz-1
    
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
                                 roi[2]-roi[0]+1,
                                 roi[3]-roi[1]+1,
                                 linewidth=1,edgecolor=cols[i %2],facecolor='none')
        ax.add_patch(rect)


    # actually queue series
    for roi in rois:
        spoolSeries(scope, maxFrames=chipROI.numberOfFrames, stack=False,
                    state={'Camera.ROI' : roi})

    # set back to original ROI
    setState(scope,state={'Camera.ROI' : curROI})
