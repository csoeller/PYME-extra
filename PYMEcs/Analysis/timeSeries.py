import numpy as np
import matplotlib.pyplot as plt

import logging
logger = logging.getLogger(__name__)

def interlaceTraces(traces):
    return np.vstack(traces).transpose().flatten()

# given a series of detected events as input
# construct a single channel trace for plotting
# this works by constructing a time point series and
# a pulse series (0 and 1 values) that goes up and down
# as dictated by the event series passed in
#
# algorithmicly it detects where series of events in consecutive
# groups begin and end and injects plotting points at those
# places
#
# the returned two arrays can be directly passed to a plot command
# to produce a "single channel"-like trace
def generateSeries(t):

    if t.shape[0] <= 1:
        return ([t[0],t[0],t[0]+1,t[0]+1],[0,1,1,0])
    
    dts = t[1:]-t[0:-1]-1
    dtg = dts[dts>0]
    nts = dtg.shape[0]

    idx, = np.where(dts>0)
    one = np.ones_like(idx)
    z = np.zeros_like(idx)

    tSer = interlaceTraces((t[idx],t[idx],t[idx+1],t[idx+1]))
    pulseSer = interlaceTraces((one,z,z,one))

    tSerAll = np.hstack(([t[0],t[0]],tSer,[t[-1],t[-1]]))
    pulseAll = np.hstack(([0,1],pulseSer,[1,0]))

    return (tSerAll, pulseAll)

# this one uses the functionality of generateSeries
# for all 'clumps' identified by the same clumpIndex
#
# it therefore requires that the 'Track single molecule trajectories"
# has been run; it must be run on the original 'Localisations' data source,
# *not* the coalesced data source
#
# by default it starts a new figure
# the plot line colour changes between clumps
# in an effectively random fashion
def plotClumpSeries(t,ci,newFig = True):
    cin = np.unique(ci)
    if newFig:
        plt.figure()

    for c in cin:
        tci = t[ci == c]
        tp, p = generateSeries(tci)
        plt.plot(tp, p)

    plt.ylim(-0.5,1.5)
