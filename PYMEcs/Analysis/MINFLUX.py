from scipy.stats import binned_statistic
from PYMEcs.IO.MINFLUX import get_stddev_property
import numpy as np
import matplotlib.pyplot as plt

def plot_stats_minflux(deltas, durations, tdiff, tdmedian, efo_or_dtovertime, times, showTimeAverages=False, dsKey=None):
    
    fig, (ax1, ax2) = plt.subplots(2, 2)
    h = ax1[0].hist(deltas,bins=40)
    dtmedian = np.median(deltas)
    ax1[0].plot([dtmedian,dtmedian],[0,h[0].max()])
    # this is time between one dye molecule and the next dye molecule being seen
    ax1[0].set_xlabel('time between traces (TBT) [s]')
    ax1[0].text(0.95, 0.8, 'median %.2f' % dtmedian, horizontalalignment='right',
             verticalalignment='bottom', transform=ax1[0].transAxes)

    h = ax1[1].hist(durations,bins=40)
    durmedian = np.median(durations)
    ax1[1].plot([durmedian,durmedian],[0,h[0].max()])
    ax1[1].set_xlabel('duration of "traces" [s]')
    ax1[1].text(0.95, 0.8, 'median %.2f' % durmedian, horizontalalignment='right',
             verticalalignment='bottom', transform=ax1[1].transAxes)

    h = ax2[0].hist(tdiff,bins=50,range=(0,0.1))
    ax2[0].plot([tdmedian,tdmedian],[0,h[0].max()])
    # these are times between repeated localisations of the same dye molecule
    ax2[0].set_xlabel('time between localisations in same trace [s]')
    ax2[0].text(0.95, 0.8, 'median %.2f' % tdmedian, horizontalalignment='right',
             verticalalignment='bottom', transform=ax2[0].transAxes)


    if showTimeAverages:
        ax2[1].plot(times,efo_or_dtovertime)
        ax2[1].set_xlabel('TBT running time average [s]')
        ax2[1].set_ylim([0, None])
    else:
        h = ax2[1].hist(1e-3*efo_or_dtovertime,bins=100,range=(0,200))
        # ax2[0].plot([tdmedian,tdmedian],[0,h[0].max()])
        ax2[1].set_xlabel('efo (photon rate kHz)')
        #ax2[0].text(0.95, 0.8, 'median %.2f' % tdmedian, horizontalalignment='right',
        #         verticalalignment='bottom', transform=ax2[0].transAxes)
    if dsKey is not None:
        plt.suptitle('Location rate analysis from datasource %s' % dsKey)
    plt.tight_layout()


# this function assumes a pandas dataframe
# the pandas frame should generally be generated via the function minflux_npy2pyme from PYMEcs.IO.MINFLUX
def analyse_locrate_pdframe(datain,use_invalid=False,showTimeAverages=True):

    if np.any(datain['vld'] < 1):
        data = datain[datain['vld'] >= 1]
        has_invalid = True
    else:
        data = datain
        has_invalid = False

    # we replace the possibly non-sequential trace ids from MINFLUX data with a set of sequential ids
    # this works better for our calculations below when we assume contiguous indices
    # for bin creation for binned_statistic calls
    uids,revids = np.unique(data['tid'].values,return_inverse=True)
    ids = np.arange(1,uids.size+1,dtype='int32')[revids]
    counts = get_stddev_property(ids,data['tid'].values,statistic='count')

    bins = np.arange(int(ids.max())+1) + 0.5
    startindex, bin_edges, binnumber = binned_statistic(ids,data.index,statistic='min', bins=bins)
    endindex, bin_edges, binnumber = binned_statistic(ids,data.index,statistic='max', bins=bins)
    
    if has_invalid and use_invalid:
        # this tries to implement the way to compute both trace durations and time between traces
        # as described in Ostersehlt, L.M. et al. (2022) ‘DNA-PAINT MINFLUX nanoscopy’, Nature Methods, 19(9), pp. 1072–1075.
        # Available at: https://doi.org/10.1038/s41592-022-01577-1.
        # TODO: still needs proper checking if this is correct. Also needs potential BOUNDS CHECKING as we currently assume
        # there are invalid localisations BOTH before the first valid locs and after the last valid locs
        durations = datain.loc[endindex+1,'tim'].values - datain.loc[startindex,'tim'].values
        deltas = datain.loc[startindex[1:],'tim'].values - datain.loc[endindex[:-1]+1,'tim'].values
    else:
        durations = data.loc[endindex,'tim'].values - data.loc[startindex,'tim'].values
        deltas = data.loc[startindex,'tim'][1:].values-data.loc[endindex,'tim'][:-1].values
    # note that we need to convert to numpy here using the values attribute, otherwise we run into an issue
    # with the 'index' mucking up how the rows are subtracted against each other
    tdiff = data['tim'].values[1:]-data['tim'].values[:-1]
    tdsmall = tdiff[tdiff <= 0.1]
    tdmedian = np.median(tdsmall)

    start_times = data.loc[startindex,'tim'][:-1].values # we use those for binning the deltas, we discard final time to match size of deltas
    
    if has_invalid and use_invalid:
        durations_proper = durations
    else:
        durations_proper = durations + tdmedian # we count one extra localisation, using the median duration
        # the extra is because we leave at least one localisation out from the total timing when we subtract ends-starts

    if showTimeAverages:
        delta_averages, bin_edges, binnumber = binned_statistic(start_times,deltas,statistic='mean', bins=50)
        delta_av_times = 0.5*(bin_edges[:-1] + bin_edges[1:]) # bin centres
        plot_stats_minflux(deltas, durations_proper, tdiff, tdmedian, delta_averages, delta_av_times, showTimeAverages=True)
    else:
        plot_stats_minflux(deltas, durations_proper, tdiff, tdmedian, data['efo'], None)


# similar version but now using a pipeline
def analyse_locrate(data,datasource='Localizations',showTimeAverages=True):
    curds = data.selectedDataSourceKey
    data.selectDataSource(datasource)
    bins = np.arange(int(data['clumpIndex'].max())+1) + 0.5
    starts, bin_edges, binnumber = binned_statistic(data['clumpIndex'],data['tim'],statistic='min', bins=bins)
    ends, bin_edges, binnumber = binned_statistic(data['clumpIndex'],data['tim'],statistic='max', bins=bins)
    durations = ends - starts
    deltas = starts[1:]-ends[:-1]
    tdiff = data['tim'][1:]-data['tim'][:-1]
    tdsmall = tdiff[tdiff <= 0.1]
    tdmedian = np.median(tdsmall)
    durations_proper = durations + tdmedian # we count one extra localisation, using the median duration
    # the extra is because we leave at least one localisation out from the total timing when we subtract ends-starts
    data.selectDataSource(curds)
    
    if showTimeAverages:
        delta_averages, bin_edges, binnumber = binned_statistic(starts[:-1],deltas,statistic='mean', bins=50)
        delta_av_times = 0.5*(bin_edges[:-1] + bin_edges[1:]) # bin centres
        plot_stats_minflux(deltas, durations_proper, tdiff, tdmedian, delta_averages, delta_av_times, showTimeAverages=True, dsKey = datasource)
    else:
        plot_stats_minflux(deltas, durations_proper, tdiff, tdmedian, data['efo'], None, dsKey = datasource)

