import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binned_statistic

from PYMEcs.IO.MINFLUX import get_stddev_property
from PYMEcs.pyme_warnings import warn


def propcheck_density_stats(ds,warning=True):
    for prop in ['clst_area','clst_vol','clst_density','clst_stdz']:
        if prop not in ds.keys():
            if warning:
                warn("required property %s not in data source" % prop)
            return False
    return True
            
def density_stats(ds,objectID='dbscanClumpID'):    
    uids, idx = np.unique(ds[objectID],return_index=True)
    area = ds['clst_area'][idx]
    vol = ds['clst_vol'][idx]
    dens = ds['clst_density'][idx]
    sz = ds['clst_stdz'][idx]

    return area, vol, dens, sz

def plot_density_stats(ds,objectID='dbscanClumpID',scatter=False):
    if not propcheck_density_stats(ds):
        return
    area, vol, dens, sz = density_stats(ds,objectID=objectID)
    fig, (ax0,ax1) = plt.subplots(2,2)
    if not scatter:
        ax0[0].boxplot(dens,labels=['Density'])
        ax0[1].boxplot(area,labels=['Area'])
        ax1[0].boxplot(vol/1e6,labels=['Volume'])
        ax1[1].boxplot(sz,labels=['Stddev Z'])
    else:
        bp_dict = ax0[0].scattered_boxplot(dens,labels=['Density'],showmeans=True)
        for line in bp_dict['means']:
            # get position data for median line
            x, y = line.get_xydata()[0] # top of median line
            # overlay median value
            ax0[0].text(x-0.25, 1.05*y, '%.0f' % y,
                        horizontalalignment='center') # draw above, centered
        for line in bp_dict['medians']:
            # get position data for median line
            x, y = line.get_xydata()[0] # top of median line
            # overlay median value
            ax0[0].text(x-0.25, 0.95*y, '%.0f' % y,
                        horizontalalignment='center',
                        verticalalignment='center') # draw above, centered
        ax0[1].scattered_boxplot(area,labels=['Area'],showmeans=True)
        ax1[0].scattered_boxplot(vol/1e6,labels=['Volume'],showmeans=True)
        ax1[1].scattered_boxplot(sz,labels=['Stddev Z'],showmeans=True)
    plt.tight_layout()

from PYMEcs.misc.matplotlib import boxswarmplot


def plot_density_stats_sns(ds,objectID='dbscanClumpID'):
    if not propcheck_density_stats(ds):
        return
    area, vol, dens, sz = density_stats(ds,objectID=objectID)
    dfdens = pd.DataFrame.from_dict(dict(density=dens))
    dfarea = pd.DataFrame.from_dict(dict(area=area))
    dfvol = pd.DataFrame.from_dict(dict(volume=vol/1e6))
    dfstdz = pd.DataFrame.from_dict(dict(std_z=sz))
    
    fig, (ax0,ax1) = plt.subplots(2,2)
    kwargs = dict(swarmsize=5,width=0.2,annotate_means=True,annotate_medians=True,swarmalpha=0.4)
    boxswarmplot(dfdens,ax=ax0[0],format="%.0f",**kwargs)
    ax0[0].set_ylim(0,1.2*dens.max())
    ax0[0].set_ylabel("#/um^2")
    boxswarmplot(dfarea,ax=ax0[1],format="%.0f",**kwargs)
    ax0[1].set_ylabel("nm^2")
    boxswarmplot(dfvol,ax=ax1[0],format="%.2f",**kwargs)
    ax1[0].set_ylabel("10^-3 um^3")
    boxswarmplot(dfstdz,ax=ax1[1],format="%.1f",**kwargs)
    ax1[1].set_ylabel("nm")
    fig.suptitle('Density stats for %d clusters' % dens.size)
    plt.tight_layout()

    # --- Save the values to csv --- (Alex B addition)
    df = pd.DataFrame({
        "Metric": ["cluster density"],
        "Mean": [np.mean(dens)],
        "Median": [np.median(dens)],
        "Unit": ["#/um^2"]})
    # --- Show the head of df in the console --- (Alex B addition)
    print(df.head())

    # --- Save as a csv file --- (Alex B addition)
    # Save the cluster size file
    #fn = os.path.splitext(ds.mdh.get('MINFLUX.Filename'))[0]
    dirpath, filename = os.path.split(ds.filename)
    fn = ds.mdh.get('MINFLUX.TimeStamp')
    df.to_csv(os.path.join(dirpath, fn + "_clusterDensity.csv"), 
              index=False, header=True)

    # save data CSV file with dbscanClustered
    df_dbscanClustered = ds.dataSources['dbscanClustered'].to_pandas()
    df_dbscanClustered.to_csv(
        os.path.join(dirpath, fn + "_dbscanClustered_data.csv"), index=False, header=True)
    # Used as a reminder for LocRate csv saving
    # print(
    #   f'\ncsv name is: {fn + "_clusterDensity.csv"}\nIf you did not load a session, csv file and figures will be saved on the desktop')
    # By default the file is saved on the Desktop, if a session file is used, it is saved in the same directory as the session file.

    return dens

def plot_stats_minflux(deltas, durations, tdintrace, efo_or_dtovertime, times,
                       showTimeAverages=False, dsKey=None, areaString=None, timestamp=None, dirPath=None):
    from scipy.stats import iqr

    # --- Create the figure (plot with 2x2 subplots) ---
    fig, (ax1, ax2) = plt.subplots(2, 2)
    
    # --- Compute some statistics for TBT plot (top-right = ax1[0]) ---
    dtmedian = np.median(deltas)
    dtmean = np.mean(deltas)
    dtiqr = iqr(deltas,rng=(10, 90)) # we are going for the 10 to 90 % range
    
    h = ax1[0].hist(deltas,bins=40,range=(0,dtmean + 2*dtiqr))
    ax1[0].plot([dtmedian,dtmedian],[0,h[0].max()])
    # this is time between one dye molecule and the next dye molecule being seen
    ax1[0].set_xlabel('time between traces (TBT) [s]')
    ax1[0].text(0.95, 0.8, 'median %.2f s' % dtmedian, horizontalalignment='right',
             verticalalignment='bottom', transform=ax1[0].transAxes)
    ax1[0].text(0.95, 0.7, '  mean %.2f s' % dtmean, horizontalalignment='right',
             verticalalignment='bottom', transform=ax1[0].transAxes)
    if areaString is not None:
        ax1[0].text(0.95, 0.6, areaString, horizontalalignment='right',
             verticalalignment='bottom', transform=ax1[0].transAxes)
    
    # --- Compute some statistics for trace duration plot (top-left = ax1[1]) ---
    durmedian = np.median(durations)
    durmean = np.mean(durations)
    duriqr = iqr(durations,rng=(10, 90))
    h = ax1[1].hist(durations,bins=40,range=(0,durmean + 2*duriqr))
    
    ax1[1].plot([durmedian,durmedian],[0,h[0].max()])
    ax1[1].set_xlabel('duration of "traces" [s]')
    ax1[1].text(0.95, 0.8, 'median %.0f ms' % (1e3*durmedian), horizontalalignment='right',
             verticalalignment='bottom', transform=ax1[1].transAxes)
    ax1[1].text(0.95, 0.7, '  mean %.0f ms' % (1e3*durmean), horizontalalignment='right',
             verticalalignment='bottom', transform=ax1[1].transAxes)
    # ax1[1].set_xlim(0,durmean + 2*duriqr) # superfluous since we are using the range keyword in hist

    # --- Compute some statistics for time between localisations in same trace (bottom-left = ax2[0]) ---
    tdintrace_ms = 1e3*tdintrace
    tdmedian = np.median(tdintrace_ms)
    tdmean = np.mean(tdintrace_ms)
    tdiqr = iqr(tdintrace_ms,rng=(10, 90))
    h = ax2[0].hist(tdintrace_ms,bins=50,range=(0,tdmean + 2*tdiqr))
    
    ax2[0].plot([tdmedian,tdmedian],[0,h[0].max()])
    # these are times between repeated localisations of the same dye molecule
    ax2[0].set_xlabel('time between localisations in same trace [ms]')
    ax2[0].text(0.95, 0.8, 'median %.0f ms' % (tdmedian), horizontalalignment='right',
             verticalalignment='bottom', transform=ax2[0].transAxes)
    ax2[0].text(0.95, 0.7, '  mean %.0f ms' % (tdmean), horizontalalignment='right',
             verticalalignment='bottom', transform=ax2[0].transAxes)
    # ax2[0].set_xlim(0,tdmean + 2*tdiqr) # superfluous since we are using the range keyword in hist

    # --- Compute the TBT running time average (bottom-right = ax2[1]) ---
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

    # --- Calculate dimensions and area of the image ---
    dimension = areaString.split(' ')[1]
    area = float(dimension.split('x')[0])*float(dimension.split('x')[1])
    
    # --- Save medians values to csv --- (Alex B addition)
    df = pd.DataFrame({
        "Metric": ["TBT (Time Between Traces)", "ROI dimension", "ROI area", "Trace Duration", "Time Between Localizations"],
        "Median": [dtmedian.round(2), dimension, area, (durmedian * 1e3).round(0), tdmedian.round(0)],  # Convert seconds -> milliseconds where needed
        "Mean": [dtmean.round(2), dimension, area, (durmean * 1e3).round(0), tdmean.round(0)],  # Convert seconds -> milliseconds where needed
        "IQR": [dtiqr.round(2), dimension, area, (duriqr * 1e3).round(0), tdiqr.round(0)],  # Convert seconds -> milliseconds where needed
        "Unit": ["s","um^2","um^2", "ms", "ms"]    })
    
    # --- Show the head of df in the console --- (Alex B addition)
    print(df.head())

    # --- Save as a csv file --- (Alex B addition)
    # Save the LocRate file
    df.to_csv(os.path.join(dirPath, timestamp + "_locRate.csv"), index=False, header=True)

    # --- Save the figure --- (Alex B addition)
    plt.savefig(os.path.join(dirPath, timestamp + '_locRate.png'),
                dpi=300, bbox_inches='tight')
    print(f'Figure and *.csv file saved to {os.path.join(dirPath, timestamp + "_locRate.png")}')
    
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
    
    dirpath, filename = os.path.split(datain.filename)
    timeStamp = datain.mdh.get('MINFLUX.TimeStamp')
    
    start_times = data.loc[startindex,'tim'][:-1].values # we use those for binning the deltas, we discard final time to match size of deltas
    
    if has_invalid and use_invalid:
        durations_proper = durations
    else:
        durations_proper = durations + tdmedian # we count one extra localisation, using the median duration
        # the extra is because we leave at least one localisation out from the total timing when we subtract ends-starts

    if showTimeAverages:
        delta_averages, bin_edges, binnumber = binned_statistic(start_times,deltas,statistic='mean', bins=50)
        delta_av_times = 0.5*(bin_edges[:-1] + bin_edges[1:]) # bin centres
        plot_stats_minflux(deltas, durations_proper, tdiff, tdmedian, 
                           delta_averages, delta_av_times, showTimeAverages=True, 
                           timestamp=timeStamp, dirPath=dirpath)
    else:
        plot_stats_minflux(deltas, durations_proper, tdiff, tdmedian, data['efo'], None, 
                           timestamp=timeStamp, dirPath=dirpath)


# similar version but now using a pipeline
def analyse_locrate(data,datasource='Localizations',showTimeAverages=True, plot=True, timestamp=None):
    curds = data.selectedDataSourceKey
    data.selectDataSource(datasource)
    bins = np.arange(int(data['clumpIndex'].max())+1) + 0.5
    counts, bin_edges, binnumber = binned_statistic(data['clumpIndex'],data['tim'],statistic='count', bins=bins)
    starts, bin_edges, binnumber = binned_statistic(data['clumpIndex'],data['tim'],statistic='min', bins=bins)
    # for some reason we seem to get empty counts, i.e. the original clumpIndices are non-consecutive
    # NOTE: investigate IO of NPY MINFLUX data why this can happen!
    starts = starts[counts > 0]
    ends, bin_edges, binnumber = binned_statistic(data['clumpIndex'],data['tim'],statistic='max', bins=bins)
    ends = ends[counts > 0]
    
    durations = ends - starts
    deltas = starts[1:]-ends[:-1]
    # now we specifically look for the deltas within a trace
    tdiff = data['tim'][1:]-data['tim'][:-1]
    tracejump = data['clumpIndex'][1:]-data['clumpIndex'][:-1] # find all positions where the trace ID changes
    tdintrace = tdiff[tracejump < 0.1] # and now we exclude all tracejump deltas
    tdmedian = np.median(tdintrace)
    durations_proper = durations + tdmedian # we count one extra localisation, using the median duration
    # the extra is because we leave at least one localisation out from the total timing when we subtract ends-starts

    lenx_um = 1e-3*(data['x'].max()-data['x'].min())
    leny_um = 1e-3*(data['y'].max()-data['y'].min())
    area_string = 'area %.1fx%.1f um^2' % (lenx_um,leny_um)
    data.selectDataSource(curds)
    
    dirpath, filename = os.path.split(data.filename)
    timeStamp = data.mdh.get('MINFLUX.TimeStamp')
    
    if plot:
        if showTimeAverages:
            delta_averages, bin_edges, binnumber = binned_statistic(starts[:-1],deltas,statistic='mean', bins=50)
            delta_av_times = 0.5*(bin_edges[:-1] + bin_edges[1:]) # bin centres
            plot_stats_minflux(deltas, durations_proper, tdintrace, delta_averages, delta_av_times,
                               showTimeAverages=True, dsKey = datasource, areaString=area_string, timestamp=timeStamp, dirPath=dirpath)
        else:
            plot_stats_minflux(deltas, durations_proper, tdintrace, data['efo'], None, dsKey = datasource, areaString=area_string, timestamp=timeStamp, dirPath=dirpath)

    return (starts,ends,deltas,durations_proper,tdintrace)
