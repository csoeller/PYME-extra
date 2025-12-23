import matplotlib.pyplot as plt
import numpy as np
import wx

import logging
logger = logging.getLogger(__file__)

from PYMEcs.pyme_warnings import warn
import PYMEcs.misc.utils as mu

def plot_errors(pipeline,ds='coalesced_nz',dsclumps='with_clumps'):
    if not ds in pipeline.dataSources:
        warn('no data source named "%s" - check recipe and ensure this is MINFLUX data' % ds)
        return
    curds = pipeline.selectedDataSourceKey
    pipeline.selectDataSource(ds)
    p = pipeline
    clumpSize = p['clumpSize']
    plt.figure()
    plt.subplot(221)
    if 'error_z' in pipeline.keys():
        plt.boxplot([p['error_x'],p['error_y'],p['error_z']],labels=['error_x','error_y','error_z'])
    else:
        plt.boxplot([p['error_x'],p['error_y']],labels=['error_x','error_y'])
    plt.ylabel('loc error - coalesced (nm)')
    pipeline.selectDataSource(dsclumps)
    plt.subplot(222)
    if 'nPhotons' in p.keys() and 'fbg' in p.keys():
        bp_dict = plt.boxplot([p['nPhotons'],p['fbg']],labels=['photons','background rate'])
        for line in bp_dict['medians']:
            # get position data for median line
            x, y = line.get_xydata()[0] # top of median line
            # overlay median value
            plt.text(x, y, '%.0f' % y,
                     horizontalalignment='right') # draw above, centered

    uids, idx = np.unique(p['clumpIndex'],return_index=True)
    plt.subplot(223)
    if 'error_z' in pipeline.keys():
        plt.boxplot([p['error_x'][idx],p['error_y'][idx],p['error_z'][idx]],
                    labels=['error_x','error_y','error_z'])
    else:
        plt.boxplot([p['error_x'][idx],p['error_y'][idx]],labels=['error_x','error_y'])
    plt.ylabel('loc error - raw (nm)')
    plt.subplot(224)
    bp_dict = plt.boxplot([clumpSize],labels=['clump size'])
    for line in bp_dict['medians']:
        # get position data for median line
        x, y = line.get_xydata()[0] # top of median line
        # overlay median value
        plt.text(x, y, '%.0f' % y,
                 horizontalalignment='right') # draw above, centered
    plt.tight_layout()
    if mu.autosave_check():
        fpath = mu.get_ds_path(p)
        plt.savefig(mu.fname_from_timestamp(fpath,p.mdh,'_locError',ext='.png'),
                    dpi=300, bbox_inches='tight')
        pipeline.selectDataSource(ds)
        df = pd.DataFrame({
            "Metric": ["Photons", "Background", "Clump Size", "Error X", "Error Y"],
            "Median": [np.median(p[key]) for key in ['nPhotons','fbg','clumpSize','error_x','error_y']],
            "Unit": ["","","","nm","nm",]
        })
        if p.mdh['MINFLUX.Is3D']: # any code needs to check for 2D vs 3D
            df.loc[df.index.max() + 1] = ["Error Z",np.median(p['error_z']),'nm']
        mu.autosave_csv(df,fpath,p.mdh,'_locError')
                                        
    pipeline.selectDataSource(curds)
    
from PYMEcs.misc.matplotlib import boxswarmplot
import pandas as pd
def _plot_clustersize_counts(cts, ctsgt1, xlabel='Cluster Size', wintitle=None, bigCfraction=None,bigcf_percluster=None, plotints=True, **kwargs):
    if 'range' in kwargs:
        enforce_xlims = True
        xlims0=kwargs['range']
        extent = xlims0[1] - xlims0[0]
        frac = 0.05
        xlims = [xlims0[0] - frac*extent, xlims0[1] + frac*extent]
    else:
        enforce_xlims = False
    fig = plt.figure()
    if (plotints):
        plotn=300
    else:
        plotn=200
    plt.subplot(plotn+21)
    h = plt.hist(cts,**kwargs,log=True)
    plt.xlabel(xlabel)
    plt.plot([np.mean(cts),np.mean(cts)],[0,h[0].max()])
    plt.plot([np.median(cts),np.median(cts)],[0,h[0].max()],'--')
    if enforce_xlims:
        plt.xlim(*xlims)
    plt.subplot(plotn+22)
    h = plt.hist(ctsgt1,**kwargs,log=True)
    plt.xlabel('%s ( > 1)' % xlabel)
    plt.plot([np.mean(ctsgt1),np.mean(ctsgt1)],[0,h[0].max()])
    plt.plot([np.median(ctsgt1),np.median(ctsgt1)],[0,h[0].max()],'--')
    if enforce_xlims:
        plt.xlim(*xlims)
    plt.subplot(plotn+23)
    dfcs = pd.DataFrame.from_dict(dict(SUclusterSize=cts))
    boxswarmplot(dfcs,format="%.1f",swarmsize=5,width=0.2,annotate_means=True,annotate_medians=True,swarmalpha=0.15,strip=True)
    plt.subplot(plotn+24)
    dfcsgt1 = pd.DataFrame.from_dict(dict(SUclusterSizeGT1=ctsgt1))
    boxswarmplot(dfcsgt1,format="%.1f",swarmsize=5,width=0.2,annotate_means=True,annotate_medians=True,swarmalpha=0.15,strip=True)
    if (plotints):
        plt.subplot(plotn+25)
        dfcs = pd.DataFrame.from_dict(dict(RyRclusterSizeInts=np.ceil(0.25*cts)))
        boxswarmplot(dfcs,format="%.1f",swarmsize=5,width=0.2,annotate_means=True,annotate_medians=True,swarmalpha=0.15,strip=True)
        plt.subplot(plotn+26)
        dfcsgt1 = pd.DataFrame.from_dict(dict(RyRclusterSizeIntsGT1=np.ceil(0.25*ctsgt1)))
        boxswarmplot(dfcsgt1,format="%.1f",swarmsize=5,width=0.2,annotate_means=True,annotate_medians=True,swarmalpha=0.15,strip=True)

    largest = cts[np.argsort(cts)][-3:]
    fraction = largest.sum(dtype=float) / cts.sum()
    msg = "3 largest Cs (%s) make up %.1f %% of SUs" %(largest,100.0*fraction)
    fig.suptitle(msg)
    
    # bp_dict = plt.boxplot([cts,ctsgt1],labels=['cluster size','clusters > 1'], showmeans=True)
    # for line in bp_dict['means']:
    #     # get position data for median line
    #     x, y = line.get_xydata()[0] # top of median line
    #     # overlay median value
    #     plt.text(x-0.25, y, '%.1f' % y,
    #              horizontalalignment='center') # draw above, centered    
    plt.tight_layout()

    if wintitle is not None:
        figtitle = "%s" % wintitle
    else:
        figtitle = ''
    if bigCfraction is not None:
        figtitle = figtitle + " bigC fraction %.1f %%" % (100*bigCfraction)
        if bigcf_percluster is not None:
            figtitle = figtitle + " per cluster %.1f %%" % (100*bigcf_percluster)
    else:
        if bigcf_percluster is not None:
            figtitle = figtitle + " bigC frac per cluster %.1f %%" % (100*bigcf_percluster)

    fig.canvas.manager.set_window_title(figtitle)

def plot_cluster_analysis(pipeline, ds='dbscanClustered',showPlot=True, return_means=False,
                          return_data=False, psu=None, bins=15, bigc_thresh=50, **kwargs):
    if not ds in pipeline.dataSources:
        warn('no data source named "%s" - check recipe and ensure this is MINFLUX data' % ds)
        return
    curds = pipeline.selectedDataSourceKey
    pipeline.selectDataSource(ds)
    p = pipeline
    uids, idx, cts = np.unique(p['dbscanClumpID'], return_index=True, return_counts=True)
    nall = p['x'].size
    if 'bigCs' in p.dataSources:
        fraction = p.dataSources['bigCs']['x'].size/float(nall)
    else:
        fraction=None
    clustersizes = p['dbscanClumpSize'][idx]
    bigc_fracpercluster = float(np.sum(clustersizes>bigc_thresh))/clustersizes.size
    ctsgt1 = cts[cts > 1.1]
    pipeline.selectDataSource(curds)
    timestamp = pipeline.mdh.get('MINFLUX.TimeStamp')
    if showPlot:
        if psu is not None:
            _plot_clustersize_counts(cts, ctsgt1,bins=bins,xlabel='# subunits',wintitle=timestamp,bigCfraction=fraction,bigcf_percluster=bigc_fracpercluster,**kwargs)
        else:
            _plot_clustersize_counts(cts, ctsgt1,bins=bins,wintitle=timestamp,bigCfraction=fraction,bigcf_percluster=bigc_fracpercluster,**kwargs)
        if psu is not None:
            _plot_clustersize_counts(cts/4.0/psu, ctsgt1/4.0/psu, xlabel='# RyRs, corrected', bins=bins,wintitle=timestamp,
                                     bigCfraction=fraction,bigcf_percluster=bigc_fracpercluster,**kwargs)
    
    csm = cts.mean()
    csgt1m = ctsgt1.mean()
    csmd = np.median(cts)
    csgt1md = np.median(ctsgt1)
    
    print("Mean cluster size: %.2f" % csm)
    print("Mean cluster size > 1: %.2f" % csgt1m)
    print("Median cluster size: %.2f" % csmd)
    print("Median cluster size > 1: %.2f" % csgt1md)

    if return_means:
        return (csm,csgt1m)

    if return_data:
        return (cts,ctsgt1)

def cluster_analysis(pipeline):
    return plot_cluster_analysis(pipeline, ds='dbscanClustered',showPlot=False,return_means=True)
    
def plot_intra_clusters_dists(pipeline, ds='dbscanClustered',bins=15,NNs=1,**kwargs):
    if not ds in pipeline.dataSources:
        warn('no data source named "%s" - check recipe and ensure this is MINFLUX data' % ds)
        return
    from scipy.spatial import KDTree
    curds = pipeline.selectedDataSourceKey
    pipeline.selectDataSource(ds)
    p = pipeline
    uids, cts = np.unique(p['dbscanClumpID'], return_counts=True)
    checkids = uids[cts > 5.0]
    dists = []
    for cid in checkids:
        coords = np.vstack([p[k][p['dbscanClumpID'] == cid] for k in ['x','y','z']]).T
        tree = KDTree(coords)
        dd, ii = tree.query(coords,k=NNs+1)
        dists.extend(list(dd[:,1:].flatten()))
    pipeline.selectDataSource(curds)
    plt.figure()
    h=plt.hist(dists,bins=bins,**kwargs)

def cornercounts(pipeline,backgroundFraction=0.0):
    curds = pipeline.selectedDataSourceKey
    allds = 'withNNdist'
    ds2 = 'group2'
    ds3 = 'group3'
    ds4 = 'group4'
    p = pipeline
    pipeline.selectDataSource(allds)
    n_all = p['x'].size
    pipeline.selectDataSource(ds2)
    n_2 = p['x'].size
    pipeline.selectDataSource(ds3)
    n_3 = p['x'].size
    pipeline.selectDataSource(ds4)
    n_4 = p['x'].size
    pipeline.selectDataSource(curds)
    # note: we count RyRs, not corners
    # double check these ideas
    ccs = np.array([(1.0-backgroundFraction)*(n_all-(n_2+n_3+n_4)), n_2/2.0, n_3/3.0, n_4/4.0])
    return ccs

def print_basiccornerstats(pipeline, ds='filtered_localizations'):
    curds = pipeline.selectedDataSourceKey
    pipeline.selectDataSource(ds)
    data = pipeline
    dx_um = (data['x'].max() - data['x'].min())/1e3
    dy_um = (data['y'].max() - data['y'].min())/1e3
    area_um2 = dx_um * dy_um
    print("x extent: %.2f um, y extent: %.2f um, area: %.1f um^2" % (dx_um,dy_um,area_um2))
    pipeline.selectDataSource('coalesced_nz')
    from scipy.stats import iqr
    z_iqr_nm = iqr(data['z'])
    z_fwhm_nm = z_iqr_nm * 2.35/(2*0.675)
    z_full_nm = 4.0 * z_iqr_nm
    print("z extent (IQR): %.1f nm, (FWHM): %.1f nm, (Full extent): %.1f nm" % (z_iqr_nm, z_fwhm_nm, z_full_nm))
    n1 = data['x'].size
    pipeline.selectDataSource('closemerged')
    n2 = data['x'].size
    print("number corners: %d, (%d closemerged)" % (n1,n2))
    print("cornerdensity %.1f corners/um^2" % (n1/area_um2))
    print("cornerdensity %.1f corners/um^2 (closemerged)" % (n2/area_um2))
    z_fwhm_ratio = z_fwhm_nm / 100.0
    print("corner volume density: %.1f corners/um^2/100nm" % (n1/area_um2/z_fwhm_ratio))
    print("corner volume density: %.1f corners/um^2/100nm (closemerged)" % (n2/area_um2/z_fwhm_ratio))
    pipeline.selectDataSource(curds)

def plot_zextent(pipeline, ds='closemerged', series_name='This series'):
    
    def set_axis_style(ax, labels):
        ax.xaxis.set_tick_params(direction='out')
        ax.xaxis.set_ticks_position('bottom')
        ax.set_xticks(np.arange(1, len(labels) + 1), labels=labels)
        ax.set_xlim(0.25, len(labels) + 0.75)
        ax.set_xlabel('Sample name')

    fig, ax2 = plt.subplots(1, 1, figsize=(9, 4))

    curds = pipeline.selectedDataSourceKey
    pipeline.selectDataSource(ds)
    zdata = pipeline['z']
    pipeline.selectDataSource(curds)
    q005,quartile1, medians, quartile3 = np.percentile(zdata, [0.5,25, 50, 75])
    zdatac = zdata - q005

    ax2.set_title('z - axis value distribution')
    parts = ax2.violinplot(
        zdatac, showmeans=True, showmedians=False,
        showextrema=False, quantiles = [0.25,0.75])

    for pc in parts['bodies']:
        pc.set_facecolor('#D43F3A')
        pc.set_edgecolor('black')
        pc.set_alpha(1)

    # set style for the axes
    labels = [series_name]
    set_axis_style(ax2, labels)

    plt.subplots_adjust(bottom=0.15, wspace=0.05)

    
from scipy.special import binom
from scipy.optimize import curve_fit

def sigpn(p):
    return pn(1,p)+pn(2,p)+pn(3,p)+pn(4,p)

def sigptot(p):
    return pn(0,p) + sigpn(p)

def pn(k,p):
    return (binom(4,k)*(np.power(p,k)*np.power((1-p),(4-k))))

def pnn(k,p):
    return (pn(k,p)/(1-pn(0,p)))

def fourcornerplot(pipeline,sigma=None,backgroundFraction=0.0,showplot=True,quiet=False):
    ccs = cornercounts(pipeline,backgroundFraction=backgroundFraction)
    ccsn = ccs/ccs.sum()
    ks = np.arange(4)+1
    popt, pcov = curve_fit(pnn, ks, ccsn,sigma=sigma)
    perr = np.sqrt(np.diag(pcov))
    p_missed = pn(0,popt[0])
    p_m_min = pn(0,popt[0]+perr[0])
    p_m_max = pn(0,popt[0]-perr[0])
    if showplot:
        plt.figure()
        ax = plt.subplot(111)
        ax.bar(ks-0.4, ccsn, width=0.4, color='b', align='center')
        ax.bar(ks, pnn(ks,popt[0]), width=0.4, color='g', align='center')
        ax.legend(['Experimental data', 'Fit'])
        plt.title("Best fit p_lab=%.3f +- %.3f" % (popt[0],perr[0]))
    if not quiet:
        print('optimal p: %.3f +- %.3f' % (popt[0],perr[0]))
        print('missed fraction: %.2f (%.2f...%.2f)' % (p_missed,p_m_min,p_m_max))
    return (popt[0],perr[0],pnn(ks,popt[0]),ccsn)


sigmaDefault = [0.15,0.05,0.03,0.03]
backgroundDefault = 0.15

def fourcornerplot_default(pipeline,sigma=sigmaDefault,backgroundFraction=backgroundDefault,showplot=True,quiet=False):
    return fourcornerplot(pipeline,sigma=sigma,backgroundFraction=backgroundFraction,showplot=showplot,quiet=False)

def subunitfit(pipeline):
    return fourcornerplot_default(pipeline,showplot=False)

def plot_tracking(pipeline,is_coalesced=False,lowess_fraction=0.05):
    p = pipeline
    has_z = 'z_nc' in pipeline.keys()
    if has_z:
        nrows = 3
    else:
        nrows = 2

    t_s = 1e-3*p['t']
    xmbm = -(p['x']-p['x_nc']) # we flip the sign from now on
    ymbm = -(p['y']-p['y_nc'])
    if has_z:
        zmbm = -(p['z']-p['z_nc'])

    if is_coalesced:
        from statsmodels.nonparametric.smoothers_lowess import lowess
        xmbms = lowess(xmbm, t_s, frac=lowess_fraction, return_sorted=False)
        ymbms = lowess(ymbm, t_s, frac=lowess_fraction, return_sorted=False)
        if has_z:
            zmbms = lowess(zmbm, t_s, frac=lowess_fraction, return_sorted=False)

    plt.figure(num='beamline monitoring corrections')
    plt.subplot(nrows,1,1)
    plt.plot(t_s,xmbm)
    if is_coalesced:
        plt.plot(t_s,xmbms)
    plt.xlabel('Time (s)')
    plt.ylabel('x-difference (nm)')
    plt.subplot(nrows,1,2)
    plt.plot(t_s,ymbm)
    if is_coalesced:
        plt.plot(t_s,ymbms)
    plt.xlabel('Time (s)')
    plt.ylabel('y-difference (nm)')
    if 'z_nc' in pipeline.keys():
        plt.subplot(nrows,1,3)
        plt.plot(t_s,zmbm)
        if is_coalesced:
            plt.plot(t_s,zmbms)
        plt.xlabel('Time (s)')
        plt.ylabel('z-difference (nm)')
    plt.tight_layout()

    if not is_coalesced:
        return # skip HF plot
    
    plt.figure(num='MBM corrections HF component')
    plt.subplot(nrows,1,1)
    plt.plot(t_s,xmbm-xmbms)
    plt.xlabel('Time (s)')
    plt.ylabel('x-difference (nm)')
    plt.grid(axis='y')
    plt.subplot(nrows,1,2)
    plt.plot(t_s,ymbm-ymbms)
    plt.grid(axis='y')
    plt.xlabel('Time (s)')
    plt.ylabel('y-difference (nm)')
    if 'z_nc' in pipeline.keys():
        plt.subplot(nrows,1,3)
        plt.plot(t_s,zmbm-zmbms)
        plt.grid(axis='y')
        plt.xlabel('Time (s)')
        plt.ylabel('z-difference (nm)')
    plt.tight_layout()


def plot_site_tracking(pipeline,fignum=None,plotSmoothingCurve=True):
    p=pipeline
    t_s = 1e-3*p['t']
    if fignum is not None:
        fig, axs = plt.subplots(2, 2,num='origami site tracks %d' % fignum)
    else:
        fig, axs = plt.subplots(2, 2)

    axs[0, 0].scatter(t_s,p['x_site_nc'],s=0.3,c='black',alpha=0.5)
    if plotSmoothingCurve:
        axs[0, 0].plot(t_s,p['x_ori']-p['x'],'r',alpha=0.4)
    axs[0, 0].set_ylim(-15,15)
    axs[0, 0].set_xlabel('t (s)')
    axs[0, 0].set_ylabel('x (nm)')
        
    axs[0, 1].scatter(t_s,p['y_site_nc'],s=0.3,c='black',alpha=0.5)
    if plotSmoothingCurve:
        axs[0, 1].plot(t_s,p['y_ori']-p['y'],'r',alpha=0.4)
    axs[0, 1].set_ylim(-15,15)
    axs[0, 1].set_xlabel('t (s)')
    axs[0, 1].set_ylabel('y (nm)')

    if p.mdh.get('MINFLUX.Is3D',False):
        axs[1, 0].scatter(t_s,p['z_site_nc'],s=0.3,c='black',alpha=0.5)
        if plotSmoothingCurve:
            axs[1, 0].plot(t_s,p['z_ori']-p['z'],'r',alpha=0.4)
        axs[1, 0].set_ylim(-15,15)
        axs[1, 0].set_xlabel('t (s)')
        axs[1, 0].set_ylabel('z (nm)')

    ax = axs[1,1]
    if plotSmoothingCurve and 'x_nc' in p.keys():
        # plot the MBM track
        ax.plot(t_s,-(p['x_ori']-p['x_nc']),alpha=0.5,label='x')
        plt.plot(t_s,-(p['y_ori']-p['y_nc']),alpha=0.5,label='y')
        if p.mdh.get('MINFLUX.Is3D',False) and 'z_nc' in p.keys():
            ax.plot(t_s,-(p['z_ori']-p['z_nc']),alpha=0.5,label='z')
        ax.set_xlabel('t (s)')
        ax.set_ylabel('MBM corr (nm)')
        ax.legend()
    else:
        axs[1, 1].plot(t_s,p['x_ori']-p['x'])
        axs[1, 1].plot(t_s,p['y_ori']-p['y'])
        if p.mdh.get('MINFLUX.Is3D',False):
            axs[1, 1].plot(t_s,p['z_ori']-p['z'])
        axs[1, 1].set_xlabel('t [s]')
        axs[1, 1].set_ylabel('orig. corr (nm)')
    plt.tight_layout()

from PYMEcs.Analysis.MINFLUX import analyse_locrate
from PYMEcs.misc.guiMsgBoxes import Error
from PYMEcs.misc.utils import unique_name
from PYMEcs.IO.MINFLUX import findmbm

from PYME.recipes.traits import HasTraits, Float, Enum, CStr, Bool, Int, List
import PYME.config

class MINFLUXSettings(HasTraits):
    withOrigamiSmoothingCurves = Bool(True,label='Plot smoothing curves',desc="if overplotting smoothing curves " +
                                      "in origami site correction analysis")
    defaultDatasourceForAnalysis = CStr('Localizations',label='default datasource for analysis',
                                        desc="the datasource key that will be used by default in the MINFLUX " +
                                        "properties functions (EFO, localisation rate, etc)") # default datasource for acquisition analysis
    defaultDatasourceCoalesced = CStr('coalesced_nz',label='default datasource for coalesced analysis',
                                        desc="the datasource key that will be used by default when a " +
                                        "coalesced data source is required")
    defaultDatasourceWithClumps = CStr('with_clumps',label='default datasource for clump analysis',
                                        desc="the datasource key that will be used by default when a " +
                                        "data source with clump info is required")
    defaultDatasourceForMBM = CStr('coalesced_nz',label='default datasource for MBM analysis and plotting',
                                        desc="the datasource key that will be used by default in the MINFLUX " +
                                        "MBM analysis") # default datasource for MBM analysis
    datasourceForClusterAnalysis = CStr(PYME.config.get('MINFLUX-clusterDS','dbscan_clustered'),label='datasource for 3D cluster analysis',
                                        desc="the datasource key that will be used to generate the 3D cluster size analysis")
    
    largeClusterThreshold = Float(50,label='Threshold for large clusters',
                                  desc='minimum number of events to classify as large cluster')
    clustercountsPlotWithInts = Bool(False,label='Include "integer" plots in cluster stats',
                                  desc='plot integer quantized cluster stats that avoid counting fractional RyR numbers')
    origamiWith_nc = Bool(False,label='add 2nd moduleset (no MBM corr)',
                          desc="if a full second module set is inserted to also analyse the origami data without any MBM corrections")
    origamiErrorLimit = Float(10.0,label='xLimit when plotting origami errors',
                              desc="sets the upper limit in x (in nm) when plotting origami site errors")


class MINFLUXSiteSettings(HasTraits):
    showPoints = Bool(True)
    plotMode = Enum(['box','violin'])
    pointsMode = Enum(['swarm','strip'])
    siteMaxNum = Int(100,label='Max number of sites for box plot',
                            desc="the maximum number of sites for which site stats boxswarmplot is generated, violinplot otherwise")
    precisionRange_nm = Float(10)

class DateString(HasTraits):
    TimeStampString = CStr('',label="Time stamp",desc='the time stamp string in format yymmdd-HHMMSS')


class MBMaxisSelection(HasTraits):
    SelectAxis = Enum(['x-y-z','x','y','z','std_x','std_y','std_z'])

class MINFLUXplottingDefaults(HasTraits):
    FontSize = Float(12)
    LineWidth = Float(1.5)

class MINFLUXanalyser():
    def __init__(self, visFr):
        self.visFr = visFr
        self.minfluxRIDs = {}
        self.origamiErrorFignum = 0
        self.origamiTrackFignum = 0
        self.analysisSettings = MINFLUXSettings()
        self.dstring = DateString()
        self.mbmAxisSelection = MBMaxisSelection()
        self.plottingDefaults = MINFLUXplottingDefaults()
        self.siteSettings = MINFLUXSiteSettings()
        
        visFr.AddMenuItem('MINFLUX', "Localisation Error analysis", self.OnErrorAnalysis)
        visFr.AddMenuItem('MINFLUX', "Cluster sizes - 3D", self.OnCluster3D)
        visFr.AddMenuItem('MINFLUX', "Cluster sizes - 2D", self.OnCluster2D)
        visFr.AddMenuItem('MINFLUX', "Analyse Localization Rate", self.OnLocalisationRate)
        visFr.AddMenuItem('MINFLUX', "EFO histogram (photon rates)", self.OnEfoAnalysis)
        visFr.AddMenuItem('MINFLUX', "plot tracking correction (if available)", self.OnTrackPlot)
        visFr.AddMenuItem('MINFLUX', "Analysis settings", self.OnMINFLUXSettings)
        visFr.AddMenuItem('MINFLUX', "Toggle MINFLUX analysis autosaving", self.OnToggleMINFLUXautosave)
        visFr.AddMenuItem('MINFLUX', "Manually create Colour panel", self.OnMINFLUXColour)

        visFr.AddMenuItem('MINFLUX>Origami', "group and analyse origami sites", self.OnOrigamiSiteRecipe)
        visFr.AddMenuItem('MINFLUX>Origami', "plot origami site correction", self.OnOrigamiSiteTrackPlot)
        visFr.AddMenuItem('MINFLUX>Origami', "plot origami error estimates", self.OnOrigamiErrorPlot)
        visFr.AddMenuItem('MINFLUX>Origami', "plot origami site stats", self.OnOrigamiSiteStats)
        visFr.AddMenuItem('MINFLUX>Origami', "add final filter for site-based corrected data", self.OnOrigamiFinalFilter)
        visFr.AddMenuItem('MINFLUX>Origami', "site settings", self.OnMINFLUXSiteSettings)

        visFr.AddMenuItem('MINFLUX>Util', "Plot temperature record matching current data series",self.OnMINFLUXplotTemperatureData)
        visFr.AddMenuItem('MINFLUX>Util', "Set MINFLUX temperature folder location", self.OnMINFLUXsetTempDataFolder)
        visFr.AddMenuItem('MINFLUX>Util', "Check if clumpIndex contiguous", self.OnClumpIndexContig)
        visFr.AddMenuItem('MINFLUX>Util', "Plot event scatter as function of position in clump", self.OnClumpScatterPosPlot)
        visFr.AddMenuItem('MINFLUX>Util', "Set plotting defaults (inc font size)", self.OnSetMINFLUXPlottingdefaults)
        visFr.AddMenuItem('MINFLUX>Util', "Estimate region size (with output filter)", self.OnEstimateMINFLUXRegionSize)
        
        visFr.AddMenuItem('MINFLUX>MBM', "Plot mean MBM info (and if present origami info)", self.OnMBMplot)
        visFr.AddMenuItem('MINFLUX>MBM', "Show MBM tracks", self.OnMBMtracks)
        visFr.AddMenuItem('MINFLUX>MBM', "Add MBM track labels to view", self.OnMBMaddTrackLabels)
        visFr.AddMenuItem('MINFLUX>MBM', "Save MBM bead trajectories to npz file", self.OnMBMSave)
        visFr.AddMenuItem('MINFLUX>MBM', "Save MBM bead settings to json file", self.OnMBMSettingsSave)
        visFr.AddMenuItem('MINFLUX>MBM', "Save MBM lowess cache", self.OnMBMLowessCacheSave)
        
        visFr.AddMenuItem('MINFLUX>RyRs', "Plot corner info", self.OnCornerplot)
        visFr.AddMenuItem('MINFLUX>RyRs', "Plot density stats", self.OnDensityStats)
        visFr.AddMenuItem('MINFLUX>RyRs', "Show cluster alpha shapes", self.OnAlphaShapes)

        visFr.AddMenuItem('MINFLUX>Zarr', "Show MBM attributes", self.OnMBMAttributes)
        visFr.AddMenuItem('MINFLUX>Zarr', "Show MFX attributes", self.OnMFXAttributes)
        visFr.AddMenuItem('MINFLUX>Zarr', "Show MFX metadata info (now in PYME metadata)", self.OnMFXInfo)
        visFr.AddMenuItem('MINFLUX>Zarr', "Convert zarr file store to zarr zip store", self.OnZarrToZipStore)
        visFr.AddMenuItem('MINFLUX>Zarr', "Run Paraflux Analysis", self.OnRunParafluxAnalysis)

        visFr.AddMenuItem('MINFLUX>Tracking', "Add traces as tracks (from clumpIndex)", self.OnAddMINFLUXTracksCI)
        visFr.AddMenuItem('MINFLUX>Tracking', "Add traces as tracks (from tid)", self.OnAddMINFLUXTracksTid)
        visFr.AddMenuItem('MINFLUX>Colour', "Plot colour stats", self.OnPlotColourStats)
        
        # this section establishes Menu entries for loading MINFLUX recipes in one click
        # these recipes should be MINFLUX processing recipes of general interest
        # and are populated from the customrecipes folder in the PYME config directories
        # code adapted from PYME.DSView.modules.recipes 
        import PYME.config
        customRecipes = PYME.config.get_custom_recipes()
        minfluxRecipes = dict((k, v) for k, v in customRecipes.items() if k.startswith('MINFLUX'))
        if len(minfluxRecipes) > 0:
            for r in minfluxRecipes:
                ID = visFr.AddMenuItem('MINFLUX>Recipes', r, self.OnLoadCustom).GetId()
                self.minfluxRIDs[ID] = minfluxRecipes[r]

    def OnEstimateMINFLUXRegionSize(self, event):
        p = self.visFr.pipeline # Get the pipeline from the GUI
        xsize = p['x'].max() - p['x'].min()
        ysize = p['y'].max() - p['y'].min()
        warn("region size is %d x % d nm (%.1f x %.1f um" % (xsize,ysize,1e-3*xsize,1e-3*ysize))
                
    def OnSetMINFLUXPlottingdefaults(self, event):
        if not self.plottingDefaults.configure_traits(kind='modal'):
            return
        from PYMEcs.misc.matplotlib import figuredefaults
        figuredefaults(fontsize=self.plottingDefaults.FontSize,linewidth=self.plottingDefaults.LineWidth)
        
    # --- Alex B provided function (to save - not yet) and plot ITR stats (Paraflux like) ---
    def OnRunParafluxAnalysis(self, event):
        from pathlib import Path

        # ======================================================================================
        # --- Select, load Zarr.zip file, convert into DataFrame, Run the analysis functions ---
        # ======================================================================================
        pipeline = self.visFr.pipeline # Get the pipeline from the GUI

        if pipeline is None:
            Error(self.visFr, "No data found. Please load a MINFLUX dataset first.")
            return
        if not pipeline.mdh['MINFLUX.Is3D']: # TODO: make paraflux analysis code 2D aware
            warn('paraflux analysis currently only implemented for 3D data, this is apparently 2D data; giving up...')
            return
        try:
            # if this is a zarr archive we should have a zarr attribute in the FitResults datasource 
            zarr_archive = pipeline.dataSources['FitResults'].zarr
        except:
            warn("data is not from a zarr archive, giving up...")
            return
        try:
            zarr_path = zarr_archive.store.path
        except AttributeError:
            warn("cannot get zarr store path from zarr object, not saving analysis data...")
            zarr_path = None
        
        # possible storage code, not yet used/implemented
        # datasources = pipeline._get_session_datasources()
        # store_path = datasources.get('FitResults')
        # store_path = Path(store_path)

        # paraflux analysis with progress dialog follows
        import PYMEcs.Analysis.Paraflux as pf
        mfx_zarrsource = pipeline.dataSources['FitResults'] # this should be a MinfluxZarrSource instance

        # check if we have a cached result
        if mfx_zarrsource._paraflux_analysis is None:
            # for example use of ProgressDialog see also
            # https://github.com/Metallicow/wxPython-Sample-Apps-and-Demos/blob/master/101_Common_Dialogs/ProgressDialog/ProgressDialog_extended.py
            progress = wx.ProgressDialog("Paraflux analysis in progress", "please wait", maximum=4,
                                         parent=self.visFr,
                                         style=wx.PD_SMOOTH
                                         | wx.PD_AUTO_HIDE)
            def upd(n):
                progress.Update(n)
                wx.Yield()

            # read all data from the zarr archive
            mfxdata = zarr_archive['mfx'][:]; upd(1)
            # processing 1st step, move data into pandas dataframe
            df_mfx, failure_map = pf.paraflux_mk_df_fm(mfxdata); upd(2)
            # Run the analysis steps
            vld_itr = pf.build_valid_df(df_mfx); upd(3)
            vld_itr = pf.compute_percentages(vld_itr)
            vld_itr = pf.analyze_failures(vld_itr, df_mfx, failure_map)
            initial_count = vld_itr['vld loc count'].iloc[0]
            vld_itr = pf.add_failure_metrics(vld_itr, initial_count); upd(4)
            mfx_zarrsource._paraflux_analysis = vld_itr
        else:
            vld_itr = mfx_zarrsource._paraflux_analysis
        
        vld_paraflux = pf.paraflux_itr_plot(vld_itr[['itr', 'passed itr %', 'CFR failure %', 'No signal % per itr pairs']])

        # here possible storage command, only if autosaving is enabled in config
        if mu.autosave_check() and zarr_path is not None:
            mu.autosave_csv(vld_itr.drop(columns='tid', errors='ignore'),
                            zarr_path,pipeline.mdh,'_iteration_stats_full')
        ### --- End of Alex B added functionality ---

    def OnClumpScatterPosPlot(self,event):
        from scipy.stats import binned_statistic
        from PYMEcs.IO.MINFLUX import get_stddev_property
        def detect_coalesced(pipeline):
            # placeholder, to be implemented
            return False

        bigClumpThreshold = 10
        pipeline = self.visFr.pipeline
        if "posInClump" not in pipeline.keys():
            warn("No posInClump info, is this MINFLUX data imported with recent PAME-extra; aborting...")
        if detect_coalesced(pipeline):
            warn("This is coalesced data, need a datasource with non-coalesced clump info, aborting...")

        bigclumps = pipeline['clumpSize'] >= bigClumpThreshold
        ids = pipeline['clumpIndex'][bigclumps]
        posIC = pipeline['posInClump'][bigclumps]
        # uids,revids = np.unique(ids,return_inverse=True)
        xIC = pipeline['x'][bigclumps]
        yIC = pipeline['y'][bigclumps]
        xctrd = get_stddev_property(ids,xIC,statistic='mean')
        yctrd = get_stddev_property(ids,yIC,statistic='mean')
        dists = np.sqrt((xIC - xctrd)**2 + (yIC - yctrd)**2)
        maxpIC = min(int(posIC.max()),500) # we do not use clumps longer than 500
        edges = -0.5+np.arange(maxpIC+2)
        idrange = (0,maxpIC)
        
        picDists, bin_edge, binno = binned_statistic(posIC, dists, statistic='mean',
                                                     bins=edges, range=idrange)
        binctrs = 0.5*(bin_edge[:-1]+bin_edge[1:])
        plt.figure()
        plt.plot(binctrs, picDists)
        plt.xlabel("position in clump")
        plt.ylabel("lateral distance from ctrd (nm)")
        plt.xlim(-5,100)
        plt.ylim(0,None)

        if pipeline.mdh.get('MINFLUX.Is3D'):
            zIC = pipeline['z'][bigclumps]
            zctrd = get_stddev_property(ids,zIC,statistic='mean')
            zdists = np.abs(zIC - zctrd)
            zpicDists, bin_edge, binno = binned_statistic(posIC, zdists, statistic='mean',
                                                          bins=edges, range=idrange)
            binctrs = 0.5*(bin_edge[:-1]+bin_edge[1:])
            plt.figure()
            plt.plot(binctrs, zpicDists)
            plt.xlabel("position in clump")
            plt.ylabel("z distance from ctrd (nm)")
            plt.xlim(-5,100)
            plt.ylim(0,None)
            
    def OnPlotColourStats(self,event):
        pipeline = self.visFr.pipeline
        chans = pipeline.colourFilter.getColourChans()
        if len(chans) < 1:
            warn("No colour channel info, check if colour filtering module is active; aborting plotting")
            return
        cols = ['m','c','y']
        if len(chans) > 3:
            raise RuntimeError("too many channels, can only deal with max 3, got %d" % len(chans))
        bins = np.linspace(0.45,0.95,int(0.5/0.005))
        fig, axs = plt.subplots(nrows=2)
        axs[0].hist(pipeline['dcr_trace'],bins=bins,color='gray',alpha=0.5,label='other')
        for chan,col in zip(chans,cols[0:len(chans)]):
            axs[0].hist(pipeline.colourFilter.get_channel_column(chan,'dcr_trace'),bins=bins,color=col,alpha=0.5,label=chan)
        axs[0].set_xlabel('dcr_trace')
        axs[0].set_ylabel('#')
        axs[0].legend(loc="upper right")
        axs[1].hist(pipeline['dcr'],bins=bins,color='gray',alpha=0.5,label='other')
        for chan,col in zip(chans,cols[0:len(chans)]):
            axs[1].hist(pipeline.colourFilter.get_channel_column(chan,'dcr'),bins=bins,color=col,alpha=0.5,label=chan)
        axs[1].set_xlabel('dcr')
        axs[1].set_ylabel('#')
        axs[1].legend(loc="upper right")
        plt.tight_layout()

    def OnMBMLowessCacheSave(self,event):
        pipeline = self.visFr.pipeline
        mod = findmbm(pipeline,return_mod=True)
        if mod is None:
            return
        mod.lowess_cachesave()

    def OnMBMAttributes(self, event):
        from  wx.lib.dialogs import ScrolledMessageDialog
        fres = self.visFr.pipeline.dataSources['FitResults']
        if 'zarr' in dir(fres):
            try:
                mbm_attrs = fres.zarr['grd']['mbm'].points.attrs['points_by_gri']
            except AttributeError:
                warn("could not access MBM attributes - do we have MBM data in zarr?")
                return
            import pprint
            mbm_attr_str = pprint.pformat(mbm_attrs,indent=4,width=120)
            with ScrolledMessageDialog(self.visFr, mbm_attr_str, "MBM attributes", size=(900,400),
                                        style=wx.RESIZE_BORDER | wx.DEFAULT_DIALOG_STYLE ) as dlg:
                dlg.ShowModal()
        else:
            warn("could not find zarr attribute - is this a MFX zarr file?")
        
    def OnMFXAttributes(self, event):
        from  wx.lib.dialogs import ScrolledMessageDialog
        fres = self.visFr.pipeline.dataSources['FitResults']
        if 'zarr' in dir(fres):
            try:
                mfx_attrs = fres.zarr['mfx'].attrs.asdict()
            except AttributeError:
                warn("could not access MFX attributes - do we have MFX data in zarr?")
                return
            import pprint
            mfx_attr_str = pprint.pformat(mfx_attrs,indent=4,width=120)
            with ScrolledMessageDialog(self.visFr, mfx_attr_str, "MFX attributes", size=(900,400),
                                        style=wx.RESIZE_BORDER | wx.DEFAULT_DIALOG_STYLE ) as dlg:
                dlg.ShowModal()
        else:
            warn("could not find zarr attribute - is this a MFX zarr file?")

    def OnMFXInfo(self, event):
        import wx.html

        fres = self.visFr.pipeline.dataSources['FitResults']
        if 'zarr' not in dir(fres):
            warn("could not find zarr attribute - is this a MFX zarr file?")
            return
        try:
            mfx_attrs = fres.zarr['mfx'].attrs.asdict()
        except AttributeError:
            warn("could not access MFX attributes - do we have MFX data in zarr?")
            return
        if '_legacy' in mfx_attrs:
            warn("legacy data detected - no useful MFX metadata in legacy data")
            return

        # if we make it to this part of the code there is some useful metadata to be looked at
        from PYMEcs.IO.MINFLUX import get_metadata_from_mfx_attrs
        md_itr_info, md_globals = get_metadata_from_mfx_attrs(mfx_attrs)

        # Create an HTML dialog
        dlg = wx.Dialog(self.visFr, title="MINFLUX Metadata Information", 
                       size=(950, 600), style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
        
        # Create HTML content
        html_window = wx.html.HtmlWindow(dlg, style=wx.html.HW_SCROLLBAR_AUTO)
        
        # Format the DataFrame as an HTML table
        html_content = "<html><body>"
        html_content += "<p><b>Note:</b> This info should now be available via the PYME metadata, please inspect your metadata tab if using the GUI.</p>"

        html_content += "<h2>MINFLUX Iteration Parameters</h2>"
        html_content += md_itr_info.to_html(
            classes='table table-striped',
            float_format=lambda x: f"{x:.2f}" if isinstance(x, float) else x
        )
        
        # Format global parameters
        html_content += "<h2>MINFLUX Global Parameters</h2>"
        html_content += "<table border=1 class='table table-striped'>"
        for key, value in sorted(md_globals.items()):
            html_content += f"<tr><td><b>{key}</b></td><td>{value}</td></tr>"
        html_content += "</table>"
        html_content += "</body></html>"
        
        html_window.SetPage(html_content)
        
        # Add OK button
        btn_sizer = wx.StdDialogButtonSizer()
        btn = wx.Button(dlg, wx.ID_OK)
        btn.SetDefault()
        btn_sizer.AddButton(btn)
        btn_sizer.Realize()
        
        # Layout
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(html_window, 1, wx.EXPAND | wx.ALL, 5)
        sizer.Add(btn_sizer, 0, wx.ALIGN_CENTER | wx.ALL, 5)
        dlg.SetSizer(sizer)
        
        dlg.ShowModal()
        dlg.Destroy()
        
    def OnZarrToZipStore(self, event):
        with wx.DirDialog(self.visFr, 'Zarr to convert ...',
                          style=wx.DD_DEFAULT_STYLE | wx.DD_DIR_MUST_EXIST) as ddialog:
            if ddialog.ShowModal() != wx.ID_OK:
                return
            fpath = ddialog.GetPath()
        from PYMEcs.misc.utils import zarrtozipstore
        from pathlib import Path
        zarr_root = Path(fpath)
        dest_dir = zarr_root.parent
        archive_name = dest_dir / zarr_root.with_suffix('.zarr').name # we make archive_name here in the calling routine so that we can check for existence etc
        
        if archive_name.with_suffix('.zarr.zip').exists():
            with wx.FileDialog(self.visFr, 'Select archive name ...',
                               wildcard='ZIP (*.zip)|*.zip',
                               defaultFile=str(archive_name.with_suffix('.zarr.zip').name),
                               defaultDir=str(archive_name.with_suffix('.zarr.zip').parent),
                               style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as fdialog:
                if fdialog.ShowModal() != wx.ID_OK:
                    return
                archive_name = Path(fdialog.GetPath())
                while archive_name.suffix in {'.zarr', '.zip'}: # this loop progressively removes the .zarr.zip suffix
                    archive_name = archive_name.with_suffix('')
                archive_name = archive_name.with_suffix('.zarr')
                # warn("got back name %s, using archive name %s" % (Path(fdialog.GetPath()).name,archive_name.name))
                
        progress = wx.ProgressDialog("converting to zarr zip store", "please wait", maximum=2,
                                     parent=self.visFr,
                                     style=wx.PD_SMOOTH | wx.PD_AUTO_HIDE
                                     )
        progress.Update(1)
        created = Path(zarrtozipstore(zarr_root,archive_name))
        progress.Update(2)
        progress.Destroy()
        
        from PYMEcs.misc.guiMsgBoxes import YesNo
        do_open = YesNo(self.visFr,("created new zip store\n\n'%s'\n\nin directory\n\n'%s'"
                        + "\n\nOpen newly created zipstore data?")
                        % (created.name,created.parent),
                        caption="Open new zarr zipstore?")

        if do_open:
            self.visFr.OpenFile(str(created))
        
    def OnDensityStats(self, event):
        from PYMEcs.Analysis.MINFLUX import plot_density_stats_sns
        plot_density_stats_sns(self.visFr.pipeline)

    def OnAlphaShapes(self, event):
        if 'cluster_shapes' not in self.visFr.pipeline.dataSources.keys():
            warn("missing data source 'cluster_shapes', will not display alpha shapes")
            return
        
        # now we add a layer to render our alpha shape polygons
        from PYME.LMVis.layers.tracks import TrackRenderLayer # NOTE: we may rename the clumpIndex variable in this layer to polyIndex or similar
        layer = TrackRenderLayer(self.visFr.pipeline, dsname='cluster_shapes', method='tracks', clump_key='polyIndex', line_width=2.0, alpha=0.5)
        self.visFr.add_layer(layer)

    def OnAddMINFLUXTracksCI(self, event):        
        # now we add a track layer to render our traces
        from PYME.LMVis.layers.tracks import TrackRenderLayer # NOTE: we may rename the clumpIndex variable in this layer to polyIndex or similar
        layer = TrackRenderLayer(self.visFr.pipeline, dsname='output', method='tracks', clump_key='clumpIndex', line_width=2.0, alpha=0.5)
        self.visFr.add_layer(layer)

    def OnAddMINFLUXTracksTid(self, event):        
        # now we add a track layer to render our traces
        from PYME.LMVis.layers.tracks import TrackRenderLayer # NOTE: we may rename the clumpIndex variable in this layer to polyIndex or similar
        layer = TrackRenderLayer(self.visFr.pipeline, dsname='output', method='tracks', clump_key='tid', line_width=2.0, alpha=0.5)
        self.visFr.add_layer(layer)

    def OnLoadCustom(self, event):
        self.visFr._recipe_manager.LoadRecipe(self.minfluxRIDs[event.GetId()])


    def OnMBMSave(self,event):
        from pathlib import Path
        pipeline = self.visFr.pipeline
        mbm = findmbm(pipeline)
        if mbm is None:
            return
        defaultFile = None
        MINFLUXts = pipeline.mdh.get('MINFLUX.TimeStamp')
        if MINFLUXts is not None:
            defaultFile = "%s__MBM-beads.npz" % MINFLUXts
        fdialog = wx.FileDialog(self.visFr, 'Save MBM beads as ...',
                                wildcard='NPZ (*.npz)|*.npz',
                                defaultFile=defaultFile,
                                style=wx.FD_SAVE)
        if fdialog.ShowModal() != wx.ID_OK:
            return

        fpath = fdialog.GetPath()
        np.savez(fpath,**mbm._raw_beads)

    def OnMBMSettingsSave(self,event):
        import json
        pipeline = self.visFr.pipeline
        mbm = findmbm(pipeline)
        if mbm is None:
            return
        mod = findmbm(pipeline,return_mod=True)
        settings = {}
        beadisgood = {}
        for bead in mbm.beadisgood:
            beadisgood[bead] = mbm.beadisgood[bead] and bead in mod._mbm_allbeads
        settings['beads'] = beadisgood
        settings['Median_window'] = mod.Median_window
        settings['Lowess_fraction'] = mod.MBM_lowess_fraction
        settings['Filename'] = mbm.name
        defaultFile = None
        MINFLUXts = pipeline.mdh.get('MINFLUX.TimeStamp')
        if MINFLUXts is not None:
            defaultFile = "%s__MBM-beads.npz-settings.json" % MINFLUXts
        fdialog = wx.FileDialog(self.visFr, 'Save MBM Settings as ...',
                                wildcard='JSON (*.json)|*.json',
                                defaultFile=defaultFile,
                                style=wx.FD_SAVE)
        if fdialog.ShowModal() != wx.ID_OK:
            return

        fpath = fdialog.GetPath()

        with open(fpath, 'w') as f:
            json.dump(settings, f, indent=4)

    def OnMBMplot(self,event):
        def drift_total(p,caxis):
            has_drift = 'driftx' in p.keys()
            has_drift_ori = 'driftx_ori' in p.keys()
            caxis_nc = "%s_nc" % caxis
            caxis_ori = "%s_ori" % caxis
            if has_drift:
                if has_drift_ori:
                    # note we have a driftaxis_ori hiding in -p[caxis_ori], so no need to explicitly add
                    drift_2ndpass = p['drift%s' % caxis] - (p[caxis_ori]-p[caxis_nc])
                    drift_1stpass = - (p[caxis_ori]-p[caxis_nc]) # again we have an implicit driftaxis_ori hiding in -p[caxis_ori]
                else:
                    drift_1stpass = p['drift%s' % caxis] - (p[caxis_ori]-p[caxis_nc])
                    drift_2ndpass = None
            return drift_1stpass, drift_2ndpass

        def plot_drift(p,ax,drift_1stpass, drift_2ndpass):
            has_drift = 'driftx' in p.keys()
            has_drift_ori = 'driftx_ori' in p.keys()
            if has_drift:
                if has_drift_ori:
                    ax.plot(t_s,drift_2ndpass, label='site-based 2nd pass')
                    ax.plot(t_s,drift_1stpass,'--', label='site-based 1st pass')
                else:
                    ax.plot(t_s,drift_1stpass, label='site-based 1st pass')

        p = self.visFr.pipeline
        if p.mdh['MINFLUX.Is3D']:
            axes = ['x','y','z']
            naxes = 3
        else:
            axes = ['x','y']
            naxes = 2
        has_drift = 'driftx' in p.keys()
        has_drift_ori = 'driftx_ori' in p.keys()
        mbm = findmbm(p,warnings=False)
        has_mbm = mbm is not None
        has_mbm2 = 'mbmx' in p.keys()

        if not has_drift and not (has_mbm or has_mbm2):
            warn("pipeline has neither drift info nor MBM info, aborting...")
        
        t_s = 1e-3*p['t']
        mbm_mean = {} # for caching
        mbm_meansm = {} # for caching
        t_sm = {}
        
        ### Fig 1 ####
        fig, axs = plt.subplots(nrows=naxes)
        for caxis, ax in zip(axes,axs):
            if has_drift:
                drift_1stpass, drift_2ndpass = drift_total(p,caxis)
                plot_drift(p,ax,drift_1stpass, drift_2ndpass)
            if has_mbm:
                mod = findmbm(p,warnings=False,return_mod=True)
                MBM_lowess_fraction = mod.MBM_lowess_fraction
                mbm_mean[caxis] = mbm.mean(caxis)
                ax.plot(mbm.t,mbm_mean[caxis],':',label='MBM mean')
                t_sm[caxis],mbm_meansm[caxis] = mod.lowess_calc(caxis) # this is now a cached version of the lowess calc!
                ax.plot(t_sm[caxis],mbm_meansm[caxis],'-.',label='MBM lowess (lf=%.2f)' % MBM_lowess_fraction)
            if has_mbm2:
                ax.plot(t_s,p['mbm%s' % caxis], label='MBM from module')
            ax.set_xlabel('time (s)')
            ax.set_ylabel('drift in %s (nm)' % caxis)
            ax.legend(loc="upper right")
            ax.set_title("Total drift (Site based plus any previous corrections)")
        fig.tight_layout()
        ### Fig 2 ####
        if has_mbm: # also plot a second figure without the non-smoothed MBM track
            fig, axs = plt.subplots(nrows=naxes)
            for caxis, ax in zip(axes,axs):
                if has_drift:
                    drift_1stpass, drift_2ndpass = drift_total(p,caxis)
                    plot_drift(p,ax,drift_1stpass, drift_2ndpass)
                if has_mbm:
                    #ax.plot(mbm.t,mbm_mean[caxis],':',label='MBM mean')
                    ax.plot(t_sm[caxis],mbm_meansm[caxis],'r-.',label='MBM lowess (lf=%.2f)' % MBM_lowess_fraction)
                ax.set_xlabel('time (s)')
                ax.set_ylabel('drift in %s (nm)' % caxis)
                ax.legend(loc="upper left")
            fig.tight_layout()
        ### Fig 3 ####
        if has_mbm: # also plot a third figure with all MBM tracks
            fig, axs = plt.subplots(nrows=naxes)
            for caxis, ax in zip(axes,axs):
                if has_drift:
                    drift_1stpass, drift_2ndpass = drift_total(p,caxis)
                    plot_drift(p,ax,drift_1stpass, drift_2ndpass)
                if has_mbm:
                    #ax.plot(mbm.t,mbm_mean[caxis],':',label='MBM mean')
                    ax.plot(t_sm[caxis],mbm_meansm[caxis],'r-.',label='MBM lowess (lf=%.2f)' % MBM_lowess_fraction)
                    mbm.plot_tracks_matplotlib(caxis,ax=ax,goodalpha=0.4)
                ax.set_xlabel('time (s)')
                ax.set_ylabel('drift in %s (nm)' % caxis)
                ax.legend(loc="upper left")
            fig.tight_layout()
        ### Fig 4 ####
        if has_mbm and has_drift: # also plot a fourth figure with a difference track for all axes
            tnew = 1e-3*p['t']
            mbmcorr = {}
            for axis in axes:
                axis_interp_msm = np.interp(tnew,t_sm[caxis],mbm_meansm[axis])          
                mbmcorr[axis] = axis_interp_msm
            fig, axs = plt.subplots(nrows=naxes)
            for caxis, ax in zip(axes,axs):
                drift_1stpass, drift_2ndpass = drift_total(p,caxis)
                if has_drift:
                    if has_drift_ori:
                        ax.plot(t_s,drift_2ndpass-mbmcorr[caxis], label='diff to site-based 2nd pass')
                    else:
                        ax.plot(t_s,drift_1stpass-mbmcorr[caxis], label='diff to site-based 1st pass')
                ax.plot([t_s.min(),t_s.max()],[0,0],'r-.')
                ax.set_xlabel('time (s)')
                ax.set_ylabel('differential drift in %s (nm)' % caxis)
                ax.legend(loc="upper left")
            fig.tight_layout()


    def OnMBMtracks(self, event):
        pipeline = self.visFr.pipeline
        mbm = findmbm(pipeline)
        if mbm is None:
            return # note that findmbm has already warned in this case
        if not self.mbmAxisSelection.configure_traits(kind='modal'):
            return
        ori_win = mbm.median_window
        mbm.median_window = 21 # go for pretty agressive smoothing
        if self.mbmAxisSelection.SelectAxis == 'x-y-z':
            fig, axes = plt.subplots(nrows=3)
            for axis,plotax in zip(['x','y','z'],axes):
                mbm.plot_tracks_matplotlib(axis,ax=plotax)
            fig.tight_layout()
        else:
            mbm.plot_tracks_matplotlib(self.mbmAxisSelection.SelectAxis)
        mbm.median_window = ori_win

    def OnMBMaddTrackLabels(self, event):
        pipeline = self.visFr.pipeline
        try:
            from PYME.LMVis.layers.labels import LabelLayer
        except:
            hasLL = False
        else:
            hasLL = True

        # note: should also add the merge module?
        # note: should also add the layer for mbm_tracks?

        if 'mbm_pos' not in pipeline.dataSources.keys():
            #warn("no datasource 'mbm_pos' which is needed for label display")
            #return
            mod = findmbm(pipeline, warnings=True, return_mod=True)
            if mod is None:
                return
            from PYME.recipes.localisations import MergeClumps
            mc = MergeClumps(pipeline.recipe,
                             inputName=mod.outputTracksCorr,
                             outputName='mbm_pos',
                             labelKey='objectID',
                             # important, otherwise we get a spurious bead labele R0
                             discardTrivial=True)
            pipeline.recipe.add_module(mc)
            pipeline.recipe.execute()
        if not hasLL:
            warn("could not load new experimental feature label layer, aborting...")
            return
        
        if hasLL:
            ll = LabelLayer(pipeline, dsname='mbm_pos', format_string='R{beadID:.0f}', cmap='grey_overflow', font_size=13, textColour='good')
            self.visFr.add_layer(ll)
            ll.update()

    def OnMINFLUXsetTempDataFolder(self, event):
        import PYME.config as config
        config_var = 'MINFLUX-temperature_folder'
        curfolder = config.get(config_var)
        if curfolder is None:
            warn("currently the MINFLUX temperature file folder is not set. Set the path in the following file dialog")
        else:
            warn("MINFLUX temperature file folder currently set to '%s'.\n\nAlter in the following dialog if needed or cancel that dialog if want to leave as is" % curfolder)
        with wx.DirDialog(self.visFr, "Choose folder containing temperature CSV files") as dialog:
            if dialog.ShowModal() == wx.ID_CANCEL:
                return
            folder = dialog.GetPath()
        if config.get(config_var) == folder:
            warn("config option '%s' already set to %s, leaving as is" % (config_var,folder))
            return # already set to this value, return

        config.update_config({config_var: folder},
                             config='user', create_backup=True)

    def OnToggleMINFLUXautosave(self, event):
        import PYME.config as config
        config_var = 'MINFLUX-autosave'

        if config.get(config_var,False):
            newval = False
        else:
            newval = True

        config.update_config({config_var: newval},
                             config='user', create_backup=False)

        warn("MINFLUX analysis autosave was set to %s" %  config.get(config_var))

        
    def OnMINFLUXplotTemperatureData(self, event):
        import PYME.config as config
        import os
        from os.path import basename
        from glob import glob

        configvar = 'MINFLUX-temperature_folder'
        folder = config.get(configvar)
        if folder is None:
            warn("Need to set Temperature file location first by setting config variable %s" % configvar)
            return
        elif not os.path.isdir(folder):
            warn(("Config variable %s is not set to a folder;\n" % (configvar)) +
                 ("needs to be a **folder** location, currently set to %s" % (folder)))
            return

        from PYMEcs.misc.utils import read_temperature_csv, set_diff, timestamp_to_datetime

        if len(self.visFr.pipeline.dataSources) == 0:
            warn("no datasources, this is probably an empty pipeline, have you loaded any data?")
            return
        
        t0 = self.visFr.pipeline.mdh.get('MINFLUX.TimeStamp')
        if t0 is None:
            warn("no MINFLUX TimeStamp in metadata, giving up")
            return
        # Convert t0 from timestamp for comparing it with timedates from csv file
        t0_dt = timestamp_to_datetime(t0)
        
        # Identify the correct temperature CSV files in the folder
        # Loop over CSVs to find matching file
        timeformat = config.get('MINFLUX-temperature_time_format', ['%d.%m.%Y %H:%M:%S',
                                                                    '%d/%m/%Y %H:%M:%S'])
        candidate_files = sorted(glob(os.path.join(folder, '*.csv')))

        for f in candidate_files:
            try:
                # print(f"\nChecking file: {basename(f)}\n")  # Debugging, Show which file is being checked

                df = read_temperature_csv(f, timeformat=timeformat)
                
                if 'datetime' not in df.columns:
                    print(f"File {f} has no 'datetime' column after parsing, skipping.")
                    continue

                logger.debug(f"successfully read and parsed temperature file {f}")

                if df['datetime'].min() <= t0_dt <= df['datetime'].max():
                    selected_file = f
                    logger.debug(f"found relevant time period in file {f}")
                    break
            except Exception as e:
                logger.debug(f"Error reading {f}: {e}")
                continue
        else:
            warn("No temperature file found that includes the MINFLUX TimeStamp")
            return

        # Read temperature data from the correct CSV file
        mtemps = read_temperature_csv(selected_file, timeformat=timeformat)
        
        set_diff(mtemps,timestamp_to_datetime(t0))
        p = self.visFr.pipeline
        range = (1e-3*p['t'].min(),1e-3*p['t'].max())
        sertemps = mtemps[mtemps['tdiff_s'].between(range[0],range[1])]
        if sertemps.empty:
            warn("no records in requested time window, is series time before or after start/end of available temperature records?\n" +
                 ("current records cover %s to %s" % (mtemps['Time'].iloc[0],mtemps['Time'].iloc[-1])) +
                 ("\nseries starts at %s" % (t0)))
            return
        else:
            # for now we make 2 subplots so that we can provide both s units and actual time
            fig, axes = plt.subplots(nrows=2, ncols=1)
            sertemps.plot('datetime','Stand',style='.-',
                                title="temperature record for series starting at %s" % t0, ax=axes[0])
            sertemps.plot('tdiff_s','Stand',style='.-', ax=axes[1])
            plt.tight_layout()
            
            fig, axes = plt.subplots(nrows=2, ncols=1)
            sertemps.plot('datetime','Box',style='.-',
                          title="temperature record for series starting at %s" % t0, ax=axes[0])
            sertemps.plot('tdiff_s','Box',style='.-', ax=axes[1])
            plt.tight_layout()

    def OnErrorAnalysis(self, event):
        plot_errors(self.visFr.pipeline,ds=self.analysisSettings.defaultDatasourceCoalesced,dsclumps=self.analysisSettings.defaultDatasourceWithClumps)

    def OnCluster3D(self, event):
        plot_cluster_analysis(self.visFr.pipeline, ds=self.analysisSettings.datasourceForClusterAnalysis,
                              bigc_thresh=self.analysisSettings.largeClusterThreshold,
                              plotints=self.analysisSettings.clustercountsPlotWithInts)

    def OnCluster2D(self, event):
        plot_cluster_analysis(self.visFr.pipeline, ds='dbscan2D')

    def OnClumpIndexContig(self, event):
        pipeline = self.visFr.pipeline
        curds = pipeline.selectedDataSourceKey
        pipeline.selectDataSource(self.analysisSettings.defaultDatasourceForAnalysis)
        if not 'clumpIndex' in pipeline.keys():
            Error(self.visFr,'no property called "clumpIndex", cannot check')
            pipeline.selectDataSource(curds)
            return
        uids = np.unique(pipeline['clumpIndex'])
        maxgap = np.max(uids[1:]-uids[:-1])
        pipeline.selectDataSource(curds)

        if maxgap > 1:
            msg = "clumpIndex not contiguous, maximal gap is %d\nCI 0..9 %s" % (maxgap,uids[0:10])
        else:
            msg = "clumpIndex is contiguous\nCI 0..9 %s" % uids[0:10]

        warn(msg)

    def OnLocalisationRate(self, event):
        pipeline = self.visFr.pipeline
        curds = pipeline.selectedDataSourceKey
        pipeline.selectDataSource(self.analysisSettings.defaultDatasourceForAnalysis)
        if not 'cfr' in pipeline.keys():
            Error(self.visFr,'no property called "cfr", likely no MINFLUX data - aborting')
            pipeline.selectDataSource(curds)
            return
        if not 'tim' in pipeline.keys():
            Error(self.visFr,'no property called "tim", you need to convert to CSV with a more recent version of PYME-Extra - aborting')
            pipeline.selectDataSource(curds)
            return
        pipeline.selectDataSource(curds)

        analyse_locrate(pipeline,datasource=self.analysisSettings.defaultDatasourceForAnalysis,showTimeAverages=True)

    def OnEfoAnalysis(self, event):
        pipeline = self.visFr.pipeline
        curds = pipeline.selectedDataSourceKey
        pipeline.selectDataSource(self.analysisSettings.defaultDatasourceForAnalysis)
        if not 'efo' in pipeline.keys():
            Error(self.visFr,'no property called "efo", likely no MINFLUX data or wrong datasource (CHECK) - aborting')
            return
        plt.figure()
        h = plt.hist(1e-3*pipeline['efo'],bins='auto',range=(0,200))
        dskey = pipeline.selectedDataSourceKey
        plt.xlabel('efo (photon rate in kHz)')
        plt.title("EFO stats, using datasource '%s'" % dskey)

        pipeline.selectDataSource(curds)

    def OnTrackPlot(self, event):
        p = self.visFr.pipeline
        curds = p.selectedDataSourceKey
        if self.analysisSettings.defaultDatasourceForMBM in p.dataSources.keys():
            # should be coalesced datasource
            p.selectDataSource(self.analysisSettings.defaultDatasourceForMBM)
            is_coalesced = 'coalesced' in self.analysisSettings.defaultDatasourceForMBM.lower()
        else:
            # try instead something that should exist
            p.selectDataSource(self.analysisSettings.defaultDatasourceForAnalysis)
            is_coalesced = 'coalesced' in self.analysisSettings.defaultDatasourceForAnalysis.lower()
        plot_tracking(p,is_coalesced,lowess_fraction=0.03)
        p.selectDataSource(curds)

    def OnOrigamiFinalFilter(self, event=None):
        from PYME.recipes.tablefilters import FilterTable
        pipeline = self.visFr.pipeline
        recipe = pipeline.recipe
        curds = pipeline.selectedDataSourceKey

        finalFiltered = unique_name('filtered_final',pipeline.dataSources.keys())

        modules = [FilterTable(recipe,inputName=curds,outputName=finalFiltered,
                               filters={'error_x' : [0,3.5],
                                        'error_x' : [0,3.5],
                                        'error_x' : [0,3.5],
                                        'efo' : [10e3,1e5],
                                        })]
        recipe.add_modules_and_execute(modules)
        pipeline.selectDataSource(finalFiltered)
        
    def OnOrigamiSiteRecipe(self, event=None):
        from PYMEcs.recipes.localisations import OrigamiSiteTrack, DBSCANClustering2
        from PYME.recipes.localisations import MergeClumps
        from PYME.recipes.tablefilters import FilterTable, Mapping
        
        pipeline = self.visFr.pipeline
        recipe = pipeline.recipe

        filters={'error_x' : [0,3.5],
                 'error_y' : [0,3.5]}
        if 'error_z' in pipeline.keys():
            filters['error_z'] = [0,3.5]
        
        preFiltered = unique_name('prefiltered',pipeline.dataSources.keys())
        corrSiteClumps = unique_name('corrected_siteclumps',pipeline.dataSources.keys())
        corrAll = unique_name('corrected_allpoints',pipeline.dataSources.keys())
        siteClumps = unique_name('siteclumps',pipeline.dataSources.keys())
        dbscanClusteredSites = unique_name('dbscanClusteredSites',pipeline.dataSources.keys())
        sites = unique_name('sites',pipeline.dataSources.keys())
        sites_c = unique_name('sites_c',pipeline.dataSources.keys())
        
        curds = pipeline.selectedDataSourceKey
        modules = [FilterTable(recipe,inputName=curds,outputName=preFiltered,
                               filters=filters),
                   DBSCANClustering2(recipe,inputName=preFiltered,outputName=dbscanClusteredSites,
                                     searchRadius = 15.0,
                                     clumpColumnName = 'siteID',
                                     sizeColumnName='siteClumpSize'),
                   FilterTable(recipe,inputName=dbscanClusteredSites,outputName=siteClumps,
                               filters={'siteClumpSize' : [3,50]}), # need a minimum clumpsize and also maximal to avoid "fused" sites
                   MergeClumps(recipe,inputName=siteClumps,outputName=sites,
                               labelKey='siteID',discardTrivial=True),
                   OrigamiSiteTrack(recipe,inputClusters=siteClumps,inputSites=sites,outputName=corrSiteClumps,smoothingBinWidthsSeconds=400,
                                    outputAllPoints=corrAll,inputAllPoints=curds,labelKey='siteID',binnedStatistic='median'), # median to play it safe
                   MergeClumps(recipe,inputName=corrSiteClumps,outputName=sites_c,
                               labelKey='siteID',discardTrivial=True)]
        recipe.add_modules_and_execute(modules)

        if self.analysisSettings.origamiWith_nc:
            preFiltered = unique_name('prefiltered_nc',pipeline.dataSources.keys())
            corrSiteClumps = unique_name('corrected_siteclumps_nc',pipeline.dataSources.keys())
            siteClumps = unique_name('siteclumps_nc',pipeline.dataSources.keys())
            dbscanClusteredSites = unique_name('dbscanClusteredSites_nc',pipeline.dataSources.keys())
            sites = unique_name('sites_nc',pipeline.dataSources.keys())
            sites_c = unique_name('sites_c_nc',pipeline.dataSources.keys())
            dbsnc = unique_name('dbs_nc',pipeline.dataSources.keys())
        
            modules = [FilterTable(recipe,inputName=curds,outputName=preFiltered,
                                   filters=filters),
                       DBSCANClustering2(recipe,inputName=preFiltered,outputName=dbscanClusteredSites,
                                         searchRadius = 15.0,
                                         clumpColumnName = 'siteID',
                                         sizeColumnName='siteClumpSize'),
                       Mapping(recipe,inputName=dbscanClusteredSites,outputName=dbsnc,
                               mappings={'x': 'x_nc', 'y': 'y_nc', 'z': 'z_nc'}),
                       FilterTable(recipe,inputName=dbsnc,outputName=siteClumps,
                                   filters={'siteClumpSize' : [3,50]}), # need a minimum clumpsize and also maximal to avoid "fused" sites
                       MergeClumps(recipe,inputName=siteClumps,outputName=sites,
                                   labelKey='siteID',discardTrivial=True),
                       OrigamiSiteTrack(recipe,inputClusters=siteClumps,inputSites=sites,outputName=corrSiteClumps,
                                        labelKey='siteID',binnedStatistic='median'),
                       MergeClumps(recipe,inputName=corrSiteClumps,outputName=sites_c,
                                   labelKey='siteID',discardTrivial=True)]
        
            recipe.add_modules_and_execute(modules)
        
        pipeline.selectDataSource(corrSiteClumps)

    def OnOrigamiSiteTrackPlot(self, event):
        p = self.visFr.pipeline
        # need to add checks if the required properties are present in the datasource!!
        plot_site_tracking(p,fignum=self.origamiTrackFignum,
                           plotSmoothingCurve=self.analysisSettings.withOrigamiSmoothingCurves)
        self.origamiTrackFignum += 1

    def OnMINFLUXSettings(self, event):
        if self.analysisSettings.configure_traits(kind='modal'):
            pass

    def OnMINFLUXSiteSettings(self, event):
        if self.siteSettings.configure_traits(kind='modal'):
            pass

    def OnOrigamiErrorPlot(self, event):
        p = self.visFr.pipeline
        # need to check if the required properties are present in the datasource
        if 'error_x_ori' not in p.keys():
            warn("property 'error_x_ori' not present, possibly not the right datasource for origami site info. Aborting...")
            return
        
        def plot_errs(ax,axisname,errkeys):
            ax.hist(p[errkeys[0]],bins='auto',alpha=0.5,density=True,label='Trace est')
            ax.hist(p[errkeys[1]],bins='auto',alpha=0.5,density=True,label='Site est')
            ax.hist(p[errkeys[2]],bins='auto',alpha=0.5,density=True,label='Site est corr')
            ax.legend()
            ax.set_xlabel('error %s (nm)' % axisname)
            ax.set_ylabel('#')
  
        fig, axs = plt.subplots(2, 2,num='origami error estimates %d' % self.origamiErrorFignum)
        plot_errs(axs[0, 0], 'x', ['error_x_ori','error_x_nc','error_x'])
        axs[0, 0].set_xlim(0,self.analysisSettings.origamiErrorLimit)
        plot_errs(axs[0, 1], 'y', ['error_y_ori','error_y_nc','error_y'])
        axs[0, 1].set_xlim(0,self.analysisSettings.origamiErrorLimit)
        if p.mdh.get('MINFLUX.Is3D'):
            plot_errs(axs[1, 0], 'z', ['error_z_ori','error_z_nc','error_z'])
            axs[1, 0].set_xlim(0,self.analysisSettings.origamiErrorLimit)
        ax = axs[1,1]
        # plot the MBM track, this way we know if we are using the _nc data or the MBM corrected data for analysis
        t_s = 1e-3*p['t']
        ax.plot(t_s,p['x_ori']-p['x_nc'],alpha=0.5,label='x')
        plt.plot(t_s,p['y_ori']-p['y_nc'],alpha=0.5,label='y')
        if 'z_nc' in p.keys():
            ax.plot(t_s,p['z_ori']-p['z_nc'],alpha=0.5,label='z')
        ax.set_xlabel('t (s)')
        ax.set_ylabel('MBM corr [nm]')
        ax.legend()
        plt.tight_layout()
        
        uids = np.unique(p['siteID']) # currently siteID is hard coded - possibly make config option
        from PYMEcs.Analysis.MINFLUX import plotsitestats
        if uids.size < self.siteSettings.siteMaxNum:
            plotsitestats(p,fignum=('origami site stats %d' % self.origamiErrorFignum))
        self.origamiErrorFignum += 1

    def OnOrigamiSiteStats(self, event):
        from PYMEcs.IO.MINFLUX import get_stddev_property
        p = self.visFr.pipeline
        plotted = False
        # need to check if the required properties are present in the datasource
        if 'error_x_ori' not in p.keys():
            warn("property 'error_x_ori' not present, possibly not the right datasource for origami site info. Aborting...")
            return
        plotmode = self.siteSettings.plotMode
        from PYMEcs.Analysis.MINFLUX import plotsitestats
        uids,idx = np.unique(p['siteID'],return_index=True) # currently siteID is hard coded - possibly make config option
        if uids.size < self.siteSettings.siteMaxNum:
            swarmsize = 3
        else:
            swarmsize = 1.5

        plotsitestats(p,fignum=('origami site stats %d' % self.origamiErrorFignum),
                      swarmsize=swarmsize,mode=plotmode,showpoints=self.siteSettings.showPoints,
                      origamiErrorLimit=self.siteSettings.precisionRange_nm,
                      strip=(self.siteSettings.pointsMode == 'strip'))
        plotted = True
        #else:
        #    warn("Number of sites (%d) > max number for plotting (%d); check settings"
        #         % (uids.size,self.analysisSettings.origamiSiteMaxNum))

        counts = get_stddev_property(p['siteID'],p['siteID'],statistic='count')
        plt.figure(num=('site visits %d' % self.origamiErrorFignum))
        ax = plt.gca()
        ctmedian = np.median(counts[idx])
        ctmean = np.mean(counts[idx])

        h = plt.hist(counts[idx],bins='auto')
        ax.plot([ctmedian,ctmedian],[0,h[0].max()])
        plt.xlabel('Number of site visits')
        plt.text(0.85, 0.8, 'median %d' % ctmedian, horizontalalignment='right',
                 verticalalignment='bottom', transform=ax.transAxes)
        plt.text(0.85, 0.7, '  mean %.1f' % ctmean, horizontalalignment='right',
                 verticalalignment='bottom', transform=ax.transAxes)
        plotted = True
        if plotted:
            self.origamiErrorFignum += 1
        
    def OnMINFLUXColour(self,event):
        from PYME.LMVis import colourPanel
        
        mw = self.visFr
        if mw.colp is None: # no colourPanel yet
            self.visFr.pipeline.selectDataSource('colour_mapped')
            mw.adding_panes=True
            mw.colp = colourPanel.colourPanel(mw, mw.pipeline, mw)
            mw.AddPage(mw.colp, caption='Colour', select=False, update=False)
            mw.adding_panes=False
        else:
            warn('Colour panel appears to already exist - not creating new colour panel')

    def OnCornerplot(self,event):
        for ds in ['withNNdist','group2','group3','group4']:
            if ds not in self.visFr.pipeline.dataSources.keys():
                warn("need datasource %s which is not present, giving up..." % ds)
                return
        fourcornerplot_default(self.visFr.pipeline)

def Plug(visFr):
    # we are trying to monkeypatch pipeline and VisGUIFrame methods to sneak MINFLUX npy IO in;
    # in future we will ask for a way to get this considered by David B for a proper hook
    # in the IO code
    from PYMEcs.IO.MINFLUX import monkeypatch_npyorzarr_io
    monkeypatch_npyorzarr_io(visFr)
        
    return MINFLUXanalyser(visFr)
