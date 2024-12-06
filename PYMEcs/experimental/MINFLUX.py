import matplotlib.pyplot as plt
import numpy as np
import wx

from PYME.warnings import warn
def plot_errors(pipeline):
    if not 'coalesced_nz' in pipeline.dataSources:
        warn('no data source named "coalesced_nz" - check recipe and ensure this is MINFLUX data')
        return
    curds = pipeline.selectedDataSourceKey
    pipeline.selectDataSource('coalesced_nz')
    p = pipeline
    clumpSize = p['clumpSize']
    plt.figure()
    plt.subplot(221)
    if 'error_z' in pipeline.keys():
        plt.boxplot([p['error_x'],p['error_y'],p['error_z']],labels=['error_x','error_y','error_z'])
    else:
        plt.boxplot([p['error_x'],p['error_y']],labels=['error_x','error_y'])
    plt.ylabel('loc error - coalesced (nm)')
    pipeline.selectDataSource('with_clumps')
    plt.subplot(222)
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
    pipeline.selectDataSource(curds)
    
from PYMEcs.misc.matplotlib import boxswarmplot
import pandas as pd
def _plot_clustersize_counts(cts, ctsgt1, xlabel='Cluster Size', wintitle=None, bigCfraction=None, **kwargs):
    fig = plt.figure()
    plt.subplot(221)
    h = plt.hist(cts,**kwargs,log=True)
    plt.xlabel(xlabel)
    plt.plot([np.mean(cts),np.mean(cts)],[0,h[0].max()])
    plt.plot([np.median(cts),np.median(cts)],[0,h[0].max()],'--')
    plt.subplot(222)
    h = plt.hist(ctsgt1,**kwargs,log=True)
    plt.xlabel('%s ( > 1)' % xlabel)
    plt.plot([np.mean(ctsgt1),np.mean(ctsgt1)],[0,h[0].max()])
    plt.plot([np.median(ctsgt1),np.median(ctsgt1)],[0,h[0].max()],'--')
    plt.subplot(223)
    dfcs = pd.DataFrame.from_dict(dict(clusterSize=cts))
    boxswarmplot(dfcs,format="%.1f",swarmsize=5,width=0.2,annotate_means=True,annotate_medians=True,swarmalpha=0.15,strip=True)
    plt.subplot(224)
    dfcsgt1 = pd.DataFrame.from_dict(dict(clusterSizeGT1=ctsgt1))
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
        
    fig.canvas.manager.set_window_title(figtitle)

def plot_cluster_analysis(pipeline, ds='dbscanClustered',showPlot=True, return_means=False, psu=None, bins=15, **kwargs):
    if not ds in pipeline.dataSources:
        warn('no data source named "%s" - check recipe and ensure this is MINFLUX data' % ds)
        return
    curds = pipeline.selectedDataSourceKey
    pipeline.selectDataSource(ds)
    p = pipeline
    uids, cts = np.unique(p['dbscanClumpID'], return_counts=True)
    nall = p['x'].size
    if 'bigCs' in p.dataSources:
        fraction = p.dataSources['bigCs']['x'].size/float(nall)
    else:
        fraction=None
    ctsgt1 = cts[cts > 1.1]
    pipeline.selectDataSource(curds)
    timestamp = pipeline.mdh.get('MINFLUX.TimeStamp')
    if showPlot:
        if psu is not None:
            _plot_clustersize_counts(cts, ctsgt1,bins=bins,xlabel='# subunits',wintitle=timestamp,bigCfraction=fraction,**kwargs)
        else:
            _plot_clustersize_counts(cts, ctsgt1,bins=bins,wintitle=timestamp,bigCfraction=fraction,**kwargs)
        if psu is not None:
            _plot_clustersize_counts(cts/4.0/psu, ctsgt1/4.0/psu, xlabel='# RyRs, corrected', bins=bins,wintitle=timestamp,
                                     bigCfraction=fraction,**kwargs)
    
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

    axs[0, 0].scatter(t_s,p['x_site_nc'],s=0.3,c='black',alpha=0.7)
    if plotSmoothingCurve:
        axs[0, 0].plot(t_s,p['x_ori']-p['x'],'r',alpha=0.4)
    axs[0, 0].set_ylim(-15,15)
    axs[0, 0].set_xlabel('t [s]')
    axs[0, 0].set_ylabel('x [nm]')
        
    axs[0, 1].scatter(t_s,p['y_site_nc'],s=0.3,c='black',alpha=0.7)
    if plotSmoothingCurve:
        axs[0, 1].plot(t_s,p['y_ori']-p['y'],'r',alpha=0.4)
    axs[0, 1].set_ylim(-15,15)
    axs[0, 1].set_xlabel('t [s]')
    axs[0, 1].set_ylabel('y [nm]')
        
    axs[1, 0].scatter(t_s,p['z_site_nc'],s=0.3,c='black',alpha=0.7)
    if plotSmoothingCurve:
        axs[1, 0].plot(t_s,p['z_ori']-p['z'],'r',alpha=0.4)
    axs[1, 0].set_ylim(-15,15)
    axs[1, 0].set_xlabel('t [s]')
    axs[1, 0].set_ylabel('z [nm]')

    ax = axs[1,1]
    if plotSmoothingCurve:
        # plot the MBM track
        ax.plot(t_s,p['x_ori']-p['x_nc'],alpha=0.5,label='x')
        plt.plot(t_s,p['y_ori']-p['y_nc'],alpha=0.5,label='y')
        if 'z_nc' in p.keys():
            ax.plot(t_s,p['z_ori']-p['z_nc'],alpha=0.5,label='z')
        ax.set_xlabel('t (s)')
        ax.set_ylabel('MBM corr [nm]')
        ax.legend()
    else:
        axs[1, 1].plot(t_s,p['x_ori']-p['x'])
        axs[1, 1].plot(t_s,p['y_ori']-p['y'])
        axs[1, 1].plot(t_s,p['z_ori']-p['z'])
        axs[1, 1].set_xlabel('t [s]')
        axs[1, 1].set_ylabel('orig. corr [nm]')
    plt.tight_layout()

from PYMEcs.Analysis.MINFLUX import analyse_locrate
from PYMEcs.misc.guiMsgBoxes import Error
from PYMEcs.misc.utils import unique_name

# we try to find an MBM collection attached
# to an MBMcorrection generated data source
# returns None if unsuccesful
def findmbm(pipeline,warnings=True,return_mod=False):
    from PYMEcs.recipes.localisations import MBMcorrection
    dsname = None
    # search/check for instance
    for mod in pipeline.recipe.modules:
        if isinstance(mod,MBMcorrection):
            dsname = mod.output
            break
    if dsname is None:
        if warnings:
            warn("we rely on mbm info present in a datasource generated by the MBMcorrection module. Can't find such a datasource, aborting...")
        return None
    mbm = pipeline.dataSources[dsname].mdh.get('Processing.MBMcorrection.mbm')
    if mbm is None:
        if warnings:
            warn("found no mbm collection in metadata of datasource 'dsname' generated by MBMcorrection module, aborting..." % dsname)
        return None
    if return_mod:
        return mod
    else:
        return mbm


def get_metadata_from_mfx_attrs(mfx_attrs):
    mfx_itrs = mfx_attrs['measurement']['threads'][0]['sequences'][0]['Itr']
    md_by_itrs = pd.DataFrame(columns=['IterationNumber','PinholeAU','ActivationLaser', 'ExcitationLaser',
                                       'ExcitationWavelength_nm', 'ExcitationPower_percent', 'ExcitationDAC',
                                       'DetectionChannel01','DetectionChannel02','BackgroundThreshold',
                                       'PhotonLimit', 'CCRLimit', 'DwellTime_ms',
                                       'PatternGeoFactor','PatternRepeat', 'PatternGeometry','Strategy'],
                              index=range(len(mfx_itrs)))
    for i, itr in enumerate(mfx_itrs):
        md_by_itrs.loc[i].IterationNumber = i
        md_by_itrs.loc[i].PinholeAU = itr['Mode']['phDiaAU']
        md_by_itrs.loc[i].ActivationLaser = itr['_activation']['laser'] if itr['_activation']['laser'] != '' else 'not used'
        md_by_itrs.loc[i].ExcitationLaser = itr['_excitation']['laser']
        md_by_itrs.loc[i].ExcitationWavelength_nm = 1e9*itr['_excitation']['wavelength']
        md_by_itrs.loc[i].ExcitationPower_percent = itr['_excitation']['power']
        md_by_itrs.loc[i].ExcitationDAC = itr['_excitation']['dac']
        md_by_itrs.loc[i].DetectionChannel01 = itr['_detection']['channels'][0]
        md_by_itrs.loc[i].DetectionChannel02 = itr['_detection']['channels'][1] if len(itr['_detection']['channels']) >1 else 'NA'
        md_by_itrs.loc[i].BackgroundThreshold = itr['bgcThreshold']
        md_by_itrs.loc[i].PhotonLimit = itr['phtLimit']
        md_by_itrs.loc[i].CCRLimit = itr['ccrLimit']
        md_by_itrs.loc[i].DwellTime_ms = 1e3*itr['patDwellTime']
        md_by_itrs.loc[i].PatternGeoFactor = itr['patGeoFactor']
        md_by_itrs.loc[i].PatternRepeat = itr['patRepeat']
        md_by_itrs.loc[i].PatternGeometry = itr['Mode']['pattern']
        md_by_itrs.loc[i].Strategy = itr['Mode']['strategy']

    return md_by_itrs

from PYME.recipes.traits import HasTraits, Float, Enum, CStr, Bool, Int, List
import PYME.config

class MINFLUXSettings(HasTraits):
    withOrigamiSmoothingCurves = Bool(True,label='Plot smoothing curves',desc="if overplotting smoothing curves " +
                                      "in origami site correction analysis")
    defaultDatasourceForAnalysis = CStr('Localizations',label='default datasource for analysis',
                                        desc="the datasource key that will be used by default in the MINFLUX " +
                                        "properties functions (EFO, localisation rate, etc)") # default datasource for acquisition analysis
    defaultDatasourceForMBM = CStr('coalesced_nz',label='default datasource for MBM analysis and plotting',
                                        desc="the datasource key that will be used by default in the MINFLUX " +
                                        "MBM analysis") # default datasource for MBM analysis
    datasourceForClusterAnalysis = CStr(PYME.config.get('MINFLUX-clusterDS','dbscan_clustered'),label='datasource for 3D cluster analysis',
                                        desc="the datasource key that will be used to generate the 3D cluster size analysis")
    
    #MBM_lowess_fraction = Float(0.03,label='lowess fraction for MBM smoothing',
    #                                    desc='lowess fraction used for smoothing of coalesced MBM data (default 0.05)')
    origamiWith_nc = Bool(False,label='add 2nd moduleset (no MBM corr)',
                          desc="if a full second module set is inserted to also analyse the origami data without any MBM corrections")

class DateString(HasTraits):
    TimeStampString = CStr('',label="Time stamp",desc='the time stamp string in format yymmdd-HHMMSS')


class MBMaxisSelection(HasTraits):
    SelectAxis = Enum(['x-y-z','x','y','z','std_x','std_y','std_z'])


class MINFLUXanalyser():
    def __init__(self, visFr):
        self.visFr = visFr
        self.minfluxRIDs = {}
        self.origamiErrorFignum = 0
        self.origamiTrackFignum = 0
        self.analysisSettings = MINFLUXSettings()
        self.dstring = DateString()
        self.mbmAxisSelection = MBMaxisSelection()
        
        visFr.AddMenuItem('MINFLUX', "Localisation Error analysis", self.OnErrorAnalysis)
        visFr.AddMenuItem('MINFLUX', "Cluster sizes - 3D", self.OnCluster3D)
        visFr.AddMenuItem('MINFLUX', "Cluster sizes - 2D", self.OnCluster2D)
        visFr.AddMenuItem('MINFLUX', "Analyse Localization Rate", self.OnLocalisationRate)
        visFr.AddMenuItem('MINFLUX', "EFO histogram (photon rates)", self.OnEfoAnalysis)
        visFr.AddMenuItem('MINFLUX', "plot tracking correction (if available)", self.OnTrackPlot)
        visFr.AddMenuItem('MINFLUX>Origami', "group and analyse origami sites", self.OnOrigamiSiteRecipe)
        visFr.AddMenuItem('MINFLUX>Origami', "plot origami site correction", self.OnOrigamiSiteTrackPlot)
        visFr.AddMenuItem('MINFLUX>Origami', "plot origami error estimates", self.OnOrigamiErrorPlot)
        visFr.AddMenuItem('MINFLUX', "Analysis settings", self.OnMINFLUXSettings)
        visFr.AddMenuItem('MINFLUX', "Manually create Colour panel", self.OnMINFLUXColour)
        visFr.AddMenuItem('MINFLUX>Util', "Plot temperature record matching current data series",self.OnMINFLUXplotTempData)
        visFr.AddMenuItem('MINFLUX>Util', "Set MINFLUX temperature file location", self.OnMINFLUXsetTempDataFile)
        visFr.AddMenuItem('MINFLUX>Util', "Check if clumpIndex contiguous", self.OnClumpIndexContig)
        visFr.AddMenuItem('MINFLUX>MBM', "Plot mean MBM info (and if present origami info)", self.OnMBMplot)
        visFr.AddMenuItem('MINFLUX>MBM', "Show MBM tracks", self.OnMBMtracks)
        visFr.AddMenuItem('MINFLUX>MBM', "Add MBM track labels to view", self.OnMBMaddTrackLabels)
        visFr.AddMenuItem('MINFLUX>MBM', "Save MBM bead trajectories to npz file", self.OnMBMSave)
        visFr.AddMenuItem('MINFLUX>MBM', "Save MBM bead settings to json file", self.OnMBMSettingsSave)
        visFr.AddMenuItem('MINFLUX>RyRs', "Plot corner info", self.OnCornerplot)
        visFr.AddMenuItem('MINFLUX>RyRs', "Plot density stats", self.OnDensityStats)
        visFr.AddMenuItem('MINFLUX>RyRs', "Show cluster alpha shapes", self.OnAlphaShapes)
        visFr.AddMenuItem('MINFLUX>Zarr', "Show MBM attributes", self.OnMBMAttributes)
        visFr.AddMenuItem('MINFLUX>Zarr', "Show MFX attributes", self.OnMFXAttributes)
        visFr.AddMenuItem('MINFLUX>Zarr', "Show MFX metadata info (experimental)", self.OnMFXInfo)
        
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
        import io
        from  wx.lib.dialogs import ScrolledMessageDialog
        fres = self.visFr.pipeline.dataSources['FitResults']
        if 'zarr' not in dir(fres):
            warn("could not find zarr attribute - is this a MFX zarr file?")
            return
        else:
            try:
                mfx_attrs = fres.zarr['mfx'].attrs.asdict()
            except AttributeError:
                warn("could not access MFX attributes - do we have MFX data in zarr?")
                return
            if '_legacy' in mfx_attrs:
                warn("legacy data detected - no useful MFX metadata in legacy data")
                return
            md_info = get_metadata_from_mfx_attrs(mfx_attrs)
            with io.StringIO() as output:
                print(md_info.to_string(show_dimensions=False,index=True,line_width=80),file=output)
                mfx_info_str = output.getvalue()
            with ScrolledMessageDialog(self.visFr, mfx_info_str, "MFX info (tentative)", size=(900,400),
                                        style=wx.RESIZE_BORDER | wx.DEFAULT_DIALOG_STYLE ) as dlg:
                dlg.ShowModal()            
    
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
        p = self.visFr.pipeline
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
        fig, axs = plt.subplots(nrows=3)
        for caxis, ax in zip(['x','y','z'],axs):
            if has_drift:
                if has_drift_ori:
                    ax.plot(t_s,p['drift%s' % caxis]+p['drift%s_ori' % caxis], label='origami 2nd pass')
                    ax.plot(t_s,p['drift%s_ori' % caxis],'--', label='origami 1st pass')
                else:
                    ax.plot(t_s,p['drift%s' % caxis], label='origami 1st pass')
            if has_mbm:
                mod = findmbm(p,warnings=False,return_mod=True)
                MBM_lowess_fraction = mod.MBM_lowess_fraction
                mbm_mean[caxis] = mbm.mean(caxis)
                ax.plot(mbm.t,mbm_mean[caxis],':',label='MBM mean')
                from statsmodels.nonparametric.smoothers_lowess import lowess
                if MBM_lowess_fraction > 1e-5:
                    mbm_meansm[caxis] = lowess(mbm_mean[caxis], mbm.t, frac=MBM_lowess_fraction,
                                               return_sorted=False)
                else:
                    mbm_meansm[caxis] = mbm_mean[caxis]
                ax.plot(mbm.t,mbm_meansm[caxis],'-.',label='MBM lowess (lf=%.2f)' % MBM_lowess_fraction)
            if has_mbm2:
                ax.plot(t_s,p['mbm%s' % caxis], label='MBM from module')
            ax.set_xlabel('time (s)')
            ax.set_ylabel('drift in %s (nm)' % caxis)
            ax.legend(loc="upper right")
        fig.tight_layout()
        if has_mbm: # also plot a second figure without the non-smoothed MBM track
            fig, axs = plt.subplots(nrows=3)
            for caxis, ax in zip(['x','y','z'],axs):
                if has_drift:
                    if has_drift_ori:
                        ax.plot(t_s,p['drift%s' % caxis]+p['drift%s_ori' % caxis], label='origami 2nd pass')
                        ax.plot(t_s,p['drift%s_ori' % caxis],'--', label='origami 1st pass')
                    else:
                        ax.plot(t_s,p['drift%s' % caxis], label='origami 1st pass')
                if has_mbm:
                    #ax.plot(mbm.t,mbm_mean[caxis],':',label='MBM mean')
                    ax.plot(mbm.t,mbm_meansm[caxis],'r-.',label='MBM lowess (lf=%.2f)' % MBM_lowess_fraction)
                ax.set_xlabel('time (s)')
                ax.set_ylabel('drift in %s (nm)' % caxis)
                ax.legend(loc="upper left")
            fig.tight_layout()
        if has_mbm: # also plot a third figure with all MBM tracks
            fig, axs = plt.subplots(nrows=3)
            for caxis, ax in zip(['x','y','z'],axs):
                if has_drift:
                    if has_drift_ori:
                        ax.plot(t_s,p['drift%s' % caxis]+p['drift%s_ori' % caxis], label='origami 2nd pass')
                        ax.plot(t_s,p['drift%s_ori' % caxis],'--', label='origami 1st pass')
                    else:
                        ax.plot(t_s,p['drift%s' % caxis], label='origami 1st pass')
                if has_mbm:
                    #ax.plot(mbm.t,mbm_mean[caxis],':',label='MBM mean')
                    ax.plot(mbm.t,mbm_meansm[caxis],'r-.',label='MBM lowess (lf=%.2f)' % MBM_lowess_fraction)
                    mbm.plot_tracks_matplotlib(caxis,ax=ax,goodalpha=0.4)
                ax.set_xlabel('time (s)')
                ax.set_ylabel('drift in %s (nm)' % caxis)
                ax.legend(loc="upper left")
            fig.tight_layout()
        if has_mbm and has_drift: # also plot a fourth figure with a difference track for all axes
            tnew = 1e-3*p['t']
            mbmcorr = {}
            for axis in ['x','y','z']:
                axis_interp_msm = np.interp(tnew,mbm.t,mbm_meansm[axis])          
                mbmcorr[axis] = axis_interp_msm
 
            fig, axs = plt.subplots(nrows=3)
            for caxis, ax in zip(['x','y','z'],axs):
                if has_drift:
                    if has_drift_ori:
                        ax.plot(t_s,p['drift%s' % caxis]+p['drift%s_ori' % caxis]-mbmcorr[caxis], label='diff to origami 2nd pass')
                    else:
                        ax.plot(t_s,p['drift%s' % caxis]-mbmcorr[caxis], label='diff to origami 1st pass')
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

    def OnMINFLUXsetTempDataFile(self, event):
        import PYME.config as config
        with wx.FileDialog(self.visFr, "Choose Temperature data file", wildcard='CSV (*.csv)|*.csv',
                           style=wx.FD_OPEN) as dialog:
            if dialog.ShowModal() == wx.ID_CANCEL:
                return
        fname = dialog.GetPath()
        
        if config.get('MINFLUX-temperature_file') == fname:
            warn("config option 'MINFLUX-temperature_file' already set to %s" % fname)
            return # already set to this value, return

        config.update_config({'MINFLUX-temperature_file': fname},
                             config='user', create_backup=True)


    def OnMINFLUXplotTempData(self, event):
        import PYME.config as config
        if config.get('MINFLUX-temperature_file') is None:
            warn("Need to set Temperature file location first")
            return
        from PYMEcs.misc.utils import read_temp_csv, set_diff, timestamp_to_datetime
        mtemps = read_temp_csv(config.get('MINFLUX-temperature_file'))
        if len(self.visFr.pipeline.dataSources) == 0:
            warn("no datasources, this is probably an empty pipeline, have you loaded any data?")
            return
        t0 = self.visFr.pipeline.mdh.get('MINFLUX.TimeStamp')
        if t0 is None:
            warn("no MINFLUX TimeStamp in metadata, giving up")
            return
        set_diff(mtemps,timestamp_to_datetime(t0))
        p = self.visFr.pipeline
        range = (1e-3*p['t'].min(),1e-3*p['t'].max())
        sertemps = mtemps[mtemps['tdiff_s'].between(range[0],range[1])]
        if sertemps.empty:
            warn("no records in requested time window, is series time before or after start/end of available temperature records?")
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
        plot_errors(self.visFr.pipeline)

    def OnCluster3D(self, event):
        plot_cluster_analysis(self.visFr.pipeline, ds=self.analysisSettings.datasourceForClusterAnalysis)

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

    def OnOrigamiSiteRecipe(self, event=None):
        from PYMEcs.recipes.localisations import OrigamiSiteTrack, DBSCANClustering2
        from PYME.recipes.localisations import MergeClumps
        from PYME.recipes.tablefilters import FilterTable, Mapping
        
        pipeline = self.visFr.pipeline
        recipe = pipeline.recipe

        preFiltered = unique_name('prefiltered',pipeline.dataSources.keys())
        corrSiteClumps = unique_name('corrected_siteclumps',pipeline.dataSources.keys())
        siteClumps = unique_name('siteclumps',pipeline.dataSources.keys())
        dbscanClusteredSites = unique_name('dbscanClusteredSites',pipeline.dataSources.keys())
        sites = unique_name('sites',pipeline.dataSources.keys())
        sites_c = unique_name('sites_c',pipeline.dataSources.keys())
        
        curds = pipeline.selectedDataSourceKey
        modules = [FilterTable(recipe,inputName=curds,outputName=preFiltered,
                               filters={'error_x' : [0,3.3],
                                        'error_z' : [0,3.3]}),
                   DBSCANClustering2(recipe,inputName=preFiltered,outputName=dbscanClusteredSites,
                                     searchRadius = 15.0,
                                     clumpColumnName = 'siteID',
                                     sizeColumnName='siteClumpSize'),
                   FilterTable(recipe,inputName=dbscanClusteredSites,outputName=siteClumps,
                               filters={'siteClumpSize' : [3,40]}),
                   MergeClumps(recipe,inputName=siteClumps,outputName=sites,
                               labelKey='siteID',discardTrivial=True),
                   OrigamiSiteTrack(recipe,inputClusters=siteClumps,inputSites=sites,outputName=corrSiteClumps,
                                    labelKey='siteID'),
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
                                   filters={'error_x' : [0,3.3],
                                            'error_z' : [0,3.3]}),
                       DBSCANClustering2(recipe,inputName=preFiltered,outputName=dbscanClusteredSites,
                                         searchRadius = 15.0,
                                         clumpColumnName = 'siteID',
                                         sizeColumnName='siteClumpSize'),
                       Mapping(recipe,inputName=dbscanClusteredSites,outputName=dbsnc,
                               mappings={'x': 'x_nc', 'y': 'y_nc', 'z': 'z_nc'}),
                       FilterTable(recipe,inputName=dbsnc,outputName=siteClumps,
                                   filters={'siteClumpSize' : [3,40]}),
                       MergeClumps(recipe,inputName=siteClumps,outputName=sites,
                                   labelKey='siteID',discardTrivial=True),
                       OrigamiSiteTrack(recipe,inputClusters=siteClumps,inputSites=sites,outputName=corrSiteClumps,
                                        labelKey='siteID'),
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

    def OnOrigamiErrorPlot(self, event):
        p = self.visFr.pipeline
        # need to add checks if the required properties are present in the datasource!!

        def plot_errs(ax,axisname,errkeys):
            ax.hist(p[errkeys[0]],bins='auto',alpha=0.5,density=True,label='Trace est')
            ax.hist(p[errkeys[1]],bins='auto',alpha=0.5,density=True,label='Site est')
            ax.hist(p[errkeys[2]],bins='auto',alpha=0.5,density=True,label='Site est corr')
            ax.legend()
            ax.set_xlabel('error %s (nm)' % axisname)
            ax.set_ylabel('#')
        
        fig, axs = plt.subplots(2, 2,num='origami error estimates %d' % self.origamiErrorFignum)
        plot_errs(axs[0, 0], 'x', ['error_x_ori','error_x_nc','error_x'])
        plot_errs(axs[0, 1], 'y', ['error_y_ori','error_y_nc','error_y'])
        plot_errs(axs[1, 0], 'z', ['error_z_ori','error_z_nc','error_z'])
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
