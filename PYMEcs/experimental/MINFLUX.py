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
    

def _plot_clustersize_counts(cts, ctsgt1, xlabel='Cluster Size', **kwargs):
    plt.figure()
    plt.subplot(221)
    h = plt.hist(cts,**kwargs)
    plt.xlabel(xlabel)
    plt.plot([np.mean(cts),np.mean(cts)],[0,h[0].max()])
    plt.plot([np.median(cts),np.median(cts)],[0,h[0].max()],'--')
    plt.subplot(222)
    h = plt.hist(ctsgt1,**kwargs)
    plt.xlabel('%s ( > 1)' % xlabel)
    plt.plot([np.mean(ctsgt1),np.mean(ctsgt1)],[0,h[0].max()])
    plt.plot([np.median(ctsgt1),np.median(ctsgt1)],[0,h[0].max()],'--')
    plt.subplot(223)
    bp_dict = plt.boxplot([cts,ctsgt1],labels=['cluster size','clusters > 1'], showmeans=True)
    for line in bp_dict['means']:
        # get position data for median line
        x, y = line.get_xydata()[0] # top of median line
        # overlay median value
        plt.text(x, y, '%.1f' % y,
                 horizontalalignment='center') # draw above, centered    
    plt.tight_layout()

def plot_cluster_analysis(pipeline, ds='dbscanClustered',showPlot=True, return_means=False, psu=None, bins=15, **kwargs):
    if not ds in pipeline.dataSources:
        warn('no data source named "%s" - check recipe and ensure this is MINFLUX data' % ds)
        return
    curds = pipeline.selectedDataSourceKey
    pipeline.selectDataSource(ds)
    p = pipeline
    uids, cts = np.unique(p['dbscanClumpID'], return_counts=True)
    ctsgt1 = cts[cts > 1.1]
    pipeline.selectDataSource(curds)
    if showPlot:
        if psu is not None:
            _plot_clustersize_counts(cts, ctsgt1,bins=bins,xlabel='# subunits',**kwargs)
        else:
            _plot_clustersize_counts(cts, ctsgt1,bins=bins,**kwargs)
        if psu is not None:
            _plot_clustersize_counts(cts/4.0/psu, ctsgt1/4.0/psu, xlabel='# RyRs, corrected', bins=bins,**kwargs)
    
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
    if showplot:
        plt.figure()
        ax = plt.subplot(111)
        ax.bar(ks-0.4, ccsn, width=0.4, color='b', align='center')
        ax.bar(ks, pnn(ks,popt[0]), width=0.4, color='g', align='center')
        ax.legend(['Experimental data', 'Fit'])
    p_missed = pn(0,popt[0])
    p_m_min = pn(0,popt[0]+perr[0])
    p_m_max = pn(0,popt[0]-perr[0])
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
    xmbm = p['x']-p['x_nc']
    ymbm = p['y']-p['y_nc']
    if has_z:
        zmbm = p['z']-p['z_nc']

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

    
from PYMEcs.Analysis.MINFLUX import analyse_locrate
from PYMEcs.misc.guiMsgBoxes import Error
from PYMEcs.misc.utils import unique_name

from PYME.recipes.traits import HasTraits, Float, Enum, CStr, Bool, Int, List

class MINFLUXSettings(HasTraits):
    withOrigamiSmoothingCurves = Bool(True,label='Plot smoothing curves',desc="if overplotting smoothing curves " +
                                      "in origami site correction analysis")
    defaultDatasourceForAnalysis = CStr('Localizations',label='default datasource for analysis',
                                        desc="the datasource key that will be used by default in the MINFLUX " +
                                        "properties functions (EFO, localisation rate, etc)") # default datasource for acquisition analysis
    defaultDatasourceForMBM = CStr('coalesced_nz',label='default datasource for MBM analysis and plotting',
                                        desc="the datasource key that will be used by default in the MINFLUX " +
                                        "MBM analysis") # default datasource for MBM analysis
    MBM_lowess_fraction = Float(0.03,label='lowess fraction for MBM smoothing',
                                        desc='lowess fraction used for smoothing of coalesced MBM data (default 0.05)')
    origamiWith_nc = Bool(False,label='add 2nd moduleset (no MBM corr)',
                          desc="if a full second module set is inserted to also analyse the origami data without any MBM corrections")

class DateString(HasTraits):
    TimeStampString = CStr('',label="Time stamp",desc='the time stamp string in format yymmdd-HHMMSS')

class MINFLUXanalyser():
    def __init__(self, visFr):
        self.visFr = visFr
        self.minfluxRIDs = {}
        self.origamiErrorFignum = 0
        self.origamiTrackFignum = 0
        self.analysisSettings = MINFLUXSettings()
        self.dstring = DateString()
        
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
        visFr.AddMenuItem('MINFLUX>MBM', "Load MBM data in npz format", self.OnMBMLoadnpz)
        visFr.AddMenuItem('MINFLUX>MBM', "Load JSON MBM bead data config", self.OnMBMLoadJSONbeads)
        
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

    def OnLoadCustom(self, event):
        self.visFr._recipe_manager.LoadRecipe(self.minfluxRIDs[event.GetId()])
        
    def OnMBMLoadnpz(self,event):
        with wx.FileDialog(self.visFr, "Select MBM data npz file", wildcard='NPZ (*.npz)|*.npz',
                           style=wx.FD_OPEN) as dialog:
            if dialog.ShowModal() == wx.ID_CANCEL:
                return
        fname = dialog.GetPath()
        from pathlib import Path
        fp = Path(fname)
        from PYMEcs.Analysis.MBMcollection import MBMCollectionDF
        mbm = MBMCollectionDF(name=fp.stem,filename=fp)

        pfp = Path(self.visFr.pipeline.filename)
        if pfp.suffix != '.npy':
            warn("likely not a MINFLUX npy dataset, extension is %s; aborting..." % pfp.suffix)
            return
        MINFLUXstem = pfp.stem
        if MINFLUXstem not in fp.stem:
            warn("different name stems in MINFLUX and MBM dataset; do they match? %s vs %s" %
                 (MINFLUXstem,fp.stem))
        self.visFr.pipeline.mbm = mbm

    def OnMBMLoadJSONbeads(self,event):
        if not 'mbm' in dir(self.visFr.pipeline):
            warn("no MBM data attached to pipeline, aborting...")
            return
        
        with wx.FileDialog(self.visFr, "Select MBM bead selection JSON file", wildcard='JSON (*.json)|*.json',
                           style=wx.FD_OPEN) as dialog:
            if dialog.ShowModal() == wx.ID_CANCEL:
                return
        fname = dialog.GetPath()
        from pathlib import Path
        fp = Path(fname)
        if not self.visFr.pipeline.mbm.name in fp.stem:
            warn("JSON file name and mbm names do not match; are you sure? Will try to load anyway but results may vary...")
        import json
        with open(fname) as f:
            mbmsettings = json.load(f)

        for bead in mbmsettings['beads']:
            self.visFr.pipeline.mbm.beadisgood[bead] =  mbmsettings['beads'][bead]
        self.visFr.pipeline.mbm.median_window = mbmsettings['Median_window']

    def OnMINFLUXsetTempDataFile(self, event):
        import PYME.config as config
        with wx.FileDialog(self.visFr, "Choose Temperature data file", wildcard='CSV (*.csv)|*.csv',
                           style=wx.FD_OPEN) as dialog:
            if dialog.ShowModal() == wx.ID_CANCEL:
                return
        fname = dialog.GetPath()
        
        if config.get('MINFLUX_temperature_file') == fname:
            warn("config option 'MINFLUX_temperature_file' already set to %s" % fname)
            return # already set to this value, return

        config.update_config({'MINFLUX_temperature_file': fname},
                             config='user', create_backup=True)


    def OnMINFLUXplotTempData(self, event):
        import PYME.config as config
        if config.get('MINFLUX_temperature_file') is None:
            warn("Need to set Temperature file location first")
            return
        from PYMEcs.misc.utils import read_temp_csv, set_diff, parse_timestamp_from_filename
        mtemps = read_temp_csv(config.get('MINFLUX_temperature_file'))
        if len(self.visFr.pipeline.dataSources) == 0:
            warn("no datasources, this is probably an empty pipeline, have you loaded any data?")
            return
        try:
            fname = self.visFr.pipeline.filename
        except AttributeError:
            warn("no filename associated with pipeline, returning")
            return
        t0 = parse_timestamp_from_filename(fname)
        if t0 is None:
            if not self.dstring.configure_traits(kind='modal'):
                return
            else:
                t0 = parse_timestamp_from_filename(self.dstring.TimeStampString)
                if t0 is None:
                    warn("entered time stamp '%s' does not parse, giving up" % self.dstring.TimeStampString)
                    return
        set_diff(mtemps,t0)
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
        plot_cluster_analysis(self.visFr.pipeline, ds='dbscanClustered')

    def OnCluster2D(self, event):
        plot_cluster_analysis(self.visFr.pipeline, ds='dbscan2D')

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
        plot_tracking(p,is_coalesced,lowess_fraction=self.analysisSettings.MBM_lowess_fraction)
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
        # also plot post correction!
        t_s = 1e-3*p['t']
        fig, axs = plt.subplots(2, 2,num='origami site tracks %d' % self.origamiTrackFignum)
        axs[0, 0].scatter(t_s,p['x_site_nc'],s=0.3,c='black',alpha=0.7)
        if self.analysisSettings.withOrigamiSmoothingCurves:
            axs[0, 0].plot(t_s,p['x_ori']-p['x'],'r',alpha=0.4)
        axs[0, 0].set_ylim(-15,15)
        axs[0, 0].set_xlabel('t [s]')
        axs[0, 0].set_ylabel('x [nm]')
        
        axs[0, 1].scatter(t_s,p['y_site_nc'],s=0.3,c='black',alpha=0.7)
        if self.analysisSettings.withOrigamiSmoothingCurves:
            axs[0, 1].plot(t_s,p['y_ori']-p['y'],'r',alpha=0.4)
        axs[0, 1].set_ylim(-15,15)
        axs[0, 1].set_xlabel('t [s]')
        axs[0, 1].set_ylabel('y [nm]')
        
        axs[1, 0].scatter(t_s,p['z_site_nc'],s=0.3,c='black',alpha=0.7)
        if self.analysisSettings.withOrigamiSmoothingCurves:
            axs[1, 0].plot(t_s,p['z_ori']-p['z'],'r',alpha=0.4)
        axs[1, 0].set_ylim(-15,15)
        axs[1, 0].set_xlabel('t [s]')
        axs[1, 0].set_ylabel('z [nm]')

        ax = axs[1,1]
        if self.analysisSettings.withOrigamiSmoothingCurves:
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


def Plug(visFr):
    # we are trying to monkeypatch pipeline and VisGUIFrame methods to sneak MINFLUX npy IO in;
    # in future we will ask for a way to get this considered by David B for a proper hook
    # in the IO code
    from PYMEcs.IO.MINFLUX import monkeypatch_npy_io
    monkeypatch_npy_io(visFr)
        
    return MINFLUXanalyser(visFr)
