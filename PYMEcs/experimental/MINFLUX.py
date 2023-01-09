import matplotlib.pyplot as plt
import numpy as np

def plot_errors(pipeline):
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
    

def plot_cluster_analysis(pipeline, ds='dbscanClustered'):
    curds = pipeline.selectedDataSourceKey
    pipeline.selectDataSource(ds)
    p = pipeline
    uids, cts = np.unique(p['dbscanClumpID'], return_counts=True)
    ctsgt1 = cts[cts > 1.1]
    pipeline.selectDataSource(curds)
    plt.figure()
    plt.subplot(221)
    h = plt.hist(cts,bins=15)
    plt.xlabel('Cluster Size')
    plt.plot([np.mean(cts),np.mean(cts)],[0,h[0].max()])
    plt.plot([np.median(cts),np.median(cts)],[0,h[0].max()],'--')
    plt.subplot(222)
    h = plt.hist(ctsgt1,bins=15)
    plt.xlabel('Cluster Size (clusters > 1)')
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
    pipeline.selectDataSource(curds)

def plot_intra_clusters_dists(pipeline, ds='dbscanClustered',bins=15):
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
        dd, ii = tree.query(coords,k=2)
        dists.extend(list(dd[:,1]))
    pipeline.selectDataSource(curds)
    plt.figure()
    h=plt.hist(dists,bins=bins)

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

def fourcornerplot(pipeline,sigma=None,backgroundFraction=0.0):
    ccs = cornercounts(pipeline,backgroundFraction=backgroundFraction)
    ccsn = ccs/ccs.sum()
    ks = np.arange(4)+1
    popt, pcov = curve_fit(pnn, ks, ccsn,sigma=sigma)
    perr = np.sqrt(np.diag(pcov))
    plt.figure()
    ax = plt.subplot(111)
    ax.bar(ks-0.4, ccsn, width=0.4, color='b', align='center')
    ax.bar(ks, pnn(ks,popt[0]), width=0.4, color='g', align='center')
    ax.legend(['Experimental data', 'Fit'])
    print('optimal p: %.3f +- %.3f' % (popt[0],perr[0]))
    p_missed = pn(0,popt[0])
    p_m_min = pn(0,popt[0]+perr[0])
    p_m_max = pn(0,popt[0]-perr[0])
    print('missed fraction: %.2f (%.2f...%.2f)' % (p_missed,p_m_min,p_m_max))
    return (popt[0],perr[0],pnn(ks,popt[0]),ccsn)


sigmaDefault = [0.1,0.1,0.02,0.02]
backgroundDefault = 0.15

def fourcornerplot_default(pipeline):
    return fourcornerplot(pipeline,sigma=sigmaDefault,backgroundFraction=backgroundDefault)

class MINFLUXanalyser():
    def __init__(self, visFr):
        self.visFr = visFr

        visFr.AddMenuItem('Experimental>MINFLUX', "Localisation Error analysis", self.OnErrorAnalysis)
        visFr.AddMenuItem('Experimental>MINFLUX', "Cluster sizes - 3D", self.OnCluster3D)
        visFr.AddMenuItem('Experimental>MINFLUX', "Cluster sizes - 2D", self.OnCluster2D)

    def OnErrorAnalysis(self, event):
        plot_errors(self.visFr.pipeline)

    def OnCluster3D(self, event):
        plot_cluster_analysis(self.visFr.pipeline, ds='dbscanClustered')

    def OnCluster2D(self, event):
        plot_cluster_analysis(self.visFr.pipeline, ds='dbscan2D')

        
def Plug(visFr):
    return MINFLUXanalyser(visFr)
