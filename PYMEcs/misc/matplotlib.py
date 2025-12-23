import matplotlib.pyplot as plt
import numpy as np
from numbers import Number

def figuredefaults(fontsize=10,linewidth=1.5):
    plt.rcParams['axes.linewidth'] = linewidth
    plt.rcParams['xtick.major.width'] = linewidth
    plt.rcParams['ytick.major.width'] = linewidth
    plt.rcParams['font.family']='arial'
    plt.rcParams['font.size'] = fontsize
    plt.rcParams['svg.fonttype'] = 'none'

def boxswarmplot(df,width=0.4,annotate_means=False,annotate_medians=False,showmeans=True,
                 meanprops=None,ax=None,swarmsize=7,swarmalpha=None,format="%.1f",
                 showpoints=True,strip=False,**kwargs):
    import seaborn as sns
    if meanprops is None:
        meanprops={'marker':'o',
                    'markerfacecolor':'white', 
                    'markeredgecolor':'black',
                    'markersize':'8'}
    flierprops = dict(marker='.', markerfacecolor='none', markersize=0, linestyle='none')
    colours = ['#72a4cdff'] * df.shape[1]
    bp = sns.boxplot(df,boxprops={'facecolor': 'none'},width=width,showmeans=showmeans,meanprops=meanprops,ax=ax,flierprops=flierprops,**kwargs)
    if showpoints:
        if strip:
            sns.stripplot(df,size=swarmsize,palette=colours,ax=ax,alpha=swarmalpha)
        else:
            sns.swarmplot(df,size=swarmsize,palette=colours,ax=ax,alpha=swarmalpha)
    if annotate_medians:
        meds = df.median().values
        for i,xtick in enumerate(bp.get_xticks()):
            bp.text(xtick-0.5*width-0.1,meds[i],format % meds[i], 
                    horizontalalignment='center',verticalalignment='center',
                    size='x-small',color='black',weight='semibold')
    if annotate_means:
        means = df.mean().values
        for i,xtick in enumerate(bp.get_xticks()):
            bp.text(xtick+0.5*width+0.1,means[i],format % means[i], 
                    horizontalalignment='center',verticalalignment='center',
                    size='x-small',color='g',weight='semibold')
    return bp

def violinswarmplot(df,width=0.8,annotate_means=False,annotate_medians=False,
                    ax=None,swarmsize=7,swarmalpha=None,format="%.1f",strip=False,
                    annotate_width=0.4,showpoints=True,**kwargs):
    import seaborn as sns
    vp = sns.violinplot(df, ax=ax, color="0.95", width=width, **kwargs)
    colours = ['#72a4cdff'] * df.shape[1]
    if showpoints:
        if strip:
            sns.stripplot(df,size=swarmsize,palette=colours,ax=ax,alpha=swarmalpha, jitter=True, zorder=1)
        else:
            sns.swarmplot(df,size=swarmsize,palette=colours,ax=ax,alpha=swarmalpha, zorder=1)
   
    if annotate_medians:
        meds = df.median().values
        for i,xtick in enumerate(vp.get_xticks()):
            vp.text(xtick-0.5*annotate_width-0.1,meds[i],format % meds[i], 
                    horizontalalignment='center',verticalalignment='center',
                    size='x-small',color='black',weight='semibold')
    if annotate_means:
        means = df.mean().values
        for i,xtick in enumerate(vp.get_xticks()):
            vp.text(xtick+0.5*annotate_width+0.1,means[i],format % means[i], 
                    horizontalalignment='center',verticalalignment='center',
                    size='x-small',color='g',weight='semibold')
    return vp


# deprecated
def _scattered_boxplot(ax, x, notch=None, sym=None, vert=None, whis=None, positions=None,
                      widths=None, patch_artist=None, bootstrap=None, usermedians=None, conf_intervals=None,
                      meanline=None, showmeans=None, showcaps=None, showbox=None,
                      showfliers="unif", hide_points_within_whiskers=False,
                      boxprops=None, labels=None, flierprops=None, medianprops=None, meanprops=None,
                      capprops=None, whiskerprops=None, manage_ticks=True, autorange=False, zorder=None, *, data=None):
    if showfliers=="classic":
        classic_fliers=True
    else:
        classic_fliers=False
    bp_dict = ax.boxplot(x, notch=notch, sym=sym, vert=vert, whis=whis, positions=positions, widths=widths,
                         patch_artist=patch_artist, bootstrap=bootstrap, usermedians=usermedians,
                         conf_intervals=conf_intervals, meanline=meanline, showmeans=showmeans, showcaps=showcaps, showbox=showbox,
                         showfliers=classic_fliers,
                         boxprops=boxprops, labels=labels, flierprops=flierprops, medianprops=medianprops, meanprops=meanprops,
                         capprops=capprops, whiskerprops=whiskerprops, manage_ticks=manage_ticks, autorange=autorange, zorder=zorder,data=data)
    N=len(x)
    datashape_message = ("List of boxplot statistics and `{0}` "
                             "values must have same the length")
    # check position
    if positions is None:
        positions = list(np.ones_like(x))
    elif len(positions) != N:
        raise ValueError(datashape_message.format("positions"))

    positions = np.array(positions)
    if len(positions) > 0 and not isinstance(positions[0], Number):
        raise TypeError("positions should be an iterable of numbers")

    # width
    if widths is None:
        widths = [np.clip(0.15 * np.ptp(positions), 0.15, 0.5)] * N
    elif np.isscalar(widths):
        widths = [widths] * N
    elif len(widths) != N:
        raise ValueError(datashape_message.format("widths"))

    if hide_points_within_whiskers:
        import matplotlib.cbook as cbook
        from matplotlib import rcParams
        if whis is None:
            whis = rcParams['boxplot.whiskers']
        if bootstrap is None:
            bootstrap = rcParams['boxplot.bootstrap']
        bxpstats = cbook.boxplot_stats(x, whis=whis, bootstrap=bootstrap,
                                       labels=labels, autorange=autorange)
    for i in range(N):
        if hide_points_within_whiskers:
            xi=bxpstats[i]['fliers']
        else:
            xi=x[i]
        if showfliers=="unif":
            jitter=np.random.uniform(-widths[i]*0.2,widths[i]*0.2,size=np.size(xi))
        elif showfliers=="normal":
            jitter=np.random.normal(loc=0.0, scale=widths[i]*0.1,size=np.size(xi))
        elif showfliers==False or showfliers=="classic":
            return
        else:
            raise NotImplementedError("showfliers='"+str(showfliers)+"' is not implemented. You can choose from 'unif', 'normal', 'classic' and False")

        ax.scatter(positions[i]+jitter,xi,alpha=0.2,marker="o", facecolors='none', edgecolors="k")

    return bp_dict

# setattr(plt.Axes, "scattered_boxplot", scattered_boxplot)
