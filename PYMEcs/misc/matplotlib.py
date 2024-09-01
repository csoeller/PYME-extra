import matplotlib.pyplot as plt
import numpy as np
from numbers import Number

def scattered_boxplot(ax, x, notch=None, sym=None, vert=None, whis=None, positions=None,
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

setattr(plt.Axes, "scattered_boxplot", scattered_boxplot)
