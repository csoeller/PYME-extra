import matplotlib.colors as mcol
import matplotlib.pyplot as plt
import numpy as np
import pylab
from PYME.misc.colormaps import cm


def r(rgb):
    return rgb[0]

def g(rgb):
    return rgb[1]

def b(rgb):
    return rgb[2]


def hot_overflow(underflowcol = 'magenta', overflowcol = 'lime', percentage=5):
    if percentage < 1:
        percentage = 1
    if percentage > 15:
        percentage = 15

    ucolrgb = mcol.hex2color(mcol.cnames[underflowcol])
    ocolrgb = mcol.hex2color(mcol.cnames[overflowcol])
    p = 0.01 * percentage
    
    # edited from hot_data in matplotlib
    hot_data = {'red':   [(0, r(ucolrgb), r(ucolrgb)),
                          (p, r(ucolrgb), 0.0416),
                          (0.365079, 1.000000, 1.000000),
                          (1.0-p, 1.0, r(ocolrgb)),
                          (1.0, r(ocolrgb), r(ocolrgb))],
                'green': [(0, g(ucolrgb), g(ucolrgb)),
                          (p, g(ucolrgb), 0.),
                          (0.365079, 0.000000, 0.000000),
                          (0.746032, 1.000000, 1.000000),
                          (1.0-p, 1.0, g(ocolrgb)),
                          (1.0, g(ocolrgb), g(ocolrgb))],
                'blue':  [(0, b(ucolrgb), b(ucolrgb)),
                          (p, b(ucolrgb), 0.),
                          (0.746032, 0.000000, 0.000000),
                          (1.0-p, 1.0, b(ocolrgb)),
                          (1.0, b(ocolrgb), b(ocolrgb))]}

    cm_hot2 = mcol.LinearSegmentedColormap('hot_overflow', hot_data)
    return cm_hot2


def grey_overflow(underflowcol = 'magenta', overflowcol = 'lime', percentage=5, greystart=0.1):
    if percentage < 1:
        percentage = 1
    if percentage > 15:
        percentage = 15

    ucolrgb = mcol.hex2color(mcol.cnames[underflowcol])
    ocolrgb = mcol.hex2color(mcol.cnames[overflowcol])
    p = 0.01 * percentage

    grey_data = {'red':   [(0, r(ucolrgb), r(ucolrgb)),
                          (p, r(ucolrgb), greystart),
                          (1.0-p, 1.0, r(ocolrgb)),
                          (1.0, r(ocolrgb), r(ocolrgb))],
                'green': [(0, g(ucolrgb), g(ucolrgb)),
                          (p, g(ucolrgb), greystart),
                          (1.0-p, 1.0, g(ocolrgb)),
                          (1.0, g(ocolrgb), g(ocolrgb))],
                'blue':  [(0, b(ucolrgb), b(ucolrgb)),
                          (p, b(ucolrgb), greystart),
                          (1.0-p, 1.0, b(ocolrgb)),
                          (1.0, b(ocolrgb), b(ocolrgb))]}

    cm_grey2 = mcol.LinearSegmentedColormap('grey_overflow', grey_data)
    return cm_grey2


def grey_overflow2(underflowcol = 'magenta', overflowcol = 'lime', percentage=5):
    if percentage < 1:
        percentage = 1
    if percentage > 30:
        percentage = 30
    clist = [(0,mcol.cnames[underflowcol]),
             # (0.01*percentage,mcol.cnames['dimgrey']),
             (0.01*percentage,(0.15,0.15,0.15)),
             (1.0-0.01*percentage,mcol.cnames['white']),
             (1.0,mcol.cnames[overflowcol])]

    return mcol.LinearSegmentedColormap.from_list('grey_overflow2',clist)

# this function needs some work
def overflow_wrap(cmapname, underflowcol = 'magenta', overflowcol = 'lime', percentage=5):
    cmap = pylab.get_cmap(cmapname)
    if percentage < 0:
        percentage = 0
    if percentage > 10:
        percentage = 10
    p = 0.01*percentage
    def mymap(c):
        ret = cmap(c)
        if c < p:
            return colors.hex2color(colors.cnames[underflowcol])
        elif c > 1.0-p:
            return colors.hex2color(colors.cnames[overflowcol])
        else:
            return cmap(c/(1.0-2*p)+p)
        
    mymap.name = cmap.name + '_overflow'
    return mymap

# FIJI inspired colour maps

# cyan hot (for pylab devide all values by 255.0)
# red: 0 0, 170 0, 255 255
# green: 0 0, 170 255, 255, 255
# blue: 0 0, 85, 255, 255 255

# cb_skyblue
# red: 0 0, 255 85
# green: 0 0, 255 180
# blue: 0 0, 255 232

# generate a 3-tupel for continuous colour and normalise to 0..1
def contc(x,y):
    return (x/255.0,y/255.0,y/255.0)

cyan_hot_d = {
    'red'   : [contc(0,0),contc(170,0),contc(255,255)],
    'green' : [contc(0,0),contc(170,255),contc(255,255)],
    'blue'  : [contc(0,0),contc(85,255),contc(255,255)]
}

cb_skyblue_d = {
    'red'   : [contc(0,0),contc(255,85)],
    'green' : [contc(0,0),contc(255,180)],
    'blue'  : [contc(0,0),contc(255,232)]
}

cb_skyblue2_d = {
    'red'   : [contc(0,0),contc(200,85),contc(255,255)],
    'green' : [contc(0,0),contc(200,180),contc(255,255)],
    'blue'  : [contc(0,0),contc(200,232),contc(255,255)]
}

cm_ch = mcol.LinearSegmentedColormap('cyan_hot',cyan_hot_d)
cm_cbskb = mcol.LinearSegmentedColormap('cb_skyblue',cb_skyblue_d)
cm_cbskb2 = mcol.LinearSegmentedColormap('cb_skyblue2',cb_skyblue2_d)

cm.update({cmap.name: cmap for cmap in [cm_ch,cm_cbskb,cm_cbskb2,hot_overflow(overflowcol='cyan',percentage=2.5)]})


############################
### if colorcet is installed
### register a few cms
############################
try:
    import colorcet as cc
    has_colorcet = True
except ImportError:
    has_colorcet = False

if has_colorcet:
    cm.update({'cet_fire':cc.m_fire,
               'cet_coolwarm': cc.m_coolwarm,
               'cet_bmw': cc.m_bmw,
               'cet_bmy': cc.m_bmy})

def main():
    import numpy as np
    a=np.outer(np.arange(0,1,0.01),np.ones(10))
    # cm1 = grey_overflow()
    # cm1 = grey_overflow2(percentage=2)
    cm1 = hot_overflow(overflowcol='cyan')
    # cm1 = overflow_wrap('jet')
    plt.imshow(a,aspect='auto',cmap=cm1,origin="lower")

if __name__ == "__main__":
    main()
