import matplotlib.colors as mcol
import matplotlib.pyplot as plt

def grey_withOandUFlow():
    clist = [(0,mcol.cnames['red']),
             (0.01,mcol.cnames['black']),
             (0.99,mcol.cnames['white']),
             (1.0,mcol.cnames['green'])]
    
    return mcol.LinearSegmentedColormap.from_list('greyOverflowGreenUnderFlowRed',clist)

def main():
    import numpy as np
    a=np.outer(np.arange(0,1,0.01),np.ones(10))
    cm1 = grey_withOandUFlow()
    plt.imshow(a,aspect='auto',cmap=cm1,origin="lower")

if __name__ == "__main__":
    main()
