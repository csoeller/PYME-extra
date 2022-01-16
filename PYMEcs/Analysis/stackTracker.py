import numpy as np
from scipy import ndimage
from numpy.fft import fftn, ifftn, fftshift, ifftshift
import matplotlib.pyplot as plt

# we use this code as a simple and preliminary way to calculate the z-factor
# this should eventually be reworked as a nicer codebase and also replace/subsume the
# offline tracker code base

# note that this may or may not work!
def zdiff(data):
    zdf = np.zeros_like(data,dtype='f')
    zdf[:,:,0:-1] = data[:,:,1:]-1.0*data[:,:,0:-1]
    zdf[:,:,-1] = zdf[:,:,-2]
    return zdf

# subsample stack from dz0 spacing to dzs sampling, dzs should be even multiple of dz0
# note: we assume the center image is the target position
def substack(stack, subsamplefactor, newszhalf=None):
        ssfac = int(subsamplefactor)
        zsz = stack.shape[2]
        zszh = zsz // 2
        if 2*zszh+1 != zsz:
                raise RuntimeError("z dimension must be odd")
        if newszhalf is not None:
                if newszhalf > zszh:
                        raise RuntimeError("new stack size must be smaller or equal old stack size")
                nstack = stack[:,:,zszh-newszhalf:zszh+newszhalf+1]
                zszh = newszhalf
        else:
                nstack = stack

        halfperiods = zszh // ssfac
        refims = nstack[:,:,zszh-halfperiods*ssfac:zszh+halfperiods*ssfac+1:ssfac]

        return (nstack, refims)
      

# this one should initialise FFTs etc for a suitable substack and return all the relevant
# items in a dict (may make class in a future version)
# return everything as a stackobject, initially just a dict with the required entries
def initialise_data(stack, subsamplefactor, vsznm, newszhalf=None):
        nstack, refimages_raw = substack(stack, subsamplefactor, newszhalf=newszhalf)
        refimages, calImages, calFTs, gradI, gradIsqr_inv, mask, X, Y = genRefData(refimages_raw)

        trackobject = {
                'stack': nstack,
                'refimages_raw': refimages_raw,
                'refimages': refimages,
                'calImages': calImages,
                'calFTs': calFTs,
                'gradI': gradI,
                'gradIsqr_inv': gradIsqr_inv,
                'mask': mask,
                'X': X,
                'Y': Y,
                'voxelsize_nm': vsznm,
                'subsamplefactor': int(subsamplefactor)}

        return trackobject


# this one takes an initialised stackobject and returns the relevant shifts by taking
# the center of the substack as the relevant level to refer to
def get_shifts_from_stackobject(trackobject):
        to = trackobject
        resx = []
        resy = []
        resz = []

        refim_zszh = to['refimages'].shape[2] // 2
        for i in range(to['stack'].shape[2]):
                image = to['stack'][:,:,i]
                driftx, drifty, driftz, cm, dm = compare(to['refimages'],
                                                         to['calImages'],
                                                         to['calFTs'],
                                                         to['gradI'],
                                                         to['gradIsqr_inv'],
                                                         refim_zszh, image,
                                                         to['mask'],
                                                         to['X'],
                                                         to['Y'],
                                                         deltaZ=to['subsamplefactor']*to['voxelsize_nm'].z)
                resx.append(driftx)
                resy.append(drifty)
                resz.append(driftz)


        dznm = np.array(resz)
        dxnm = to['voxelsize_nm'].x*np.array(resx)
        dynm = to['voxelsize_nm'].y*np.array(resy)

        return (dxnm,dynm,dznm)

def fit_and_plot_zf(dxnm,dynm,dznm,trackobject):

        to = trackobject
        zszh = dznm.shape[0] // 2
        dzexp = dznm[zszh-4:zszh+5]
        dztheo = (np.arange(dznm.shape[0])-zszh)*to['voxelsize_nm'].z

        x = dztheo[zszh-3:zszh+4]
        y = dznm[zszh-3:zszh+4]
        m, b = np.polyfit(x,y,1)

        plt.figure()
        plt.plot(dztheo[zszh-4:zszh+5],dznm[zszh-4:zszh+5],'-o')
        plt.plot(dztheo[zszh-4:zszh+5],m*dztheo[zszh-4:zszh+5],'--')

        zfactor = 1.0 / m
        font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 14,
        }
        plt.text(-150, 50, 'Z-factor = %.2f' % zfactor, fontdict=font)



def genRefData(refimages_raw, bdry=10, useSimplediff=False):
        X, Y = np.mgrid[0.0:refimages_raw.shape[0], 0.0:refimages_raw.shape[1]]
        X -= np.ceil(refimages_raw.shape[0]*0.5)
        Y -= np.ceil(refimages_raw.shape[1]*0.5)

        mask = np.ones_like(refimages_raw[:,:,0])
        mask[:bdry, :] = 0
        mask[-bdry:, :] = 0
        mask[:, :bdry] = 0
        mask[:,-bdry:] = 0

        refimages = np.zeros_like(refimages_raw,dtype='f')
        calImages = np.zeros_like(refimages_raw,dtype='f')
        calFTs = np.zeros_like(refimages_raw, dtype='complex64')       

        for i in range(refimages_raw.shape[2]):
            d = refimages_raw[:,:,i] # should we work with offset corrected images or irrelevant?
            ref = d/d.mean() - 1

            refimages[:,:,i] = ref
                
            calFTs[:,:,i] = ifftn(ref)
            calImages[:,:,i] = ref*mask

        if useSimplediff:
                gradI = zdiff(calImages).reshape(-1, calImages.shape[2])
        else:
                gradI = np.gradient(calImages)[2].reshape(-1, calImages.shape[2])
        gradIsqr_inv = np.hstack([1./np.dot(gradI[:,i], gradI[:,i]) for i in range(calImages.shape[2])])

        return refimages, calImages, calFTs, gradI, gradIsqr_inv, mask, X, Y


def compare(refImages, calImages, calFTs, gradI, gradIsqr_inv, posInd, image, mask, X, Y, deltaZ = 0.2):        
    d = 1.0*image
    dm = d/d.mean() - 1

    FA = calFTs[:,:,posInd]
    refA = calImages[:,:,posInd]

    gradIA = gradI[:,posInd]
    gradIsqr_invA = gradIsqr_inv[posInd]

    C = ifftshift(np.abs(ifftn(fftn(dm)*FA)))
        
    Cm = C.max()    
        
    Cp = np.maximum(C - 0.5*Cm, 0)
    Cpsum = Cp.sum()
        
    dx = (X*Cp).sum()/Cpsum
    dy = (Y*Cp).sum()/Cpsum
        
    ds = ndimage.shift(dm, [-dx, -dy])*mask
        
    ds_A = (ds - refA)
        
    #calculate z offset between actual position and calibration position
    dz = deltaZ*np.dot(ds_A.ravel(), gradIA)*gradIsqr_invA

    return dx, dy, dz, Cm, dm
