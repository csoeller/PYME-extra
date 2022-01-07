# -*- coding: utf-8 -*-

import numpy as np
from scipy import ndimage

def genRef(refimage, normalised=True):
        X, Y = np.mgrid[0.0:refimage.shape[0], 0.0:refimage.shape[1]]
        X -= refimage.shape[0]/2
        Y -= refimage.shape[1]/2
        bdry = 20
        mask = np.ones_like(refimage[:,:,0],dtype='float64')
        mask[:bdry, :] = 0
        mask[-bdry:, :] = 0
        mask[:, :bdry] = 0
        mask[:,-bdry:] = 0

        maskg = ndimage.gaussian_filter(mask.astype('float'),5)

        calImages = np.zeros(refimage.shape[:2] + (21,))
        calFTs = np.zeros(refimage.shape[:2] + (21,), dtype='complex64')
        #refimage3D = refimage.squeeze()

        for i in range(refimage.shape[2]):
                if not normalised:
                        d = refimage[:,:,i]
                        ref = d/d.mean() - 1
                else:
                        ref = refimage[:,:,i]
                
                calFTs[:,:,i] = np.fft.ifftn(ref)
                calImages[:,:,i] = ref*maskg

        dz = np.gradient(calImages)[2].reshape(-1, 21)
        dzn = np.hstack([1./np.dot(dz[:,i], dz[:,i]) for i in range(21)])

        return calImages, calFTs, dz, dzn, maskg, X, Y

def compare(calImages, calFTs, dz, dzn, posInd, image, mask, X, Y, normalised=False, deltaZ = 0.2):        
        d = 1.0*image

        if not normalised:
                dm = d/d.mean() - 1
        else:
                dm = d
                
        FA = calFTs[:,:,posInd]
        refA = calImages[:,:,posInd] 

        ddz = dz[:,posInd]
        dznn = dzn[posInd]

        C = np.fft.ifftshift(np.abs(np.fft.ifftn(np.fft.fftn(dm)*FA)))
        #C = ifftshift(np.abs(ifftn(fftn(A)*ifftn(B))))
        Cm = C.max()    
        
        Cp = np.maximum(C - 0.5*Cm, 0)
        Cpsum = Cp.sum()
        
        dx = (X*Cp).sum()/Cpsum
        dy = (Y*Cp).sum()/Cpsum

        ds = ndimage.shift(dm, [-dx, -dy])*mask
        ds_A = (ds - refA)

        dzz = deltaZ*np.dot(ds_A.ravel(), ddz)*dznn

        return dx, dy, dzz, Cm, dm
