import sys
import argparse

def main():
    # options parsing
    op = argparse.ArgumentParser(description='convert old style PSFs to new style.')
    op.add_argument('psfname', metavar='psfname', default=None,
                    help='filename of old style PSF')
    op.add_argument('outname', metavar='outname', default=None,
                    help='filename for converted PSF')
    op.add_argument('-d', '--debug', action='store_true')
    
    args = op.parse_args()

    try:
        import cPickle as pickle
    except ImportError:
        import pickle
    from PYME.IO.image import ImageStack
    
    # we need to trick pickle into thinking that
    #    PYME.Acquire.MetaDataHandler is still available
    import PYME.IO.MetaDataHandler as MD
    sys.modules['PYME.Acquire.MetaDataHandler'] = MD
    with open(args.psfname,'rb') as inpfi:
        psf = pickle.load(inpfi)

    if args.debug:
        print "psf: %s, %s" % (repr(psf[0].shape), repr(psf[0].dtype))
        print psf[1]

    del sys.modules['PYME.Acquire.MetaDataHandler']

    psfdata = psf[0]
    psfvox = psf[1]
    mdh = MD.NestedClassMDHandler()
    mdh['voxelsize'] = psfvox
    import os
    filename, file_extension = os.path.splitext(args.outname)
    if file_extension.startswith('.tif'):
        mdh['ImageType'] = 'PSF'
    # print(mdh)
    ImageStack(psfdata,mdh=mdh).Save(filename=args.outname)

    
if __name__ == "__main__":
    main()
