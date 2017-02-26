import sys
import argparse

def main():
    # options parsing
    op = argparse.ArgumentParser(description='convert old style PSFs to new style.')
    op.add_argument('psfname', metavar='psfname', default=None,
                    help='filename of old style PSF')
    op.add_argument('outname', metavar='outname', default=None,
                    help='filename for converted PSF')
    args = op.parse_args()

    try:
        import cPickle as pickle
    except ImportError:
        import pickle

    import PYME.IO.MetaDataHandler as MD
    sys.modules['PYME.Acquire.MetaDataHandler'] = MD
    with open(args.psfname,'rb') as inpfi:
        psf = pickle.load(inpfi)

    del sys.modules['PYME.Acquire.MetaDataHandler']
    with open(args.outname,'wb') as outfi:
        pickle.dump(psf,outfi)

if __name__ == "__main__":
    main()
