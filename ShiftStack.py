import argparse
from astropy.coordinates import SkyCoord
from shift_stack.ImgSeq import ImgSeq
from shift_stack.FitsFile import FitsFile

def _parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "source",
        help="Specify the path that contains FITS files.",
    )
    parser.add_argument(
        "speed_ra",
        help=""
    )
    parser.add_argument(
        "speed_dec",
        help=""
    )
    parser.add_argument(
        "width",
        help=""
    )
    parser.add_argument(
        "height",
        help=""
    )
    parser.add_argument(
        "Output",
        help="Specify the path that you wanna store the output FITS file."
    )
    return parser.parse_args()

if __name__ == '__main__':
    args = _parse_args()
    coord = SkyCoord()
    imseq = ImgSeq(args.source)
    imseq.shift_stack()

