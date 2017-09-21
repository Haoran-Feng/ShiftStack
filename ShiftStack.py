import argparse
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from shift_stack.ImgSeq import ImgSeq
from shift_stack.FitsFile import FitsFile

def _parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "source",
        help="Specify the path that contains FITS files.",
    )
    parser.add_argument(
        "Center_ra",
        type=str,
        help=""
    )
    parser.add_argument(
        "Center_dec",
        type=str,
        help=""
    )
    parser.add_argument(
        "speed_ra",
        type=float,
        help=""
    )
    parser.add_argument(
        "speed_dec",
        type=float,
        help=""
    )
    parser.add_argument(
        "width",
        type=int,
        help=""
    )
    parser.add_argument(
        "height",
        type=int,
        help=""
    )
    parser.add_argument(
        "Output",
        help="Specify the name that you wanna store the output FITS file."
    )
    return parser.parse_args()

if __name__ == '__main__':
    args = _parse_args()
    coord = SkyCoord(args.Center_ra, args.Center_dec, frame='icrs', unit=(u.hourangle, u.deg))
    imseq = ImgSeq(args.source)
    result = imseq.shift_stack(coord, target_speed_ra=args.speed_ra, target_speed_dec=args.speed_dec,
                               region_width=args.width, region_height=args.height, if_remove_field_star=True,
                               threshold=1.8)
    fits.writeto(args.Output, result)

