import argparse
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from shift_stack.ImgSeq import ImgSeq
from shift_stack.FitsFile import FitsFile

import os

def _parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "source",
        metavar="FitsFilesPath",
        type=str,
        help="Specify the path that contains FITS files.",
        default=os.getcwd()
    )

    parser.add_argument(
        "--coord_pixel",
        dest='center_coord_pixel',
        nargs=3,
        metavar=("FitsName", "X", "Y"),
        type=str,
        default=None,
        help="Specify a Fits File and stacking around one of its pixels."
    )

    parser.add_argument(
        "--coord",
        dest='center_coord',
        nargs=2,
        metavar=('(RA in hh:mm:ss.xxx)', '(DEC in dd:mm:ss.xxx)'),
        type=str,
        default=None,
        help="The center of stacking image."
    )

    parser.add_argument(
        "--speed",
        dest='speed',
        nargs=2,
        metavar=('(RA in arcsec/hr)', "(DEC in arcsec/hr)"),
        type=float,
        help="Speed of this stacking. Default to be 0, 0.",
        default=(0, 0)
    )
    parser.add_argument(
        "--list",
        action="store_true",
        dest="if_list",
        help="Flag for listing FITS image information."
    )
    parser.add_argument(
        "--number",
        dest="fits_collection_to_stack",
        nargs=2,
        type=int,
        metavar=("START", "END"),
        help="Use --list to check the FITS list before selecting images to stack"

    )
    parser.add_argument(
        "--size",
        type=int,
        nargs=2,
        dest='size',
        metavar=('Height', 'Width'),
        help="Height and Width of the output image.",
        default=(1000, 1000)
    )

    parser.add_argument(
        "--threshold",
        nargs=1,
        dest='threshold',
        metavar='Threshold.',
        default=1.8
    )

    parser.add_argument(
        "--nomask",
        action='store_false',
        dest='mask',
        help='Flag for removing field stars.'
    )
    parser.add_argument(
        "Output",
        help="Specify the name that you wanna store the output FITS file."
    )
    return parser.parse_args()

'''
def list_fits_files(imseq: ImgSeq):
    print("No.   Cover                         Time            ")
    print("----------------------------------------------------")
    for item, index in enumerate(imseq.fits_objects):
        corner_row = [item.height-1, item.height-1, 0, 0]
        corner_col = [0, item.width-1, 0, item.width-1]
        wcs_world_coord = item.wcs.wcs_pix2world()
        
        
        left_down = SkyCoord(item.wcs.wcs_pix2world(0, 0, True), )
        print("{0}  {1} {2} {3}".format(index, ))
        print("     {0} {1}")
'''

if __name__ == '__main__':
    args = _parse_args()

    if args.if_list:
        imseq = ImgSeq(args.source)

    if args.center_coord_pixel is None:
        coord = SkyCoord(args.center_coord[0], args.center_coord[1], frame='icrs', unit=(u.hourangle, u.deg))
    else:
        basis_fits = FitsFile(args.source, args.center_coord_pixel[0])
        coord = basis_fits.wcs.wcs_pix2world(int(args.center_coord_pixel[1]), int(args.center_coord_pixel[2]),True)
        pix_coord = basis_fits.wcs.wcs_world2pix(coord[0], coord[1], True)
        coord = SkyCoord(coord[0], coord[1], frame='icrs', unit=(u.deg, u.deg))

    imseq = ImgSeq(args.source)
    result = imseq.shift_stack(coord, target_speed_ra=args.speed[0], target_speed_dec=args.speed[1],
                               region_width=args.size[1], region_height=args.size[0], if_remove_field_star=args.mask,
                               threshold=args.threshold)
    fits.writeto(args.Output, result, overwrite=True)

