from shift_stack.FitsFile import FitsFile
from astropy.coordinates import SkyCoord
import os
from enum import Enum
import numpy as np
import astropy.units as u
from astropy.utils.console import ProgressBar



class ImgSeq:
    def __init__(self, path: str):
        """
        ImgSeq constructor, read the given path and ready all fits files.
        :param path: Path that contains all the fits files.
        """
        if not path.endswith('/'):
            path += '/'
        if os.path.isdir(path):
            file_list = []
            for file in os.listdir(path):
                if file.endswith('_SIP.fits'):
                    pass
                elif file.endswith('.fits'):
                    file_list.append(file)
        else:
            print("%s not a directory!" % path)
            raise SystemExit

        # Make FitsFile objects for each fits file in the path
        self.fits_objects = [FitsFile(path, item) for item in file_list]

        # Sort FitsFile objects by their time
        self.fits_objects.sort(key=lambda item: item.timestamp)

        # Get the width and height of these fits file.
        self.max_width = self.fits_objects[0].width
        self.max_height = self.fits_objects[0].height

    def get_possible_region(self):
        pass

    def shift_stack(self, origin_coord: SkyCoord,
                          target_speed_ra: float=0, # arcsec/hr
                          target_speed_dec: float=0,# arcsec/hr
                          region_width: int=1000,
                          region_height: int=1000,
                          if_remove_field_star: bool = False,
                          threshold: float=50) -> np.ndarray:
        self.fits_time_gap = \
            list([(self.fits_objects[i + 1].timestamp - self.fits_objects[i].timestamp) for i in range(len(self.fits_objects)-1)])
        print("%s FITS images to stack, with speed_ra=%sarcsec/hr, speed_dec=%sarcsec/hr"
              % (len(self.fits_objects), target_speed_ra, target_speed_dec))

        self.fits_objects[0].get_data(origin_coord, region_width, region_height)
        if if_remove_field_star:
            self.fits_objects[0].set_mask(threshold=threshold, r=5.5)
        self.img_data = self.fits_objects[0].data
        for index, gap in enumerate(self.fits_time_gap):
            print("FITS No.%s stacking..." % (index + 2), end='')
            delta_ra = gap * (target_speed_ra / 3600) / 3600
            delta_dec = gap * (target_speed_dec / 3600) / 3600
            new_origin = SkyCoord(origin_coord.ra.deg + delta_ra,
                                  origin_coord.dec.deg + delta_dec, unit='deg')
            origin_coord = new_origin
            self.fits_objects[index + 1].get_data(new_origin, region_width, region_height)
            if if_remove_field_star:
                self.fits_objects[index + 1].set_mask(threshold=threshold, r=5.5)
            self.img_data += self.fits_objects[index + 1].data
            print('Done')

        return self.img_data.T



