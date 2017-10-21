from shift_stack.FitsFile import FitsFile
from astropy.coordinates import SkyCoord
import os
from enum import Enum
import numpy as np
import astropy.units as u
from tqdm import tqdm
import skimage.morphology as sm


class ImgSeq:
    """
    Class ImgSeq (Image Sequence) manipulates all the image files and
    the processing procedure.
    """
    def __init__(self, path: str,
                 list_all_fits: bool = False,
                 origin_coord: SkyCoord = None,
                 ROI_height = None,
                 ROI_width  = None,
                 target_speed_ra: float = None,
                 target_speed_dec: float= None,
                 index_begin = None,
                 index_end   = None,
                 if_mask = True,
                 dilation_radius = 2
                 ):

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

        if list_all_fits:
            print("No.   Cover                         Time            ")
            print("--------------------------------------------------------------------------------------")
            for index, item in enumerate(self.fits_objects):
                corner_row = [item.height - 1, item.height - 1, 0, 0]
                corner_col = [0, item.width - 1, 0, item.width - 1]
                wcs_world_coord = item.wcs.wcs_pix2world(corner_col, corner_row, True)
                coord = SkyCoord(wcs_world_coord[0], wcs_world_coord[1], unit='deg')
                print("{0}     {1} {2} {3}".format(index, coord[0].to_string('hmsdms'),
                                                   coord[1].to_string('hmsdms'), item.time))
                print("       {0} {1}".format(coord[2].to_string('hmsdms'), coord[3].to_string('hmsdms')))

            print("--------------------------------------------------------------------------------------")
        else:
            # Member variables
            self.origin_coord = origin_coord
            self.ROI_height = ROI_height
            self.ROI_width = ROI_width
            self.target_speed_ra = target_speed_ra
            self.target_speed_dec = target_speed_dec
            self.if_mask = if_mask
            self.dilation_radius = dilation_radius
            max_fits_number = self.fits_objects.__len__()

            if index_begin is None:
                index_begin = 0
                index_end   = self.fits_objects.__len__() - 1
            else:
                index_begin = index_begin if 0 <= index_begin < max_fits_number else 0
                index_end = index_end if index_begin <= index_end < max_fits_number else max_fits_number

            self.fits_objects = self.fits_objects[index_begin: index_end+1]
            # print("Using fits image {0}->{1}".format(index_begin, index_end+1))
            self.actual_fits_number = index_end + 1 - index_begin
            time_gap = [self.fits_objects[i].timestamp - self.fits_objects[0].timestamp
                             for i in range(len(self.fits_objects))]
            time_gap = np.array(time_gap)
            coord_ra = self.origin_coord.ra.deg + time_gap * (self.target_speed_ra / 3600) / 3600
            coord_dec = self.origin_coord.dec.deg + time_gap * (self.target_speed_dec / 3600) / 3600
            self.coordinates = SkyCoord(coord_ra, coord_dec, unit='deg')

    def read_into_memory(self):
        """
        Read the data we need into memory (Optional, cost more resource but faster)
        """
        # Ready memory
        self.img_data_list = np.zeros((self.actual_fits_number, self.ROI_height, self.ROI_width),
                                      dtype=np.uint16)
        self.globalback = np.zeros(self.actual_fits_number, dtype=np.uint16)
        self.sigmback   = np.zeros(self.actual_fits_number, dtype=np.float32)
        with tqdm(total=self.actual_fits_number, ncols=90, desc='Read data into memory') as pbar:
            for index, item in enumerate(self.fits_objects):
                self.img_data_list[index] = item.get_data(self.coordinates[index],
                                                          self.ROI_height, self.ROI_width)
                self.globalback[index] = np.uint16(item.back)
                self.sigmback[index] = item.sigma
                pbar.update()

    def remove_background(self):
        with tqdm(total=self.actual_fits_number, ncols=90, desc="Removing image background") as pbar:
            for index, item in enumerate(self.img_data_list):
                back = self.globalback[index]
                item[item <= back] = 0
                item[item > back] -= back
                pbar.update()

    def generate_mask(self):
        self.static_sum = np.zeros((self.ROI_height, self.ROI_width), dtype=np.uint16)
        with tqdm(total=self.actual_fits_number, ncols=90, desc="Generating mask...") as pbar:
            for item in self.img_data_list:
                self.static_sum += item
                pbar.update()
            # bw_threshold = 50 * self.actual_fits_number * self.sigmback.mean()
            bw_threshold = 40 * (65535 // self.actual_fits_number)
            mask = np.array(self.static_sum, dtype=np.uint16)
            mask[mask <= bw_threshold] = 0
            mask[mask > bw_threshold]  = 1
            self.mask = np.array(mask, dtype=np.bool)
            sm.binary_dilation(mask, selem=sm.disk(self.dilation_radius, dtype=np.bool), out=self.mask)

    def shift_stack_in_memory(self)->np.ndarray:

        # ADD some notations here, something like speed_ra speed_dec

        out_image = np.zeros((self.ROI_height, self.ROI_width), dtype=np.uint32)
        with tqdm(total=self.actual_fits_number, ncols=90, desc="Stacking images") as pbar:
            for item in self.img_data_list:
                if self.if_mask:
                    item[self.mask] = 0
                out_image += item
                pbar.update()
        return out_image






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

        self.fits_objects[0].get_data(origin_coord, region_height, region_width)
        if if_remove_field_star:
            self.fits_objects[0].set_mask(threshold=threshold)
        with tqdm(total=len(self.fits_objects), ncols=90) as pbar:
            self.img_data = self.fits_objects[0].data
            pbar.update(1)
            for index, gap in enumerate(self.fits_time_gap):
                # print("FITS No.%s stacking..." % (index + 2))
                delta_ra = gap * (target_speed_ra / 3600) / 3600
                delta_dec = gap * (target_speed_dec / 3600) / 3600
                new_origin = SkyCoord(origin_coord.ra.deg + delta_ra,
                                      origin_coord.dec.deg + delta_dec, unit='deg')
                origin_coord = new_origin
                self.fits_objects[index + 1].get_data(new_origin, region_height, region_width)
                if if_remove_field_star:
                    self.fits_objects[index + 1].set_mask(threshold=threshold)
                self.img_data += self.fits_objects[index + 1].data
                pbar.update(1)

        return self.img_data
