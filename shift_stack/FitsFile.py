from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
import time
import subprocess
import os
import numpy as np
import sep


def exec_pv2sip(file_full_name: str) -> bool:
    try:
        p = subprocess.Popen(['wcs-pv2sip', file_full_name, file_full_name.replace(".fits", "_SIP.fits")])
        p.wait()
    except:
        print("Error converting PV2SIP!")
        return False
    else:
        return True




class FitsFile:
    def __init__(self, file_path: str, file_name: str):
        """
        Constructor of class FitsFile.
        :param file_path: Directory that contains the FITS file.
        :param file_name: FITS filename.
        """
        if not file_path.endswith("/"):
            file_path += "/"
        self.file_full_name = file_path + file_name
        self.file_sip_name = self.file_full_name.replace(".fits", "_SIP.fits")
        try:
            file = fits.open(self.file_full_name)
        except FileNotFoundError:
            print("File %s not found in %s" % (file_name, file_path))
            raise SystemExit
        else:
            self.header = file[0].header
            file.close()
            self.time = self.header['date-wrt']  # Get UTC time.
            self.timestamp = time.mktime(time.strptime(self.time[0:19], "%Y-%m-%dT%H:%M:%S"))
            self.exposure_time = self.header['exptime']  # Get exposure time.
            self.width = self.header['naxis1']
            self.height= self.header['naxis2']
            self.sigma = self.header['sigmback']
            self.back  = self.header['meanback']

            if 'PV1_5' in self.header.keys():
                # SCAMP header
                if not os.path.exists(self.file_sip_name):
                    exec_pv2sip(self.file_full_name)
                file_sip = fits.open(self.file_sip_name)
                self.wcs = WCS(file_sip[0].header)
                self.wcs.fix()
                file_sip.close()
            elif 'A_0_0' in self.header.keys():
                # SIP header
                file_sip = fits.open(self.file_full_name)
                self.wcs = WCS(file_sip[0].header)
                self.wcs.fix()
                file_sip.close()

    def get_data(self, origin_coord: SkyCoord=None, height: int=4000, width: int=4000) -> np.ndarray:
        """
        Get image data from FITS file related with current FitsFile instance.
        :return: Image Data Matrix
        """
        # file = fits.open(self.file_full_name)
        self.roi_height = height
        self.roi_width = width
        if origin_coord is None:
            origin_pix = (int(self.height / 2 - height / 2), int(self.width / 2 - width / 2))
        else:
            origin_pix_reverse = self.wcs.wcs_world2pix(origin_coord.ra.deg, origin_coord.dec.deg, True)
            origin_pix_reverse = [round(float(item)) for item in origin_pix_reverse]
            origin_pix = [origin_pix_reverse[1], origin_pix_reverse[0]]

        if not width/2 < origin_pix[0] < self.width - width // 2 and height/2 < origin_pix[1] < self.height - height//2:
            print("ROI out of boundary in fits %s" % self.file_full_name)
            return np.zeros((height, width), dtype=np.uint16)

        try:
            file = fits.open(self.file_full_name)
            self.orgin_height = origin_pix[0]
            self.orgin_width = origin_pix[1]
            data = file[0].data[self.orgin_height - height//2:self.orgin_height + height//2,
                                     self.orgin_width - width//2:self.orgin_width + width//2]
            data = np.array(data, dtype=np.uint16)
            file.close()

        except TypeError:
            print('file {0} has been damaged...'.format(self.file_full_name))
            data = np.zeros((height, width), dtype=np.uint16)

        return data

    """

    def set_mask(self, threshold: float):
       
        # Using sep to detect bright stars
        threshold = float(threshold[0])
        self.data = self.data.copy(order='C')
        self.data = np.array(self.data, dtype=np.int32)
        try:
            bkg = sep.Background(self.data)
            # print('backgroundvalue', bkg.globalback, 'globalrms', bkg.globalrms, sep=' ')

        except IndexError:
            self.data = np.zeros((self.roi_height, self.roi_width), dtype=np.int32)
        else:
            self.data[self.mask_function(bkg, threshold)] = 0

    def mask_function(self, bkg, threshold: float) -> np.ndarray:
        # data[data < 0] = 0
        # data[data > threshold * bkg.globalrms] = 0
        bkg.subfrom(self.data)
        objs = sep.extract(self.data, threshold)
        mask = np.zeros((self.roi_height, self.roi_height), dtype=bool)
        sep.mask_ellipse(mask, objs['x'], objs['y'], objs['a'], objs['b'], objs['theta'], r= 5.0)
        return mask

    """

