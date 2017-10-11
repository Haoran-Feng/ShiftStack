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
            self.height = self.header['naxis2']

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
            # print("While stacking {0}, ra:{1} dec:{2} pix X:{3} pix Y:{4}".format(self.file_full_name,
            #                                                                       origin_coord.ra.deg,
            #                                                                       origin_coord.dec.deg,
            #                                                                       origin_pix[1],
            #                                                                       origin_pix[0]))

        if not width/2 < origin_pix[0] < self.width - width // 2 and height/2 < origin_pix[1] < self.height - height//2:
            print("ROI out of boundary in fits %s" % self.file_full_name)
            self.data = np.zeros((height, width), dtype=np.int32)
            return self.data

        try:
            file = fits.open(self.file_full_name)
            self.orgin_height = origin_pix[0]
            self.orgin_width = origin_pix[1]
            self.data = file[0].data[self.orgin_height - height//2:self.orgin_height + height//2,
                                     self.orgin_width - width//2:self.orgin_width + width//2]
            self.data = np.array(self.data, dtype=np.int32)
            file.close()

        except TypeError:
            print('file {0} has been damaged...'.format(self.file_full_name))
            self.data = np.zeros((width, height), dtype=np.int32)

        return self.data

    def set_mask(self, threshold: float):
        """
        Set mask on bright stars.
        :return:
        """
        # Using sep to detect bright stars
        self.data = self.data.copy(order='C')
        self.data = np.array(self.data, dtype=np.int32)
        try:
            bkg = sep.Background(self.data)
            print('backgroundvalue', bkg.globalback, 'globalrms', bkg.globalrms, sep=' ')
        except IndexError:
            self.data = np.zeros((self.roi_width, self.roi_height), dtype=np.int32)
        else:
            self.data = self.data - bkg
            self.data[self.data < 0] = 0
            self.data[self.data > threshold * bkg.globalrms] = 0
