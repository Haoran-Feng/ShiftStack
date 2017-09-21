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

    def get_data(self, origin_coord: SkyCoord=None, width: int=4000, height: int=4000) -> np.ndarray:
        """
        Get image data from FITS file related with current FitsFile instance.
        :return: Image Data Matrix
        """
        # file = fits.open(self.file_full_name)
        self.roi_width = width
        self.roi_height = height
        if origin_coord is None:
            origin_pix = (int(self.width / 2 - width / 2), int(self.height / 2 - height / 2))
        else:
            origin_pix = self.wcs.wcs_world2pix(origin_coord.ra.deg, origin_coord.dec.deg, True)
            origin_pix = [round(float(item)) for item in origin_pix]

        if not width/2 < origin_pix[0] < self.width - width // 2 and height/2 < origin_pix[1] < self.height - height//2:
            print("ROI out of boundary in fits %s" % self.file_full_name)
            self.data = np.zeros((width, height), dtype=np.int32)
            return self.data

        try:
            file = fits.open(self.file_full_name)
            self.orgin_x = origin_pix[0]
            self.orgin_y = origin_pix[1]
            self.data = file[0].data.T[self.orgin_x - width//2:self.orgin_x + width//2,
                                       self.orgin_y - height//2:self.orgin_y + height//2]
            self.data = np.array(self.data, dtype=np.int32)
            file.close()

        except TypeError:
            print('file {0} has been damaged...'.format(self.file_full_name))
            self.data = np.zeros((width, height), dtype=np.int32)

        return self.data

    def set_mask(self, threshold: float=2, r: float=1.5):
        """
        Set mask on bright stars.
        :return:
        """
        # Using sep to detect bright stars
        self.data = self.data.copy(order='C')
        self.data = np.array(self.data, dtype=np.int32)
        try:
            bkg = sep.Background(self.data)
        except IndexError:
            self.data = np.zeros((self.roi_width, self.roi_height), dtype=np.int32)
        else:
            self.data = self.data - bkg
            objects = sep.extract(self.data, threshold, err=bkg.globalrms)
            self.objects_list = [item for item in objects]
            self.objects_list.sort(key=lambda item:item['peak'], reverse=False)
            # Using sep to mask bright stars
            self.mask = np.zeros(self.data.shape, dtype=bool)
            sep.mask_ellipse(self.mask, objects['x'], objects['y'], objects['a'], objects['b'], objects['theta'], r=r)
            self.data[self.mask] = 0
        # fits.writeto(self.file_full_name + "_mid2.fits", self.data)

