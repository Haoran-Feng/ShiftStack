from astropy.io import fits
import sep
import numpy as np

if __name__ == '__main__':
    f = fits.getdata("L20160204_03666_090227+1657_90s_SR_123476.fits")
    f = np.array(f, dtype=np.int32)
    background = sep.Background(f)
    print(background.globalrms)
    print(background.globalback)
    data = f - background
#    objects = sep.extract(data, thresh=1.8, err=background.globalrms)
#    array = np.zeros(data.shape, dtype=bool)
#    sep.mask_ellipse(array, objects['x'], objects['y'], objects['a'], objects['b'], objects['theta'], r=6)
#    data[array] = 0
    data[data > 2 * background.globalrms] = data.min()

    fits.writeto('mask_test_result.fits', data, overwrite=True)