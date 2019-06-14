"""
================================================
Finding backside Active regions using STEREO 304
================================================

How you can to find the brightest regions in an AIA image and
count the approximate number of regions of interest using ndimage.
"""
# sphinx_gallery_thumbnail_number = 2

from scipy import ndimage
import matplotlib.pyplot as plt

import astropy.units as u

import sunpy.map
from sunpy.net import Fido, attrs as a


###############################################################################
# We start with the sample data
stereo = (a.vso.Source('STEREO_A') &
          a.Instrument('EUVI') & a.Wavelength(304*u.Angstrom) &
          a.Time('2012-06-12T10:00:00', '2012-06-12T10:10:00'))
result = Fido.search(stereo)
downloaded_file = Fido.fetch(result)
print(downloaded_file)

m = sunpy.map.Map(downloaded_file)

##############################################################################
# First we make a mask, which tells us which regions are bright. We
# choose the criterion that the data should be at least 5% of the maximum
# value. Pixels with intensity values greater than this are included in the
# mask, while all other pixels are excluded.
pixel_area = m.wcs.wcs.cdelt[0] * m.wcs.wcs.cdelt[1] * u.arcsec ** 2
normed_data = m.data.copy() / (m.exposure_time * pixel_area)
normed_map = sunpy.map.Map(normed_data, m.meta)

normed_data[normed_data.value < 500] = 0
data = ndimage.gaussian_filter(normed_data.value, 20)
#data[data < 500] = 0

filtered_map = sunpy.map.Map(data, m.meta)


filtered_map.peek()
labels, n = ndimage.label(data)

plt.figure()
ax = plt.subplot(projection=aiamap)
plt.imshow(m.data)
plt.contour(labels)
plt.figtext(0.3, 0.2, 'Number of regions = {}'.format(n), color='white')

plt.show()
