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
from astropy.visualization.stretch import SinhStretch, LogStretch
from astropy.visualization import ImageNormalize, MinMaxInterval


###############################################################################
# We start with the sample data
stereo = (a.vso.Source('STEREO_A') &
          a.Instrument('EUVI') & a.Wavelength(304*u.Angstrom) &
          a.Time('2012-06-12T10:00:00', '2012-06-12T10:10:00'))
#result = Fido.search(stereo)
#downloaded_file = Fido.fetch(result)
downloaded_file = '/Users/schriste/sunpy/data/secchi_l0_a_img_euvi_20120612_20120612_100615_n4eua.fts'
print(downloaded_file)

m = sunpy.map.Map(downloaded_file)#.rotate()

##############################################################################
# First we make a mask, which tells us which regions are bright. We
# choose the criterion that the data should be at least 5% of the maximum
# value. Pixels with intensity values greater than this are included in the
# mask, while all other pixels are excluded.
pixel_area = m.scale.axis2 * m.scale.axis1
normed_data = m.data.copy() / (m.exposure_time * pixel_area)
normed_map = sunpy.map.Map(normed_data, m.meta)

normed_data[normed_data.value < 500] = 0
data = ndimage.gaussian_filter(normed_data.value, 9)

filtered_map = sunpy.map.Map(data, m.meta)

#filtered_map.peek()
labels, n = ndimage.label(data)

com = []
for i in range(n):
    ndimage.measurements.center_of_mass(labels == i)
    com.append(ndimage.measurements.center_of_mass(labels == i))

plt.figure()
ax = plt.subplot(projection=m)
# Create interval object
interval = MinMaxInterval()
vmin, vmax = interval.get_limits(m.data)
norm = ImageNormalize(vmin=vmin, vmax=vmax*0.8, stretch=LogStretch())

ax.imshow(m.data, norm=norm, cmap=plt.cm.bone)
ax.contour(labels, colors='white')
plt.figtext(0.3, 0.2, 'Number of regions = {}'.format(n))

plt.show()
