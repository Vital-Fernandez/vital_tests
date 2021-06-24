import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style, quantity_support
plt.style.use(astropy_mpl_style)
quantity_support()

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

target_coords = SkyCoord(274.088176, -3.404173, unit="deg")
print(target_coords)

las_campanas = EarthLocation(lat=29.0182*u.deg, lon=70.6915*u.deg, height=390*u.m)
utcoffset = -4*u.hour  # Chilean time
time = Time('2021-6-10 23:00:00') - utcoffset

orion_AltAz = target_coords.transform_to(AltAz(obstime=time, location=las_campanas))

midnight = Time('2012-7-13 00:00:00') - utcoffset
delta_midnight = np.linspace(-2, 10, 100)*u.hour
frame_July13night = AltAz(obstime=midnight+delta_midnight,
                          location=las_campanas)
orion_AltAz_night_period = target_coords.transform_to(frame_July13night)
orion_Airmass_night_period = orion_AltAz_night_period.secz

plt.plot(delta_midnight, orion_Airmass_night_period)
# plt.xlim(-2, 10)
# plt.ylim(1, 4)
plt.xlabel('Hours from EDT Midnight')
plt.ylabel('Airmass [Zenith angle]')
plt.show()