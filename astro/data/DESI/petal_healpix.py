import numpy as np
from astropy.table import Table, join
import matplotlib.pyplot as plt

## Making the matplotlib plots look nicer
settings = {
    'font.size':14,
    'axes.linewidth':2.0,
    'xtick.major.size':6.0,
    'xtick.minor.size':4.0,
    'xtick.major.width':2.0,
    'xtick.minor.width':1.5,
    'xtick.direction':'in',
    'xtick.minor.visible':True,
    'xtick.top':True,
    'ytick.major.size':6.0,
    'ytick.minor.size':4.0,
    'ytick.major.width':2.0,
    'ytick.minor.width':1.5,
    'ytick.direction':'in',
    'ytick.minor.visible':True,
    'ytick.right':True
}

plt.rcParams.update(**settings)

specprod = 'fuji'    # Internal name for the EDR
specprod_dir = '.'
print(specprod_dir)

# Tile file from SURVEY='sv3' and PROGRAM='dark'
tile_file = f'{specprod_dir}/ztile-sv3-dark-cumulative.fits'

t_tile = Table.read(tile_file, hdu=1)
ra_tile = t_tile['TARGET_RA']
dec_tile = t_tile['TARGET_DEC']

# We will use Tile #30 for the example below
selec_tile = t_tile['TILEID'] == 30

# print(f'N rows in sv3-dark tile file = {len(t_tile)}')
# print(f'N rows on tile #30 = {len(t_tile[selec_tile])}')

# Healpix file from SURVEY='sv3' and PROGRAM='dark'
zpix_file = f'{specprod_dir}/zpix-sv3-dark.fits'
t_zpix_sv3 = Table.read(zpix_file, hdu=1)

## Rosette 1 tiles (PROGRAM='dark')
is_ros1_tile = ((t_tile['TILEID'] >= 28) & (t_tile['TILEID']<=53)) | (t_tile['TILEID']==445)
ros1_tile = t_tile['TARGETID','TILEID'][is_ros1_tile]

## Print the list of TILEIDs from sv3-dark
print(np.unique(ros1_tile['TILEID']))

## Keep only zpix entries that are in Rosette 1, 'sv3', 'dark'
t_zpix = join(ros1_tile, t_zpix_sv3, join_type='inner', keys='TARGETID', metadata_conflicts='silent')

print(f'Resulting number of rows after joining with Rosette 1 dark tiles = {len(t_zpix)}')

## Print the list of TILEIDs from sv3-dark
print(np.unique(ros1_tile['TILEID']))

## Keep only zpix entries that are in Rosette 1, 'sv3', 'dark'
t_zpix = join(ros1_tile, t_zpix_sv3, join_type='inner', keys='TARGETID', metadata_conflicts='silent')

print(f'Resulting number of rows after joining with Rosette 1 dark tiles = {len(t_zpix)}')

## RA, Dec for the matched zpix table
ra_pix = t_zpix['TARGET_RA']
dec_pix = t_zpix['TARGET_DEC']

# List if unique HEALPIX values in the Rosette 1 area
hpixs = np.unique(t_zpix['HEALPIX'])

print(f"N(healpix) = {len(np.unique(t_zpix['HEALPIX']))}")
print(" ")
print(hpixs)

## Define a circle to label the petals around the field-of-view

# central poition (xc, yc) with a radius r
xc = 179.7+0.05  #offset to accommodate font size
yc = 0.-0.05     #offset to accommodate font size
r = 1.75         #deg
theta = np.arange(270,-90,-36)*np.pi/180.

# define arrays of 10 positions around the circle
xcirc = xc + r*np.cos(theta)
ycirc = yc + r*np.sin(theta)

## Pick a list of colors for the HEALPixels
colors = ['gold','thistle','tab:purple','tab:green','tab:orange','mediumblue','mediumturquoise','crimson','k','tab:pink','tab:grey',\
          'tab:cyan','tab:red','k','tab:blue','limegreen','tab:olive','slateblue','tab:brown']

## Two panel figure
fig, (ax1, ax2) = plt.subplots(1,2, figsize=(16.2,8), sharex=True, sharey=True)
plt.subplots_adjust(wspace=0.05)

# Loop over petal_loc from 0 to 9 (10 petals)
for petal in range(10):
    rpetal = selec_tile&(t_tile['PETAL_LOC']==petal)
    ax1.text(xcirc[petal],ycirc[petal],f'{petal}',color='dimgrey')
    ax1.scatter(ra_tile[rpetal],dec_tile[rpetal],s=3,label=f'Petal={petal}')

# Loop over the healpixels covering Rosette 1
for i,hpix in enumerate(hpixs):
    rhpx = (t_zpix['HEALPIX']==hpix)
    ax2.scatter(ra_pix[rhpx],dec_pix[rhpx],s=1,color=colors[i])

# Titles for each panel
ax1.set_title('Tile 30 color-coded per PETAL_LOC')
ax2.set_title('Rosette 1 color-coded per HEALPIX')

# Common plotting ranges for both panels
plt.xlim(181.6,177.6)
plt.ylim(-1.95,1.95)

# Axis labels
ax1.set_xlabel('RA [deg]')
ax1.set_ylabel('Dec [deg]')
ax2.set_xlabel('RA [deg]')
plt.show()