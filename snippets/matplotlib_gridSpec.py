import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename

filename = get_pkg_data_filename('galactic_center/gc_msx_e.fits')

hdu = fits.open(filename)[0]
wcs = WCS(hdu.header)
fig = plt.figure()
ax = fig.add_subplot(projection=None)

for i in range(10):
    ax.clear()
    ax.imshow(hdu.data, vmin=-2.e-5, vmax=2.e-4, origin='lower')
    ax.update({'xlabel': 'Galactic Longitude', 'ylabel': 'Galactic Latitude'})
    fig.canvas.draw()
    plt.pause(10)
plt.close(fig)




# import matplotlib.pyplot as plt
# import numpy as np
#
# from astropy.wcs import WCS
# from astropy.io import fits
# from astropy.utils.data import get_pkg_data_filename
#
#
# class ClickToDrawPoints(object):
#
#     def __init__(self, ax):
#         self.ax = ax
#         self.fig = ax.figure
#         self.xy = []
#         self.points = ax.scatter([], [], s=200, color='red', picker=20)
#         self.fig.canvas.mpl_connect('button_press_event', self.on_click)
#
#     def on_click(self, event):
#         if event.inaxes is None:
#             return
#         self.xy.append([event.xdata, event.ydata])
#         self.points.set_offsets(self.xy)
#         self.ax.draw_artist(self.points)
#         self.fig.canvas.blit(self.ax.bbox)
#
#     def show(self):
#         plt.show()
#
# def main():
#
#     filename = get_pkg_data_filename('galactic_center/gc_msx_e.fits')
#
#     hdu = fits.open(filename)[0]
#     wcs = WCS(hdu.header)
#
#     fig = plt.figure()
#     ax = fig.add_subplot(projection=wcs)
#     ax.imshow(hdu.data, vmin=-2.e-5, vmax=2.e-4, origin='lower')
#     ax.update({'xlabel': 'Galactic Longitude', 'ylabel': 'Galactic Latitude'})
#
#     ClickToDrawPoints(ax).show()
#
# main()



# plt.subplot(projection=wcs)
# plt.imshow(hdu.data, vmin=-2.e-5, vmax=2.e-4, origin='lower')
# plt.grid(color='white', ls='solid')
# plt.xlabel('Galactic Longitude')
# plt.ylabel('Galactic Latitude')
#
# plt.show()

# from matplotlib import pyplot as plt, gridspec
#
# fig = plt.figure()
# gs = gridspec.GridSpec(nrows=1,
#                        ncols=2,
#                        figure=fig,
#                        width_ratios=[1, 2],
#                        height_ratios=[1])
# ax1 = fig.add_subplot(gs[0])
# ax1.text(0.5, 0.5, 'ax1: gs[0]', fontsize=12, fontweight="bold", va="center", ha="center")
# ax2 = fig.add_subplot(gs[1])
# ax2.text(0.5, 0.5, 'ax2: gs[1]', fontsize=12, fontweight="bold", va="center", ha="center")
# plt.show()

# fig = plt.figure(figsize=(7,7))
# gs = gridspec.GridSpec(nrows=3,
#                        ncols=3,
#                        figure=fig,
#                        width_ratios= [1, 1, 1],
#                        height_ratios=[1, 1, 1],
#                        wspace=0.3,
#                        hspace=0.3)
# ax1 = fig.add_subplot(gs[0, 0])
# ax1.text(0.5, 0.5, 'ax1: gs[0, 0]', fontsize=12, fontweight="bold", va="center", ha="center")  # adding text to ax1
# ax2 = fig.add_subplot(gs[0, 1:3])
# ax2.text(0.5, 0.5, 'ax2: gs[0, 1:3]', fontsize=12, fontweight="bold", va="center", ha="center")
# ax3 = fig.add_subplot(gs[1:3, 0:2])
# ax3.text(0.5, 0.5, 'ax3: gs[1:3, 0:2]', fontsize=12, fontweight="bold", va="center", ha="center")
# ax4 = fig.add_subplot(gs[1:3, 2])
# ax4.text(0.5, 0.5, 'ax4: gs[1:3, 2]', fontsize=12, fontweight="bold", va="center", ha="center")
# plt.show()