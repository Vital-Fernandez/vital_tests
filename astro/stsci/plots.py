import lime
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
from regions import CircleSkyRegion, PointSkyRegion
from astropy.visualization import ZScaleInterval, ImageNormalize
import astropy.units as u
import matplotlib.patheffects as pe
import matplotlib.gridspec as gridspec

from matplotlib import pyplot as plt, colors, cm, rc_context
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
# Optional: pretty sexagesimal output
from astropy.coordinates import SkyCoord
import astropy.units as u

def connect_wcs_click(ax, wcs):
    """
    Attach a mouse-click callback to a WCSAxes plot that prints
    pixel and sky coordinates (RA, Dec) in the terminal.
    """

    def onclick(event):
        # Only react to clicks inside this axes
        if event.inaxes != ax:
            return

        if event.button != 3:
            return

        x = event.xdata
        y = event.ydata

        if x is None or y is None:
            return

        # Convert pixel -> world (RA, Dec)
        ra, dec = wcs.pixel_to_world_values(x, y)

        # print(f"Pixel: x={x:.2f}, y={y:.2f}")
        # print(f"RA:    {ra:.8f} deg")
        # print(f"DEC:   {dec:.8f} deg")


        c = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')
        print(f'Coordinate: ["{c.ra.to_string(unit=u.hourangle, precision=3, pad=True)}", "{c.dec.to_string(precision=3, alwayssign=True, pad=True)}"]')

        # print("Sexagesimal:")
        # print("  RA :", c.ra.to_string(unit=u.hourangle, sep=':'))
        # print("  DEC:", c.dec.to_string(sep=':'))
        # print("-" * 40)

    # Connect the callback
    cid = ax.figure.canvas.mpl_connect('button_press_event', onclick)
    return cid


def cos_image_plotter(imdata, wcs, object_label, aper_dict=None, disp_axis=False, orientation_axis=False, dec_sum=None,
                      comps_dict=None, instr_im=None, mask_points=None, overlay=None, title=None):

    with (rc_context(lime.theme.fig_defaults({"figure.dpi": 300, "figure.figsize": [4, 4], "legend.fontsize": 3, "legend.fontsize" : 7}))):

        # Generate figure
        fig = plt.figure()

        # Generate axis
        if dec_sum is None:
            ax = [fig.add_subplot(projection=wcs, slices=('x', 'y'))]
        else:
            gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.05, figure=fig)
            ax = [None, None]
            ax[0] = fig.add_subplot(gs[0], projection=wcs, slices=('x', 'y'))
            ax[1] = fig.add_subplot(gs[1], sharex=ax[0])

        # Image
        if instr_im == 'COS':
            ax[0].imshow(imdata, interpolation="none")
        else:
            norm = ImageNormalize(imdata, interval=ZScaleInterval())
            ax[0].imshow(imdata, interpolation="none", norm=norm)

        # Pixel sum
        if dec_sum is not None:
            x_dec_sum = np.arange(dec_sum.size)
            ax[1].step(x_dec_sum, dec_sum, color="cyan")
            ax[1].set_xlim(ax[0].get_xlim())

            if comps_dict is not None:
                x_comps = np.linspace(x_dec_sum.min(), x_dec_sum.max(), comps_dict['c_'].size)
                ax[1].plot(x_comps, comps_dict['c_'], "--", lw=1.2)
                ax[1].plot(x_comps, comps_dict['g1_']+comps_dict['c_'], "--", lw=1.2, label=f'FWHM = {comps_dict['g1_fwhm']:.2f}')
                ax[1].legend()

        # Apertures
        if aper_dict is not None:

            cmap = plt.get_cmap("plasma_r")
            color_list = cmap(np.linspace(0, 1, len(aper_dict)))

            for k, items in enumerate(aper_dict.items()):
                aper_obj, aper_data = items
                center = SkyCoord(ra=aper_data['ra'], dec=aper_data['dec'], frame="icrs")
                sky_region = CircleSkyRegion(center=center, radius=1.25 * u.arcsec)  # , radius=aper_data['radius'])
                pixel_region = sky_region.to_pixel(wcs)
                pixel_region.plot(ax=ax[0], linewidth=0.75, label=aper_obj, linestyle=':', color=color_list[k])

                # Dispersion direction
                if disp_axis:
                    pa_disp = (aper_data['pa'] - 90 * u.deg) % (360 * u.deg)
                    arrow_len = 5.0 * u.arcsec  # adjust length
                    end_disp = center.directional_offset_by(pa_disp, arrow_len.to(u.deg))
                    x0, y0 = wcs.world_to_pixel(center)
                    x1, y1 = wcs.world_to_pixel(end_disp)
                    ax[0].annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", color="orange", lw=0.5))
                    ax[0].text(x1, y1, "disp", color="orange", fontsize=4, ha="left", va="bottom")

                # ---------- NE arrows with origin in *axes fractions* ----------
                if orientation_axis:
                    if k == 0:
                        ne_origin_frac = (0.90, 0.80)
                        px0 = ne_origin_frac[0] * (imdata.shape[1] - 1)
                        py0 = ne_origin_frac[1] * (imdata.shape[0] - 1)
                        p0_world = wcs.pixel_to_world(px0, py0)

                        # North (PA=0°), East (PA=90°) in the "east of north" convention
                        pN_world = p0_world.directional_offset_by(0 * u.deg, 0.001 * u.deg)
                        pE_world = p0_world.directional_offset_by(90 * u.deg, 0.001 * u.deg)

                        xN, yN = wcs.world_to_pixel(pN_world)
                        xE, yE = wcs.world_to_pixel(pE_world)

                        # draw arrows (thick + outlined for visibility)
                        outline = [pe.Stroke(linewidth=4.5, foreground="black"), pe.Normal()]

                        ax[0].annotate("", xy=(xN, yN), xytext=(px0, py0),
                                    arrowprops=dict(arrowstyle="-|>", color="white", lw=1.5),
                                    annotation_clip=False, zorder=10)
                        ax[0].annotate("", xy=(xE, yE), xytext=(px0, py0),
                                    arrowprops=dict(arrowstyle="-|>", color="white", lw=1.5),
                                    annotation_clip=False, zorder=10)

                        # labels near tips
                        txtN = ax[0].text(xN, yN, "N", color="white", fontsize=8, ha="center", va="bottom", zorder=11)
                        txtE = ax[0].text(xE, yE, "E", color="white", fontsize=8, ha="left", va="center", zorder=11)


        # Plot mirror mask
        if mask_points is not None:
            c1 = SkyCoord(mask_points[0][0], mask_points[0][1])
            c2 = SkyCoord(mask_points[1][0], mask_points[1][1])
            ax[0].plot([c1.ra.deg, c2.ra.deg], [c1.dec.deg, c2.dec.deg], linestyle='-', transform=ax[0].get_transform('icrs'))


        if overlay is not None:
            ax[0].imshow(overlay, origin="lower")

            # if aper_dict is not None:
            #
            #     # Get the size of the image
            #     ny, nx = imdata.shape
            #     yy, xx = np.mgrid[:ny, :nx]
            #
            #     # Use the favored apperture or use the first
            #     if apert_fav is None:
            #         apert_fav = list(aper_dict.keys())[0]
            #     aper_data = aper_dict[apert_fav]
            #
            #     # Generate the circular mask
            #     sky_region = CircleSkyRegion(center=SkyCoord(ra=aper_data['ra'], dec=aper_data['dec'], frame="icrs"),
            #                                  radius=1.25 * u.arcsec)
            #     circle_pix = sky_region.to_pixel(wcs=wcs)  # -> CirclePixelRegion
            #     circle_mask = circle_pix.to_mask(mode="center").to_image((ny, nx)).astype(bool)
            #
            #     # Generate the line mask
            #     x1, y1 = wcs.world_to_pixel(c1)
            #     x2, y2 = wcs.world_to_pixel(c2)
            #     sign = (x2 - x1) * (yy - y1) - (y2 - y1) * (xx - x1)
            #
            #     # Combine them and create an overlay
            #     side_mask = (sign < 0) & circle_mask
            #     overlay = np.zeros((*side_mask.shape, 4))
            #     overlay[..., 3] = 0.0  # fully transparent everywhere
            #     overlay[side_mask] = [0, 0, 0, 1]  # black & opaque where mask is True
            #     ax[0].imshow(overlay, origin="lower")

        # Labels
        if title is None:
            title = f'{object_label} ' + ('' if aper_dict is None or len(aper_dict) == 0 else f"({len(aper_dict)} data sets)")

        ax[0].update({'ylabel': "Dec", 'xlabel': "RA", 'title': title})

        # Grids
        ax[0].grid(True, linestyle='--', linewidth=0.25, alpha=1, color='orange')

        # Remove repeated entries
        handles, labels = ax[0].get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax[0].legend(by_label.values(), by_label.keys(), loc='lower center')

        # Display figure
        plt.tight_layout()

        # Connection to print the coordinates
        cid = connect_wcs_click(ax[0], wcs)

        if dec_sum is not None:
            plt.setp(ax[0].get_xticklabels(), visible=False)

        plt.show()

    return