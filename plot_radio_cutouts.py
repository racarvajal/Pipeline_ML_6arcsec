#!/usr/bin/env python

# Plot cutouts from
# predicted new radio AGNs
# in radio wavelenghts

import glob
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.visualization.wcsaxes import add_beam, add_scalebar
from astropy.cosmology import Planck18 as cosmo
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patheffects as mpe
from matplotlib import font_manager as fm
from tqdm import tqdm
import global_variables as gv

mpl.rcdefaults()
plt.rcParams['text.usetex'] = True

depth_S82      = (52 * u.uJy).to(u.Jy)
depth_HETDEX   = (71 * u.uJy).to(u.Jy)

pe1            = [mpe.withStroke(linewidth=1.0, foreground='black', alpha=1)]
pe4            = [mpe.withStroke(linewidth=2.1, foreground='white', alpha=1)]
pe5            = [mpe.withStroke(linewidth=0.5, foreground='white', alpha=0.5)]
pe6            = [mpe.withStroke(linewidth=2.1, foreground='black', alpha=1)]

parquet_HETDEX    = gv.preds_path + 'HETDEX_full_prediction.parquet'
parquet_S82       = gv.preds_path + 'S82_full_prediction.parquet'
candidates_HETDEX = pd.read_parquet(parquet_HETDEX, columns=['RA_ICRS', 'DE_ICRS', 'pred_Z_rAGN', 'ID'])
candidates_S82    = pd.read_parquet(parquet_S82, columns=['RA_ICRS', 'DE_ICRS', 'pred_Z_rAGN', 'ID'])

file_names = glob.glob('cutouts_fits/*/*')

for fits_file in tqdm(file_names, total=len(file_names)):
    id_name     = int(fits_file.split('_')[-1].split('.')[0])
    with fits.open(fits_file) as hdu_l:
        file_header = hdu_l[0].header
        file_header['TIMESYS'] = 'utc'
        file_image  = hdu_l[0].data
    file_wcs    = WCS(file_header).celestial
    # update header
    file_header.update(file_wcs.to_header())

    percents = np.percentile(file_image, q=[68.26, 99.9937])

    if 'S82' in fits_file:
        noise_sigma   = depth_S82
        candidates_df = candidates_S82
    if 'HETDEX' in fits_file:
        noise_sigma   = depth_HETDEX
        candidates_df = candidates_HETDEX
    cont_levels = [2, 3, 4, 5]  # [1, 2, 3, 4, 5]
    cont_widths = [1.45, 1.30, 1.15, 1.00]  # [1.60, 1.45, 1.30, 1.15, 1.00]
    cont_alpha  = 1.0

    # norm_img = mcolors.TwoSlopeNorm(vcenter=0.0, vmin=percents[0], vmax=percents[1])
    # norm_img = mcolors.Normalize(vmin=percents[0], vmax=percents[1])
    norm_img = mcolors.PowerNorm(vmin=percents[0], vmax=percents[1], gamma=0.25)
    # norm_img = mcolors.LogNorm(vmin=percents[0])
    # norm_img = mcolors.SymLogNorm(linthresh=0.03, linscale=1, base=10, vmin=-percents[1], vmax=percents[1])

    fig  = plt.figure(figsize=(6.5, 6.5))
    axs  = fig.add_subplot(111, projection=file_wcs, slices=('x', 'y'))
    axs.imshow(file_image, origin='lower', norm=norm_img, zorder=-1, cmap='plasma', interpolation='nearest', aspect='equal')

    snr_data       = file_image / noise_sigma.value
    snr_max        = np.nanmin([np.nanmax(snr_data), 20])
    snr_max_level  = np.floor(snr_max)
    snr_levels     = np.arange(2, snr_max_level + 1)
    snr_linewidths = np.logspace(np.log10(1.75), np.log10(0.25), len(snr_levels))
    contours   = axs.contour(snr_data, levels=snr_levels, colors='k',
                                alpha=cont_alpha, linewidths=snr_linewidths, zorder=1)
    contours.set(path_effects=pe5)
    base_xlim = axs.get_xlim()
    base_ylim = axs.get_ylim()

    centre_coord = file_wcs.pixel_to_world_values([file_header['NAXIS1'] / 2], [file_header['NAXIS2'] / 2])
    # axs.plot(centre_coord[0], centre_coord[1], marker='+', ms=180, color='white', transform=axs.get_transform('icrs'), alpha=0.75, mew=2)   # not needed with offset coordinates
    
    add_beam(axs, major=file_header['BMAJ'] * u.deg, minor=file_header['BMIN'] * u.deg, angle=file_header['BPA'] * u.deg, frame=False, facecolor='white', lw=2.5, edgecolor='k')
    
    # Add scale bar in physical units
    length_bar_arcsec = 15 * u.arcsec
    distance_ang_z    = cosmo.angular_diameter_distance(candidates_df.loc[id_name, 'pred_Z_rAGN'])
    length_bar_kpc    = (length_bar_arcsec * distance_ang_z).to(u.kpc, u.dimensionless_angles())
    label_scale_bar_arcs  = rf'$\mathbf{{{length_bar_arcsec.value:.0f} ~ arcsec}}$'
    label_scale_bar_kpc   = rf'$\mathbf{{{length_bar_kpc.value:.1f} ~ kpc}}$'
    label_scale_bar_joint = label_scale_bar_arcs + '\n' + label_scale_bar_kpc
    bar_size_vertical = 0.75 if 'HETDEX' in fits_file else 1.75
    tick_prop = fm.FontProperties(size=22)
    # tick_prop = fm.FontProperties('outline')
    add_scalebar(axs, length_bar_arcsec, label=label_scale_bar_joint, color='white', fontproperties=tick_prop, size_vertical=bar_size_vertical, path_effects=pe1, fill_bar=True, label_top=True)

    axs.annotate(text=f'$\mathbf{{ID\!: {id_name:>09,d}}}$'.replace(',', '\\,'), xy=(0.016, 0.89),
                     xycoords='axes fraction', fontsize=30, ha='left', va='bottom', color='whitesmoke', bbox=dict(boxstyle='round, pad=0.1',
                     fc='gray', ec='k', lw=2), fontweight='bold')

    axs.annotate(text=f'$\mathbf{{Pred. ~ z\!: {candidates_df.loc[id_name, "pred_Z_rAGN"]:.2f}}}$',
                     xy=(0.016, 0.79), xycoords='axes fraction', fontsize=30, ha='left', va='bottom',
                     color='whitesmoke',
                     bbox=dict(boxstyle='round, pad=0.1', fc='gray', ec='k', lw=2), fontweight='bold')

    # Use offset coordinates
    centre_coord = SkyCoord(*centre_coord, unit=u.deg)
    ra           = axs.coords['ra']
    ra.set_auto_axislabel(False)
    dec          = axs.coords['dec']
    ra.set_coord_type('longitude', 180 * u.deg)
    overlay      = axs.get_coords_overlay(centre_coord.skyoffset_frame())
    ra.set_ticklabel_visible(False)
    dec.set_ticklabel_visible(False)
    ra.set_ticks_visible(False)
    dec.set_ticks_visible(False)
    lon          = overlay['lon']
    lon.set_coord_type('longitude', 180 * u.deg)
    lon.set_format_unit(u.arcsec)
    # lon.set_ticklabel(rotation=45)
    lon.set_ticks_position('b')
    lon.set_ticklabel_position('b')
    lon.set_axislabel_position('b')
    lon.tick_params(which='major', labelsize=28, direction='in')
    lon.tick_params(which='major', length=8, width=1.5, direction='in')
    lat          = overlay['lat']
    lat.set_format_unit(u.arcsec)
    lat.set_ticklabel()
    lat.set_ticks_position('l')
    lat.set_ticklabel_position('l')
    lat.set_axislabel_position('l')
    lat.tick_params(which='major', labelsize=28, direction='in')
    lat.tick_params(which='major', length=8, width=1.5, direction='in')

    # Format axes
    # axs.set_xlabel('$\mathrm{R.A.}$', fontsize=22)
    # axs.set_ylabel('$\mathrm{Declination}$', fontsize=22)
    lon.set_axislabel('$\Delta \mathrm{R.A.}$', fontsize=24)
    lat.set_axislabel('$\Delta \mathrm{Dec.}$', fontsize=24)
    # axs.tick_params(axis='both', which='major', labelsize=14, direction='in')
    # axs.tick_params(which='major', length=8, width=1.5, direction='in')
    plt.setp(axs.spines.values(), linewidth=3.5)
    plt.setp(axs.spines.values(), linewidth=3.5)
    plt.setp(axs.get_yticklabels(), visible=False)
    overlay.grid(True, zorder=1, lw=0.50)
    plt.tight_layout()
    if 'S82' in fits_file:
        folder_save = 'S82'
    if 'HETDEX' in fits_file:
        folder_save = 'HETDEX'
    plt.savefig(f'cutouts_plots/{folder_save}/cutout_{folder_save}_{id_name:08d}.pdf', bbox_inches='tight')
    # plt.show()
    plt.clf()
    plt.close()
