#!/usr/bin/env python

# Plot cutouts from
# predicted new radio AGNs
# in radio wavelenghts

import glob
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.visualization.wcsaxes import add_beam, add_scalebar, SphericalCircle
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import font_manager as fm

mpl.rcdefaults()
plt.rcParams['text.usetex'] = True


file_names = glob.glob('cutouts_fits/*/*')

for fits_file in file_names:
    id_name     = int(fits_file.split('_')[-1].split('.')[0])
    hdu_0       = fits.open(fits_file)[0]
    file_header = hdu_0.header
    file_image  = hdu_0.data
    file_wcs    = WCS(file_header).celestial

    percents = np.percentile(file_image, q=[68.26, 99.9937])

    # norm_img = mcolors.TwoSlopeNorm(vcenter=0.0, vmin=percents[0], vmax=percents[1])
    # norm_img = mcolors.Normalize(vmin=percents[0], vmax=percents[1])
    norm_img = mcolors.PowerNorm(vmin=percents[0], vmax=percents[1], gamma=0.25)
    # norm_img = mcolors.LogNorm(vmin=percents[0])
    # norm_img = mcolors.SymLogNorm(linthresh=0.03, linscale=1, base=10, vmin=-percents[1], vmax=percents[1])

    fig  = plt.figure(figsize=(6.5, 6.5))
    axs  = fig.add_subplot(111, projection=file_wcs)
    plt.imshow(file_image, origin='lower', norm=norm_img, zorder=-1, cmap='plasma')

    centre_coord = file_wcs.pixel_to_world_values([file_header['NAXIS1'] / 2], [file_header['NAXIS2'] / 2])
    axs.plot(centre_coord[0], centre_coord[1], marker='+', ms=180, color='red', transform=axs.get_transform('icrs'), alpha=0.5)
    
    add_beam(axs, major=file_header['BMAJ'] * u.deg, minor=file_header['BMIN'] * u.deg, angle=file_header['BPA'] * u.deg, frame=False, facecolor='white', lw=2.5, edgecolor='k')
    
    bar_size_vertical = 0.75 if 'HETDEX' in fits_file else 1.75
    tick_prop = fm.FontProperties(size=20)
    add_scalebar(axs, 15 * u.arcsec, label='$15\,\mathrm{arcsec}$', color='white', fontproperties=tick_prop, size_vertical=bar_size_vertical)

    axs.annotate(text=f'$\mathrm{{ID}}:$ ${id_name:09,d}$'.replace(',', '$\,$'), xy=(0.016, 0.89),
                     xycoords='axes fraction', fontsize=30,ha='left', va='bottom', color='whitesmoke', bbox=dict(boxstyle='round,pad=0.1',
                      fc='gray', ec='k', lw=2))


    axs.set_xlabel('$\mathrm{R.A.}$', fontsize=22)
    axs.set_ylabel('$\mathrm{Declination}$', fontsize=22)
    axs.tick_params(axis='both', which='major', labelsize=14, direction='in')
    axs.tick_params(which='major', length=8, width=1.5, direction='in')
    plt.setp(axs.spines.values(), linewidth=3.5)
    plt.setp(axs.spines.values(), linewidth=3.5)
    plt.setp(axs.get_yticklabels(), visible=False)
    plt.tight_layout()
    if 'S82' in fits_file:
        folder_save = 'S82'
    if 'HETDEX' in fits_file:
        folder_save = 'HETDEX'
    plt.savefig(f'cutouts_plots/{folder_save}/cutout_{folder_save}_{id_name:08d}.pdf', bbox_inches='tight')
    # plt.show()
    plt.clf()
