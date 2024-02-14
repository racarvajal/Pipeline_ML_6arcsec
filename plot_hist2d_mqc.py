#!/usr/bin/env python

# Plot 2D histogram
# of sources in the
# Million Quasar Catalog (MQC)

from os import path, rename
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
from matplotlib import ticker
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patheffects as mpe
import matplotlib.patches as mpatches
from matplotlib.ticker import ScalarFormatter
from astropy.visualization import LogStretch, PowerStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
import colorcet as cc
import cmasher as cmr
import skyproj
import global_variables as gv

# Path effects for labels and plots.
pe1            = [mpe.Stroke(linewidth=5.0, foreground='black'),
                  mpe.Stroke(foreground='white', alpha=1),
                  mpe.Normal()]
pe2            = [mpe.Stroke(linewidth=2.0, foreground='white'),
                  mpe.Normal()]
pe3            = [mpe.Stroke(linewidth=3.1, foreground='black'),
                  mpe.Normal()]

mpl.rcdefaults()
plt.rcParams['text.usetex'] = True
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16
plt.rcParams['path.effects'] = pe2


# temporal path for mqc data (mac)
mqc_version = '7_4d' # '8', '7_4d', '7_4c', etc.
tmp_path  = '/Users/rcarvajal/Documents/Astro/Science/Code/Catalogs/'
file_name = f'milliquas_v{mqc_version}.fits'

with fits.open(tmp_path + file_name) as hdul:
    data = hdul[1].data

data_ra  = data['RA']
data_dec = data['DEC']

ax_props = {'axes.linewidth':3.5, 'axes.labelsize':20}
fig  = plt.figure(figsize=(9,5), constrained_layout=True)
ax1  = fig.add_subplot(111)
sp = skyproj.HammerSkyproj(ax=ax1, rcparams=ax_props, lon_0=0)

norm_sources  = ImageNormalize(stretch=PowerStretch(0.30))
sp.draw_hpxbin(data_ra, data_dec, cmap='cet_bmy', nside=128, xsize=400)
sp.draw_colorbar(pad=0.05,
                label=r'$\mathrm{Number ~ of ~ sources ~ per ~ bin}$',
                fontsize=16, fraction=0.025)
sp.set_xlabel(r'$\mathrm{Right ~ Ascension}$')
sp.set_ylabel(r'$\mathrm{Declination}$')
sp.ax.tick_params(axis='both', which='major', labelsize=30)
# plt.tight_layout()
# plt.savefig(gv.plots_path + f'MQC_v{mqc_version}_source_density_map.pdf', bbox_inches='tight')
plt.show()
 