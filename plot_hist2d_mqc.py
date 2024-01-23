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
import skyproj
import global_variables as gv

mpl.rcdefaults()
plt.rcParams['text.usetex'] = True

# temporal path for mqc data (mac)
tmp_path  = '/Users/rcarvajal/Documents/Astro/Science/Code/Catalogs/'
file_name = 'milliquas_v8.fits'

with fits.open(tmp_path + file_name) as hdul:
    data = hdul[1].data

data_ra  = data['RA']
data_dec = data['DEC']

ax_props = {'xtick.labelsize':14, 'ytick.labelsize':14}
fig  = plt.figure(figsize=(9,5))
ax1  = fig.add_subplot(111)
sp = skyproj.HammerSkyproj(ax=ax1, rcparams=ax_props, lon_0=0)

norm_sources  = ImageNormalize(stretch=PowerStretch(0.30))
sp.draw_hpxbin(data_ra, data_dec, cmap='cet_gray_r', nside=128, xsize=400)
sp.draw_colorbar(pad=0.05)

# plt.tight_layout()
# plt.savefig(gv.plots_path + 'MQC_source_density_map.pdf')
plt.show()
 