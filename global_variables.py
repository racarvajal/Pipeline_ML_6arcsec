#!/usr/bin/env python

# File with most used
# variables in this project.
# File paths, file names, etc.

# Paths
cat_path           = '../../Catalogs/'  # relative path to the same directory
plots_path         = 'plots/'
models_path        = 'models/'
preds_path         = 'pred_rAGN/'
moc_path           = 'moc_files/'
tmp_shap_path      = 'tmp_shap/'
indices_path       = 'subsets_indices/'

# Boolean flags
use_5sigma         = True  # use files with 5-sigma magnitude imputation

# File versions
mqc_version        = '7_4d'  # '7_2'

# Catalogue file names
fits_HETDEX         = 'CatWISE2020_LoTSS_area_large_LoTSS_6arcsec_PS1_ALLWISE_2MASS_MQC_74D_SDSS_DR16_1_1arcsec.fits'  # 15136878 objects (1.5e7)
fits_S82            = 'CatWISE2020_VLAS82_area_large_VLAS82_6arcsec_PS1_ALLWISE_2MASS_MQC_74D_SDSS_DR16_1_1arcsec.fits'  # 3590306 objects
fits_COSMOS         = ''
fits_S82_Ananna_17  = 'CatWISE2020_VLAS82_area_large_VLAS82_6arcsec_PS1_ALLWISE_2MASS_MQC_74D_SDSS_DR16_1_1arcsec_Ananna_17.fits'  # 804 objects
# Non imputed h5 files
file_non_imp_HETDEX = fits_HETDEX.replace('.fits', '_non_imp.h5')
file_non_imp_S82    = fits_S82.replace('.fits', '_non_imp.h5')
file_non_imp_COSMOS = fits_COSMOS.replace('.fits', '_non_imp.h5')
# Imputed h5 files
file_HETDEX         = fits_HETDEX.replace('.fits', '_imp.h5')  # 15136878 objects (1.5e7)
file_S82            = fits_S82.replace('.fits', '_imp.h5')
file_COSMOS         = fits_COSMOS.replace('.fits', '_imp.h5')
file_S82_Ananna_17  = fits_S82_Ananna_17.replace('.fits', '_imp.h5')

# Fields properties
# Areas (deg2)
area_HETDEX         = 424
area_S82            = 92
area_COSMOS         = 4  # Not real value. Placeholder

# Model names with train, test, calibration, and validation sub-sets
star_model         = 'classification_star_no_star_mar_15_2023'
AGN_gal_model      = 'classification_AGN_galaxy_mar_16_2023'
radio_model        = 'classification_LOFAR_detect_mar_17_2023'
full_z_model       = 'regression_z_mar_18_2023'
high_z_model       = 'regression_high_z_mar_19_2023'

# Models for galaxies (test)
radio_galaxies_model   = 'classification_LOFAR_detect_galaxies_mar_25_2023'
z_radio_galaxies_model = 'regression_z_radio_galaxies_mar_26_2023'

# Calibrated models
cal_str_model      = 'cal_' + star_model    + '.joblib'
cal_AGN_gal_model  = 'cal_' + AGN_gal_model + '.joblib'
cal_radio_model    = 'cal_' + radio_model   + '.joblib'

cal_radio_gals_model = 'cal_' + radio_galaxies_model + '.joblib'

# Seeds
seed               = 42

# Thresholds
# Beta for beta-scores
beta_F             = 1.1  # beta positive real value
# Naive values
naive_star_thresh  = 0.5
naive_AGN_thresh   = 0.5
naive_radio_thresh = 0.5

naive_radio_gals_thresh = 0.5
# Values obtained with train, test, calibration, and validation sub-sets
# PR-optimised models (with train+test sub-set)
star_thresh        = 0.1873511777 # old value
AGN_thresh         = 0.5000115951
radio_thresh       = 0.9911368526

radio_gals_thresh  = 0.5910151612

# Calibrated and PR-optimised models (with calibration sub-set)
cal_str_thresh     = 0.6007345636412931 # old value
cal_AGN_thresh     = 0.34895396724527294
cal_radio_thresh   = 0.24489520601404396

cal_radio_gals_thresh = 0.19957242325544386
# High redshift limit
high_z_limit       = 2.0  # 3.6

# Colours and colormaps
# cmap_bands         = 'cmr.pride'
# cmap_shap          = 'cmr.guppy'  # cmr.pride, cet_CET_R3 cmr.wildfire cmr.guppy
# cmap_conf_matr     = 'cet_dimgray_r'
# cmap_z_plots       = 'cet_linear_kryw_5_100_c64_r'

cmap_bands         = 'cmr.rainforest'
cmap_shap          = 'cmr.guppy_r'  # cmr.pride, cet_CET_R3 cmr.wildfire cmr.guppy cmr.ember
cmap_conf_matr     = 'cmr.neutral_r'
cmap_z_plots       = 'cmr.fall_r'
cmap_dens_plots    = 'cmr.neutral_r'
cmap_hists         = 'cmr.fusion'