#!/usr/bin/env python

# Code to create catalogue files from
# HETDEX, Stripe 82, and COSMOS data.
# It converts fluxes to magnitudes,
# imputes missing data and
# creates new features out of
# original quantities (colours,
# ratios, flags, etc.).

import gc
from itertools import combinations
import numpy as np
import pandas as pd
from astropy.table import Table
from astropy import units as u
import global_variables as gv

# --------------------
# Generic helpers
# --------------------

def fix_dtypes_table(initial_tab):
    for col in initial_tab.colnames:
        if initial_tab[col].dtype.name == 'float64':
            initial_tab[col] = initial_tab[col].astype(np.float32)
        elif 'float' in initial_tab[col].dtype.name:
            initial_tab[col].fill_value = np.nan
        elif initial_tab[col].dtype.name == 'int64':
            initial_tab[col] = initial_tab[col].astype(np.int32)
        elif 'bytes' in initial_tab[col].dtype.name:
            initial_tab[col] = initial_tab[col].astype(str)
    # Special case
    if 'QPCT' in initial_tab.colnames:
        initial_tab['QPCT'] = initial_tab['QPCT'].astype(np.int32)
    if 'TYPE' in initial_tab.colnames:
        initial_tab['TYPE'] = initial_tab['TYPE'].filled('')
    return initial_tab

def create_band_count(mags_df, magnitude_cols, feat_name):
    na_bool = 1 - mags_df.loc[:, magnitude_cols].isna().astype(int)
    return pd.DataFrame({feat_name: na_bool.sum(axis=1).astype(np.int16)})

def create_MQC_filter(initial_tab, AGN_types):
    filters_array = [
        np.array(np.char.find(initial_tab['TYPE'].data, agn_type) != -1)
        for agn_type in AGN_types
    ]
    return np.bitwise_or.reduce(filters_array)

def create_AGN_gal_flags(initial_tab, imputed_df, AGN_types, mqc_version):
    filt_NLAGN = create_MQC_filter(initial_tab, AGN_types[mqc_version])
    is_str     = (np.array(initial_tab['spCl'] == 'STAR  ')).astype(np.int8)
    tmp_AGN_0  = np.array(initial_tab['Z'] > 0)
    tmp_AGN_1  = np.array((initial_tab['Z'] * 10) % 1 == 0)  # non spec z only
    tmp_AGN_2  = np.array(
        [(st.startswith('B') or st.startswith('R') or
          st.startswith('X') or st.startswith('2'))
         for st in initial_tab['TYPE']]
    )
    is_SDSS_QSO = (np.array(initial_tab['spCl'] == 'QSO   ')).astype(np.int8)
    is_AGN      = (tmp_AGN_0 & ~(tmp_AGN_1 & tmp_AGN_2) & filt_NLAGN).astype(np.int8)
    is_SDSS_gal = (np.array(initial_tab['spCl'] == 'GALAXY')).astype(np.int8)
    is_gal      = (is_SDSS_gal & ~is_AGN).astype(np.int8)

    imputed_df['is_str']      = is_str
    imputed_df['is_SDSS_QSO'] = is_SDSS_QSO
    imputed_df['is_AGN']      = is_AGN
    imputed_df['is_SDSS_gal'] = is_SDSS_gal
    imputed_df['is_gal']      = is_gal
    return imputed_df

def create_radio_detect(imputed_df, initial_tab, radio_cols):
    filters_array = [
        (np.array(initial_tab[radio_col] > 0) & np.isfinite(initial_tab[radio_col]))
        for radio_col in radio_cols
    ]
    or_in_arrays = np.bitwise_or.reduce(filters_array)
    imputed_df['radio_detect'] = or_in_arrays.astype(np.int8)
    for radio_col in radio_cols:
        name = radio_col.split('_')[-1] + '_detect'
        imputed_df[name] = (
            (np.array(initial_tab[radio_col] > 0) & np.isfinite(initial_tab[radio_col]))
        ).astype(np.int8)
    return imputed_df

def create_imputation_count(mags_df, magnitude_cols, magnitude_limits, feat_name):
    imputation_bool_df = pd.DataFrame()
    for mag in magnitude_cols:
        imputation_bool_df[mag] = np.array(
            mags_df.loc[:, mag] == magnitude_limits[mag]
        ).astype(np.int8)
    return pd.DataFrame({feat_name: imputation_bool_df.sum(axis=1).astype(np.int16)})

def create_colours(df, mag_list, mag_names_short_ver):
    for mags_pair in combinations(mag_list, 2):
        colour_name = mag_names_short_ver[mags_pair[0]] + '_' + mag_names_short_ver[mags_pair[1]]
        df[colour_name] = df[mags_pair[0]] - df[mags_pair[1]]
    return df

def save_hdf_and_parquet(df, h5_path, parquet_path, key='df', engine='fastparquet'):
    df.to_hdf(h5_path, key=key)
    df.to_parquet(parquet_path, engine=engine)

# --------------------
# Global config
# --------------------

run_HETDEX_flag = True
run_S82_flag    = False
run_COSMOS_flag = False

run_S82_full    = True

save_HETDEX_flag = False
save_S82_flag    = False
save_COSMOS_flag = False

all_vega_cols  = ['W1mproPM', 'W2mproPM', 'W1mag', 'W2mag', 'W3mag', 'W4mag',
                  'Jmag', 'Hmag', 'Kmag',
                  'e_W1mproPM', 'e_W2mproPM', 'e_W1mag', 'e_W2mag',
                  'e_W3mag', 'e_W4mag', 'e_Jmag', 'e_Hmag', 'e_Kmag']
vega_cols      = ['W1mproPM', 'W2mproPM', 'W1mag', 'W2mag', 'W3mag', 'W4mag',
                  'Jmag', 'Hmag', 'Kmag']
vega_shift     = {'W1mproPM': 2.699, 'W2mproPM': 3.339, 'W1mag': 2.699,
                  'W2mag': 3.339, 'W3mag': 5.174, 'W4mag': 6.620,
                  'Jmag': 0.910, 'Hmag': 1.390, 'Kmag': 1.850}

mag_cols_for_colours = ['gmag', 'rmag', 'imag', 'zmag', 'ymag',
                        'Jmag', 'Hmag', 'Kmag',
                        'W1mproPM', 'W2mproPM', 'W3mag', 'W4mag']
mag_names_short      = {'gmag': 'g', 'rmag': 'r', 'imag': 'i', 'zmag': 'z',
                        'ymag': 'y', 'Jmag': 'J', 'Hmag': 'H', 'Kmag': 'K',
                        'W1mproPM': 'W1', 'W2mproPM': 'W2',
                        'W3mag': 'W3', 'W4mag': 'W4'}

AGN_types_list = {'7_4d': ['Q', 'A', 'B', 'L', 'K', 'N', 'R', 'X', '2']}

# Your 5-sigma limits dicts (same as original)
mag_cols_lim_5sigma = {'W1mproPM': 20.13, 'W2mproPM': 19.81, 'Sint_LOFAR_AB': 17.52, 'Total_flux_VLASS_AB': 15.21,
                       'TotalFlux_LoLSS_AB': 12.91, 'Stotal_TGSS_AB': 11.18, 'Fint_VLAS82_AB': 17.86, 
                       'Flux_COSMOSVLA3_AB': 21.25, 'W1mag': 19.6, 'W2mag': 19.34, 'W3mag': 16.67, 'W4mag': 14.62, 
                       'gmag': 23.3, 'rmag': 23.2, 'imag': 23.1, 'zmag': 22.3, 'ymag': 21.4, 'FUVmag': 20.0, 
                       'NUVmag': 21.0, 'FEP': 57.9, 'Jmag': 17.45, 'Hmag': 17.24, 'Kmag': 16.59}   # keep your original values here
flx_cols_lim_5sigma = {col: (mag_cols_lim_5sigma[col + '_AB'] * u.mag(u.AB)).to(u.mJy).value for col in ['Sint_LOFAR', 'Total_flux_VLASS',
                                                                                                 'TotalFlux_LoLSS', 'Stotal_TGSS', 
                                                                                                 'Fint_VLAS82', 'Flux_COSMOSVLA3']}   # same as original

for key in mag_cols_lim_5sigma:
    mag_cols_lim_5sigma[key] = np.float32(mag_cols_lim_5sigma[key])
mag_cols_lim = {'5sigma': mag_cols_lim_5sigma}


def process_hetdex():
    if not run_HETDEX_flag:
        return

    print('-' * 40)
    print('Working with HETDEX data')
    print('Reading files')
    H = Table.read(gv.cat_path + gv.fits_HETDEX, format='fits')

    print('Fixing dtypes')
    H = fix_dtypes_table(H)
    for col in ['Speak_LOFAR', 'rms_LOFAR']:
        if col in H.colnames:
            try:
                H[col] = H[col].filled(np.nan)
            except Exception:
                pass

    id_cols = ['objID', 'RA_ICRS', 'DE_ICRS', 'Name',
               'RA_MILLI', 'DEC_MILLI', 'TYPE', 'Z', 'zsp', 'spCl']

    # Convert fluxes to magnitudes
    print('Convert fluxes to magnitudes')
    mJy_cols = [c for c in H.colnames
                if getattr(H[c], 'unit', None) == 'mJy'
                and not (c.startswith('e') or c.startswith('E'))]

    for col in mJy_cols:
        try:
            H[col] = H[col].filled(np.nan)
        except Exception:
            pass

    for col in mJy_cols:
        H[col + '_AB'] = H[col].to(u.mag(u.AB))

    # Transform Vega to AB
    print('Transforming Vega to AB')
    for col in vega_cols:
        if col in H.colnames:
            H[col] += vega_shift[col]

    for col in H.colnames:
        if getattr(H[col], 'unit', None) == u.mag:
            H[col].unit = u.mag(u.AB)
            H[col]      = u.Magnitude(H[col])

    magnitude_cols = [c for c in H.colnames
                      if getattr(H[c], 'unit', None) == u.mag(u.AB)
                      and not (c.startswith('e') or c.startswith('E')
                               or c.endswith('MILLI') or c.endswith('SDSS'))]

    radio_cols = ['Sint_LOFAR'] if 'Sint_LOFAR' in H.colnames else []

    # Build pandas frames once
    clean_df = H[id_cols].to_pandas()
    mags_df  = H[magnitude_cols].to_pandas()
    flxs_df  = H[radio_cols].to_pandas() if radio_cols else pd.DataFrame(index=mags_df.index)

    # Imputed container starts with radio + flags
    imputed_df = pd.DataFrame(index=mags_df.index)

    print('Creating flags for X-ray and radio detections')
    if radio_cols:
        imputed_df = create_radio_detect(imputed_df, H, radio_cols)

    # radio measurements (no copy)
    for col in radio_cols:
        imputed_df[col]        = H[col]
        ab_col = col + '_AB'
        if ab_col in H.colnames:
            imputed_df[ab_col] = H[ab_col]
    for opt_col in ['Speak_LOFAR', 'rms_LOFAR']:
        if opt_col in H.colnames:
            imputed_df[opt_col] = H[opt_col]

    print('Creating flags for AGN/Galaxy/Star classification')
    imputed_df = create_AGN_gal_flags(H, imputed_df, AGN_types_list, gv.mqc_version)

    # Drop H except for columns already converted
    del H
    gc.collect()

    # Remove columns with high nullity (use mags_df)
    print('Removing columns with high nullity')
    kept_mags = []
    for col in magnitude_cols:
        filt_temp = np.isfinite(mags_df[col].values)
        if np.sum(~filt_temp) > int(np.ceil(1.0 * len(mags_df[col]))):
            print(f'column: {col}\t-\t n_bad: {np.sum(~filt_temp)}\tREMOVED')
        else:
            kept_mags.append(col)
    magnitude_cols = kept_mags
    mags_df = mags_df[magnitude_cols]

    # Derived features
    print('Creating new features:')
    print('Creating counter of valid measurements')
    band_count_df = create_band_count(mags_df, mag_cols_for_colours, 'band_num')

    print('Imputing values')
    imputed_phot_df = mags_df.copy()
    non_imputed_df  = mags_df.copy()

    for col in magnitude_cols:
        lim = mag_cols_lim['5sigma'][col]
        tmp = imputed_phot_df[col].fillna(np.float32(lim))
        tmp = tmp.mask(tmp > lim, lim)
        imputed_phot_df[col] = tmp

    for col in radio_cols:
        flx_lim = flx_cols_lim_5sigma[col]
        tmp = flxs_df[col].fillna(np.float32(flx_lim))
        tmp = tmp.mask(tmp < flx_lim, flx_lim)
        imputed_df[col] = tmp
        non_imputed_df[col] = flxs_df[col]

    print('Creating colours (non-imputed & imputed)')
    non_imp_full = pd.concat([imputed_df, non_imputed_df], axis=1)
    imp_full     = pd.concat([imputed_df, imputed_phot_df], axis=1)

    non_imp_full = create_colours(non_imp_full, mag_cols_for_colours, mag_names_short)
    imp_full     = create_colours(imp_full,     mag_cols_for_colours, mag_names_short)

    print('Creating counter of imputed measurements')
    imputed_count_df = create_imputation_count(
        imputed_phot_df, mag_cols_for_colours, mag_cols_lim['5sigma'], 'num_imputed'
    )

    if save_HETDEX_flag:
        print('Joining and saving HETDEX tables')

        # NON-imputed
        cat_non_imp = pd.concat(
            [clean_df, band_count_df, non_imp_full],
            axis=1
        )
        save_hdf_and_parquet(
            cat_non_imp,
            gv.cat_path + gv.file_non_imp_HETDEX,
            gv.cat_path + gv.file_non_imp_HETDEX.replace('.h5', '.parquet')
        )
        save_hdf_and_parquet(
            cat_non_imp.loc[:, mag_cols_for_colours],
            gv.preds_path + 'HETDEX_mags_non_imputed.h5',
            gv.preds_path + 'HETDEX_mags_non_imputed.parquet'
        )
        del cat_non_imp
        gc.collect()

        # Imputed
        cat_imp = pd.concat(
            [clean_df, band_count_df, imputed_count_df, imp_full],
            axis=1
        )
        save_hdf_and_parquet(
            cat_imp,
            gv.cat_path + gv.file_HETDEX,
            gv.cat_path + gv.file_HETDEX.replace('.h5', '.parquet')
        )
        save_hdf_and_parquet(
            cat_imp.loc[:, mag_cols_for_colours],
            gv.preds_path + 'HETDEX_mags_imputed.h5',
            gv.preds_path + 'HETDEX_mags_imputed.parquet'
        )
        del cat_imp

    # clean up
    del clean_df, band_count_df, imputed_phot_df, non_imputed_df
    del mags_df, flxs_df, imputed_df, non_imp_full, imp_full, imputed_count_df
    gc.collect()

if __name__ == '__main__':
    if run_HETDEX_flag:
        process_hetdex()
    if run_S82_flag:
        # implement process_s82() similarly, then call it here
        pass
    if run_COSMOS_flag:
        # implement process_cosmos() similarly, then call it here
        pass
    print('EOF')
