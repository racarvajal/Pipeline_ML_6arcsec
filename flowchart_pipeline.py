#!/usr/bin/env python

import schemdraw
from schemdraw import flow
import schemdraw.elements as elm
import matplotlib as mpl
import matplotlib.pyplot as plt
import global_variables as gv

mpl.rcdefaults()
plt.rcParams['text.usetex'] = True

save_plot_flag  = False
draw_new_branch = True

with schemdraw.Drawing(show=False, margin=-0.66) as d:
    init           = flow.Start(h=1.5, w=2.75).at((0, 0)).label(r'$\mathrm{SOURCE}$' '\n' r'$\mathrm{FROM}$' '\n' r'$\mathrm{CATALOGUE}$')
    AGN_gal_model  = flow.Decision(h=2.25, w=3.75, W=r'$\mathrm{Predicted}$' '\n' r'$\mathrm{as ~ AGN}$', E=r'$\mathrm{Predicted}$' '\n' r'$\mathrm{as ~ SFG}$').at((0, -2.5)).label(r'$\mathrm{AGN/SFG}$' '\n' r'$\mathrm{CLASSIFICATION}$' '\n' r'$\mathrm{MODEL}$')
    elm.Wire('-', arrow='->').at(init.S).to(AGN_gal_model.N)
    rAGN_model     = flow.Decision(h=2.25, w=4, S=r'$\mathrm{Predicted}$' '\n' r'$\mathrm{as ~ radio}$', E=r'$\mathrm{Predicted}$' '\n' r'$\mathrm{as ~ no ~ radio}$').at((-5, -4.5)).label(r'$\mathrm{AGN ~ RADIO}$' '\n' r'$\mathrm{DETECTION}$' '\n' r'$\mathrm{MODEL}$')
    rGal_model     = flow.Decision(h=2.25, w=4, S=r'$\mathrm{Predicted}$' '\n' r'$\mathrm{as ~ radio}$', W=r'$\mathrm{Predicted}$' '\n' r'$\mathrm{as ~ no ~ radio}$').at((5, -4.5)).label(r'$\mathrm{SFG ~ RADIO}$' '\n' r'$\mathrm{DETECTION}$' '\n' r'$\mathrm{MODEL}$')
    elm.Wire('-|', arrow='->').at(AGN_gal_model.W).to(rAGN_model.N)
    elm.Wire('-|', arrow='->').at(AGN_gal_model.E).to(rGal_model.N)
    discarded      = flow.StateEnd(r=1.30).at((0, -6.5)).label(r'$\mathrm{NON{-}RADIO}$' '\n' r'$\mathrm{SOURCE}$', ofst=(0.0, -0.06))
    z_rAGN_model   = flow.Box(h=1.5, w=2.5).at((-5, -7.75)).label(r'$\mathrm{rAGN ~ Z}$' '\n' r'$\mathrm{PREDICTION}$' '\n' r'$\mathrm{MODEL}$')
    z_rGal_model   = flow.Box(h=1.5, w=2.5).at((5, -7.75)).label(r'$\mathrm{rSFG ~ Z}$' '\n' r'$\mathrm{PREDICTION}$' '\n' r'$\mathrm{MODEL}$')
    elm.Wire('-', arrow='->').at(rAGN_model.S).to(z_rAGN_model.N)
    elm.Wire('-', arrow='->').at(rGal_model.S).to(z_rGal_model.N)
    elm.Wire('-|', arrow='->').at(rAGN_model.E).to(discarded.N)
    elm.Wire('-|', arrow='->').at(rGal_model.W).to(discarded.N)
    final_rAGN     = flow.Box(h=1.5, w=2.5).at((-5, -10.25)).label(r'$\mathrm{PREDICTED}$' '\n' r'$\mathrm{rAGN}$' '\n' r'$\mathrm{W/REDSHIFT}$')
    final_rGal     = flow.Box(h=1.5, w=2.5).at((5, -10.25)).label(r'$\mathrm{PREDICTED}$' '\n' r'$\mathrm{rSFG}$' '\n' r'$\mathrm{W/REDSHIFT}$')
    if draw_new_branch:
        new_branch_box = elm.EncircleBox([rGal_model, z_rGal_model, final_rGal], padx=-3.48, pady=.3).linestyle('--').linewidth(3).label(r'$\mathrm{THIS ~ WORK}$', ofst=(1.1, 0.15), rotate=0).at((5, -8.2))
        rGal_model.style(color='black', ls='-')
        new_branch_box.color('#3498DB')#.fill('#EBF5FB')
    elm.Wire('-', arrow='->').at(z_rAGN_model.S).to(final_rAGN.N)
    elm.Wire('-', arrow='->').at(z_rGal_model.S).to(final_rGal.N)
    final_compile  = flow.Box(h=1.5, w=2.5).at((0, -11.75)).label(r'$\mathrm{COMPILE}$' '\n' r'$\mathrm{PREDICTIONS}$')
    elm.Wire('|-', arrow='->').at(final_rAGN.S).to(final_compile.W)
    elm.Wire('|-', arrow='->').at(final_rGal.S).to(final_compile.E)
    final_state    = flow.StateEnd(r=1.30).at((1.30, -15.25)).label(r'$\mathrm{FINAL}$' '\n' r'$\mathrm{PREDICTED}$' '\n' r'$\mathrm{SOURCES}$')
    elm.Wire('-', arrow='->').at(final_compile.S).to(final_state.N)

    final_ghost    = flow.Start(h=1.0, w=2.5).at((0, -18.0))
    d.draw(show=True)
    if save_plot_flag:
        name_file = 'flowchart_pipeline_extended_ltx.pdf'
        if draw_new_branch:
            name_file = 'flowchart_pipeline_extended_ltx_highlight.pdf'
        d.save(gv.plots_path + name_file)
print('EOF')
