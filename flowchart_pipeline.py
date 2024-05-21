#!/usr/bin/env python

import schemdraw
from schemdraw import flow
import schemdraw.elements as elm
import matplotlib as mpl
import matplotlib.pyplot as plt
import global_variables as gv

mpl.rcdefaults()
plt.rcParams['text.usetex'] = True

save_plot_flag = False

with schemdraw.Drawing(show=False, margin=-0.66) as d:
    init           = flow.Start(h=1.5, w=2.75).at((0, 0)).label('SOURCE\nFROM\nCATALOGUE')
    AGN_gal_model  = flow.Decision(h=2.25, w=3.75, W='Predicted\nas AGN', E='Predicted\nas SFG').at((0, -2.5)).label('AGN/SFG\nCLASSIFICATION\nMODEL')
    elm.Wire('-', arrow='->').at(init.S).to(AGN_gal_model.N)
    rAGN_model     = flow.Decision(h=2.25, w=4, S='Predicted\nas radio', E='Predicted\nas no radio').at((-5, -4.5)).label('AGN RADIO\nDETECTION\nMODEL')
    rGal_model     = flow.Decision(h=2.25, w=4, S='Predicted\nas radio', W='Predicted\nas no radio').at((5, -4.5)).label('SFG RADIO\nDETECTION\nMODEL')
    elm.Wire('-|', arrow='->').at(AGN_gal_model.W).to(rAGN_model.N)
    elm.Wire('-|', arrow='->').at(AGN_gal_model.E).to(rGal_model.N)
    discarded      = flow.StateEnd(r=1.25).at((0, -6.5)).label('NON-RADIO\nSOURCE')
    z_rAGN_model   = flow.Box(h=1.5, w=2.5).at((-5, -7.75)).label('RAGN Z\nPREDICTION\nMODEL')
    z_rGal_model   = flow.Box(h=1.5, w=2.5).at((5, -7.75)).label('RSFG Z\nPREDICTION\nMODEL')
    elm.Wire('-', arrow='->').at(rAGN_model.S).to(z_rAGN_model.N)
    elm.Wire('-', arrow='->').at(rGal_model.S).to(z_rGal_model.N)
    elm.Wire('-|', arrow='->').at(rAGN_model.E).to(discarded.N)
    elm.Wire('-|', arrow='->').at(rGal_model.W).to(discarded.N)
    final_rAGN     = flow.Box(h=1.5, w=2.5).at((-5, -10.25)).label('PREDICTED\nRADIO AGN\nW/REDSHIFT')
    final_rGal     = flow.Box(h=1.5, w=2.5).at((5, -10.25)).label('PREDICTED\nRADIO SFG\nW/REDSHIFT')
    elm.Wire('-', arrow='->').at(z_rAGN_model.S).to(final_rAGN.N)
    elm.Wire('-', arrow='->').at(z_rGal_model.S).to(final_rGal.N)
    final_compile  = flow.Box(h=1.5, w=2.5).at((0, -11.75)).label('COMPILE\nPREDICTIONS')
    elm.Wire('|-', arrow='->').at(final_rAGN.S).to(final_compile.W)
    elm.Wire('|-', arrow='->').at(final_rGal.S).to(final_compile.E)
    final_state    = flow.StateEnd(r=1.25).at((1.25, -15.25)).label('FINAL\nPREDICTED\nSOURCES')
    elm.Wire('-', arrow='->').at(final_compile.S).to(final_state.N)

    final_ghost    = flow.Start(h=1.0, w=2.5).at((0, -18.0))
    # d.draw(show=True)
    if save_plot_flag:
        d.save(gv.plots_path + 'flowchart_pipeline_extended.pdf')
print('EOF')