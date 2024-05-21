#!/usr/bin/env python

import schemdraw
from schemdraw import flow
import schemdraw.elements as elm
import matplotlib as mpl
import matplotlib.pyplot as plt
import global_variables as gv

mpl.rcdefaults()
plt.rcParams['text.usetex'] = False

save_plot_flag = True

with schemdraw.Drawing(show=False) as d:
    # d.config(unit=.75)
    training_data   = flow.Start(h=1.5, w=2.25).at((0, 0)).label(r'$\mathrm{TRAINING}$' + '\n' + r'$\mathrm{DATA}$')

    base_model_1    = flow.Box(h=1, w=2.75).at((5, 2.0)).label(r'$\mathrm{BASE ~ MODEL ~ 1}$')
    base_model_2    = flow.Box(h=1, w=2.75).at((5, 0.5)).label(r'$\mathrm{BASE ~ MODEL ~ 2}$')
    base_model_3    = flow.Box(h=1, w=2.75).at((5, -1.0)).label(r'$\mathrm{BASE ~ MODEL ~ 3}$')
    elm.Dot(radius=0.05).at((5, -2.25))
    elm.Dot(radius=0.05).at((5, -2.50))
    elm.Dot(radius=0.05).at((5, -2.75))
    base_model_n    = flow.Box(h=1, w=2.75).at((5, -3.0)).label(r'$\mathrm{BASE ~ MODEL ~ N}$')

    dot_a = elm.Dot().at((1.75, -0.75))

    elm.Wire('-', arrow='').at(training_data.E).to(dot_a.center)

    elm.Wire('|-', arrow='->').at(dot_a.center).to(base_model_1.W)
    elm.Wire('|-', arrow='->').at(dot_a.center).to(base_model_2.W)
    elm.Wire('|-', arrow='->').at(dot_a.center).to(base_model_3.W)
    elm.Wire('|-', arrow='->').at(dot_a.center).to(base_model_n.W)

    prediction_1    = flow.Data(h=1, w=2.75).at((8.0, 1.5)).label(r'$\mathrm{PREDICTION ~ 1}$')
    prediction_2    = flow.Data(h=1, w=2.75).at((8.0, 0.0)).label(r'$\mathrm{PREDICTION ~ 2}$')
    prediction_3    = flow.Data(h=1, w=2.75).at((8.0, -1.5)).label(r'$\mathrm{PREDICTION ~ 3}$')
    elm.Dot(radius=0.05).at((9.25, -2.25))
    elm.Dot(radius=0.05).at((9.25, -2.50))
    elm.Dot(radius=0.05).at((9.25, -2.75))
    prediction_n    = flow.Data(h=1, w=2.75).at((8.0, -3.5)).label(r'$\mathrm{PREDICTION ~ N}$')

    elm.Wire('-', arrow='->').at(base_model_1.E).to(prediction_1.W)
    elm.Wire('-', arrow='->').at(base_model_2.E).to(prediction_2.W)
    elm.Wire('-', arrow='->').at(base_model_3.E).to(prediction_3.W)
    elm.Wire('-', arrow='->').at(base_model_n.E).to(prediction_n.W)

    dot_b = elm.Dot().at((12.5, -0.75))

    elm.Wire('-|', arrow='').at(prediction_1.E).to(dot_b.center)
    elm.Wire('-|', arrow='').at(prediction_2.E).to(dot_b.center)
    elm.Wire('-|', arrow='').at(prediction_3.E).to(dot_b.center)
    elm.Wire('-|', arrow='').at(prediction_n.E).to(dot_b.center)

    dot_c = elm.Dot(radius=0).at((7.25, -5))
    elm.Wire('|-', arrow='->').at(dot_a.center).to(dot_c.center)
    elm.Wire('-|', arrow='').at(dot_c.center).to(dot_b.center)

    meta_model      = flow.Box(h=1, w=2.75).at((15, -1.25)).label(r'$\mathrm{META ~ MODEL ~ 1}$')

    elm.Wire('-', arrow='->').at(dot_b.center).to(meta_model.W)

    final_prediction = flow.Start(h=1, w=2.75).at((18.5, -0.25)).label(r'$\mathrm{FINAL}$' + '\n' + r'$\mathrm{PREDICTION}$')

    elm.Wire('-', arrow='->').at(meta_model.E).to(final_prediction.W)

    final_ghost = flow.Start(h=1.0, w=2.5).at((0, -28.5))
    d.draw(show=True)
    if save_plot_flag:
        d.save(gv.plots_path + 'flowchart_model_stacking.pdf')
print('EOF')
