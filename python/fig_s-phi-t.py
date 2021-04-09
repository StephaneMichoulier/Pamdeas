#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys

# Arg passed by plot-eden
mass_or_size = str(sys.argv[1])

path1 = "../build/"
path2 = "../build/"

# PACED data
radius0 = [10]#,10, 20, 50, 75, 100, 150, 200, 250, 300]

# r-s plot
plt.clf()
fig, axes = plt.subplots(1, 1, sharex = "col", sharey = "row", gridspec_kw = {"hspace": 0,
            "wspace": 0, "left": 0.11, "right": 0.99, "bottom": 0.1, "top": .975},
            figsize = (7., 4.8), dpi = 300)

for radius in radius0:
    file = path1 + "outputpor_" + mass_or_size + "_" + str(radius) + ".out"
    with open(file, "r") as f:
        next(f) # skip first line
        t =[]
        phi = []
        s = []
        for ligne in f:
            ligne = ligne.split()
            t.append(float(ligne[0]))
            phi.append(float(ligne[3]))
            s.append(float(ligne[4]))
    line = axes.scatter(s, phi, c = t, s = 0.5, norm = colors.LogNorm(1.0, max(t)), cmap="rainbow")
    line.set_label(r"$R_0=$" + str(radius) + " au")

clb=fig.colorbar(line, pad = 0.005)
clb.set_label(r"$t$ (yr)")

# Axes and labels
axes.set(xlim = (5e-8, 5e6), xscale = "log", ylim = (1e-5, 1), yscale = "log",
             xlabel = r"$s$ (m)", ylabel = r"$\phi$")
axes.set_xticks([1e-6, 1e-4, 1e-2, 1, 1e2, 1e4, 1e6])
axes.xaxis.set_minor_locator(plt.LogLocator(base = 10.0, numticks = 40))
axes.xaxis.set_minor_formatter(plt.NullFormatter())
axes.tick_params(direction = "in", which = "both", right = True, top = True)
axes.xaxis.set_label_coords(0.5, -0.05)
axes.yaxis.set_label_coords(-0.08, 0.5)
axes.legend(fontsize = "small", loc = "lower left")

plt.savefig(path2 + "fig_" + mass_or_size + "_s-phi-t.png")

print("fig_" + mass_or_size + "_s-phi-t.png done")
