#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys

# Arg passed by plot-eden
mass_or_size = str(sys.argv[1])

# Path
path1 = "../build/"
path2 = "../build/"

# PACED data
radius0 = [5, 10, 20, 50, 75, 100, 150, 200, 250, 300]

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
        r = []
        s = []
        for ligne in f:
            ligne = ligne.split()
            t.append(float(ligne[0]))
            r.append(float(ligne[1]))
            s.append(float(ligne[4]))
    line = axes.scatter(t, r, c = s, s = 0.5, norm = colors.LogNorm(1e-7, 1e4), cmap = "inferno")
    line.set_label(r"$R_0=$" + str(radius) + " au")

clb=fig.colorbar(line, pad = 0.005)
clb.set_label(r"$s$ (m)")

# Axes and labels
axes.set(xlim = (1, 1e8), xscale = "log", ylim = (2, 320), yscale = "log",
             xlabel = r"$t$ (yr)", ylabel = r"$R$ (AU)")
axes.xaxis.set_minor_locator(plt.LogLocator(base = 10.0, numticks = 40))
axes.xaxis.set_minor_formatter(plt.NullFormatter())
axes.tick_params(direction = "in", which = "both", right = True, top = True)
axes.xaxis.set_label_coords(0.5, -0.05)
axes.yaxis.set_label_coords(-0.08, 0.5)
axes.legend(fontsize = "small", loc = "lower left")

plt.savefig(path2 + "fig_" + mass_or_size + "_t-r-s.png")

print("fig_" + mass_or_size + "_t-r-s.png done")

