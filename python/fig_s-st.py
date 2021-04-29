#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import sys

# Arg passed by plot-eden
mass_or_size = str(sys.argv[1])

# Path
path1 = "../build/"
path2 = "../build/"

# PACED data
radius0 = [10, 15, 25, 50, 75, 100, 150, 200, 250, 300]
col = ["k", "g", "b", "tab:orange", "tab:purple", "tab:brown", "c", "tab:pink", "m", "y", 
       "tab:blue", "tab:green", "tab:red", "tab:gray", "tab:olive", "tab:cyan"]

# r-s plot
plt.clf()
fig, axes = plt.subplots(1, 1, sharex = "col", sharey = "row", gridspec_kw = {"hspace": 0,
            "wspace": 0, "left": 0.11, "right": 0.99, "bottom": 0.1, "top": .975},
            figsize = (6.4, 4.8), dpi = 300)

for radius, col in zip(radius0, col):
    file = path1 + "outputpor_" + mass_or_size + "_" + str(radius) + ".out"
    with open(file, "r") as f:
        next(f) # skip first line
        s = []
        st = []
        for ligne in f:
            ligne = ligne.split()
            s.append(float(ligne[4]))
            st.append(float(ligne[5]))
    line = axes.scatter(s, st, color = col, s = 0.1)
    line.set_label(r"$R_0=$" + str(radius) + " au")

# Axes and labels
axes.set(xlim = (5e-8, 5e6), xscale = "log", ylim = (1e-7, 1e5), yscale = "log",
             xlabel = r"$s$ (m)", ylabel = r"$St$")
axes.set_xticks([1e-6, 1e-4, 1e-2, 1, 1e2, 1e4, 1e6])
axes.xaxis.set_minor_locator(plt.LogLocator(base = 10.0, numticks = 40))
axes.xaxis.set_minor_formatter(plt.NullFormatter())
axes.tick_params(direction = "in", which = "both", right = True, top = True)
axes.xaxis.set_label_coords(0.5, -0.05)
axes.yaxis.set_label_coords(-0.08, 0.5)
axes.legend(fontsize = "small", loc = "lower left")

plt.axhline(y=1,lw=0.5, color="black")

plt.savefig(path2 + "fig_" + mass_or_size + "_s-st.png")

print("fig_" + mass_or_size + "_s-st.png done")
