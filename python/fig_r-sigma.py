#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib as mpl

# Path
path1 = "../build/"
path2 = "../build/"

# r-sigma plot
plt.clf()
fig, (axe1, axe2, axe3) = plt.subplots(1, 3, gridspec_kw = {"hspace": 0, "wspace": 0.3}, figsize = (15, 4.8), dpi = 300)

file = path1 + "disc_profiles.out"
with open(file, "r") as f:
    next(f) # skip first line
    r = []
    sigma = []
    dustfrac = []
    P = []
    for ligne in f:
        ligne = ligne.split()
        r.append(float(ligne[0]))
        sigma.append(float(ligne[3]))
        dustfrac.append(float(ligne[5]))
        P.append(float(ligne[6]))
axe1.scatter(r, sigma, s = 0.5)
axe2.scatter(r, P, s = 0.5)
axe3.scatter(r, dustfrac, s = 0.5)

# Plot 1, Axes and labels
axe1.set(xlim = (1, 350), ylim=(5e-1, 1e4), xscale = "log", yscale = "log",
             xlabel = r"$R$ (AU)", ylabel = r"$\Sigma$ (kg/m²)")
axe1.set_xticks([3, 10, 50, 100, 150, 200, 300])
axe1.xaxis.set_minor_locator(plt.LogLocator(base = 10.0, numticks = 40))
axe1.xaxis.set_minor_formatter(plt.NullFormatter())
axe1.tick_params(direction = "in", which = "both", right = True, top = True)
axe1.xaxis.set_label_coords(0.5, -0.05)
axe1.yaxis.set_label_coords(-0.11, 0.5)

# Plot 2, Axes and labels
axe2.set(xlim = (1, 350), ylim = (1e-9, 1e-1), xscale = "log", yscale = "log",
             xlabel = r"$R$ (AU)", ylabel = r"$P$ (Pa)")
axe2.set_xticks([3, 10, 50, 100, 150, 200, 300])
axe2.xaxis.set_minor_locator(plt.LogLocator(base = 10.0, numticks = 40))
axe2.xaxis.set_minor_formatter(plt.NullFormatter())
axe2.tick_params(direction = "in", which = "both", right = True, top = True)
axe2.xaxis.set_label_coords(0.5, -0.05)
axe2.yaxis.set_label_coords(-0.11, 0.5)

# Plot 3, Axes and labels
axe3.set(xlim = (1, 350), ylim = (5e-3, 0.5),xscale = "log", yscale = "log",
             xlabel = r"$R$ (AU)", ylabel = r"$Dustfrac$")
axe3.set_xticks([3, 10, 50, 100, 150, 200, 300])
axe3.xaxis.set_minor_locator(plt.LogLocator(base = 10.0, numticks = 40))
axe3.xaxis.set_minor_formatter(plt.NullFormatter())
axe3.tick_params(direction = "in", which = "both", right = True, top = True)
axe3.xaxis.set_label_coords(0.5, -0.05)
axe3.yaxis.set_label_coords(-0.11, 0.5)

plt.savefig(path2 + "fig_r-sigma.png")

print("fig_r-sigma.png done")
