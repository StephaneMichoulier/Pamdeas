#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib as mpl
import sys

# Arg passed by plot-eden
mass_or_size = str(sys.argv[1])

# Path
path1 = "../build/"
path2 = "../build/"

# PACED data
r0 = [5,10, 20, 50, 75, 100, 150, 200, 250, 300]

density = 917.0

# r-s plot
plt.clf()
fig, axes = plt.subplots(1, 1, sharex="col", sharey="row", gridspec_kw={"hspace": 0,
            "wspace": 0, "left": 0.11, "right": 0.99, "bottom": 0.1, "top": .975},
            figsize=(7., 4.8), dpi=300)

for rad in r0:
    file = path1+"outputpor_"+mass_or_size+"_"+str(rad)+".out"
    with open(file, "r") as f:
        next(f) # skip first line
        t =[]
        s = []
        r = []
        phi=[]
        aero=[]
        for ligne in f:
            ligne = ligne.split()
            t.append(float(ligne[0]))
            s.append(float(ligne[4]))
            r.append(float(ligne[1]))
            phi.append(float(ligne[3]))
    aero=np.multiply(s,phi)
    line = axes.scatter(r, aero*density, c=t, s=0.5, norm=colors.LogNorm(1,1e5),cmap="rainbow")
    line.set_label(r"$R_0=$"+str(rad)+" au")

clb=fig.colorbar(line, pad=0.005)
clb.set_label(r"$t$ (yr)")

# Axes and labels
axes.set(xlim=(1,310), ylim=(5e-8,5e6), xscale="log",yscale="log",
             xlabel=r"$R$ (AU)", ylabel=r"$\rho_s \phi s$ (kg/m²)")
#axes.set_xticks()
axes.xaxis.set_minor_locator(plt.LogLocator(base=10.0, numticks=40))
axes.xaxis.set_minor_formatter(plt.NullFormatter())
axes.tick_params(direction="in", which="both", right=True, top=True)
axes.xaxis.set_label_coords(0.5, -0.05)
axes.yaxis.set_label_coords(-0.08, 0.5)
axes.legend(fontsize="small", loc="upper right")

plt.savefig(path2+"fig_"+mass_or_size+"_r-rhosphis-t.png")

print("fig_"+mass_or_size+"_r-rhosphis-t.py done")
