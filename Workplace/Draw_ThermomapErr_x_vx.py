# -*- coding: utf-8 -*-
"""
Created on Fri May  8 18:01:18 2026

@author: Mike
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df = pd.read_csv("ThermomapErr_x_vx.csv", index_col=False, header=None)
theor_err = np.zeros((4, 30))
SNN_err = np.zeros((4, 30))
no_prediction = np.zeros((4, 30))
for i, row in df.iterrows():
    no_prediction[int(row[1])][int(row[0])] = row[2]
    theor_err[int(row[1])][int(row[0])] = row[3]
    SNN_err[int(row[1])][int(row[0])] = row[4]

vmin = min(np.min(theor_err), np.min(SNN_err))
vmax = max(np.max(theor_err), np.max(SNN_err))

fig, axs = plt.subplots(3)
norm = mpl.colors.Normalize(vmin = vmin, vmax = vmax)
im = axs[0].imshow(no_prediction, cmap = 'seismic')
axs[0].set_title("part of endless loops")
fig.colorbar(im, ax = axs[0])
im = axs[1].imshow(theor_err, cmap = 'seismic', norm=norm)
axs[1].set_title("theoretical error")
fig.colorbar(im, ax = axs[1])
im = axs[2].imshow(SNN_err, cmap = 'seismic', norm=norm)
axs[2].set_title("SNN error")
fig.colorbar(im, ax = axs[2])
plt.show()
