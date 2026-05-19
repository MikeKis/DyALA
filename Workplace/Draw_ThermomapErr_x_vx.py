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
dfRF = pd.read_csv("ThermomapErr_x_vx_RF.csv", index_col=False, header=None)
theor_err = np.zeros((4, 30))
SNN_err = np.zeros((4, 30))
RF_err = np.zeros((4, 30))
no_prediction = np.zeros((4, 30))
no_prediction_RF = np.zeros((4, 30))
for i, row in df.iterrows():
    no_prediction[int(row[1])][int(row[0])] = row[2]
    theor_err[int(row[1])][int(row[0])] = row[3]
    SNN_err[int(row[1])][int(row[0])] = row[4]
for i, row in dfRF.iterrows():
    no_prediction_RF[int(row[1])][int(row[0])] = row[2]
    RF_err[int(row[1])][int(row[0])] = row[3]

vmin = min(np.min(theor_err), np.min(SNN_err), np.min(RF_err))
vmax = max(np.max(theor_err), np.max(SNN_err), np.max(RF_err))
vmaxnopred = max(np.max(no_prediction), np.max(no_prediction_RF))
norm = mpl.colors.Normalize(vmin = vmin, vmax = vmax)
normnopred = mpl.colors.Normalize(vmin = 0, vmax = vmaxnopred)

fig, axs = plt.subplots(5)
im = axs[0].imshow(no_prediction, cmap = 'seismic', norm=normnopred)
axs[0].set_title("part of endless loops (SNN)")
fig.colorbar(im, ax = axs[0])
im = axs[1].imshow(no_prediction_RF, cmap = 'seismic', norm=normnopred)
axs[1].set_title("part of endless loops (RF)")
fig.colorbar(im, ax = axs[1])
im = axs[2].imshow(theor_err, cmap = 'seismic', norm=norm)
axs[2].set_title("theoretical error")
fig.colorbar(im, ax = axs[2])
im = axs[3].imshow(SNN_err, cmap = 'seismic', norm=norm)
axs[3].set_title("SNN error")
fig.colorbar(im, ax = axs[3])
im = axs[4].imshow(RF_err, cmap = 'seismic', norm=norm)
axs[4].set_title("RF error")
fig.colorbar(im, ax = axs[4])
plt.show()
