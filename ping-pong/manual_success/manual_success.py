# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 16:07:53 2024

@author: Mike
"""

import numpy as np
from sklearn.ensemble import RandomForestClassifier
import matplotlib.pyplot as plt
import random

spike_file = "inpstatic.txt"
target_file = "rewstastates.txt"
ntactsperstate = 10
nlearningstates = 120000

nVelocityZones = 9
nSpatialZones = 30
RacketSize = 0.18

duration1 = 10

# Read receptor protocol file, aggregate spike counts (in ntactsperstate tacts), read target labels (0 - no class label)

def read_data(filX, strY, ncols, nrows):
    X = np.zeros((nrows, ncols))
    Y = np.zeros(nrows)
    for i in range(nrows):
        for j in range(ntactsperstate):
            s = filX.readline()
            if len(s) < ncols:
                print("Too few records in ", spike_file)
                exit(-1)
            for k in range(ncols):
                if s[k] == '@':
                    X[i][k] += 1
        s = strY[i]
        try:
            Y[i] = int(s) + 1
        except:
            Y[i] = 0
    return X, Y

with open(spike_file) as filspikes, open(target_file) as filtargets:
    targets = filtargets.readlines()
    s = filspikes.readline()
    ninputs = len(s) - 1
    filspikes.seek(0)
    X, Y = read_data(filspikes, targets, ninputs, len(targets))

########################################################################################
# Determine median of velocity bins

minSpotPassageTime_ms = 200
maxSpotPassageTime_ms = 700

def rMakeBallVelocity():
    while True:
        i1 = random.randint(minSpotPassageTime_ms, maxSpotPassageTime_ms - 1);
        i2 = random.randint(0, maxSpotPassageTime_ms - 1);
        if i2 <= i1:
            break
    return 1. / i1;
    
vr_samples = np.zeros(9000)
for i in range(len(vr_samples)):
    rBallVelocity = rMakeBallVelocity()
    rBallMovementDirection = np.random.random() * np.pi / 2
    vr_samples[i] = rBallVelocity * np.sin(rBallMovementDirection)
np.sort(vr_samples)
vr_VelocityZoneMedian = np.zeros((nVelocityZones - 1) // 2)
for i in range(len(vr_VelocityZoneMedian)):
    vr_VelocityZoneMedian[i] = vr_samples[int((i + 1) * 2 * len(vr_samples) / 9)]

########################################################################################
# For states with x < 0 and vx < 0 determine states which lead to success - manual_pred.  

manual_pred = np.zeros(len(X))
for i in range(len(X)):
    ind = 0
    d = 0
    x = -0.5 + 0.5 / nSpatialZones
    n = 0;
    for indx in range(nSpatialZones):
        d += X[i][ind] * x
        n += X[i][ind]
        x += 1 / nSpatialZones
        ind += 1
    x = d / n if n > 0 else -1
    d = 0
    y = -0.5 + 0.5 / nSpatialZones
    n = 0;
    for indx in range(nSpatialZones):
        d += X[i][ind] * y
        n += X[i][ind]
        y += 1 / nSpatialZones
        ind += 1
    y = d / n if n > 0 else -1
    d = 0
    n = 0;
    for ivx in range(nVelocityZones):
        d += X[i][ind] * (0 if ivx == nVelocityZones // 2 else -vr_VelocityZoneMedian[nVelocityZones // 2 - ivx - 1] if ivx < nVelocityZones // 2 else vr_VelocityZoneMedian[ivx - nVelocityZones // 2 - 1])
        n += X[i][ind]
        ind += 1
    vx = d / n if n > 0 else 0
    d = 0
    n = 0;
    for ivy in range(nVelocityZones):
        d += X[i][ind] * (0 if ivy == nVelocityZones // 2 else -vr_VelocityZoneMedian[nVelocityZones // 2 - ivy - 1] if ivy < nVelocityZones // 2 else vr_VelocityZoneMedian[ivy - nVelocityZones // 2 - 1])
        n += X[i][ind]
        ind += 1
    vy = d / n if n > 0 else 0
    d = 0
    ry = -0.5 + 0.5 / nSpatialZones
    n = 0;
    for indx in range(nSpatialZones):
        d += X[i][ind] * ry
        n += X[i][ind]
        ry += 1 / nSpatialZones
        ind += 1
    ry = d / n if n > 0 else -1
    dx = x + 0.5
    if vx < 0 and x < 0:
        yint = y - dx * vy / vx + 0.5
        yint -= np.floor(yint / 2) * 2
        ry += 0.5
        if ry - RacketSize / 2 < yint < ry + RacketSize / 2:
            manual_pred[i] = 1
        else:
            ry = 2 - ry
            if ry - RacketSize / 2 < yint < ry + RacketSize / 2:
                manual_pred[i] = 1

########################################################################################
# Determine theoretical limit for good state determination on the basis of the receptor signal

real = 0
predicted = 0
correcty_predicted = 0
for i in range(len(X)):
    if Y[i] > 0:
        if Y[i] == 2:    
            real += 1
        if manual_pred[i] == 1:
            predicted += 1
            if Y[i] == 2:
                correcty_predicted += 1
rPrecision = correcty_predicted / predicted if predicted > 0 else  0.
rRecall = correcty_predicted / real if real > 0 else 0.
F =  2 / (1 / rPrecision + 1 / rRecall) if rPrecision * rRecall > 0 else 0.
print("THEORETICAL: Precision=", rPrecision, " Recall=", rRecall, " F=", F)

########################################################################################
    
clf = RandomForestClassifier()

def cla(learning_time):
    Xtrain = X[0:learning_time]
    Ytrain = Y[0:learning_time]
    Xtest = X[learning_time:]
    Ytest = Y[learning_time:]
    xx = []
    yy = []
    for i in range(learning_time):
        if Ytrain[i] > 0:
            yy.append(Ytrain[i] - 1)
            xx.append(Xtrain[i])
    clf.fit(xx, yy)
    xx = []
    yy = []
    for i in range(len(Ytest)):
        if Ytest[i] > 0:
            yy.append(Ytest[i] - 1)
            xx.append(Xtest[i])
    pred = clf.predict(xx)
    real = 0
    predicted = 0
    correcty_predicted = 0
    for i in range(len(yy)):
        p = int(pred[i])
        if yy[i] == 1:
            real += 1
        if p == 1:
            predicted += 1
            if p == yy[i]:
                correcty_predicted += 1
    rPrecision = correcty_predicted / predicted if predicted > 0 else  0.
    rRecall = correcty_predicted / real if real > 0 else 0.
    F =  2 / (1 / rPrecision + 1 / rRecall) if rPrecision * rRecall > 0 else 0.
    print("done ", learning_time)
    return rPrecision, rRecall, F

pre = np.zeros(9)
rec = np.zeros(9)
f = np.zeros(9)
tacts = range(20000, 200000, 20000)
res = []
for i in range(len(tacts)):
    pre[i], rec[i], f[i] = cla(tacts[i])  
fig, ax = plt.subplots()
ax.plot(tacts, pre, color = "red")
ax.plot(tacts, rec, color = "green")
ax.plot(tacts, f, color = "blue")
plt.show()
