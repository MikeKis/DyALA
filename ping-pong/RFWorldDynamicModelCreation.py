# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 15:38:02 2026

@author: Mike
"""

import csv
import random
from sklearn.ensemble import RandomForestClassifier
import numpy as np
import joblib

filein = "ping_pong_state_discrete.csv"
fileout = "ping_pong_state_changes.csv"
filemodel = 'RFpingpongdynamics.joblib'
RFmaxdepth = 30
minPeriod = 10
nPeriods = 4
nStates = [30, 30, 9, 9, 30]
RacketMovementPeriod = 17.74622

print("Reading " + filein + "...")
offset = [0]
for i in range(len(nStates)):
    offset.append(offset[-1] + nStates[i] * nPeriods)
tactcnt = [0 for i in nStates]
LastState = [-1 for i in nStates]
LastCompleteState = [0 for i in range(offset[-1])]
NewCompleteStates = []
nTargets = sum(nStates[:-1]) * nPeriods
with open(filein, newline = '') as filin:
    csr = csv.reader(filin)
    for row in csr:
        tact = int(row[0])
        CurrentState = [int(i) for i in row[1:]]
        for i in range(len(nStates)):
           tactcnt[i] = 1 if CurrentState[i] != LastState[i] else tactcnt[i] + 1
        CurrentCompleteState = [0 for i in range(offset[-1])]
        for i in range(len(nStates)):
            t = minPeriod
            T = tactcnt[i]
            LongestPeriod = True
            for j in range(nPeriods - 1):
                T -= t
                if T < 0:
                    LongestPeriod = False
                    break
                t *= 2
            if LongestPeriod:
                j = nPeriods - 1
            CurrentCompleteState[offset[i] + CurrentState[i] * nPeriods + j] = 1
        if CurrentCompleteState[:nTargets] != LastCompleteState[:nTargets]:
            LastCompleteState = CurrentCompleteState
            NewCompleteStates.append([tact, CurrentCompleteState])
        LastState = CurrentState
    
ntactsperCompleteState = NewCompleteStates[-1][0] / len(NewCompleteStates)
RacketMovementPeriod = int(RacketMovementPeriod / ntactsperCompleteState)
        
with open(fileout, "w") as filout:
    for i in NewCompleteStates:
        filout.write(str(i[0]))
        for j in i[1]:
            filout.write(f",{j}")
        filout.write('\n')

OnlyCompleteStates = [i[1] for i in NewCompleteStates]

def getRF(ind):
    target = [i[ind] for i in OnlyCompleteStates[1:]]
    if target == [0] * len(target):
        return None
    ret = RandomForestClassifier(max_depth = RFmaxdepth)
    ret.fit(OnlyCompleteStates[:-1], target)
    return ret

print("Creating predictors...")
Predictors = []
for i in range(nTargets):
    Predictors.append(getRF(i))
    print(f"Predictor #{i} created")
    
joblib.dump(Predictors, filemodel)    
print("RF world dynamics model saved.")