# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 15:38:02 2026

@author: Mike
"""

import joblib
import random
import numpy as np
from ReadStateChanges import ReadStateChanges

wd = "/home/mike/E/DyALA/Workplace/"
filemodel = 'RFpingpongdynamics.joblib'
filein = "ping_pong_state_discrete.csv"
fileout = "RFPredictonResults.csv"
nTests = 3000

nPeriods = 4
minPeriod = 10
RacketMovementPeriod = 17.74622
nStates = [30, 30, 9, 9, 30]
offset = [0]
for i in range(len(nStates)):
    offset.append(offset[-1] + nStates[i] * nPeriods)

print("Loading RF model...")
Predictors = joblib.load(wd + filemodel)    

NewCompleteStates, nTargets = ReadStateChanges(wd + filein)
ntactsperCompleteState = NewCompleteStates[-1][0] / len(NewCompleteStates)
RacketMovementPeriod = int(RacketMovementPeriod / ntactsperCompleteState)
OnlyCompleteStates = [i[1] for i in NewCompleteStates]
    
def SymbolicState(CompleteState):
    ret = []
    for i in offset[:-1]:
        ind = i
        while CompleteState[ind] == 0:
            ind += 1
        v = (ind - i) // nPeriods
        t = ind - i - v * nPeriods
        ret.append([v, t])
    return ret

def NextCompleteState(Predictors, CompleteState, SymSta):
    global tactcnt, tactcntRacket
    ExpiredStates = [i for i in range(len(SymSta) - 1) if tactcnt[i] > minPeriod * 2 ** SymSta[i][1]]
    AllowedTransfer = [True for i in CompleteState[:nTargets]]
    if len(ExpiredStates) > 0:
        for i in range(len(nStates) - 1):
            if not i in ExpiredStates:
                AllowedTransfer[offset[i]:offset[i + 1]] =[False] * (offset[i + 1] - offset[i])
    pre = [0 if i == None else i.predict_proba([CompleteState])[0][np.where(i.classes_ == 1)[0][0]] for i in Predictors]
    pmax = 0
    for i in range(nTargets):
        if CompleteState[i] == 0 and AllowedTransfer[i] and pre[i] > pmax:
            pmax = pre[i]
            Changedto = i
    if pmax == 0:
        print("No next state!")
        exit(0)
    for i in tactcnt[:-1]:
        i += ntactsperCompleteState
    dim = 0
    while offset[dim] <= Changedto:
        dim += 1
    dim -= 1
    tactcnt[dim] = 0
    NewCompleteState = CompleteState[:offset[dim]] + [0] * (nStates[dim] * nPeriods) + CompleteState[offset[dim + 1]:offset[-2]] + [0] * (nStates[-1] * nPeriods)
    NewCompleteState[Changedto] = 1
    OldRacketPosition = SymSta[-1][0]
    if tactcntRacket > RacketMovementPeriod:
        tactcntRacket -= RacketMovementPeriod
        SymSta[-1][0] += -1 if random.random() < 0.5 else 0
        SymSta[-1][0] = max(0, min(nStates[-1] - 1, SymSta[-1][0]))
    tactcntRacket += ntactsperCompleteState
    if SymSta[-1][0] != OldRacketPosition:
       tactcnt[-1] = 0 
    else:
       tactcnt[-1] += ntactsperCompleteState
    t = minPeriod
    T = tactcnt[-1]
    LongestPeriod = True
    for j in range(nPeriods - 1):
        T -= t
        if T < 0:
            LongestPeriod = False
            break
        t *= 2
    if LongestPeriod:
        j = nPeriods - 1
    NewCompleteState[offset[-2] + SymSta[-1][0] * nPeriods + j] = 1
    return NewCompleteState, SymbolicState(NewCompleteState)

with open(wd + fileout, "wt", buffering=1) as filout:
    filout.write("real,pred\n");    
    res = []
    for i in range(nTests):
        start = int(random.random() * len(OnlyCompleteStates))
        StartingState = OnlyCompleteStates[start]
        SymSta = SymbolicState(StartingState)
        if SymSta[0][0] != 0:
            j = start + 1
            while j < len(OnlyCompleteStates) and SymbolicState(OnlyCompleteStates[j])[0][0] != 0:
                j += 1
            if j < len(OnlyCompleteStates):
                yreal = SymbolicState(OnlyCompleteStates[j])[1][0]
                tactcnt = [0 for i in nStates]
                tactcntRacket = 0
                k = 0
                Visited = set()
                Cycled = False
                while SymSta[0][0] != 0 and k < 3000:
                    print(f"test {i} step {k}")
                    StartingState, SymSta = NextCompleteState(Predictors, StartingState, SymSta)
                    t = tuple(tuple(item) for item in SymSta)
                    if t in Visited:
                        Cycled = True
                        break
                    k += 1
                    Visited.add(t)
                if k == 3000 or Cycled:
                    print("no prediction")
                    ypred = -1
                else: 
                    print(f"steps to prediction: {k}")
                    ypred = SymSta[1][0]
                res.append([yreal, ypred])
                filout.write(f"{yreal},{ypred}\n");    
        
    print(res)