# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 11:52:26 2025

@author: Kiselev_Mi
"""

fileori = "WorldDynamics.snn.csv"
fileres = "WorldDynamics.upper.snn.csv"

PresynapticId = -2

with open(fileori, 'rt') as filin, open(fileres, 'wt') as filout:
    lstr = filin.readlines()
    lstr = [s.strip() for s in lstr]
    i = lstr.index("$PopulationProperties")
    k = i + 3
    Sections = []
    while lstr[k][0] != '$':
        Sections.append([lstr[k], 0])
        k += 1
    m = Sections.index(["MEMORYT", 0])
    for j in range(len(Sections)):
        k = lstr.index("$Neurons", k)
        k += 3
        Sections[j][1] = int(lstr[k].split(',')[0])
    indMEMORYT = [Sections[m][1], Sections[m + 1][1]]
    del lstr[i+3:i+7]
    i = lstr.index("$Neurons", i)
    k = i
    for j in range(4):
        k = lstr.index("$Neurons", k + 1)
    del lstr[i:k]
    m = i
    to_delete = []
    while True:
        try:
            i = lstr.index("$Links", i)
        except ValueError:
            break
        i += 3
        while lstr[i][0] != '$':
            if lstr[i].split(',')[PresynapticId][0] == '-':
                to_delete.append(i)    
            i += 1
    lstr = [lstr[l] for l in range(len(lstr)) if not l in to_delete]
    while True:
        try:
            m = lstr.index("$Links", m)
        except ValueError:
            break
        m += 3
        while lstr[m][0] != '$':
            lstrrow = lstr[m].split(',')
            src = int(lstrrow[PresynapticId]) - 1
            if indMEMORYT[0] <= src < indMEMORYT[1]:
                src = -(src - indMEMORYT[0]) - 1
                lstrrow[PresynapticId] = str(src)
                lstr[m] = ','.join(lstrrow)
            m += 1
    for s in lstr:
        filout.write(s + '\n')
