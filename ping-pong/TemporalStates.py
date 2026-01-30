# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 15:38:02 2026

@author: Mike
"""

import csv
filein = "ping_pong_state_discrete.csv"
minPeriod = 10
nPeriods = 4
nStates = [30, 30, 9, 9, 30]

with open(filein, newline = '') as filin:
    csr = csv.reader(filin)
    for row in csr:

