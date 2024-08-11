#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 17:04:20 2024

@author: mikhail
"""

period_duration = [100, 130, 200]
reward_file = "rewstatic.txt"
state_file = "rewstastates.txt"
state_duration = 10

perbeg = [sum(period_duration[:i+1]) for i in range(len(period_duration))]

with open(reward_file, "rt") as filpro:
    pro = filpro.readlines()
    rewards = [i for i,x in enumerate(pro) if x[0] == '@']
    state_beg = []
    for i in rewards:
        for j in range(len(perbeg)):
            state_beg.append([i - perbeg[-j - 1], j])
            
with open(state_file, "wt") as filstate:
    js = 0
    jr = 0
    current_state = -1
    for i in range(0, len(pro), state_duration):
        if js < len(state_beg) and i + state_duration / 2 >= state_beg[js][0]:
            current_state = state_beg[js][1]
            js += 1
        elif jr < len(rewards) and  i > rewards[jr]:
            current_state = -1
            jr += 1
        if current_state == -1:
            filstate.write('-\n')
        else:
            filstate.write("%d\n" % (current_state,))
                
    
    
    