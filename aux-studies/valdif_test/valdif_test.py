from numpy import random
from numpy import pi
from numpy import cos
import numpy as np
import matplotlib.pyplot as plt

random.seed(30)

test_len = 2000000
N = 100

def gentest(fmin, fmax, Tp):
    random.seed(30)
    v = [(fmin + 0.5 * (1 - cos(2 * pi * i / Tp)) * (fmax - fmin)) / N for i in range(test_len)]
    ret = [0 for i in range(test_len)]
    for i in range(test_len):
        for j in range(N):
            if random.random() < v[i]:
                ret[i] += 1
    return ret

test = gentest(0, 30, 1000)

aW = [0.01, 0.015, 0.023, 0.039, 0.081]
aTminus = [10, 12, 17, 32, 131]
tau = 10
rel = 1 - 1 / tau
active = [False for i in range(len(aW))]
au = [0 for i in range(len(aW))]
decplus = [[] for i in range(len(aW))]
decminus = [[] for i in range(len(aW))]
tactlastspike = [-1000 for i in range(len(aW))]
for i in range(test_len):
    for j in range(len(aW)):
        au[j] *= rel
        au[j] += test[i] * aW[j]
        if au[j] > 1:
            au[j] = 0
            if not active[j] and i - tactlastspike[j] <= 3:
                decplus[j].append(i)
                active[j] = True
            tactlastspike[j] = i
        elif active[j] and i - tactlastspike[j] > aTminus[j]:
            decminus[j].append(i)
            active[j] = False

def show_dec(tactbeg, tactend):
    x = range(tactbeg, tactend)
    y = test[tactbeg:tactend]
    c = ['k' for i in range(len(x))]
    alpha = [1 for i in range(len(x))]
    for i in range(len(decplus)):
        for tact in decplus[i]:
            if tactbeg <= tact < tactend:
                reltact = tact - tactbeg
                if c[reltact] != 'k':
                    print('conflict at ', tact)
                c[reltact] = 'r'
                alpha[reltact] = 1 - i / len(decplus)
        for tact in decminus[i]:
            if tactbeg <= tact < tactend:
                reltact = tact - tactbeg
                if c[reltact] != 'k':
                    print('conflict at ', tact)
                c[reltact] = 'b'
                alpha[reltact] = 1 - i / len(decplus)
    fig, ax = plt.subplots()
    ax.scatter(x, y, s=10, alpha=alpha, c=c)
    plt.show()

show_dec(1000, 2000)
