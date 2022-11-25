import csv
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import numpy as np
import math

nSpatialZones = 30
nVelocityZones = 9
nRelPos = 5

nGoalLevels = 4

nActions = 2
nNeuronsperAction = 100

mpl.rc('font', size=10)
coo = np.arange(0, 30)

state = []
with open('ping_pong_state.csv', newline = '') as fil:
    csr = csv.reader(fil)
    for row in csr:
        state.append([])
        for v in row:
            if len(v) > 0:
                state[-1].append(float(v))
        state[-1].append((-0.5 - state[-1][0]) * (-0.5 - state[-1][0]) + (state[-1][1] - state[-1][4]) * (state[-1][1] - state[-1][4]))

neuLspikes = []
with open('spikes.3.lst', newline = '') as fil:
    csr = csv.reader(fil)
    for row in csr:
        if len(neuLspikes) == nActions * nNeuronsperAction:
            break
        neuLspikes.append([]);
        Up = len(neuLspikes) > nNeuronsperAction
        for v in row:
            tact = int(v)
            if not Up:

            neuLspikes.append([tact])




RecField = [[]]
maxs = np.zeros((5, 5))
with open('Bayes.csv', newline = '') as fil:
    csr = csv.reader(fil)
    for row in csr:
        if len(row) == 0:
            RecField.append([])
        else:
            RecField[-1].append([np.zeros(nSpatialZones), np.zeros(nSpatialZones), np.zeros(nVelocityZones), np.zeros(nVelocityZones), np.zeros(nSpatialZones), np.zeros((nRelPos, nRelPos))])
            i = 0
            for v in row:
                ind = i
                if ind < nSpatialZones:
                    RecField[-1][-1][0][ind] = int(v)
                    if len(RecField[-1]) < nGoalLevels:
                        maxs[len(RecField) - 1][0] = max(maxs[len(RecField) - 1][0], int(v))
                ind -= nSpatialZones
                if 0 <= ind < nSpatialZones:
                    RecField[-1][-1][1][ind] = int(v)
                    if len(RecField[-1]) < nGoalLevels:
                        maxs[len(RecField) - 1][1] = max(maxs[len(RecField) - 1][1], int(v))
                ind -= nSpatialZones
                if 0 <= ind < nVelocityZones:
                    RecField[-1][-1][2][ind] = int(v)
                    if len(RecField[-1]) < nGoalLevels:
                        maxs[len(RecField) - 1][2] = max(maxs[len(RecField) - 1][2], int(v))
                ind -= nVelocityZones
                if 0 <= ind < nVelocityZones:
                    RecField[-1][-1][3][ind] = int(v)
                    if len(RecField[-1]) < nGoalLevels:
                        maxs[len(RecField) - 1][3] = max(maxs[len(RecField) - 1][3], int(v))
                ind -= nVelocityZones
                if 0 <= ind < nSpatialZones:
                    RecField[-1][-1][4][ind] = int(v)
                    if len(RecField[-1]) < nGoalLevels:
                        maxs[len(RecField) - 1][1] = max(maxs[len(RecField) - 1][1], int(v))
                ind -= nSpatialZones
                if 0 <= ind:
                    RecField[-1][-1][5][int(ind / nRelPos)][ind % nRelPos] = int(v)
                    if len(RecField[-1]) < nGoalLevels:
                        maxs[len(RecField) - 1][4] = max(maxs[len(RecField) - 1][4], int(v))
                i += 1

fig, ax = plt.subplots()
axes_hist = [[] for i in range(len(RecField[0]))]
i = 1
for j in range(len(RecField[0])):
    for k in range(len(RecField[0][0]) - 1):
        axes_hist[j].append(plt.subplot(len(RecField[0]), len(RecField[0][0]) - 1, i))
        i += 1

def draw_Bayes(tact):
    indtact = int(tact)
    i = 0
    for j in range(len(RecField[0])):
        for k in range(len(RecField[0][0]) - 1):
            axes_hist[j][k].clear()
        if j < len(RecField[0]) - 1:
            axes_hist[j][0].set_ylim(0, maxs[indtact][0])
        axes_hist[j][0].plot(RecField[indtact][j][0])
        if j == 0:
            axes_hist[j][0].set_title('x')
        if j < len(RecField[0]) - 1:
            axes_hist[j][1].set_ylim(0, maxs[indtact][1])
        axes_hist[j][1].plot(coo, RecField[indtact][j][1], coo, RecField[indtact][j][4])
        if j == 0:
            axes_hist[j][1].set_title('y')
        if j < len(RecField[0]) - 1:
            axes_hist[j][2].set_ylim(0, maxs[indtact][2])
        axes_hist[j][2].plot(RecField[indtact][j][2])
        if j == 0:
            axes_hist[j][2].set_title('vx')
        if j < len(RecField[0]) - 1:
            axes_hist[j][3].set_ylim(0, maxs[indtact][3])
        axes_hist[j][3].plot(RecField[indtact][j][3])
        if j == 0:
            axes_hist[j][3].set_title('vy')
        if j < len(RecField[0]) - 1:
            axes_hist[j][4].imshow(RecField[indtact][j][5], cmap = 'viridis', vmin = 0, vmax = maxs[indtact][4])
        else:
            axes_hist[j][4].imshow(RecField[indtact][j][5], cmap = 'viridis')

draw_Bayes(0)

axsli = plt.axes([0.25, 0.03, 0.65, 0.03])
sli = Slider(axsli, 'tact', 0., len(RecField) - 1, valinit = 0, valstep = 1, valfmt = "%d")
axsli.xaxis.set_visible(True)
axsli.set_xticks(range(len(RecField) - 1))

def update_sli_Bayes(val):
    tac = sli.val
    draw_Bayes(tac)

sli.on_changed(update_sli_Bayes)

plt.get_current_fig_manager().canvas.set_window_title('Bayes model')
plt.show()

rews = []
with open('rews.csv', newline = '') as fil:
    csr = csv.reader(fil)
    for row in csr:
        rews.append([int(row[0]), int(row[1])])

ns = np.zeros((5, 4))
nscor = np.zeros((5, 4))
nscora = np.zeros((5, 4))
nr = [[] for i in range(nGoalLevels)]
for p in rews:
    if p[0] >= 1000000:
        break
    ns[int(p[0] / 200000)][p[1]] += 1
    if state[p[0]][-1] < state[p[0] - 1][-1]:
        nscor[int(p[0] / 200000)][p[1]] += 1
        if state[p[0]][4] != state[p[0] - 1][4]:
            nscora[int(p[0] / 200000)][p[1]] += 1
    if p[0] >= 800000:
        nr[p[1]].append(math.sqrt(state[p[0]][-1]))
data_cum = ns.cumsum(axis=1)

plt.style.use('bmh')
fig, ax = plt.subplots()
means = []
for i in nr:
    ax.hist(i, histtype="stepfilled", bins=25, alpha=0.8, density=True)
    means.append(np.mean(i))
ax.set_title("rews")
plt.show()

print(means)

category_colors = ['b', 'g', 'r', 'c']

fig, ax = plt.subplots(figsize=(9.2, 5))
ax.invert_yaxis()
ax.xaxis.set_visible(False)
ax.set_xlim(0, np.sum(ns, axis=1).max())

for i in range(4):
    widths = ns[:, i]
    starts = data_cum[:, i] - widths
    color = category_colors[i]
    rects = ax.barh(range(5), widths, left=starts, height=0.5, label=i, color=color)

    #r, g, b, _ = color
    #text_color = 'white' if r * g * b < 0.5 else 'darkgrey'
    #ax.bar_label(rects, label_type='center', color=text_color)

plt.show()

fig, axs = plt.subplots(1, 1, figsize = (10, 3))
for j in range(4):
    x = range(5)
    y = [nscor[t][j] / ns[t][j] for t in x]
    axs.plot(x, y, 'r', label = "c%d" % (j,), linewidth = 3.5 - j)
    y = [nscora[t][j] / ns[t][j] for t in x]
    axs.plot(x, y, 'b', label = "ca%d" % (j,), linewidth = 3.5 - j)
leg = axs.legend(loc = 'best', ncol = 2, mode = "expand", shadow = True, fancybox = True)
leg.get_frame().set_alpha(0.5)
plt.title('cor rewards')
plt.show()

