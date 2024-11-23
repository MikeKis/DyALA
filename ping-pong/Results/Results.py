import csv
import matplotlib as mpl
from collections import namedtuple
from matplotlib.widgets import Slider
import numpy as np
# import statistics
# import bisect
# import math
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

file = "monitoring.102.csv"

SaturatedResourceThreshold = 30
norm = mpl.colors.TwoSlopeNorm(vmin = -SaturatedResourceThreshold, vcenter = 0, vmax = SaturatedResourceThreshold)
normdif = mpl.colors.TwoSlopeNorm(vmin = -0.3, vcenter = 0, vmax = 0.3)

nActions = 2
nSpatialZones = 30
nVelocityZones = 10
nRelPos = 5
coo = np.arange(0, nSpatialZones)

indrecinp = 2

#--------------------------------------------------------------------------------------------
# PARAMETERS OF INTERMEDIATE GOAL BLOCK
# We assume that group includes 1 neuron
nRecognizingGroupsLevel0 = 1
LevelNGroupMultiplier = 1
nLevels = 7
#--------------------------------------------------------------------------------------------

nNeuronsonLevel = [int(nRecognizingGroupsLevel0 * LevelNGroupMultiplier ** i) for i in range(nLevels)]
nPrimaryStateRecognizers = sum(nNeuronsonLevel)

def get_neuron_level(ind):
    i = ind
    lev = 0
    while i >= 0:
        i -= nNeuronsonLevel[lev]
        lev += 1
    return lev - 1

sec = []
CharTime = []
ThresholdExcessIncrementMultiplier = []
AbsRefT = []
WINC = []
Inh = []
minW = []
maxW = []
SectionName = []

secint = []
SectionIntensity = namedtuple('SectionIntensity', ['sec', 'tact', 'relfre'])
lin = []
Link = namedtuple('Link', 'tact,neu,type,indsyn,Is_pla,delay,src,effw,basew,W')
tact = []
neusec = []
migrations = []
ReceptorBounds = []

neuint = [0 for i in range(nPrimaryStateRecognizers)]
stability = [[] for i in range(nPrimaryStateRecognizers)]
intensity = [[] for i in range(nPrimaryStateRecognizers)]
neubaseL = -1

# It is guaranteed that all records are ordered by tact

with open(file, newline = '') as fil:
    bHeader = True
    csr = csv.reader(fil)
    for row in csr:
        if bHeader:
            bHeader = False
        elif row[0] == "secsta":
            if row[12] == "L":
                indsecL = len(SectionName)
            sec.append(float(row[1]))
            CharTime.append(float(row[2]))
            ThresholdExcessIncrementMultiplier.append(row[5])
            AbsRefT.append(float(row[6]))
            WINC.append(float(row[7]))
            Inh.append(float(row[8]))
            minW.append(float(row[9]))
            maxW.append(float(row[10]))
            SectionName.append(row[12])
        elif row[0] == "rec":
            ReceptorBounds.append([int(row[1]), int(row[2])])
        elif row[0] == "secint":
            if float(row[2]) >= 0:
                secint.append(SectionIntensity(float(row[2]), float(row[1]), float(row[3])))
        elif row[0] == "lin":
            lin.append(Link(float(row[1]), int(row[2]), float(row[3]), float(row[4]), float(row[5]) != 0, float(row[6]), int(row[7]), float(row[8]), float(row[9]), float(row[10])))
        elif row[0] == "neu->sec":
            tac = int(row[1])
            neu = int(row[2])
            s = int(row[3])
            
            if s == indsecL and neubaseL == -1:
                neubaseL = neu

            # It is assumed that the 1st neu->sec section is ordered ny neu

            if len(tact) == 0 or tac != tact[-1]:
                tact.append(tac)
                neusec.append([] if tac == 0 else [-1 for i in range(len(neusec[0]))])
                migrations.append(0)
            if tac == 0:
                neusec[-1].append(s)
            else:
                neusec[-1][neu] = s
        elif row[0] == "neudyn":
            migrations[-1] += 1
        elif row[0] == "neu":
            neu = int(row[2])
            if neusec[-1][neu] == indsecL:
                neuint[neu - neubaseL] += int(row[3])
                stability[neu - neubaseL].append(float(row[13]))
                intensity[neu - neubaseL].append(int(row[3]))

c = list(mcolors.TABLEAU_COLORS.values())
cmap = plt.colormaps["plasma"]
fig, axs = plt.subplots()
l = 0
for j in range(nLevels):
    for k in range(nNeuronsonLevel[j]):
        axs.plot(stability[l], color = cmap((nLevels - j) / nLevels))
        l += 1
fig.suptitle('Neuron stability', fontsize=16)
plt.show()
fig, axs = plt.subplots()
l = 0
for j in range(nLevels):
    for k in range(nNeuronsonLevel[j]):
        axs.plot(intensity[l], color = cmap((nLevels - j) / nLevels))
        l += 1
fig.suptitle('Neuron intensity', fontsize=16)
plt.show()

effw = [t.effw if t.type == 0 and t.Is_pla else t.W for t in lin]

i = 0
deffw_sum = []
deffw_level = []
sat_sup_dist = []
tactNo = -1
guess_shift = -1

NSuppressedLinks = 0
NSaturatedLinks = 0
NPlasticLinks = 0

RecField = []

print("Link data processing...")
cnt = 0
lastneu = -1
while i <= len(lin):
    end = i == len(lin)
    if i == 0 or end or lin[i].tact != lin[i - 1].tact:
        guess_shift = cnt
        cnt = 0
    if end or lin[i].Is_pla:
        if end or lin[i].neu != lastneu:
            if len(sat_sup_dist) > 0:
                sat_sup_dist[-1][0].append(NSuppressedLinks)
                sat_sup_dist[-1][1].append(NSaturatedLinks)
            if not end:
                lastneu = lin[i].neu
                NSuppressedLinks = 0
                NSaturatedLinks = 0
        if not end:
            if tactNo < 0 or lin[i].tact != tact[tactNo]:
                sat_sup_dist.append([[], []])
                RecField.append([[np.zeros(nSpatialZones), np.zeros(nSpatialZones), np.zeros(nVelocityZones), np.zeros(nVelocityZones), np.zeros(nSpatialZones), np.zeros((nRelPos, nRelPos))] for i in range(nPrimaryStateRecognizers)])
                tactNo += 1
                print('Now tact %d is processed' % tact[tactNo])
            if tactNo > 0:
                if len(deffw_sum) < tactNo:
                    deffw_sum.append(0)
                    deffw_level.append([0 for i in range(nLevels)])
                if i - guess_shift < 0:
                    guess_shift = i
                if lin[i - guess_shift].neu < lin[i].neu:
                    while lin[i - guess_shift].neu < lin[i].neu:
                        guess_shift -= 1
                elif lin[i - guess_shift].neu > lin[i].neu:
                    while lin[i - guess_shift].neu > lin[i].neu:
                        guess_shift += 1
                if lin[i - guess_shift].src > lin[i].src:
                    while lin[i - guess_shift].src > lin[i].src and lin[i - guess_shift].neu == lin[i].neu:
                        guess_shift -= 1
                elif lin[i - guess_shift].src < lin[i].src:
                    while lin[i - guess_shift].src < lin[i].src and lin[i - guess_shift].neu == lin[i].neu:
                        guess_shift += 1
                if lin[i - guess_shift].src == lin[i].src and lin[i - guess_shift].neu == lin[i].neu:
                    deffw_sum[-1] += abs(effw[i] - effw[i - guess_shift])
                else:
                    deffw_sum[-1] += abs(effw[i])

                if neusec[tactNo][lin[i].neu] == indsecL and lin[i].src < 0 and ReceptorBounds[indrecinp][0] <= -lin[i].src - 1 < ReceptorBounds[indrecinp][1]:
                    if lin[i - guess_shift].src == lin[i].src and lin[i - guess_shift].neu == lin[i].neu:
                        deffw_level[-1][get_neuron_level(lin[i].neu - neubaseL)] += abs(effw[i] - effw[i - guess_shift])
                    else:
                        deffw_level[-1][get_neuron_level(lin[i].neu - neubaseL)] += abs(effw[i])
                    ind = -lin[i].src - 1 - ReceptorBounds[indrecinp][0]
                    if ind < nSpatialZones:
                        RecField[-1][lin[i].neu - neubaseL][0][ind] = lin[i].W
                    ind -= nSpatialZones
                    if 0 <= ind < nSpatialZones:
                        RecField[-1][lin[i].neu - neubaseL][1][ind] = lin[i].W
                    ind -= nSpatialZones
                    if 0 <= ind < nVelocityZones:
                        RecField[-1][lin[i].neu - neubaseL][2][ind] = lin[i].W
                    ind -= nVelocityZones
                    if 0 <= ind < nVelocityZones:
                        RecField[-1][lin[i].neu - neubaseL][3][ind] = lin[i].W
                    ind -= nVelocityZones
                    if 0 <= ind < nSpatialZones:
                        RecField[-1][lin[i].neu - neubaseL][4][ind] = lin[i].W
                    ind -= nSpatialZones
                    if 0 <= ind:
                        RecField[-1][lin[i].neu - neubaseL][5][int(ind / nRelPos)][ind % nRelPos] = lin[i].W

            if lin[i].W < 0:
                NSuppressedLinks += 1
            elif lin[i].W > 30:
                NSaturatedLinks += 1
            if tactNo == 0:
                NPlasticLinks += 1
    cnt += 1
    i += 1

deffwT = np.array(deffw_level).T
fig, axs = plt.subplots()
l = 0
for j in range(nLevels):
    axs.plot(range(0, len(deffwT[j]) * 200, 200), deffwT[j], color = cmap((nLevels - j) / nLevels))
axs.set_xlabel('Time (sec)')
fig.suptitle('Synaptic weight redistribution volume', fontsize=16)
plt.show()

mpl.rc('font', size=10)

def draw_rec_fields(tac, rec_fields):
    indtact = tact.index(tac)
    i = 0
    if indtact == 0 or not DrawDifference:
        for neu in range(len(rec_fields[0])):
            for k in axes_hist[neu]:
                k.clear()
            axes_hist[neu][0].set_ylim(-30, 30)
            axes_hist[neu][0].plot(rec_fields[indtact][neu][0])
            if neu == 0:
                axes_hist[neu][0].set_title('x')
            axes_hist[neu][1].set_ylim(-30, 30)
            axes_hist[neu][1].plot(coo, rec_fields[indtact][neu][1], coo, rec_fields[indtact][neu][4])
            if neu == 0:
                axes_hist[neu][1].set_title('y')
            axes_hist[neu][2].set_ylim(-30, 30)
            axes_hist[neu][2].plot(rec_fields[indtact][neu][2])
            if neu == 0:
                axes_hist[neu][2].set_title('vx')
            axes_hist[neu][3].set_ylim(-30, 30)
            axes_hist[neu][3].plot(rec_fields[indtact][neu][3])
            if neu == 0:
                axes_hist[neu][3].set_title('vy')
            axes_hist[neu][4].imshow(rec_fields[indtact][neu][5], cmap = 'seismic', norm = norm)
            axes_hist[neu][4].set(xticks=[], yticks=[])
        fig.suptitle('L receptive fields for level {}'.format(m), fontsize=10)
    else:
        for neu in range(len(rec_fields[0])):
            for k in axes_hist[neu]:
                k.clear()
            axes_hist[neu][0].set_ylim(-0.3, 0.3)
            axes_hist[neu][0].plot(rec_fields[indtact][neu][0] - rec_fields[indtact - 1][neu][0])
            if neu == 0:
                axes_hist[neu][0].set_title('x')
            axes_hist[neu][1].set_ylim(-0.3, 0.3)
            axes_hist[neu][1].plot(coo, rec_fields[indtact][neu][1] - rec_fields[indtact - 1][neu][1], coo, rec_fields[indtact][neu][4] - rec_fields[indtact - 1][neu][4])
            if neu == 0:
                axes_hist[neu][1].set_title('y')
            axes_hist[neu][2].set_ylim(-0.3, 0.3)
            axes_hist[neu][2].plot(rec_fields[indtact][neu][2] - rec_fields[indtact - 1][neu][2])
            if neu == 0:
                axes_hist[neu][2].set_title('vx')
            axes_hist[neu][3].set_ylim(-0.3, 0.3)
            axes_hist[neu][3].plot(rec_fields[indtact][neu][3] - rec_fields[indtact - 1][neu][3])
            if neu == 0:
                axes_hist[neu][3].set_title('vy')
            axes_hist[neu][4].imshow(rec_fields[indtact][neu][5] - rec_fields[indtact - 1][neu][5], cmap = 'seismic', norm = normdif)
            axes_hist[neu][4].set(xticks=[], yticks=[])
        fig.suptitle('L receptive field dynamics for level {}'.format(m), fontsize=10)

n = 0
DrawDifference = False
for m in range(nLevels):
    rf =[a[n : n + nNeuronsonLevel[m]] for a in RecField]
    fig, ax = plt.subplots()
    axes_hist = [[] for i in range(nNeuronsonLevel[m])]
    i = 1
    for neu in range(nNeuronsonLevel[m]):
        for k in range(len(RecField[0][0]) - 1):
            axes_hist[neu].append(plt.subplot(nNeuronsonLevel[m], len(RecField[0][0]) - 1, i))
            i += 1

    draw_rec_fields(0, rf)

    axsli = plt.axes([0.25, 0.03, 0.65, 0.03])
    sli = Slider(axsli, 'tact', 0., tact[-1], valinit = 0, valstep = tact[1], valfmt = "%d")
    axsli.xaxis.set_visible(True)
    axsli.set_xticks(tact)

    def update_sli(val):
        tac = sli.val
        draw_rec_fields(tac, rf)

    sli.on_changed(update_sli)

    plt.show()
    n += nNeuronsonLevel[m]


ReceptorSectionBoundaries = [1,2,135]   # It is assumed that first 4 sections are primary and secondary evaluation.
indLA = [616, 668]
SectionNames = ["L", "REWGATE", "LevelMeasurement", "$$$Reward", "MEMLEVEL", "$$$Punishment", "LACT", "GATEREW", "GATEPUN", "GATEREWINT", "GATEPUNINT", "FINALGATEACT"]

nPrimaryStateRecognizers = 298

sec = []
CharTime = []
ThrDecayT = []
ThresholdExcessIncrementMultiplier = []
AbsRefT = []
WINC = []
RelWDec = []
Inh = []

secint = []
SectionIntensity = namedtuple('SectionIntensity', ['sec', 'tact', 'relfre'])
lin = []
Link = namedtuple('Link', 'tact,neu,type,indsyn,minW,maxW,Is_pla,delay,src,effw,basew,W')
tact = []
neusec = []
migrations = []
neuintA = [0 for i in range(indLA[1] - indLA[0])]
stabilityA = [[] for i in range(indLA[1] - indLA[0])]

maxtact = 2000000

# It is guaranteed that all records are ordered by tact

with open(file, newline = '') as fil:
    bHeader = True
    csr = csv.reader(fil)
    for row in csr:
        if bHeader:
            bHeader = False
        elif row[0] == "secsta":
            sec.append(float(row[1]))
            CharTime.append(float(row[2]))
            ThrDecayT.append(float(row[3]))
            ThresholdExcessIncrementMultiplier.append(row[4])
            AbsRefT.append(float(row[5]))
            WINC.append(float(row[6]))
            RelWDec.append(float(row[7]))
            Inh.append(float(row[8]))
        elif row[0] == "secint":
            if float(row[2]) >= 0:
                secint.append(SectionIntensity(float(row[2]), float(row[1]), float(row[3])))
        elif row[0] == "lin":
            lin.append(Link(float(row[1]), int(row[2]), float(row[3]), float(row[4]), float(row[5]), float(row[6]), float(row[7]) != 0, float(row[8]), int(row[9]), float(row[10]), float(row[11]), float(row[12])))
        elif row[0] == "neu->sec":
            tac = int(row[1])
            if tac >= maxtact:
                break
            neu = int(row[2])
            s = int(row[3])

            # It is assumed that the 1st neu->sec section is ordered ny neu

            if len(tact) == 0 or tac != tact[-1]:
                tact.append(tac)
                neusec.append([] if tac == 0 else [-1 for i in range(len(neusec[0]))])
                migrations.append(0)
            if tac == 0:
                neusec[-1].append(s)
            else:
                neusec[-1][neu] = s
        elif row[0] == "neudyn":
            migrations[-1] += 1
        elif row[0] == "neu":
            neu = int(row[2])
            if indLA[0] <= neu < indLA[1]:
                neuintA[neu - indLA[0]] += int(row[3])
                stabilityA[neu - indLA[0]].append(float(row[13]))

nNeuronsperAction = int((indLA[1] - indLA[0]) / nActions)

ActionRepresentatives = []
for j in range(nActions):
    ActionRepresentatives.append(np.argsort(np.array(neuintA[j * nNeuronsperAction : (j + 1) * nNeuronsperAction]))[-10:] + j * nNeuronsperAction)

nSectionsperNetwork = len(sec)
print('total number of sections ', len(sec))
print('sections per network ', nSectionsperNetwork)
nSectionsperNetwork = int(nSectionsperNetwork)

SecInt = {}
fig, axs = plt.subplots(1, 1, figsize = (10, 3))
i = 0
for j in range(nSectionsperNetwork):
    s = sec[i]
    x = [t.tact for t in secint if t.sec == s]
    y = [t.relfre for t in secint if t.sec == s]
    SecInt[s] = statistics.mean(y)
    axs.plot(x, y, label = SectionNames[j], linewidth = 1)
    i += 1
leg = axs.legend(loc = 'best', ncol = 2, mode = "expand", shadow = True, fancybox = True)
leg.get_frame().set_alpha(0.5)
plt.title('Section relative activity')
plt.show()

fig, ax = plt.subplots(2, 1)
for j in range(nNeuronsperAction):
    ax[0].plot(stabilityA[j])
for j in range(nNeuronsperAction):
    ax[1].plot(stabilityA[nNeuronsperAction + j])
plt.title('Neuron stability')
plt.show()

print("mean section intensity:")
print(SecInt)

effw = [t.effw if t.type == 0 else t.W for t in lin]

i = 0
deffw_sum = {}
sat_sup_dist = []
tactNo = -1
guess_shift = -1
maxNLinksofThisType = {}

NLinksofThisType = {}
NSuppressedLinksofThisType = {}
NSaturatedLinksofThisType = {}

RecField = []
RecFieldA = []

print("Link data processing...")
cnt = 0
lastneu = -1
while i <= len(lin):
    end = i == len(lin)
    if i == 0 or end or lin[i].tact != lin[i - 1].tact:
        guess_shift = cnt
        cnt = 0
    if end or lin[i].Is_pla:
        if end or lin[i].neu != lastneu:
            for str, l in NLinksofThisType.items():
                if l > maxNLinksofThisType.setdefault(str, 0):
                    maxNLinksofThisType[str] = l
                if l > sat_sup_dist[-1].setdefault(str, [0, [], []])[0]:
                    sat_sup_dist[-1][str][0] = l
                sat_sup_dist[-1][str][1].append(NSuppressedLinksofThisType.get(str, 0))
                sat_sup_dist[-1][str][2].append(NSaturatedLinksofThisType.get(str, 0))
            if not end:
                lastneu = lin[i].neu
                NLinksofThisType.clear();
                NSuppressedLinksofThisType.clear();
                NSaturatedLinksofThisType.clear();
        if not end:
            if tactNo < 0 or lin[i].tact != tact[tactNo]:
                sat_sup_dist.append({})
                RecFieldA.append([[np.zeros(nSpatialZones), np.zeros(nSpatialZones), np.zeros(nVelocityZones), np.zeros(nVelocityZones), np.zeros(nSpatialZones), np.zeros((nRelPos, nRelPos))] for i in range(indLA[1] - indLA[0])])
                tactNo += 1
                print('Now tact %d is processed' % tact[tactNo])
            strsecfrom = "R%d" % (bisect.bisect_right(ReceptorSectionBoundaries, lin[i].src),) if lin[i].src >= 0 else "%d" % (neusec[tactNo][-1 - lin[i].src],)
            secto = neusec[tactNo][lin[i].neu]
            strLink = strsecfrom + '->' + "%d" % (secto,)
            if tactNo == 0 and deffw_sum.get(strLink) == None:
                deffw_sum[strLink] = []
            if tactNo > 0:
                if len(deffw_sum[strLink]) < tactNo:
                    deffw_sum[strLink].append(0)
                if i - guess_shift < 0:
                    guess_shift = i
                if lin[i - guess_shift].neu < lin[i].neu:
                    while lin[i - guess_shift].neu < lin[i].neu:
                        guess_shift -= 1
                elif lin[i - guess_shift].neu > lin[i].neu:
                    while lin[i - guess_shift].neu > lin[i].neu:
                        guess_shift += 1
                if lin[i - guess_shift].src > lin[i].src:
                    while lin[i - guess_shift].src > lin[i].src and lin[i - guess_shift].neu == lin[i].neu:
                        guess_shift -= 1
                elif lin[i - guess_shift].src < lin[i].src:
                    while lin[i - guess_shift].src < lin[i].src and lin[i - guess_shift].neu == lin[i].neu:
                        guess_shift += 1
                if lin[i - guess_shift].src == lin[i].src and lin[i - guess_shift].neu == lin[i].neu:
                    deffw_sum[strLink][-1] += abs(effw[i] - effw[i - guess_shift])
                else:
                    deffw_sum[strLink][-1] += abs(effw[i])

                if  indLA[0] <= lin[i].neu < indLA[1] and ReceptorSectionBoundaries[1] <= lin[i].src < ReceptorSectionBoundaries[2]:
                    ind = lin[i].src - ReceptorSectionBoundaries[1]
                    if ind < nSpatialZones:
                        RecFieldA[-1][lin[i].neu - indLA[0]][0][ind] = lin[i].W
                    ind -= nSpatialZones
                    if 0 <= ind < nSpatialZones:
                        RecFieldA[-1][lin[i].neu - indLA[0]][1][ind] = lin[i].W
                    ind -= nSpatialZones
                    if 0 <= ind < nVelocityZones:
                        RecFieldA[-1][lin[i].neu - indLA[0]][2][ind] = lin[i].W
                    ind -= nVelocityZones
                    if 0 <= ind < nVelocityZones:
                        RecFieldA[-1][lin[i].neu - indLA[0]][3][ind] = lin[i].W
                    ind -= nVelocityZones
                    if 0 <= ind < nSpatialZones:
                        RecFieldA[-1][lin[i].neu - indLA[0]][4][ind] = lin[i].W
                    ind -= nSpatialZones
                    if 0 <= ind:
                        RecFieldA[-1][lin[i].neu - indLA[0]][5][int(ind / nRelPos)][ind % nRelPos] = lin[i].W

            if strsecfrom[0] != 'R':
                strsecfrom = "%d" % (neusec[tactNo][-1 - lin[i].src] % nSectionsperNetwork,)
            secto = secto % nSectionsperNetwork
            strLink = strsecfrom + '->' + "%d" % (secto,)
            if lin[i].W < 0:
                NSuppressedLinksofThisType.setdefault(strLink, 0)
                NSuppressedLinksofThisType[strLink] += 1
            elif lin[i].W > 30:
                NSaturatedLinksofThisType.setdefault(strLink, 0)
                NSaturatedLinksofThisType[strLink] += 1
            NLinksofThisType.setdefault(strLink, 0)
            NLinksofThisType[strLink] += 1
    cnt += 1
    i += 1

rfmaxA = []

for rf in RecFieldA:
    rfm = []
    for j in range(nActions):
        z = rf[j * nNeuronsperAction : (j + 1) * nNeuronsperAction]
        y = [[a[k] for a in z] for k in range(len(z[0]))]
        ww = [np.amax(k, axis=0) for k in y]
        rfm.append(ww)
    rfmaxA.append(rfm)

print("Maximum number (in the whole history) of links between sections")
print(maxNLinksofThisType)

nPostsynapticNeuronsforLinkType = [0 for i in range(len(sat_sup_dist[0]))]
for dic in sat_sup_dist:
    i = 0
    for str, l in dic.items():
        nPostsynapticNeuronsforLinkType[i] = max(nPostsynapticNeuronsforLinkType[i], len(l[1]))
        i += 1

print("Maximum number (in the whole history) of postsynaptic neurons for each link type")
print(nPostsynapticNeuronsforLinkType)

fig, axs = plt.subplots(1, 1, figsize = (10, 3))
i = 0
for d in deffw_sum.keys():
    sectar = int(d.split('>')[1])
    if i * nSectionsperNetwork <= sectar < (i + 1) * nSectionsperNetwork:
        if max(deffw_sum[d]) > 0:
            axs.plot(tact[1:], deffw_sum[d], label = d, linewidth = 1)
        else:
            for s in sat_sup_dist:
                s.pop(d, 0)
leg = axs.legend(loc = 'best', ncol = 2, mode = "expand", shadow = True, fancybox = True)
leg.get_frame().set_alpha(0.5)
i += 1
plt.title('Weight change dynamics')
plt.show()

mpl.rc('font', size=10)
coo = np.arange(0, 30)
norm = mpl.colors.TwoSlopeNorm(vmin = -30, vcenter = 0, vmax = 30)

def draw_rec_fields(tac, reps, rec_fields):
    indtact = tact.index(tac)
    i = 0
    for j in range(len(reps)):
        for k in range(len(rec_fields[0][0]) - 1):
            axes_hist[j][k].clear()
        axes_hist[j][0].set_ylim(-30, 30)
        axes_hist[j][0].plot(rec_fields[indtact][reps[j]][0])
        if j == 0:
            axes_hist[j][0].set_title('x')
        axes_hist[j][1].set_ylim(-30, 30)
        axes_hist[j][1].plot(coo, rec_fields[indtact][reps[j]][1], coo, rec_fields[indtact][reps[j]][4])
        if j == 0:
            axes_hist[j][1].set_title('y')
        axes_hist[j][2].set_ylim(-30, 30)
        axes_hist[j][2].plot(rec_fields[indtact][reps[j]][2])
        if j == 0:
            axes_hist[j][2].set_title('vx')
        axes_hist[j][3].set_ylim(-30, 30)
        axes_hist[j][3].plot(rec_fields[indtact][reps[j]][3])
        if j == 0:
            axes_hist[j][3].set_title('vy')
        axes_hist[j][4].imshow(rec_fields[indtact][reps[j]][5], cmap = 'seismic', norm = norm)

curact = 'down'

for reps in ActionRepresentatives:

    fig, ax = plt.subplots()
    axes_hist = [[] for i in range(10)]
    i = 1
    for j in range(10):
        for k in range(len(RecFieldA[0][0]) - 1):
            axes_hist[j].append(plt.subplot(10, len(RecFieldA[0][0]) - 1, i))
            i += 1

    draw_rec_fields(0, reps, RecFieldA)

    axsli = plt.axes([0.25, 0.03, 0.65, 0.03])
    sli = Slider(axsli, 'tact', 0., tact[-1], valinit = 0, valstep = tact[1], valfmt = "%d")
    axsli.xaxis.set_visible(True)
    axsli.set_xticks(tact)

    def update_sli(val):
        tac = sli.val
        draw_rec_fields(tac, reps, RecFieldA)

    sli.on_changed(update_sli)

    plt.get_current_fig_manager().canvas.setWindowTitle('L receptive fields (action %s)' % (curact,))
    plt.show()
    curact = 'up'

fig, ax = plt.subplots()
axes_hist = [[] for i in range(nActions)]
i = 1
for j in range(nActions):
    for k in range(len(RecFieldA[0][0]) - 1):
        axes_hist[j].append(plt.subplot(nActions, len(RecFieldA[0][0]) - 1, i))
        i += 1
plt.title('L receptive fields')

def draw_rec_fields_all_actions(tac):
    indtact = tact.index(tac)
    i = 0
    for j in range(nActions):
        for k in range(len(RecFieldA[0][0]) - 1):
            axes_hist[j][k].clear()
        axes_hist[j][0].set_ylim(-30, 30)
        axes_hist[j][0].plot(rfmaxA[indtact][j][0])
        if j == 0:
            axes_hist[j][0].set_title('x')
        axes_hist[j][1].set_ylim(-30, 30)
        axes_hist[j][1].plot(coo, rfmaxA[indtact][j][1], coo, rfmaxA[indtact][j][4])
        if j == 0:
            axes_hist[j][1].set_title('y')
        axes_hist[j][2].set_ylim(-30, 30)
        axes_hist[j][2].plot(rfmaxA[indtact][j][2])
        if j == 0:
            axes_hist[j][2].set_title('vx')
        axes_hist[j][3].set_ylim(-30, 30)
        axes_hist[j][3].plot(rfmaxA[indtact][j][3])
        if j == 0:
            axes_hist[j][3].set_title('vy')
        axes_hist[j][4].imshow(rfmaxA[indtact][j][5], cmap = 'seismic', norm = norm)

draw_rec_fields_all_actions(0)

axsli = plt.axes([0.25, 0.03, 0.65, 0.03])
sli = Slider(axsli, 'tact', 0., tact[-1], valinit = 0, valstep = tact[1], valfmt = "%d")
axsli.xaxis.set_visible(True)
axsli.set_xticks(tact)

def update_sli_all_actions(val):
    tac = sli.val
    draw_rec_fields_all_actions(tac)

sli.on_changed(update_sli_all_actions)

plt.get_current_fig_manager().canvas.setWindowTitle('L - max weights for actions')
plt.show()
