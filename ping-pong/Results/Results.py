from readmon import readmon

experiment = 10000

file = "monitoring.%d.csv" % (experiment,)
fileprotxt = "spikes.%d.txt" % (experiment,)

SaturatedResourceThreshold = 30

mon = readmon(file)

sec = mon.sec
CharTime = mon.CharTime
ThresholdExcessIncrementMultiplier = mon.ThresholdExcessIncrementMultiplier
AbsRefT = mon.AbsRefT
WINC = mon.WINC
Inh = mon.Inh
minW = mon.minW
maxW = mon.maxW
SectionName = mon.SectionName
secint = mon.secint
lin = mon.lin
tact = mon.tact
neusec = mon.neusec
migrations = mon.migrations
ReceptorBounds = mon.ReceptorBounds
neus = mon.Meaning

linlast = [[a.neu, a.src, a.effw] for a in lin if a.tact == lin[-1].tact and a.type == 3]
neulast = [[a.neu, a.meaning] for a in neus if a.tact == lin[-1].tact]

meaning = {}
for i in neulast:
    meaning[i[0]] = i[1]

with open("links_fin.csv", "wt") as filout:
    for i in neulast:
        for j in linlast:
            if j[0] == i[0]:
                s = i[1] + (",%f," % (j[2],)) + (meaning[j[1] - 1] if j[1] > 0 else "R%d" % (-j[1] - 1,)) + '\n'
                filout.write(s)
