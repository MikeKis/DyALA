import csv
import matplotlib
import matplotlib.pyplot as plt

log = "9.log"
state = "ping_pong_state.csv"

effwtot = 0
tactrew = []
valrew = []
dwytot = [0 for i in range(30)]
with open(log, newline = '') as fil:
    csr = csv.reader(fil, delimiter = '\t')
    for row in csr:
        tact = int(row[0])
        type = row[1]
        if type == "SPIK":
            neu = int(row[2])
            if neu == 1:
                lastspike = tact
        # elif type == "EFFW":
        #     neu = int(row[2])
        #     val = int(row[4])
        #     if neu == 1:
        #         effwtot += val
        elif type == "REWA":
            val = int(row[2])
            neu = int(row[3])
            pre = row[5][5:]
            if neu == 1 and pre[0] == 'y':
                y = int(pre[1:])
                dwytot[y] += val
                if y == 29 and val > 0:
                    tactrew.append(lastspike)
                    valrew.append(val)

print(effwtot)
print(dwytot)

tac = []
zone = []
y = []
vy = []
ry = []
indcur = 0
with open(state, newline = '') as fil:
    csr = csv.reader(fil)
    for row in csr:
        tact = int(row[0])
        if indcur < len(tactrew) and tact >= tactrew[indcur]:
            tac.append(tact)
            zone.append(int(row[1]))
            y.append(float(row[3]))
            vy.append(float(row[5]))
            ry.append(float(row[6]))
            indcur += 1

# cmap = matplotlib.cm.plasma
fig, ax = plt.subplots()
ax.scatter(y, ry, c=valrew)
sm = plt.cm.ScalarMappable(norm=plt.Normalize(vmin = min(valrew), vmax=max(valrew)))
fig.colorbar(sm)
ax.set_xlabel("ball y")
ax.set_ylabel("racket y")
plt.show()

