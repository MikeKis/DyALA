import csv

logfile = '9.log'
nactions = 2

previousXvel = 0
previousX = 0
successes = []
failures = []
tactLEV0entry = -1
with open('ping_pong_state.csv', newline = '') as fil:
    csr = csv.reader(fil)
    for row in csr:
        currentXvel = float(row[4])
        currentX = float(row[2])
        if currentXvel == -previousXvel and currentX < -0.4:
            successes.append([int(row[0]), tactLEV0entry])
            tactLEV0entry = -1
        elif currentX > previousX + 0.4:
            failures.append([int(row[0]), tactLEV0entry])
            tactLEV0entry = -1
        previousXvel = currentXvel
        previousX = currentX
        if row[1] == '0' and tactLEV0entry == -1:
            tactLEV0entry = int(row[0])

last_actions = []
actions = []
suppressed_actions = []

with open(logfile, newline = '') as fil:
    csr = csv.reader(fil, delimiter='\t')
    for row in csr:
        cur_tact = int(row[0])
        if row[1] == "SPIK":
            if row[3][:2] == "LA":
                last_actions.append([cur_tact, int(row[3][2:3])]) 
            elif row[3][:12] == "FINALGATEACT":
                action = int(row[3][12:13])
                last_actions.remove([cur_tact - 1, action])
                actions.append([cur_tact, action])
        suppressed = list(filter(lambda x: x[0] < cur_tact - 1, last_actions))
        for a in suppressed:
            suppressed_actions.append(a)
            last_actions.remove(a)

suc = []
for a in successes:
    if a[1] >= 0:
        nact = sum(True for i in actions if a[1] <= i[0] < a[0])
        nsup = sum(True for i in suppressed_actions if a[1] <= i[0] < a[0])
        suc.append([a[0], nact, nsup])

fai = []
for a in failures:
    if a[1] >= 0:
        nact = sum(True for i in actions if a[1] <= i[0] < a[0])
        nsup = sum(True for i in suppressed_actions if a[1] <= i[0] < a[0])
        fai.append([a[0], nact, nsup])

print(fai)