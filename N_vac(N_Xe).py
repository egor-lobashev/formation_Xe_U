import numpy as np

f = open('sia_Xe_vacs', 'r')

C_min = 0
C_max = 100500

vacs_on_Xe = dict()
for s in f:
    if s == '\n':
        continue

    l = list(map(int, s[:-1].split(' ')))
    if l[1] == 0:
        continue
    
    C = l[0]
    if (C >= C_min) and (C <= C_max):
        if l[1] in vacs_on_Xe.keys():
            vacs_on_Xe[l[1]].append(l[2])
        else:
            vacs_on_Xe[l[1]] = [l[2]]

def R(vacs):
    a = 3.556
    v_0 = a**3 / 2
    return (3*vacs*v_0 /4/np.pi) ** (1/3)

print('Xe_all vX_all')
for k in vacs_on_Xe.keys():
    if len(vacs_on_Xe[k]) != 1:
        # print(k, np.mean(vacs_on_Xe[k]), np.std(vacs_on_Xe[k]) / (len(vacs_on_Xe[k]) - 1)**0.5, sep = ', ')
        print(k, np.mean(vacs_on_Xe[k])/k, sep = ' ')
