import numpy as np
import copy
import sys

a = 3.556
L = 20*a

f = open('D.txt', 'r')
D = eval(f.readline())
f.close()

D_a = []
for i in range(121):
    D_a.append(D[i] if i<28 else 20/i)

def N_Xe_over_N_vac(N_Xe):
    return 0.364 + 1.2 * np.exp(-0.63 * N_Xe)

def C_vac(sizes):
    return np.sum(sizes / N_Xe_over_N_vac(sizes))

def R(N_Xe):
    return a * (3/(8*np.pi) * N_Xe/N_Xe_over_N_vac(N_Xe)) ** (1/3)

bubbles = []
steps = 100
C_sia_0 = 104.7
N_max = 120

matrix = R(np.array(tuple(range(1, N_max+1))))[:, np.newaxis] + R(np.array(tuple(range(1, N_max+1)))) + 0
six_Dt = 6 * 1/steps * np.array(D_a)

sizes = np.ones(127, dtype = 'i')
# py_C_0 = np.array(np.round(np.array([1557, 96, 111, 41, 26, 9, 9, 1, 2, 1, 1]) / 20), dtype = 'i')
# sizes = []
# for i in range(len(py_C_0)):
#     sizes += [1+i] * py_C_0[i]
# sizes = np.array(sizes)

bubbles = np.random.rand(3*len(sizes)) * L
bubbles = bubbles.reshape((len(sizes), 3))

C_sia = C_sia_0
C_vac_0 = C_vac(sizes)

for t in range(100 * steps + 1):
    N = bubbles.shape[0]

    new_bubbles = []
    new_sizes = []
    joined = np.zeros(N, dtype = 'bool')

    for i in range(N):
        for j in range(i+1, N):
            dr = bubbles[i] - bubbles[j]
            dr %= L
            for k in range(3):
                dr[k] = dr[k] if dr[k] <= L/2 else dr[k] - L
            dr = np.sum(dr**2) ** 0.5

            if not (joined[i] or joined[j]) and dr <= matrix[sizes[i]-1, sizes[j]-1]:
                # print( bubbles[i],  bubbles[j], dr, matrix[sizes[i]-1, sizes[j]-1])
                joined[i] = True
                joined[j] = True
                
                new_sizes.append(sizes[i] + sizes[j] if sizes[i] + sizes[j] <= N_max else 120)

                dr_bub = (bubbles[j]-bubbles[i])%L
                for k in range(3):
                    dr_bub[k] = dr_bub[k] if dr_bub[k] <= L/2 else dr_bub[k] - L
                new_bubbles.append(bubbles[i] + dr_bub * sizes[i]/new_sizes[-1])

                break
        if not joined[i]:
            new_sizes.append(sizes[i])
            new_bubbles.append(bubbles[i])

    bubbles = np.array(new_bubbles)
    sizes = np.array(np.round(new_sizes), dtype = 'i')
    C_sia = C_sia_0 + C_vac(sizes) - C_vac_0
    lamb = sqrt(six_Dt * (C_sia/L**3) / (250/(20*a)**3))
    N = bubbles.shape[0]
                
    ### move
    phi = np.random.rand(N) * 2 * np.pi
    cos_theta = (np.random.rand(N) - 0.5) * 2
    sin_theta = np.sqrt(1 - cos_theta**2)

    bubbles[:, 0] += lamb[sizes] * sin_theta * np.cos(phi)
    bubbles[:, 1] += lamb[sizes] * sin_theta * np.sin(phi)
    bubbles[:, 2] += lamb[sizes] * cos_theta

    ### print
    if t in np.array((1, 25, 50, 75, 100)) * steps:
        f = open(sys.argv[1], 'a')
        print(list(sizes), file = f)
        f.close()
    # print('xp_', t, ' = ', list(sizes), sep = '')

    # if (sizes == new_sizes).all():
    #     dr = bubbles - new_bubbles
    #     dr %= L
    #     for k in range(3):
    #         for l in range(len(dr)):
    #             dr[l, k] = dr[l, k] if dr[l, k] <= L/2 else dr[l, k] - L
    #     dr = np.sum(dr**2, axis = 1) ** 0.5
    #     print(dr)
    #     print('\n')