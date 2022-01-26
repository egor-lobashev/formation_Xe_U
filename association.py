import numpy as np
import sys
import matplotlib.pyplot as plt
from itertools import combinations

a = 3.556
L = 20*a

f = open('C.txt', 'r')
C = eval(f.readline())
f.close()

f = open('D_over_C_sia.txt', 'r')
D_over_C_sia_pure = eval(f.readline())
f.close()

D_over_C_sia = []
for n in range(121):
    # D_a.append(D[n] if n<1 else 20/n)
    D_over_C_sia.append(D_over_C_sia_pure[n] if n<30 else 58000/n) # 20/n

f = open('sia.txt', 'r')
sia = eval(f.readline())
f.close()

def distance_between_bubbles(step_1, step_2):
    ans = 0
    for i in range(1, 4):
        delta_xi = (step_1[i] - step_2[i]) % L
        ans += min((delta_xi - L) ** 2, delta_xi ** 2)
    return ans**0.5

def association_radii(superdict, radii_will_associate, radii_all, n_min, n_max):
    for key in superdict.keys():
        if key[0] > n_max or key[0] < n_min:
            continue

        last_step = superdict[key][0][-1][0]
        associated_with = {key[1]}

        for v in superdict.values():
            if v[0][0][0] == last_step + 1 and \
                    key[1] in v[1]:
                associated_with = set(v[1])
                break
        
        for v in superdict.values():
            # if v[1][0] in associated_with:
            #     continue
            associated = (v[1][0] in associated_with)
            size = len(v[1])
            
            for step_1 in v[0]:
                for step_2 in superdict[key][0]:
                    if step_1[0] != step_2[0]:
                        continue

                    if size not in radii_all.keys():
                        radii_all[size] = []
                        
                    radii_all[size].append(distance_between_bubbles(step_1, step_2))

                    if v[1][0] in associated_with:
                        if size not in radii_will_associate.keys():
                            radii_will_associate[size] = []
                            
                        radii_will_associate[size].append(distance_between_bubbles(step_1, step_2))

def association_count_calculate(superdict, association_count):
    for k in superdict.keys():
        last_timestep = superdict[k][0][-1][0]

        if last_timestep == 100:
            continue

        IDs_after_association = set()

        for k1 in superdict.keys():
            if superdict[k1][0][0][0] == last_timestep + 1 and k[1] in superdict[k1][1]:
                IDs_after_association = set(superdict[k1][1])
                break

        association_members = []

        for k1 in superdict.keys():
            if (k1[1] in IDs_after_association) and (superdict[k1][0][-1][0] == last_timestep):
                association_members.append(k1[0])
        
        association_members.sort()
        association_members = list(association_members)

        N = len(association_members)
        new_key = tuple(association_members + [last_timestep])

        if N == 2:
            if new_key not in association_count.keys():
                association_count[new_key] = 0
            
            association_count[new_key] += 0.5

        elif N > 2:
            # comb = np.math.factorial(N) / np.math.factorial(N-2) / 2

            for c in combinations(association_members, 2):
                new_key = (c[0], c[1], last_timestep)

                if new_key not in association_count.keys():
                    association_count[new_key] = 0
                
                association_count[new_key] += 1/(N-1)/N

def calculate_all_start_data(association_count, n_min=0, n_max=0):
    for filename in sys.argv[1:]:
        if filename == '-':
            break
            
        f = open(filename, 'r')
        superdict = eval(f.readline())
        f.close()
        
        association_count_calculate(superdict, association_count)

        if n_max != 0:
            association_radii(superdict, radii_will_associate, radii_all, n_min, n_max)

def interesting_associations(range1, range2, time_range):
    ans = 0
    for n in range(range1[0], range1[1] + 1):
        for m in range(range2[0], range2[1] + 1):
            for t in range(time_range[0], time_range[1] + 1):
                key = (min(n, m), max(n, m), t)
                if key in association_count.keys():
                    ans += association_count[key]

    return ans

def equal_or_not_intersecting(range1, range2): # not tested
    set1 = set(range(range1[0], range1[1] + 1))
    set2 = set(range(range2[0], range2[1] + 1))
    return len(set1 | set2) in (len(set1) + len(set2), len(set1))

def Cn_Cm_Cs_list(range1, range2):
    if not equal_or_not_intersecting(range1, range2):
        return 'ranges are not equal or intersecting'

    ans = []

    for t in range(0, 101):
        ans.append(0)
        for n in range(range1[0], range1[1]+1):
            if (t, n) not in C.keys():
                continue

            for m in range(range2[0], range2[1]+1):
                if (t, m) not in C.keys():
                    continue
                
                ans[t] += C[(t,n)]*C[(t,m)]
        ans[t] *= sia[t]
    
    return ans

def factor_before_Cn_Cm_Cs(range1, range2):
    Cn_Cm_Cs = Cn_Cm_Cs_list(range1, range2)
    Cn_Cm_Cs_integral = np.sum(Cn_Cm_Cs)

    last_C = interesting_associations(range1, range2, [0, 100])

    f = last_C / Cn_Cm_Cs_integral
    return f, last_C

def print____t__associations__Cn_Cm_Cs____integrals_from_0_to_t(range1, range2):
    Cn_Cm_Cs = Cn_Cm_Cs_list(range1, range2)

    for t in range(101):
        time_range = [0, t]
        print(t, interesting_associations(range1, range2, time_range),
            np.sum(Cn_Cm_Cs[t]) * f, sep = ', ')

def R(N_Xe):
    return 2.44*N_Xe**(1/3) - 0.74*np.exp(-0.18*N_Xe)
# 2.44*x^{1/3} - 0.74*exp(-0.18*x)

def R_plus(n):
    return 2

def f_theor(n, m):
    return 4*np.pi * (D_over_C_sia[n] + D_over_C_sia[m]) * (R(n) + R(m) + R_plus(n) + R_plus(m))
    # return 4*np.pi * (D_a[n] + D_a[m]) / (250/L**3) * (R(n) + R(m) + R_plus(n) + R_plus(m))

def print_matix():
    print('factor = [', end = '')
    for n in range(1, 121):
        for m in range(1, 121):
            if n == 1 and m == 1:
                f = 1191345.54999
            elif n == m:
                f /= 2
            else:
                f = f_theor(n, m)

            # if n > 1 or m > 1:
            #     f *= 1.2
            # if n > 5 or m > 5:
            #     f *= 1.5
            # if n > 8 or m > 8:
            #     f *= 2.2/1.5
            print(f, end = ' ')
        print(';', end = '')
    print(']', end = '')

def print__f__and__f_theor(n1, n2_range):
    for n in n2_range:
        f_t = f_theor(n1, n)

        f, N = factor_before_Cn_Cm_Cs([n1, n1], [n, n])
        f *= L**6 * 20
        print(n, f, f_t, sep = ', ')

radii_will_associate = dict()
radii_all = dict()
association_count = dict()

calculate_all_start_data(association_count, 1, 1)

# print__f__and__f_theor(1, range(1, 10))
print_matix()