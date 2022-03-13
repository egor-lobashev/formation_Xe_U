import numpy as np

f = open('t_size', 'r')

t_min = 99
t_max = 99
ans = dict()
for s in f:
    if s == '\n':
        continue

    l = list(map(int, s[:-1].split(' ')))
    if l[1] == 0:
        continue
    
    t = l[0]
    if (t >= t_min) and (t <= t_max):
        if l[1] in ans.keys():
            ans[l[1]] += 1
        else:
            ans[l[1]] = 1

f.close()

for k in ans.keys():
    # print(', '.join(map(str, [k] * ans[k])), end = ', ')
    print('[', k, ']*', ans[k], sep = '', end = '+')
print()