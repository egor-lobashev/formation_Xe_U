import numpy as np
import sys

def only_good_timesteps(trajectory, timestep_min, timestep_max):
    for step in range(len(trajectory)-1,-1,-1):
        if trajectory[step][0] < timestep_min or trajectory[step][0] > timestep_max:
            trajectory.pop(step)

concentrations = dict()

for filename in sys.argv[1:]:
    f = open(filename, 'r')
    d = eval(f.readline())
    f.close()
    
    for key in d.keys():
        only_good_timesteps(d[key][0], 0, 100)
        for tr_step in d[key][0]:
            if (tr_step[0], key[0]) not in concentrations.keys():
                concentrations[(tr_step[0], key[0])] = 0
            concentrations[(tr_step[0], key[0])] += 1
print(concentrations)