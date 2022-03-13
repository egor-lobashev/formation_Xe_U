import numpy as np
import sys

def unwrap_bubble_trajectory(trajectory, cell_size):
    for i in range(len(trajectory) - 1):
        dr = np.array(trajectory[i+1][1:]) - np.array(trajectory[i][1:])
        dr = dr % cell_size
        dr = np.min((dr, abs(dr - cell_size)), axis = 0) * np.sign(cell_size/2 - dr)
        for j in range(3):
            trajectory[i+1][j+1] = trajectory[i][j+1] + dr[j]
        
def dr2(r1, r2):
    dr = np.array(r2) - np.array(r1)
    return np.dot(dr, dr)

def plot_data_for_a_trajectory(trajectory, plot_data, C_sia):
    for i in range(1, len(trajectory)):
        if len(plot_data) <= i:
            plot_data.append([])
        for j in range(len(trajectory) - i):
            plot_data[i].append(dr2(trajectory[j][1:], trajectory[j+i][1:]) / C_sia[j])

def only_good_timesteps(trajectory, timestep_min, timestep_max):
    for step in range(len(trajectory)-1,-1,-1):
        if trajectory[step][0] < timestep_min or trajectory[step][0] > timestep_max:
            trajectory.pop(step)

# def print_length_distribution():
#     distribution = dict()

#   for filename in sys.argv[1:]:
#         f = open(filename, 'r')
#         superdict = eval(f.readline())
#         f.close()

#         for key in superdict.keys():
#             if superdict[key][0][-1][0] == 100:
#                 continue
#             if key[0] not in distribution.keys():
#                 distribution[key[0]] = []
#             distribution[key[0]].append(len(superdict[key][0]))
#     for k in distribution.keys():
#         print(k, np.mean(distribution[k]), len(distribution[k]), sep = ', ')

def D(bubble_size, timestep_min, timestep_max, cell_size, C_sia, only_single_steps=False):
    trajectories_with_size = []

    for filename in sys.argv[1:]:
        f = open(filename, 'r')
        superdict = eval(f.readline())
        f.close()
        
        for key in superdict.keys():
            if key[0] == bubble_size and len(superdict[key][0]) > 1:
                trajectory = superdict[key][0]
                if trajectory[-1][0] - trajectory[0][0] != len(trajectory) - 1:
                    continue
                trajectories_with_size.append(trajectory)

    plot_data = [[0]]

    for trajectory in trajectories_with_size:
        only_good_timesteps(trajectory, timestep_min, timestep_max)
        if len(trajectory) == 0:
            continue
        unwrap_bubble_trajectory(trajectory, cell_size)
        plot_data_for_a_trajectory(trajectory, plot_data, C_sia)


    # for i in range(1, len(plot_data)-1):
    #     data = plot_data[i]
    #     print(i, np.mean(data), np.std(data)/(len(plot_data[i])/i)**0.5, sep = ', ')
    ###### print(np.max(plot_data[1])**0.5)
    if len(plot_data) == 1:
        return '?'
    else:
        # return np.mean(plot_data[1])/6, np.std(plot_data[1])/6/len(plot_data[1])**0.5
        return list(map(np.mean, (plot_data))) #np.std(plot_data[1])/6/len(plot_data[1])**0.5

    """
    # total_data = 0
    dr2_over_dt = []

    max_length_to_consider = min(1, len(plot_data)-1)  if only_single_steps else len(plot_data)-1
    for i in range(1, max_length_to_consider + 1):
        # total_data += len(plot_data[i])
        dr2_over_dt.extend(plot_data[i])
        # print(i, np.mean(plot_data[i]), len(plot_data[i]), sep = ', ')
    

    if len(dr2_over_dt) == 0:
        return '?'
    
    # dr2_over_dt /= total_data

    print(np.mean(dr2_over_dt)/6, np.std(dr2_over_dt)/6/len(dr2_over_dt), sep = ', ')
    return np.mean(dr2_over_dt)/6
    """

# print_length_distribution()

timestep_min = 0 # 40
timestep_max = 1e10

f = open('C_sia_approx.txt', 'r')
C_sia = eval(f.readline())
f.close()

C_sia = [1]*1000

a = 3.556
cr_cells = 10
cell_size = a * cr_cells

# for bubble_size in range(1, 120):
#     D_1 = D(bubble_size, timestep_min, timestep_max, cell_size, C_sia, True)
#     if D_1 != '?':
#         print('', D_1[0], D_1[1], sep = ' ')

print(D(10, timestep_min, timestep_max, cell_size, C_sia, True))
# \left(x_{1},y_{1}+z_{1}t\right)
