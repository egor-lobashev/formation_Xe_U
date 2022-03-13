from ovito.io import *
from ovito.modifiers import *
from ovito.pipeline import *
from ovito.io import export_file
import numpy as np
import copy
import sys

def set_modifiers(pipeline):
    com_shift = AffineTransformationModifier()
    com_shift.operate_on.remove('cell')
    wrap = WrapPeriodicImagesModifier()
    WS = WignerSeitzAnalysisModifier()
    WS.reference = FileSource()
    WS.reference.load('initial_20.atom')
    WS.affine_mapping = ReferenceConfigurationModifier.AffineMapping.ToReference
    WS.per_type_occupancies = True
    not_void = ExpressionSelectionModifier(expression = 'Occupancy.1 >= 1')
    highlight_void = DeleteSelectedModifier()
    clusters = ClusterAnalysisModifier()
    clusters.sort_by_size = True
    clusters.cutoff = 4
    this_cluster = ExpressionSelectionModifier()

    pipeline.modifiers.append(com_shift) # 0
    pipeline.modifiers.append(wrap) # 1
    pipeline.modifiers.append(WS) # 2
    pipeline.modifiers.append(not_void) # 3
    pipeline.modifiers.append(highlight_void) # 4
    pipeline.modifiers.append(clusters) # 5
    pipeline.modifiers.append(this_cluster) # 6

def number_of_U_and_Xe(pipeline):
    first_before = pipeline.modifiers[2].output_displaced
    second_before = pipeline.modifiers[6].expression
    third_before = pipeline.modifiers[3].expression
    pipeline.modifiers[3].expression = '1==0'

    pipeline.modifiers[2].output_displaced = True
    pipeline.modifiers[6].expression = 'ParticleType == 3'
    N_Xe = pipeline.compute(0).attributes['ExpressionSelection.count.2']

    pipeline.modifiers[6].expression = 'ParticleType == 1'
    N_U = pipeline.compute(0).attributes['ExpressionSelection.count.2']

    pipeline.modifiers[2].output_displaced = first_before
    pipeline.modifiers[6].expression = second_before
    pipeline.modifiers[3].expression = third_before

    return N_U, N_Xe

def set_parameters(filename):
    pipeline = import_file(filename, columns =
        ["Particle Identifier", "Particle Type", "Position.X", "Position.Y", "Position.Z"])

    set_modifiers(pipeline)
    N_U, N_Xe = number_of_U_and_Xe(pipeline)
    N_atoms = N_U + N_Xe
    sia_0 = N_atoms - 16000
    a = 3.554# pipeline.compute(0).cell[0][0] / 20

    return pipeline, a, N_atoms, N_Xe, sia_0

def com_displacement(atoms, a):
    m =  np.array(atoms.particle_types)*np.ones((3, 1)) == 3
    pos = np.ma.array(atoms['Position'], mask = m)

    pos1 = pos % (a/2)
    pos2 = (pos + a/4) % (a/2)

    com1 = pos1.mean(axis = 0)
    com2 = pos2.mean(axis = 0) - a/4

    std1 = (pos1**2).mean(axis = 0) - pos1.mean(axis = 0)**2
    std2 = (pos2**2).mean(axis = 0) - pos2.mean(axis = 0)**2

    com = [0]*3
    for j in range(3):
        com[j] = com2[j] if std1[j] > std2[j] else com1[j]
    
    return com

def set_shift(pipeline, timestep):
    pipeline.modifiers[2].output_displaced = True
    pipeline.modifiers[4].enabled = False
    pipeline.modifiers[0].enabled = False

    shift = com_displacement(pipeline.compute(timestep).particles, a)

    pipeline.modifiers[2].output_displaced = False
    pipeline.modifiers[4].enabled = True
    pipeline.modifiers[0].enabled = True

    vacs_min = 100500
    i1i2 = (0, 0)
    for i1 in (0, 1):
        for i2 in (0, 1):
            pipeline.modifiers[0].transformation = [[1, 0, 0, -shift[0] + (a/2)*i1],
                                                    [0, 1, 0, -shift[1] + (a/2)*i2],
                                                    [0, 0, 1, -shift[2]]]

            data = pipeline.compute(timestep)
            vacs = data.attributes['WignerSeitz.vacancy_count']
            if vacs < vacs_min:
                vacs_min = vacs
                i1i2 = (i1, i2)

    i1, i2 = i1i2
    pipeline.modifiers[0].transformation = [[1, 0, 0, -shift[0] + (a/2)*i1],
                                            [0, 1, 0, -shift[1] + (a/2)*i2],
                                            [0, 0, 1, -shift[2]]]
    # print(timestep, pipeline.modifiers[0].transformation[:,3], sep = ', ')

def t_size(pipeline):
    ans = ''
    for timestep in range(0, pipeline.source.num_frames):
        if timestep % 30 == 0:
            set_shift(pipeline, timestep)

        pipeline.modifiers[2].output_displaced = True
        data = pipeline.compute(timestep)
        clusters_list = list(data.tables['clusters'].xy()[:,1])

        for cluster in range(1, len(clusters_list) + 1):
            pipeline.modifiers[6].expression = 'Cluster == ' + str(cluster)
            data = pipeline.compute(timestep)
            N_Xe = data.attributes['ExpressionSelection.count.2']

            ans += str(timestep) + ' ' + str(N_Xe) + '\n'
    return ans

def Sia(pipeline, timestep):
    data = pipeline.compute(timestep)
    clusters_list = list(data.tables['clusters'].xy()[:,1])
    
    return N_atoms - (16000 + N_Xe) + np.sum(clusters_list)
    # return N_atoms - (16000 + N_Xe) + np.sum(clusters_list)

def Singles(pipeline, timestep):
    data = pipeline.compute(timestep)
    clusters_list = list(data.tables['clusters'].xy()[:,1])

    return clusters_list.count(1)

def Singles_Xe(pipeline, timestep, singles):
    data = pipeline.compute(timestep)
    clusters_list = list(data.tables['clusters'].xy()[:,1])

    pipeline.modifiers[6].expression = \
        'Occupancy.3 == 1 && Cluster >= ' + str(len(clusters_list) - singles)
    data = pipeline.compute(timestep)
    return data.attributes['ExpressionSelection.count.2']

def Vacs_in_cluster(pipeline, timestep, cluster):
    pipeline.modifiers[6].expression = 'Cluster == ' + str(cluster)
    data = pipeline.compute(timestep)
    return data.attributes['ExpressionSelection.count.2']

def Xenon_in_cluster(pipeline, timestep, cluster):
    pipeline.modifiers[6].expression = 'Occupancy.3 == 1 && Cluster == ' + str(cluster)
    data = pipeline.compute(timestep)
    return data.attributes['ExpressionSelection.count.2']

def print__sia_Xe_vacs(pipeline, timestep):
    data = pipeline.compute(timestep)
    clusters_list = list(data.tables['clusters'].xy()[:,1])
    ans = ''

    for cluster in range(1, len(clusters_list) + 1):
        vacs_in_cluster = Vacs_in_cluster(pipeline, timestep, cluster)
        Xe_in_cluster = Xenon_in_cluster(pipeline, timestep, cluster)

        ans += ' '.join(list(map(str, (sia, Xe_in_cluster, vacs_in_cluster)))) + '\n'
    print(ans)

def Sia_last(pipeline):
    last_timestep = pipeline.source.num_frames - 1
    set_shift(pipeline, last_timestep)

    sia_last = Sia(pipeline, last_timestep)
    return sia_last

ave_sia = []

for filename in sys.argv[1:]:
    print(filename)
    pipeline, a, N_atoms, N_Xe, sia_0 = set_parameters(filename)
    # print(t_size(pipeline))

    if len(ave_sia) == 0:
        ave_sia = [0]*pipeline.source.num_frames
    shift_period = 3

    sia_last = Sia_last(pipeline)

    set_shift(pipeline, 0)
    for timestep in range(0, pipeline.source.num_frames):
        if timestep % shift_period == shift_period//2:
            set_shift(pipeline, timestep)

        singles = Singles(pipeline, timestep)
        # singles_Xe = Singles_Xe(pipeline, timestep, singles)
        sia = Sia(pipeline, timestep)
        ave_sia[timestep] += sia

        # print__sia_Xe_vacs(pipeline, timestep)

print('t_' + sys.argv[1][-10:-5] + ', v_' + sys.argv[1][-10:-5])
for i in range(101):
    print(i, ave_sia[i]/(len(sys.argv) - 1), sep = ' ')

### (104.714285714, 118.428571429, 128.380952381, 138.428571429, 146.619047619, 155.952380952, 162.80952381, 165.80952381, 172.666666667, 177.857142857, 178.904761905, 184.857142857, 185.476190476, 191.047619048, 193.619047619, 196.0, 202.476190476, 202.666666667, 203.571428571, 208.571428571, 211.380952381, 214.666666667, 212.761904762, 215.761904762, 215.619047619, 220.285714286, 221.952380952, 223.952380952, 225.523809524, 227.80952381, 218.666666667, 221.238095238, 224.428571429, 224.619047619, 226.19047619, 225.333333333, 228.142857143, 227.380952381, 227.428571429, 229.476190476, 233.238095238, 231.476190476, 233.238095238, 236.238095238, 234.285714286, 237.428571429, 238.666666667, 241.714285714, 238.380952381, 242.571428571, 239.380952381, 239.952380952, 243.095238095, 243.666666667, 247.0, 246.714285714, 246.428571429, 247.952380952, 247.095238095, 250.476190476, 238.904761905, 237.238095238, 237.80952381, 243.19047619, 244.142857143, 240.904761905, 242.428571429, 244.428571429, 243.619047619, 243.714285714, 246.047619048, 247.0, 247.523809524, 245.0, 248.476190476, 250.761904762, 254.904761905, 255.952380952, 256.047619048, 256.571428571, 259.333333333, 257.904761905, 262.0, 268.285714286, 263.095238095, 265.857142857, 272.80952381, 276.285714286, 281.0, 282.857142857, 249.476190476, 249.333333333, 247.333333333, 245.857142857, 250.238095238, 249.904761905, 251.047619048, 249.523809524, 251.571428571, 250.80952381, 252.80952381)
### 104.3-2.394*x+28*x**0.5732

# print(ave_sia)
# print(np.mean(ave_sia))