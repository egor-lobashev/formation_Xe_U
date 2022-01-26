from ovito.io import *
from ovito.modifiers import *
from ovito.pipeline import *
from ovito.io import export_file
import numpy as np
import copy
import sys

filename = sys.argv[1]
pipeline = import_file(filename, sort_particles = True, columns =
    ["Particle Identifier", "Particle Type", "Position.X", "Position.Y", "Position.Z"])

def set_modifiers(pipeline):
    wrap = WrapPeriodicImagesModifier()
    not_Xe = ExpressionSelectionModifier(expression = 'ParticleType != 3')
    only_Xe = DeleteSelectedModifier()
    clusters = ClusterAnalysisModifier()
    clusters.sort_by_size = True
    clusters.unwrap_particles = True
    clusters.cutoff = 4.5
    this_cluster = ExpressionSelectionModifier()
    invert = InvertSelectionModifier()
    delete = DeleteSelectedModifier()

    pipeline.modifiers.append(wrap) # 0
    pipeline.modifiers.append(not_Xe) # 1
    pipeline.modifiers.append(only_Xe) # 2
    pipeline.modifiers.append(clusters) # 3
    pipeline.modifiers.append(this_cluster) # 4
    pipeline.modifiers.append(invert) # 5
    pipeline.modifiers.append(delete) # 6

    pipeline.modifiers[5].enabled = False
    pipeline.modifiers[6].enabled = False

def number_of_Xe(pipeline):
    before = pipeline.modifiers[4].expression
    pipeline.modifiers[4].expression = 'ParticleType == 3'

    ans = pipeline.compute(0).attributes['ExpressionSelection.count.2']

    pipeline.modifiers[4].expression = before
    return ans

set_modifiers(pipeline)
N_Xe = number_of_Xe(pipeline)
N_atoms = pipeline.compute(0).particles.count
extra_atoms_0 = N_atoms - 16000
a = pipeline.compute(0).cell[0][0] / 20

#####################
def center_of_mass(pipeline, timestep):
    before_first = pipeline.modifiers[5].enabled
    before_second = pipeline.modifiers[6].enabled

    pipeline.modifiers[5].enabled = True
    pipeline.modifiers[6].enabled = True

    positions = pipeline.compute(timestep).particles.positions
    com = np.mean(positions, axis = 0)

    pipeline.modifiers[5].enabled = before_first
    pipeline.modifiers[6].enabled = before_second

    return com

superdict = dict()
for timestep in range(0, pipeline.source.num_frames):
    data = pipeline.compute(timestep)
    clusters_list = list(data.tables['clusters'].xy()[:,1])

    for cluster in range(1, len(clusters_list) + 1):
        pipeline.modifiers[4].expression = 'Cluster == ' + str(cluster)
        selection = pipeline.compute(timestep).particles.selection
        selected = np.where(selection == 1)[0]

        key = (len(selected), selected[0])
        if key not in superdict.keys():
            superdict[key] = ([], list(selected))

        com = center_of_mass(pipeline, timestep)

        superdict[key][0].append([timestep] + list(com))
print(superdict)
