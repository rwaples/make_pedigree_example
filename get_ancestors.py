import pandas as pd
import numpy as np
import graph_tool.all as gt
import gzip


pedigree_path = str(snakemake.input.pedigree)
out_path = str(snakemake.output.batch)

max_gen_updown = int(snakemake.params.max_gen_updown)
n = int(snakemake.wildcards.n)
nsplit = int(snakemake.params.nsplit)

# read in the pedigree
# pedigree should include lines for each ind listed as parent
# missing parents have id of -1 
# we need three columns, id of the target ind, and the ids of the two parents
pedigree = pd.read_hdf(pedigree_path)
pedigree = pedigree[['pid', 'parent1', 'parent2']].copy()

idx_of_id = dict(zip(
    pedigree['pid'].values,
    pedigree.index.values
    )
)
idx_of_id['0'] = -1
idx_of_id[0] = -1
idx_of_id[-1] = -1

# graph nodes with be indexed / label from 0 and counting up
# we assign each ind an index, and assign them to node with that index
parent1_of_idx = dict(zip(pedigree.index.values,
                       [idx_of_id[x] for x in pedigree['parent1']]))
parent2_of_idx = dict(zip(pedigree.index.values,
                       [idx_of_id[x] for x in pedigree['parent2']]))

pedigree['parent1_idx'] =  [idx_of_id[x] for x in pedigree['parent1']]
pedigree['parent2_idx'] =  [idx_of_id[x] for x in pedigree['parent2']]


# edges from child -> parent
z = [(x.Index, x.parent1_idx) for x in pedigree.itertuples() if x.parent1_idx >= 0]
y = [(x.Index, x.parent2_idx) for x in pedigree.itertuples() if x.parent2_idx >= 0]

# make the directed graph with edges going from the child to the parent
# up to two edges for each row in the pedigree
graph = gt.Graph(y+z, directed=True)

# here I have split up for multiprocessing, but not really necessary 
splits = np.array_split(np.arange(graph.num_vertices()), nsplit)

anc = list()
# for each node (ind)
for i in splits[n]:
	# distance to all other nodes
    dist_map = graph.new_vp("int16_t", val=np.iinfo(np.int16).max)
    dist = gt.shortest_distance(
        graph,
        source=i,
        target = None,
        dist_map = dist_map,
        return_reached = False
    )
    arr = dist.a
	
	# record distances above 0 and below the threshold (e.g. 12)
    ancestors = np.where(
        np.logical_and(
            arr < max_gen_updown, arr > 0
        )
    )[0]
    
    for a in ancestors:
        anc.append((i, a, arr[a]))

# output is one row per child-ancestor relationship, also record the distance in generations
anc = pd.DataFrame(anc)
if len(anc)>0:
    anc.columns = ['child', 'ancestor', 'ngen']
    anc.to_csv(out_path, sep = '\t', index=None, header=None)
else:
    with gzip.open(out_path, "w") as fh:
        pass
