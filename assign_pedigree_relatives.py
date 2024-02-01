import pandas as pd
import numpy as np


# input
matched_ancestors = str(snakemake.input.matched_ancestors)
# output
assigned_relatives = str(snakemake.output.assigned_relatives)
relative_table = str(snakemake.output.relative_table)


dtypes = {
    'descendant_x': np.int32,
    'descendant_y': np.int32,
    'gsep_x': np.int8,
    'gsep_y': np.int8,
    'ancestor': np.int32,
 }
matched_ancestors = pd.read_csv(matched_ancestors, sep = '\t', dtype=dtypes)

# find most recent common ancestor(s)
matched_ancestors['gsep_t'] = matched_ancestors['gsep_x'] + matched_ancestors['gsep_y']
mrca = matched_ancestors.groupby(['descendant_x', 'descendant_y'])[['gsep_t']].min().reset_index()
mrca.columns = ['descendant_x', 'descendant_y', 'gsep_min']
# only retain relationship at least this recent
minrels = mrca.merge(matched_ancestors).query('gsep_t == gsep_min')
# calculate kinship due to each shared ancestor
# relatedness coefficient is 2x kinship coeffienct 
minrels['kinship'] = np.power(0.5, minrels['gsep_t']+1)

# find total kinship based on all shared ancestors
# sum of each shared ancestor
totalrel = minrels.groupby(
    ['descendant_x', 'descendant_y'])[['kinship', 'gsep_x','gsep_y']].agg(
    {'kinship': [sum, 'size'], 
     'gsep_x': lambda x: int(np.mean(x)), 
     'gsep_y': lambda x: int(np.mean(x))})

totalrel = totalrel.reset_index()
totalrel.columns = ['descendant_x', 'descendant_y', 'kinship_sum', 'n_shared_anc', 'gsep_x', 'gsep_y']
del minrels
del mrca
del matched_ancestors


# x is the number of generations up
# y is the number of generations down
# z is the number of inheritance paths
# assume x is >= y

rel_of_xyz = {
    # direct ancestors
    (1, 0, 1): 'PO',
    (2, 0, 1): '1G',
    (3, 0, 1): '2G',
    (4, 0, 1): '3G',
    (5, 0, 1): '4G',
    # siblings
    (1, 1, 2): 'FS',
    (1, 1, 1): 'HS',
    # avuncular
    (2, 1, 2): 'Av',
    (2, 1, 1): 'HAv',
    (3, 1, 2): '1GAv',
    (3, 1, 1): '1GHAv',
    (4, 1, 2): '2GAv',
    (4, 1, 1): '2GHAv',
    (5, 1, 2): '3GAv',
    (5, 1, 1): '3GHAv',
    # 1st cousins
    (2, 2, 2): '1C', # also twice first half cousins?
    (2, 2, 1): 'H1C',
    (3, 2, 2): '1C1R',
    (3, 2, 1): 'H1C1R',
    (4, 2, 2): '1C2R',
    (4, 2, 1): 'H1C2R',
    (5, 2, 2): '1C3R',
    (5, 2, 1): 'H1C3R',    
    # 2nd cousins 
    (3, 3, 2): '2C',
    (3, 3, 1): 'H2C',
    (4, 3, 2): '2C1R',
    (4, 3, 1): 'H2C1R',
    (5, 3, 2): '2C2R',
    (5, 3, 1): 'H2C2R',   
    # 3rd cousins
    (4, 4, 2): '3C',
    (4, 4, 1): 'H3C',
    (5, 4, 2): '3C1R',
    (5, 4, 1): 'H3C1R',
    # 4th cousins
    (5, 5, 2): '4C',
    (5, 5, 1): 'H4C',
    
    # extended relatives
    (2, 2, 4): '2x1C',
    (2, 2, 3): '1.5x1C',
    (3, 2, 4): '2x1C1R',
    (3, 2, 3): '1.5x1C1R',
    (4, 2, 4): '2x1C2R',
    (4, 2, 3): '1.5x1C2R',
    
    (3, 3, 4): '2x2C',
    (3, 3, 3): '1.5x2C',
    (4, 3, 4): '2x2C1R',
    (4, 3, 3): '1.5x2C1R',
}

rel_cats = pd.DataFrame({
    'gens_up': [x[0] for x in rel_of_xyz.keys()], 
    'gen_down': [x[1] for x in rel_of_xyz.keys()],
    'nshared_ancestors': [x[2] for x in rel_of_xyz.keys()],
    'rel': [x for x in rel_of_xyz.values()],
    'kinship': [x[2] * 0.5**(x[0]+x[1]+1) for x in rel_of_xyz.keys()],
})

# write these out to a file for reference
rel_cats.to_csv(relative_table, sep='\t',index=None)

# assign relative categories
rels = [] 
for row in totalrel[['gsep_x','gsep_y', 'n_shared_anc']].itertuples(index = False):
    #ensure the order matches 
    rel_tup = tuple(row) if (row[0]>=row[1]) else tuple((row[1], row[0], row[2]))
    rels.append(rel_of_xyz.get(rel_tup, 'undef'))

totalrel['rel'] = rels
del rels 

# change dtypes to take less memory
totalrel[['descendant_x','descendant_y']] = totalrel[['descendant_x','descendant_y']].astype(np.int32)
totalrel[['n_shared_anc','gsep_x', 'gsep_y']] = totalrel[['n_shared_anc','gsep_x', 'gsep_y']].astype(np.int8)
totalrel[['kinship_sum']] = totalrel[['kinship_sum']].astype(np.float32)

# write data frame out to compressed hdf file
totalrel.to_hdf(
    assigned_relatives, 
    mode='w', 
    complevel=3,
    complib='blosc:zstd',
    key='df'
)
