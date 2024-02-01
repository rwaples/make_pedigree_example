import pandas as pd
import numpy as np
import os


# input
ancestors_path = str(snakemake.input.ancestors)
ancestors_batch = [str(x) for x in snakemake.input.splits]
# output
OUTFILE = str(snakemake.output.relatives)

# we want to do an inner merge of anc_df with itself, on the ancestor column
# merging can be done in parts
# left data is kept whole, right data is split as needed
# read the entire data set
anc_df = pd.read_csv(ancestors_path, sep='\t', header=None)
anc_df.columns = ['descendant', 'ancestor', 'gsep']


for f in ancestors_batch:
    try:
    	# read in the parts
        merge_in=pd.read_csv(f, sep ='\t', header=None)
        merge_in.columns = ['descendant', 'ancestor', 'gsep']
        m = pd.merge(
            merge_in, anc_df, 
            on = 'ancestor', 
            how = 'inner',
            suffixes = ['_x', '_y'],
        )
        m = m.query('descendant_x != descendant_y')  # get rid of descendant matches
        if os.path.isfile(OUTFILE):
            header = None
        else:
            header = True
        # columns give the two related individuals, number of generations up and down to their shared ancestor, and the ancestor
        # append to the output
        m = m[['descendant_x', 'descendant_y', 'gsep_x', 'gsep_y', 'ancestor']]
        m.to_csv(OUTFILE, mode='a', index=None, sep='\t', header=header)
    except pd.errors.EmptyDataError:
        pass

print('done with merging')


# above code does not record direct descendants
# record them out to 5 generations
direct_ancestors = anc_df.copy()
rel_of_gsep = {1:'P0', 2: '1G', 3:'2G', 4:'3G', 5:'4G'}
direct_ancestors['rel'] = direct_ancestors['gsep'].map(rel_of_gsep)
rel_types = pd.CategoricalDtype(categories=['P0', '1G', '2G', '3G', '4G'], ordered=True)
direct_ancestors['rel'].astype(rel_types)
direct_ancestors = direct_ancestors.sort_values(['descendant', 'ancestor'])

direct_ancestors.columns = ['descendant_x', 'descendant_y', 'gsep_x', 'rel']
direct_ancestors['gsep_y'] = 0
direct_ancestors['ancestor'] = direct_ancestors['descendant_y']
direct_ancestors = direct_ancestors[['descendant_x', 'descendant_y', 'gsep_x', 'gsep_y', 'ancestor']]

direct_ancestors.to_csv(OUTFILE, mode='a', index=None, sep='\t', header=None)
