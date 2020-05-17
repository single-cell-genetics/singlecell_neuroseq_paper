import pandas as pd
import random
import yaml
import os
random.seed(10)

dataset = 'pool1_13_noddd_D30'

# set up output file directory
#outfile_template = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_neuroseq/data/dataset_configuration_files/'+dataset+'_permutations/{}.yaml'
outfile_template = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_neuroseq/data/data_processed/'+dataset+'_permutations/{perm_name}/{perm_name}.yaml'

# load the configuration that will be subset in permutations
base_config_file = '/nfs/leia/research/stegle/dseaton/hipsci/singlecell_neuroseq/data/dataset_configuration_files/{}.yaml'.format(dataset)


with open(base_config_file, 'r') as f:
    base_config_dict = yaml.load(f)

pool_list = base_config_dict['sample_selection']['pool_id']

n_permutations = 40
n_subsample = 20

permutation_list = ['perm_{}'.format(x) for x in range(1,n_permutations+1)]
permutation_dict = dict([(x,[]) for x in permutation_list])


# to ensure each pool is in an equal number of permutations, loop over pools
for pool_id in pool_list:
    subsample = random.sample(permutation_list, n_subsample)
    for perm in subsample:
        permutation_dict[perm].append(pool_id)


for perm in permutation_dict.keys():
    filename = outfile_template.format(perm_name = perm)
    outdir = os.path.dirname(filename)
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
    config_dict = base_config_dict
    config_dict['sample_selection']['pool_id'] = permutation_dict[perm]
    with open(filename, 'w') as f:
        yaml.dump(config_dict, f)
