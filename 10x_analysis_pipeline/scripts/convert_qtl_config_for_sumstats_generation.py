import qtl_config_utils
import sys

in_config_file = sys.argv[1]
out_config_file = sys.argv[2]

config = qtl_config_utils.load_config(in_config_file)

# calculating marginal results for every SNP - no permutations, and a larger window size
config['w'] = 500000
config['np'] = 0
config['output_summary_stats'] = True

qtl_config_utils.write_config(config, out_config_file)
