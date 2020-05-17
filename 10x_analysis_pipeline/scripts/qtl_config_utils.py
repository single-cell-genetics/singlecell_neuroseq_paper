import yaml

def write_config(config_dict, filename):
    with open(filename, 'w') as f:
        yaml.dump(config_dict, f)

def load_config(filename):
    with open(filename, 'r') as f:
        config_dict = yaml.load(f)
    return config_dict

