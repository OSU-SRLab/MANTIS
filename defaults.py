# @file defaults.py
# @author Esko Kautto (esko.kautto@osumc.edu)
# @updated 2016-07-06
#
# Contains helpers for using the three-tier priority loading of configuration
# settings that the main MANTIS script utilizes.

import os

"""
Helper class that stores data needed for the parameter loading.
"""
class Parameter(object):
    def __init__(self, key=None, arg_key=None, required=False, default=None):
        self.__key = key
        self.__arg_key = arg_key
        if arg_key is None:
            self.__arg_key = key
        self.__required = required
        self.__default = default
        # end .__init__()

    def key(self):
        return self.__key

    def arg_key(self):
        return self.__arg_key

    def required(self):
        return self.__required

    def default(self):
        return self.__default

    # end Parameter class definition.


"""
Loads any configured values from the config file and
returns them as a dict object.
"""
def load_config_file(filepath):
    config = {}
    with open(filepath, 'Ur') as f:
        for line in f:
            line = line.strip()
            if len(line) and line[0] != '#':
                if '=' not in line:
                    print('Invalid line in config file!')
                    print(line)
                    exit(1)
                line = line.split('=', 2)
                key = line[0].strip()
                value = line[1].strip()
                config[key] = value
        f.close()
    return config
    # end load_config_file()

"""
Performs the three-tier priority loading of config settings,
where command line arguments get the highest priority, followed
by config file settings, and finally followed by the default
values for the parameters.
"""
def load_settings(params, config_file_path, args):
    settings = {}

    config_values = {}
    if os.path.isfile(config_file_path):
        config_values = load_config_file(config_file_path)
    
    for param in params:
        key = param.key()
        required = param.required()
        arg_key = param.arg_key()

        if hasattr(args, arg_key) and getattr(args, arg_key) is not None:
            # Command line args take precedence.
            settings[key] = getattr(args, arg_key)
        elif key in config_values:
            # Configuration file takes second priorty.
            settings[key] = config_values[key]
        elif param.default() is not None:
            # Default values take third priority.
            settings[key] = param.default()

        if required and (key not in settings or settings[key] is None):
            # required value was not set
            print('Error: Required parameter \'{0}\' not set!'.format(key))
            exit(1)

    return settings
    # end load_settings()