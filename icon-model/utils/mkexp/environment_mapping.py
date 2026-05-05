'''
Mapping from build info to environment name

$Id$
'''

from fnmatch import fnmatch

from configobj import ConfigObj
from feedback import debug

class Map:

    def __init__(self, map_config_name):
        self.map_config = ConfigObj(map_config_name, interpolation=False)

    def get_environment(self, pre_config):
        '''
        Return the name of the first environment that matches the given config

        Maps contain a list of sections, one for each environment.
        Their keys contain ``fnmatch`` ("glob") patterns that are matched
        against the value of their respective keys in the given system config.
        If a key contains a list of patterns, only one needs to match.
        The first section where all keys match is selected and returned,
        meaning that reordering of map sections may change the result.
        If no match is found, ``None`` is returned.
        '''
        for env in self.map_config.sections:
            debug(f"{env=}")
            match = True
            for key in self.map_config[env]:
                debug(f"{key=}")
                if key not in pre_config:
                    debug("not in pre_config")
                    match = False
                    break
                else:
                    globs = self.map_config[env][key]
                    # Always treat values as sequence
                    if not isinstance(globs, (list, tuple)):
                        globs = [globs]
                    debug(f"{globs=}")
                    debug(f"{pre_config[key]=}")
                    # Find match in sequence
                    if not any(
                        [fnmatch(pre_config[key], glob) for glob in globs]
                    ):
                        debug("no match found")
                        match = False
                        break
            # Quit after first successful match
            if match:
                return env
        return None
