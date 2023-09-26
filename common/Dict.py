
class Dict(dict):
    '''
    For test, simulate command line parsing by add get values method through '.' to dict class
    '''
    __setattr__ = dict.__setitem__
    __getattr__ = dict.__getitem__
