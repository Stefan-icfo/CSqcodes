from cs_experiment import *

class CS_meta(CSExperiment()):
    def __init__(self, name: str, address: str, **kwargs):
        super().__init__(name, address, **kwargs)