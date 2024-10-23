from instruments import Triton
import numpy as np

def relevant_T(cutoff=2.5):
    if Triton.T5()>cutoff:
        return Triton.T5()
    else: 
        return Triton.MC()