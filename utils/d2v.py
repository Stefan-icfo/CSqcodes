import math
def d2v(dBm):
    vpk = 10**((dBm-10)/20)
    v=vpk/math.sqrt(2)
    return v