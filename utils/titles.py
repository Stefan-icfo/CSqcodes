from qcodes.instrument.parameter import Parameter
from qcodes.instrument.sweep_values import SweepFixedValues
from experiment_functions.multi_sweep import SweepMultiParam
from si_prefix import si_format

def build_simple_title(custom_prefix, build_title_dict):
    '''
    SweepMultiParam : is a swept object for multigate (only BILT) sweeps
    SweepFixedValues : is a swept object.
    else : for regular parameters

    e.g. 
    bilt.sweep_multi_channel([...]) is a SweepMultiParam
    zurich.oscs.oscs1.freq.sweep(start=1, stop=10, step=1) is a SweepFixedValues
    triton.T4 is a single regular parameter.
    
    '''
    out_string = custom_prefix
    for key in build_title_dict:

        if isinstance(build_title_dict[key], SweepMultiParam):
            sweep_object = build_title_dict[key]
            out_string += key + ":" + \
                f"{sweep_object.parameter.label}{sweep_object.parameter.unit},"
        elif isinstance(build_title_dict[key], SweepFixedValues):
            sweep_object = build_title_dict[key]
            out_string += key + ":" + \
                f"{sweep_object.parameter.instrument.name_parts[-1]}({si_format(sweep_object[0])},{si_format(sweep_object[-1])}){sweep_object.parameter.unit},"
        else:
            out_string += key + ":" + \
                f"{si_format(build_title_dict[key].get())}" + \
                f"{build_title_dict[key].unit},"

    return out_string[:-1]
