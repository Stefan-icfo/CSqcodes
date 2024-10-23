def sample_name(prefix:str='run_',specs:dict={},comments:str='v1'):
    out_string = prefix
    for key in specs:
        out_string += "{:.{}f}_".format(specs[key],5) + key + '_'
#        out_string += "{:.{}f}_".format(specs[key],5) + "{}_".format(specs[key])+ key + '_'

    out_string += comments
    print(comments)
    return out_string
