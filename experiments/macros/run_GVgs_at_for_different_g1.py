from instruments import exp,qdac
for g1V in [0.7,0.9,1.1,1.3,1.5,1.7,1.9,0,-0.2,-0.4,-0.6,-0.8,-1,-1.2]:
    qdac.ramp_multi_ch_fast([qdac.ch01,qdac.ch06],[g1V,-2.5])
    exp.GVG_fun(pre_ramping_required=True)
