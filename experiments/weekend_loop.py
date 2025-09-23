

V_sit=exp.sit_at_max_Isens()
max_sens_vgo,max_sens_Vgcs,max_sens=exp.linesweep_parallel_LFsens(costum_prefix="g23halfwaydotg1adjust",start_vgi=V_sit-45e-3,stop_vgi=V_sit+15e-3,start_vgi_scan=V_sit)

qdac.ramp_multi_ch_slowly([1,2,3,4,5,6],[-4,0.3,0.3,0.3,0.3,0.8],step_size=5e-2,ramp_speed=5e-3) #half g2 g3 dot


vg1,vgcs,sensmax=exp.linesweep_parallel_LFsens(costum_prefix='g1') #g1 dot
qdac.ch01.dc_constant_V(vg1)
time.sleep(10)
f_m=exp.find_mech_mode(start_drive=20e-3,freq_range=[115e6,130e6],start_step_pitch=2e3,found_range=5e6,find_sitpos=True)
zurich.set_mixdown(f_m)
for n in range(2):
    %run experiments/zurich_experiments/spectrum_0825.py
%run experiments/run_swept_driven_vs_sitpos_ampsweep_g1_190925.py
%run experiments/run_thermomech_vs_sitpos_g3_18_09_25.py
%run experiments/run_thermomech_vs_sitpos_g3_18_09_25.py


qdac.ramp_multi_ch_slowly([1,2,3,4,5,6],[-4,0.3,0.3,0.3,0.3,0.95],step_size=5e-2,ramp_speed=5e-3)
vg1,vgcs,sensmax=exp.linesweep_parallel_LFsens(costum_prefix='g1')

#if time for intervention: g1g2 dot 
#qdac.ramp_multi_ch_slowly([1,2,3,4,5,6],[-2,-2,0.3,0.3,0.3,0.8],step_size=5e-2,ramp_speed=5e-3)

###################

##################

qdac.ch07.dc_constant_V(140e-6)
qdac.set_gates_to_metadata_config(1146) #halfway between g3 and g4
exp.sit_at_max_Isens(side="left")
Vg1,Vgcs,sens=exp.linesweep_parallel_LFsens(costum_prefix='g2_dot_adjust_g1',start_vgo=qdac.ch01.dc_constant_V()-0.05,stop_vgo=qdac.ch01.dc_constant_V()-0.05,start_vgi_scan=qdac.ch06.dc_constant_V())
qdac.ch01.dc_constant_V(Vg1)
time.sleep(10)
exp.sit_at_max_Isens(side="left")
f_max=exp.find_mech_mode(start_drive=20e-3,freq_range=[144e6,146e6],start_step_pitch=1e3,found_range=1e6,find_sitpos=False)
zurich.set_mixdown(f_max)
time.sleep(50)
for n in range(50):
    %run experiments/zurich_experiments/spectrum_0825.py

qdac.set_gates_to_metadata_config(1172) #3/4 g3 and 1/4 g4
exp.sit_at_max_Isens(side="left")
Vg1,Vgcs,sens=exp.linesweep_parallel_LFsens(costum_prefix='g2_dot_adjust_g1',start_vgo=qdac.ch01.dc_constant_V()-0.05,stop_vgo=qdac.ch01.dc_constant_V()-0.05,start_vgi_scan=qdac.ch06.dc_constant_V())
qdac.ch01.dc_constant_V(Vg1)
time.sleep(10)
exp.sit_at_max_Isens(side="left")
f_max=exp.find_mech_mode(start_drive=20e-3,freq_range=[144e6,146e6],start_step_pitch=1e3,found_range=1e6,find_sitpos=False)
zurich.set_mixdown(f_max)
time.sleep(50)
for n in range(50):
    %run experiments/zurich_experiments/spectrum_0825.py

qdac.set_gates_to_metadata_config(1025)  #1/4 g3 and 3/4 g4
exp.sit_at_max_Isens(side="left")
Vg1,Vgcs,sens=exp.linesweep_parallel_LFsens(costum_prefix='g2_dot_adjust_g1',start_vgo=qdac.ch01.dc_constant_V()-0.05,stop_vgo=qdac.ch01.dc_constant_V()-0.05,start_vgi_scan=qdac.ch06.dc_constant_V())
qdac.ch01.dc_constant_V(Vg1)
time.sleep(10)
exp.sit_at_max_Isens(side="left")
f_max=exp.find_mech_mode(start_drive=20e-3,freq_range=[147e6,149e6],start_step_pitch=1e3,found_range=1e6,find_sitpos=False)
zurich.set_mixdown(f_max)
#zurich.set_mixdown(148.105e6)
time.sleep(50)
for n in range(50):
    %run experiments/zurich_experiments/spectrum_0825.py

qdac.ramp_multi_ch_slowly([1,2,3,4,5,6],[0.5,0.3,-4.7,0.3,0.3,0.9],step_size=5e-2,ramp_speed=5e-3) #back to g3 dot
qdac.set_gates_to_json_config("g3dot_18_09_25")
Vg1,Vgcs,sens=exp.linesweep_parallel_LFsens(costum_prefix='g2_dot_adjust_g1',start_vgo=qdac.ch01.dc_constant_V()-0.05,stop_vgo=qdac.ch01.dc_constant_V()-0.05,start_vgi_scan=qdac.ch06.dc_constant_V())
qdac.ch01.dc_constant_V(Vg1)
time.sleep(10)
exp.sit_at_max_Isens(side="left")
f_max=exp.find_mech_mode(start_drive=20e-3,freq_range=[149e6,151e6],start_step_pitch=1e3,found_range=0.5e6,find_sitpos=False)
zurich.set_mixdown(f_max)
time.sleep(50)
for n in range(50):
    %run experiments/zurich_experiments/spectrum_0825.py

qdac.ramp_multi_ch_slowly([1,2,3,4,5,6],[0.5,0.3,-3.5,0.3,0.3,0.9],step_size=5e-2,ramp_speed=5e-3) #reduce hole nr
Vg1,Vgcs,sens=exp.linesweep_parallel_LFsens(costum_prefix='g2_dot_adjust_g1',start_vgo=qdac.ch01.dc_constant_V()-0.09,stop_vgo=qdac.ch01.dc_constant_V()+0.2,start_vgi_scan=qdac.ch06.dc_constant_V())
qdac.ramp_multi_ch_slowly([1],[Vg1],step_size=5e-2,ramp_speed=5e-3)
time.sleep(10)
exp.sit_at_max_Isens(side="left")
f_max=exp.find_mech_mode(start_drive=20e-3,freq_range=[135e6,160e6],start_step_pitch=1e3,found_range=0.5e6,find_sitpos=False)
zurich.set_mixdown(f_max)
time.sleep(50)
for n in range(50):
    %run experiments/zurich_experiments/spectrum_0825.py

qdac.ramp_multi_ch_slowly([1,2,3,4,5,6],[0.5,0.3,-2.3,0.3,0.3,0.9],step_size=5e-2,ramp_speed=5e-3) #reduce hole nr
Vg1,Vgcs,sens=exp.linesweep_parallel_LFsens(costum_prefix='g2_dot_adjust_g1',start_vgo=qdac.ch01.dc_constant_V()-0.09,stop_vgo=qdac.ch01.dc_constant_V()+0.2,start_vgi_scan=qdac.ch06.dc_constant_V())
qdac.ramp_multi_ch_slowly([1],[Vg1],step_size=5e-2,ramp_speed=5e-3)
time.sleep(10)
exp.sit_at_max_Isens(side="left")
f_max=exp.find_mech_mode(start_drive=20e-3,freq_range=[135e6,160e6],start_step_pitch=1e3,found_range=0.5e6,find_sitpos=False)
zurich.set_mixdown(f_max)
time.sleep(50)
for n in range(50):
    %run experiments/zurich_experiments/spectrum_0825.py

qdac.ramp_multi_ch_slowly([1,2,3,4,5,6],[0.5,0.3,-1.1,0.3,0.3,0.9],step_size=5e-2,ramp_speed=5e-3)
Vg1,Vgcs,sens=exp.linesweep_parallel_LFsens(costum_prefix='g2_dot_adjust_g1',start_vgo=qdac.ch01.dc_constant_V()-0.09,stop_vgo=qdac.ch01.dc_constant_V()+0.2,start_vgi_scan=qdac.ch06.dc_constant_V())
qdac.ramp_multi_ch_slowly([1],[Vg1],step_size=5e-2,ramp_speed=5e-3)
time.sleep(10)
exp.sit_at_max_Isens(side="left")
f_max=exp.find_mech_mode(start_drive=20e-3,freq_range=[135e6,160e6],start_step_pitch=1e3,found_range=0.5e6,find_sitpos=False)
zurich.set_mixdown(f_max)
time.sleep(50)
for n in range(50):
    %run experiments/zurich_experiments/spectrum_0825.py



qdac.ramp_multi_ch_slowly([1,2,3,4,5,6],[-4,0.3,0.3,0.3,0.3,0.95],step_size=5e-2,ramp_speed=5e-3)
vg1,vgcs,sensmax=exp.linesweep_parallel_LFsens(costum_prefix='g1')