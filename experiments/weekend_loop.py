import time

from experiments.zurich_experiments.spectrum_0925B import run_thermomech_temp_meas

exp.GVG_fun(step_num=8*200)
exp.sit_at_max_Isens(side="left")
zurich.set_mixdown(2.5e6)
time.sleep(100)
run_thermomech_temp_meas(reps_nodrive=50,exp_name="for_sensitivity",take_time_resolved_spectrum=False,avg_num=1,background_id=None)

zurich.set_mixdown(2.6e6)
time.sleep(100)
run_thermomech_temp_meas(reps_nodrive=50,exp_name="for_sensitivity",take_time_resolved_spectrum=False,avg_num=1,background_id=None)

zurich.set_mixdown(3.5e6)
time.sleep(100)
run_thermomech_temp_meas(reps_nodrive=50,exp_name="for_sensitivity",take_time_resolved_spectrum=False,avg_num=1,background_id=None)

zurich.set_mixdown(12.5e6)
time.sleep(100)
run_thermomech_temp_meas(reps_nodrive=50,exp_name="for_sensitivity",take_time_resolved_spectrum=False,avg_num=1,background_id=None)

pos_list1 = [0.412826, 0.563126, 0.693387, 0.828657, 0.973948, 1.11924, 1.26453, 1.40481, 1.5501, 1.69539, 1.84068, 1.98597, 2.12124, 2.26152, 2.4018, 2.53707, 2.67234, 2.81263, 2.95291, 3.08818, 3.22345, 3.35872, 3.49399, 3.61924, 3.74449, 3.87976]#midpoints of g2 linesweep on 1710
pos_list2 =[0.452906, 0.598196, 0.723447, 0.866232, 1.00902, 1.15681, 1.2996, 1.43988, 1.58768, 1.73046, 1.87826, 2.02104, 2.15381, 2.2991, 2.43437, 2.57214, 2.70491, 2.8502, 2.98547, 3.12325, 3.25601, 3.39379, 3.52655, 3.6493, 3.77705, 3.91483]#closer to top asymmetry=3
pos_list2 =[ 0.723447, 0.866232, 1.00902, 1.15681, 1.2996, 1.43988, 1.58768, 1.73046, 1.87826, 2.02104, 2.15381, 2.2991, 2.43437, 2.57214, 2.70491, 2.8502, 2.98547, 3.12325, 3.25601, 3.39379, 3.52655, 3.6493, 3.77705, 3.91483]#closer to top asymmetry=3
exp.go_through_gate_pos(pos_list=pos_list2)#high points
exp.go_through_gate_pos(pos_list=pos_list1)#midpoints


exp.linesweep_parallel_LFsens_extended(costum_prefix='try_move_15_e_to_right',find_startpos=True,
main_gate=qdac.ch02.dc_constant_V,
aux_gates=[qdac.ch03.dc_constant_V,qdac.ch01.dc_constant_V],
increments=[-1.02419,-0.277561])

%run experiments/DQD_charge_stability_viaCS_standardD.py


