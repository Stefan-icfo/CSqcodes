
run_id=812
from dataprocessing.processing_fkt import *
from dataprocessing.extract_fkts import *
gate2,gate4,data_2d=extract_2d(run_id,
                  data_2d_name="signal_shift_Vxn_deriv",
                  setpoints1_name='QDAC_ch02_dc_constant_V',  # gate 2
                   setpoints2_name='QDAC_ch04_dc_constant_V',  # gate 4
                  plot=False)

ICT_points(gate2, gate4, data_2d, run_id=run_id, threshold_std=4.2, plot=True)