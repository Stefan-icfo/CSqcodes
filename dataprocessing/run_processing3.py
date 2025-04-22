#qdac.set_gates_to_json_config('singledot_thermomech_march25')
#time.sleep(500)
qdac.ch06.dc_constant_V(-1.296)
for m in range(20):
  time.sleep(10)
  zurich.set_mixdown(166.397e6)
  %run experiments/take_spectra_and_fit_noslope.py
  zurich.set_mixdown(exp.max_thermomech_freq)
  for n in range(25):
    %run experiments/ringupdown_v4.py
    %run experiments/takedemodtimetrace.py
  
  current_CVS=qdac.ch06.dc_constant_V()
  qdac.ch06.dc_constant_V(current_CVS+1e-3)

  