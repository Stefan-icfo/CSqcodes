# run_GVG_experiment.py

from cs_experiment import CSExperiment

# 1) Create an experiment instance, reading defaults from experiment_parameters.py
expt = CSExperiment()

# 2) Do ad-hoc overrides up top:
#  change  sweep range:

expt.start_vg_cs = -1.6
expt.stop_vg_cs = -1.4
expt.step_num_cs = 200*30

# 3) Define a function to run GVG
def run_GVG():
    """
    We can call the GVG_fun method from the expt instance.
    Adjust additional arguments if needed.
    """
    # Example: run the measurement, get data
    result = expt.GVG_fun(run=True, return_data=True, return_only_Vg_and_G=True)
    
    if result is not None:
        Vg, G = result
        #print("Ran GVG_fun. Vg:", Vg)
        #print("Conductances G:", G)

if __name__ == "__main__":
    run_GVG()
