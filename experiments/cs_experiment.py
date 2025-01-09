# cs_experiment.py

import numpy as np
import time
from tqdm import tqdm

# Import your parameters
import experiment_parameters as params

# Example placeholders for instruments
from instruments import station, zurich, qdac, Triton
from qcodes.dataset import Measurement, new_experiment

# Utility functions
from utils.sample_name import sample_name
from utils.zi_uhfli_GVg_setup import zi_uhfli_GVg_setup
from utils.d2v import d2v
from utils.v2d import v2d
from utils.rms2pk import rms2pk
from utils.CS_utils import save_metadata_var, get_var_name


class CSExperiment:
    """
    A class that copies defaults from experiment_parameters.py
    and provides measurement methods (GVG_fun, other_fun, etc.).
    """

    def __init__(self):
        # Copy from experiment_parameters
        self.device_name = params.device_name
        self.tc = params.tc
        self.attn_dB_source = params.attn_dB_source
        self.source_amplitude_instrumentlevel_GVg = params.source_amplitude_instrumentlevel_GVg
        self.mix_down_f = params.mix_down_f
        self.x_avg = params.x_avg
        self.y_avg = params.y_avg
        self.start_vg_cs = params.start_vg_cs
        self.stop_vg_cs = params.stop_vg_cs
        self.step_num_cs = params.step_num_cs

    def GVG_fun(
        self,
        #run=False,
        save_in_database=True,
        return_data=False,
        return_only_Vg_and_G=True,
        reverse=False,
        pre_ramping_required=False,
    ):
        """
        Example measurement that uses self.xxx (from experiment_parameters).
        Ad-hoc overrides can be done by directly changing self.xxx in the run file.
        """
        #if not run:
        #    print("GVG_fun: run=False, skipping measurement.")
        #    return

        tc = self.tc
        vsd_dB = self.attn_dB_source
        amp_lvl = self.source_amplitude_instrumentlevel_GVg
        f_mix = self.mix_down_f
        x_avg = self.x_avg
        y_avg = self.y_avg
        start_vg = self.start_vg_cs
        stop_vg = self.stop_vg_cs
        step_num = self.step_num_cs

        # Instrument references
        gate = qdac.ch06
        freq = zurich.oscs.oscs0.freq
        source_amplitude_param = zurich.sigouts.sigouts0.amplitudes.amplitudes0.value
        measured_parameter = zurich.demods.demods0.sample

        # Initialize hardware
        freq(f_mix)
        source_amplitude_param(amp_lvl)
        if pre_ramping_required:
            print(f"Pre-ramping gate to {start_vg}")
            qdac.ramp_multi_ch_slowly(channels=[gate], final_vgs=[start_vg])
        gate.ramp_ch(start_vg)

        vsdac = d2v(v2d(np.sqrt(1/2) * amp_lvl) - vsd_dB) / 10
        vgdc_sweep = gate.dc_constant_V.sweep(start=start_vg, stop=stop_vg, num=step_num)

        device_name = self.device_name
        prefix_name = 'Conductance_rf_'
        exp_name = f"vsac@inst={amp_lvl*1e3} mV"
        exp_dict = dict(mV=vsdac * 1000)
        postfix = (
            f"_g1={round(qdac.ch01.dc_constant_V(), 2)},"
            f"g2={round(qdac.ch02.dc_constant_V(), 2)},"
            f"g3={round(qdac.ch03.dc_constant_V(), 2)},"
            f"g4={round(qdac.ch04.dc_constant_V(), 2)},"
            f"g5={round(qdac.ch05.dc_constant_V(), 2)}"
        )
        gate.label = 'cs_gate'
        exp_name = sample_name(prefix_name, exp_dict, postfix)

        # Prepare data arrays
        Glist, Vlist, Ilist, Phaselist, Rlist = [], [], [], [], []

        if save_in_database:
            experiment = new_experiment(name=exp_name, sample_name=device_name)
            meas = Measurement(exp=experiment)
            meas.register_parameter(vgdc_sweep.parameter)
            meas.register_custom_parameter('G', unit='S', setpoints=[vgdc_sweep.parameter])
            meas.register_custom_parameter('V_r', unit='V', setpoints=[vgdc_sweep.parameter])
            meas.register_custom_parameter('Phase', unit='rad', setpoints=[vgdc_sweep.parameter])
            meas.register_custom_parameter('R', unit='Ohm', setpoints=[vgdc_sweep.parameter])
            meas.register_custom_parameter('I', unit='A', setpoints=[vgdc_sweep.parameter])

            with meas.run() as datasaver:
                varnames = [str(name) for name in [tc, vsd_dB, amp_lvl, x_avg, y_avg]]
                save_metadata_var(datasaver.dataset, varnames, [tc, vsd_dB, amp_lvl, x_avg, y_avg])

                for vgdc_value in tqdm(vgdc_sweep, desc='Gate voltage Sweep'):
                    gate.ramp_ch(vgdc_value)
                    time.sleep(1.1 * tc)

                    # Some measurement
                    #_ = measured_parameter()
                    theta_calc, v_r_calc, I, G = zurich.phase_voltage_current_conductance_compensate(vsdac)
                    R = 1 / G

                    datasaver.add_result(
                        ('R', R),
                        ('G', G),
                        ('V_r', v_r_calc),
                        ('Phase', theta_calc),
                        ('I', I),
                        (vgdc_sweep.parameter, vgdc_value)
                    )

                    if return_data:
                        Glist.append(G)
                        Vlist.append(v_r_calc)
                        Ilist.append(I)
                        Phaselist.append(theta_calc)
                        Rlist.append(R)
        else:
            # Not saving in DB
            for vgdc_value in tqdm(vgdc_sweep, desc='Gate voltage Sweep'):
                gate.ramp_ch(vgdc_value)
                time.sleep(1.1 * tc)

                _ = measured_parameter()
                theta_calc, v_r_calc, I, G = zurich.phase_voltage_current_conductance_compensate(vsdac)
                R = 1 / G

                if return_data:
                    Glist.append(G)
                    Vlist.append(v_r_calc)
                    Ilist.append(I)
                    Phaselist.append(theta_calc)
                    Rlist.append(R)

        if return_data and reverse:
            Glist.reverse()
            Vlist.reverse()
            Ilist.reverse()
            Phaselist.reverse()
            Rlist.reverse()

        if return_data:
            Vglist = list(vgdc_sweep)
            if return_only_Vg_and_G:
                return np.array(Vglist), np.array(Glist)
            else:
                return (
                    np.array(Vglist),
                    np.array(Glist),
                    np.array(Vlist),
                    np.array(Ilist),
                    np.array(Phaselist),
                    np.array(Rlist),
                )

    