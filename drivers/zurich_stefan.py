from os import name
import zhinst.qcodes as ziqc
from experiment_functions.multi_sweep import SweepMultiParam
from qcodes.instrument.parameter import DelegateParameter, ScaledParameter
from qcodes.instrument.channel import InstrumentChannel, ChannelList
import numpy as np


class MyZurich(ziqc.UHFLI):
    '''
    [use_cache_conductance]

    This parameter (bool) of source.conductance
    that controls whether to use the last returned value of the 
    "source.voltage" function. e.g.

    measured_parameter = zurich.source.voltage()
    measured_parameter2 = zurich.source.conductance
    zurich.source.conductance.use_cache_conductance = True

    Will not read the voltage from the Zurich twice, whereas the 
    following will:
    
    measured_parameter = zurich.source.voltage()
    measured_parameter2 = zurich.source.conductance
    zurich.source.conductance.use_cache_conductance = False

    Note here that the order is important. Using the cached value 
    relies on it being the right one.

    '''
    def __init__(self, manual_instr, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.manual_instr = manual_instr

        '''
        Convention:

        Measurement always at frequency of osc0 on demod0.
        Source frequency set on osc1 and output on amplitudes4.
        Gate output always on amplitudes6 (modulation SB C-M)
        
        '''

        # ---------- SOURCE ----------------

        source = InstrumentChannel(parent=self, name='source')
        source.add_parameter(
            name='freq',
            parameter_class=DelegateParameter,
            source=self.oscs.oscs1.freq,
        )
        source.add_parameter(
            name='amplitude',
            parameter_class=DelegateParameter,
            source=self.sigouts.sigouts0.amplitudes.amplitudes4.value,
        )
        # source.add_parameter(
        #     name='voltage',
        #     parameter_class=DelegateParameter,
        #     source=self.demods.demods0.sample,
        #     label='voltage',
        #     snapshot_exclude=True
        # )

        source.add_parameter(
            name='demod_sample',
            parameter_class=DelegateParameter,
            source=self.demods.demods0.sample,
            label='demod_sample',
            snapshot_exclude=True
        )
        
        source.add_parameter(
            name='DC_voltage',
            parameter_class=DelegateParameter,
            source=self.sigouts.sigouts0.offset,
        )               

        
        source.add_parameter('source_power_at_CNT',
                             label='Source power at CNT',
                             unit='dBm',
                             docstring='Source power at CNT',
                             get_cmd=self._get_source_power_at_CNT,
                             set_cmd=False,
                             )

        self.add_parameter('time_constant',
                            parameter_class=DelegateParameter,
                            source=self.demods.demods0.timeconstant,)

        self.add_submodule('source', source)
        self.source.amplitude_vrms = ScaledParameter(output=self.source.amplitude, division=np.sqrt(2),
                                            label='Vrms_out')

        # ---------- GATE ----------------
        gate = InstrumentChannel(parent=self, name='gate')
        gate.add_parameter('amplitude',
                            parameter_class=DelegateParameter,
                            source=self.sigouts.sigouts1.amplitudes.amplitudes6.value,
        )


        # ---------- Compensation ----------------
        compensation = InstrumentChannel(parent=self, name='compensation')
        compensation.add_parameter('amplitude',
                            parameter_class=DelegateParameter,
                            source=self.sigouts.sigouts1.amplitudes.amplitudes7.value,
        )

        self.add_submodule('compensation', compensation)

        
    def _get_R(self):
        
        sample = self.source.demod_sample()

        X = sample['x'][0]
        Y = sample['y'][0]

        R = np.sqrt(X**2+Y**2)
        return R

    def _get_demod_complex(self): 

        sample = self.source.demod_sample()  

        X = sample['x'][0]
        Y = sample['y'][0]
    
        demod_complex = X + 1j*Y

        return demod_complex

    def _get_conductance(self):
        if self.source.conductance.use_cache_conductance:
            voltage = self.source.demod_complex.cache()
        else:
            voltage = self.source.demod_complex()

        V_in = self.source.amplitude()/np.sqrt(2)

        Pw_source = (V_in**2)/50             # [Vrms] to [W]
        P_source = 10*np.log10(Pw_source*1e3)  # [W] to [dBm]
        Psd_ac = P_source + self.manual_instr.source_line_attn()        # Substract attenuation
        # [dBm] to [Vrms] (factor 2 is due to impedence mismatch)
        Vsd_ac = 2*(1/np.sqrt(2))*10**((Psd_ac-10)/20)

        gain_RT = 200
        gain_HEMT = 5.64
        Z_tot = 7521

        # #G calculation
        I = voltage/(gain_RT*gain_HEMT*Z_tot)
        G = 1/((Vsd_ac/I)-Z_tot)  # ACTUALLY WE THINK (Chris and Victor) that there is an error here!!! Error related to subtracting a real value instead of an absolute one. See _get_conductance_Ithaco fucn for details.  

        return G

    def _get_conductance_Ithaco(self):
        ''' ZI output and input should be set to HZ. NOT 50 ohm.
        Measured value ("voltage") is in general complex such that we can see mag and phase of the conductance.
        '''
        if self.source.conductance.use_cache_conductance:
            voltage = self.source.demod_complex.cache()
        else:
            voltage = self.source.demod_complex()

        V_in = self.source.amplitude_vrms()  # Output Vrms voltage of the ZI. 

        attenuation = 1/101  # currently we have no attentuation, but we might add a voltage divider or whatever. 
        Vsd_ac = V_in * attenuation

        R_extra = 3e3 + 10e3 # extra resistance in the source DC line between output of the voltage source and the Ithaco current amp. e.g. on PCB could have 2 x 100 kOhm SMD resistors. 
        Gain_Ithaco_VperA = 1.00e10  # uncalibrated.  But basically 1.05e6 in the 10^6 setting. 
        Gain_preamp = 1  # uncalibrated
        Gain = Gain_Ithaco_VperA * Gain_preamp

        fac1 = Vsd_ac / voltage * Gain

        R_CNT = (abs(fac1) - R_extra)*np.exp(1j*np.angle(fac1))
        G_CNT = 1 / R_CNT
        return G_CNT

    def _get_resistance_Ithaco(self):
        ''' ZI output and input should be set to HZ. NOT 50 ohm.
        Measured value ("voltage") is in general complex such that we can see mag and phase of the conductance.
        '''
        if self.source.conductance.use_cache_conductance:
            voltage = self.source.demod_complex()
        else:
            voltage = self.source.demod_complex()

        V_in = self.source.amplitude_vrms()  # Output Vrms voltage of the ZI. 

        attenuation = 1/101  # currently we a voltage divider
        Vsd_ac = V_in * attenuation

        # R_extra = 218e3 + 10e3 # extra resistance in the source DC line between output of the voltage source and the Ithaco current amp. e.g. on PCB could have 2 x 100 kOhm SMD resistors. 
        R_extra = 3e3 + 10e3 # extra resistance in the source DC line between output of the voltage source and the Ithaco current amp. e.g. on PCB could have 2 x 100 kOhm SMD resistors. 
        Gain_Ithaco_VperA = 1.00e10  # uncalibrated.  But basically 1.05e6 in the 10^6 setting. 
        Gain_preamp = 1  # uncalibrated
        Gain = Gain_Ithaco_VperA * Gain_preamp

        fac1 = Vsd_ac / voltage * Gain

        R_CNT = (abs(fac1) - R_extra)*np.exp(1j*np.angle(fac1))
        G_CNT = 1 / R_CNT
        return R_CNT

    def _get_source_power_at_CNT(self):
        V_out_of_ZI_pk = self.source.amplitude()/np.sqrt(2)

        Pw_source_out_of_ZI = (V_out_of_ZI_pk**2)/50             # [Vrms] to [W]
        P_source_out_of_ZI = 10*np.log10(Pw_source_out_of_ZI*1e3)  # [W] to [dBm]
        Psd_ac_dBm = P_source_out_of_ZI + self.manual_instr.source_line_attn()        # Substract attenuation
        return Psd_ac_dBm

    def _get_gate_power_at_CNT(self):
        V_out_of_ZI_pk = self.gate.amplitude()/np.sqrt(2)

        Pw_gate_out_of_ZI = (V_out_of_ZI_pk**2)/50             # [Vrms] to [W]
        P_gate_out_of_ZI = 10*np.log10(Pw_gate_out_of_ZI*1e3)  # [W] to [dBm]
        Psd_ac_dBm = P_gate_out_of_ZI + self.manual_instr.gate_line_attn()        # Substract attenuation
        return Psd_ac_dBm


    def sweep_multi_output(self, output_list):
        '''
        Channels list will be the list of dictionaries with the form:
        {
            'channel': zurich.oscs.oscs0.freq,
            'sweep_values':{
                'start': x,
                'stop': y,
                'step': z, #or it could be 'num': 2
            }

        }
        '''
        return SweepMultiParam(output_list)
