from os import name
import zhinst.qcodes as ziqc
from qcodes.instrument_drivers.zurich_instruments.uhfli import UHFLI
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

        Measurement always at frequency of osc0 on demod3.
        Source frequency set on osc1 and output on amplitudes4.
        Gate output always on amplitudes6 (modulation SB C-M)
        
        '''

        # ---------- SOURCE ----------------

        source = InstrumentChannel(parent=self, name='source')
        source.add_parameter(
            name='freq',
            parameter_class=DelegateParameter,
            source=self.oscs.oscs0.freq,
        )
        source.add_parameter(
            name='amplitude',
            parameter_class=DelegateParameter,
            source=self.sigouts.sigouts0.amplitudes.amplitudes4.value,
        )
        source.add_parameter(
            name='voltage',
            parameter_class=DelegateParameter,
            source=self.demods.demods2.sample,
            label='voltage',
            snapshot_exclude=True
        )
        # source.add_parameter('demod3_complex',
        #                    label='Demod Complex Value',
        #                    unit='V',
        #                    docstring='X and Y value of the demodulator 1 in complex form',
        #                    get_cmd=self._get_complex,
        #                    )
        source.add_parameter('conductance',
                             label='conductance',
                             unit='S',
                             docstring='Conductance of the channel',
                             get_cmd=self._get_conductance,
                             get_parser=complex,
                             set_cmd=False,
                             )
        source.add_parameter('source_power_at_CNT',
                             label='Source power at CNT',
                             unit='dBm',
                             docstring='Source power at CNT',
                             get_cmd=self._get_source_power_at_CNT,
                             set_cmd=False,
                             )

        self.add_submodule('source', source)
        self.source.amplitude_vrms = ScaledParameter(output=self.source.amplitude, division=np.sqrt(2),
                                            label='Vrms_out')
        self.source.conductance.use_cache_conductance = False



        # ---------- GATE ----------------
        gate = InstrumentChannel(parent=self, name='gate')
        gate.add_parameter('amplitude',
                            parameter_class=DelegateParameter,
                            source=self.sigouts.sigouts1.amplitudes.amplitudes6.value,
        )
        #gate.add_parameter('gate_power_at_CNT',
        #                    label='Gate power at CNT',
        #                    unit='dBm',
        #                    docstring='Gate power at CNT',
        #                    get_cmd=self._get_gate_power_at_CNT,
        #                    set_cmd=False,
        #                    )
        self.add_submodule('gate', gate)


        # ---------- Compensation ----------------
        compensation = InstrumentChannel(parent=self, name='compensation')
        compensation.add_parameter('amplitude',
                            parameter_class=DelegateParameter,
                            source=self.sigouts.sigouts1.amplitudes.amplitudes7.value,
        )

        self.add_submodule('compensation', compensation)
    
    # def _get_complex(self):
    #     tmp = self.self.demods.demods0.sample()
    #     return tmp['x'][0]+1j*tmp['y'][0]
        
    def _get_conductance(self):
        if self.source.conductance.use_cache_conductance:
            voltage = self.source.voltage.cache()
        else:
            voltage = self.source.voltage()

        V_in = self.source.amplitude()/np.sqrt(2)

        Pw_source = (V_in**2)/50             # [Vrms] to [W]
        P_source = 10*np.log10(Pw_source*1e3)  # [W] to [dBm]
        Psd_ac = P_source + self.manual_instr.vsd_attn()        # Substract attenuation
        # [dBm] to [Vrms] (factor 2 is due to impedence mismatch)
        Vsd_ac = 2*(1/np.sqrt(2))*10**((Psd_ac-10)/20)

        gain_RT = 200
        gain_HEMT = 5.64
        Z_tot = 7521

        # #G calculation
        I = voltage/(gain_RT*gain_HEMT*Z_tot)
        G = 1/((Vsd_ac/I)-Z_tot)

        return G

    
    def _get_source_power_at_CNT(self):
        V_out_of_ZI_pk = self.source.amplitude()/np.sqrt(2)

        Pw_source_out_of_ZI = (V_out_of_ZI_pk**2)/50             # [Vrms] to [W]
        P_source_out_of_ZI = 10*np.log10(Pw_source_out_of_ZI*1e3)  # [W] to [dBm]
        Psd_ac_dBm = P_source_out_of_ZI + self.manual_instr.vsd_attn()        # Substract attenuation
        return Psd_ac_dBm

    def _get_gate_power_at_CNT(self):
        V_out_of_ZI_pk = self.gate.amplitude()/np.sqrt(2)

        Pw_gate_out_of_ZI = (V_out_of_ZI_pk**2)/50             # [Vrms] to [W]
        P_gate_out_of_ZI = 10*np.log10(Pw_gate_out_of_ZI*1e3)  # [W] to [dBm]
        Psd_ac_dBm = P_gate_out_of_ZI + self.manual_instr.vg_attn()        # Substract attenuation
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
