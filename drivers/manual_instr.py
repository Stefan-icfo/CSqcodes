from qcodes.instrument.base import Instrument
from qcodes.instrument.parameter import ManualParameter

class Manual_Inst(Instrument):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.add_parameter(
            name='source_v_divider',
            parameter_class=ManualParameter,
            label='divider of the source voltage',
        )

        self.add_parameter(
            name='source_line_attn',
            parameter_class=ManualParameter,
            label='attenuation of the source line',
            initial_value=-57,
            unit='dB'
        )

        self.add_parameter(
            name='source_power',
            parameter_class=ManualParameter,
            label='Power set on the zurich to the source line',
            unit='Vrms'
        )
        
        self.add_parameter(
            name='gate_line_attn',
            parameter_class=ManualParameter,
            label='attenuation of the gate line',
            initial_value=-36,
            unit='dB'
        )
        
        self.add_parameter(
            name='v_final',
            parameter_class=ManualParameter,
            label='Final voltage set',
            unit='V'
            )
        
        self.add_parameter(
            name='v_initial',
            parameter_class=ManualParameter,
            label='Initial voltage set',
            unit='V'
            )

        self.add_parameter(
            name='gain_RT',
            parameter_class=ManualParameter,
            label='Gain of the Room temperature amplifier',
            initial_value=200,
            unit=None,
            )

        self.add_parameter(
            name='gain_HEMT',
            parameter_class=ManualParameter,
            label='Gain of the HEMT',
            initial_value=5.64,
            unit=None,
            )
                    
        self.add_parameter(
            name='Z_total',
            parameter_class=ManualParameter,
            label='Total impedance of the circuit',
            initial_value=7521,
            unit='ohm',
            )

        self.add_parameter(
            name='output_attenuator',
            parameter_class=ManualParameter,
            label='Attenuation at the otuput of the VNA',
            initial_value=0,
            unit='dB'
        )
