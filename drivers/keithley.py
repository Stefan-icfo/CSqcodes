from qcodes.instrument_drivers.tektronix.Keithley_2400 import Keithley_2400
from qcodes.instrument.parameter import DelegateParameter, ScaledParameter

class MyKeithley(Keithley_2400):
    

    def __init__(self, name: str, address: str, **kwargs):
        super().__init__(name, address, **kwargs)

        self.add_parameter('R_CNT',
                        label='R_CNT',
                        unit='Ohm',
                        docstring='CNT inferred resistance',
                        get_cmd=self._get_R_CNT,
                        set_cmd=False,
                        )

        self.add_parameter('curr_bkg',
                        label='curr_bkg',
                        unit='A',
                        docstring='current bkg',
                        initial_value=0,
                        )
        self.R_CNT.use_cache_conductance = False

    def _get_R_CNT(self):
        if self.R_CNT.use_cache_conductance:
            voltage = self.volt.cache()
        else:
            voltage = self.volt()
        I_measure = self.curr() - self.curr_bkg()
        V_source = voltage
        R_CNT_Ohm = (V_source/(I_measure)) # resistance in Ohm
        return R_CNT_Ohm