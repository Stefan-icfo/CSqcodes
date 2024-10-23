from typing import Any
import numpy as np
from qcodes import VisaInstrument
from qcodes.utils.helpers import create_on_off_val_mapping
from qcodes.utils.validators import Enum, Strings


class Keithley_2400(VisaInstrument):
    """
    QCoDeS driver for the Keithley 2400 voltage source.
    """
    def __init__(
            self,
            name: str,
            address: str,
            **kwargs: Any):
        super().__init__(name, address, terminator='\n', **kwargs)

        self.add_parameter('rangev',
                           get_cmd='SENS:VOLT:RANG?',
                           get_parser=float,
                           set_cmd='SOUR:VOLT:RANG {:f}',
                           label='Voltage range')

        self.add_parameter('rangei',
                           get_cmd='SENS:CURR:RANG?',
                           get_parser=float,
                           set_cmd='SOUR:CURR:RANG {:f}',
                           label='Current range')

        self.add_parameter('compliancev',
                           get_cmd='SENS:VOLT:PROT?',
                           get_parser=float,
                           set_cmd='SENS:VOLT:PROT {:f}',
                           label='Voltage Compliance')

        self.add_parameter('compliancei',
                           get_cmd='SENS:CURR:PROT?',
                           get_parser=float,
                           set_cmd='SENS:CURR:PROT {:f}',
                           label='Current Compliance')

        self.add_parameter('volt',
                           get_cmd=self._get_read_output_protected,
                           get_parser=self._volt_parser,
                           set_cmd=':SOUR:VOLT:LEV {:.8f}',
                           label='Voltage',
                           unit='V',
                           docstring="Sets voltage in 'VOLT' mode. "
                                     "Get returns measured voltage if "
                                     "sensing 'VOLT' otherwise it returns "
                                     "setpoint value. "
                                     "Note that it is an error to read voltage with "
                                     "output off")

        self.add_parameter('curr',
                           get_cmd=self._get_read_output_protected,
                           get_parser=self._curr_parser,
                           set_cmd=':SOUR:CURR:LEV {:.8f}',
                           label='Current',
                           unit='A',
                           docstring="Sets current in 'CURR' mode. "
                                     "Get returns measured current if "
                                     "sensing 'CURR' otherwise it returns "
                                     "setpoint value. "
                                     "Note that it is an error to read current with "
                                     "output off")

        self.add_parameter('mode',
                           vals=Enum('VOLT', 'CURR'),
                           get_cmd=':SOUR:FUNC?',
                           set_cmd=self._set_mode_and_sense,
                           label='Mode')

        self.add_parameter('sense',
                           vals=Strings(),
                           get_cmd=':SENS:FUNC?',
                           set_cmd=':SENS:FUNC "{:s}"',
                           label='Sense mode')

        self.add_parameter(
            'output',
            set_cmd=':OUTP:STAT {}',
            get_cmd=':OUTP:STAT?',
            val_mapping=create_on_off_val_mapping(on_val="1", off_val="0")
        )

        self.add_parameter('nplcv',
                           get_cmd='SENS:VOLT:NPLC?',
                           get_parser=float,
                           set_cmd='SENS:VOLT:NPLC {:f}',
                           label='Voltage integration time')

        self.add_parameter('nplci',
                           get_cmd='SENS:CURR:NPLC?',
                           get_parser=float,
                           set_cmd='SENS:CURR:NPLC {:f}',
                           label='Current integration time')

        self.add_parameter('resistance',
                           get_cmd=self._get_read_output_protected,
                           get_parser=self._resistance_parser,
                           label='Resistance',
                           unit='Ohm',
                           docstring="Measure resistance from current and voltage "
                                     "Note that it is an error to read current "
                                     "and voltage with output off")

        self.add_parameter('T_cernox',
                            get_cmd=self._get_temperature_cernox,
                            get_parser=float,
                            label='Temperature Cernox',
                            unit='K')

        # self.calibration = np.genfromtxt("C:\\Users\\LAB-nanooptomechanic\\Documents\\Param\\Python_Programs\\tblg_measurement\\utils\\DataCernox_X170036.txt",skip_header = 2)
        
        self.write(':TRIG:COUN 1;:FORM:ELEM VOLT,CURR')
        # This line sends 2 commands to the instrument:
        # ":TRIG:COUN 1" sets the trigger count to 1 so that each READ? returns
        # only 1 measurement result.
        # ":FORM:ELEM VOLT,CURR" sets the output string formatting of the the
        # Keithley 2400 to return "{voltage}, {current}".
        # Default value on instrument reset is "VOLT, CURR, RES, TIME, STATUS";
        # however, resistance, status, and time are unused in this driver and
        # so are omitted.
        # These commands do not reset the instrument but do the minimal amount
        # to ensure that voltage and current parameters can be read from the
        # instrument, in the event that output formatting of the instrument was
        # previously changed to some other unknown state.
        self.connect_message()

    def _get_temperature_cernox(self):
        
        # T = self.calibration[:,0]
        # R = self.calibration[:,1]
        # T_cernox = T[np.argmin(np.abs(R-self.volt()/self.curr()))]
        
        R = self.volt()/self.curr()
        Z = np.log10(R)

        # From 1.4K to 14K 
        if R > 1034 and R <= 2.395e4:
            
            ZL = 2.96562273769
            ZU = ZU = 4.55345569955

            k = ((Z-ZL) - (ZU-Z)) / (ZU-ZL)
            
            coeff = [5.433503,
                    -6.259940, 
                    2.819944,
                    -1.060638,
                    0.335777, 
                    -0.085914, 
                    0.014611, 
                    0.000012, 
                    -0.001543, 
                    0.000629]

            T_cernox = 0
            
            for index, val in enumerate(coeff):
                T_cernox += val * np.cos(index*np.arccos(k)) 

            return T_cernox

        # From 14K to 80K 
        elif R > 231.1 and R <= 1034:

            ZL = 2.3169667156
            ZU = 3.0727151027

            k = ((Z-ZL) - (ZU-Z)) / (ZU-ZL)

            coeff = [42.225148, 
                    -37.834442,
                    8.633893,
                    -1.157072,
                    0.129661,
                    -0.007318,
                    -0.005425]

            T_cernox = 0
            
            for index, val in enumerate(coeff):
                T_cernox += val * np.cos(index*np.arccos(k)) 

            return T_cernox

        # From 80K to 325K 
        elif R <= 231.1 and R >= 60.68:

            ZL = 1.77653011834
            ZU = 2.41652337684

            k = ((Z-ZL) - (ZU-Z)) / (ZU-ZL)
            
            coeff = [176.380194,
                    -126.645073,
                    23.123875,
                    -3.406444,
                    0.645149,
                    -0.128011,
                    0.017116,
                    -0.003354,
                    0.001902]

            T_cernox = 0
            
            for index, val in enumerate(coeff):
                T_cernox += val * np.cos(index*np.arccos(k)) 

            return T_cernox
    
    def _get_read_output_protected(self) -> str:
        """
        This wrapper function around ":READ?" exists because calling
        ":READ?" on an instrument with output disabled is an error.
        So first we check that output is on and if not we return
        nan for volt, curr etc.
        """
        output = self.output.get_latest()
        if output is None:
            # if get_latest returns None we have
            # to ask the instrument for the status of output
            output = self.output.get()

        if output == 1:
            msg = self.ask(':READ?')
        else:
            raise RuntimeError("Cannot perform read with output off")
        return msg

    def _set_mode_and_sense(self, msg: str) -> None:
        # This helps set the correct read out curr/volt
        if msg == 'VOLT':
            self.sense('CURR')
        elif msg == 'CURR':
            self.sense('VOLT')
        else:
            raise AttributeError('Mode does not exist')
        self.write(f':SOUR:FUNC {msg:s}')

    def reset(self) -> None:
        """
        Reset the instrument. When the instrument is reset, it performs the
        following actions.

            Returns the SourceMeter to the GPIB default conditions.

            Cancels all pending commands.

            Cancels all previously send `*OPC` and `*OPC?`
        """
        self.write(':*RST')

    def _volt_parser(self, msg: str) -> float:
        fields = [float(x) for x in msg.split(',')]
        return fields[0]

    def _curr_parser(self, msg: str) -> float:
        fields = [float(x) for x in msg.split(',')]
        return fields[1]

    def _resistance_parser(self, msg: str) -> float:
        fields = [float(x) for x in msg.split(',')]
        res = fields[0] / fields[1]
        return res
