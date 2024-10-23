import numpy as np
from functools import partial
from typing import Optional, Any, Tuple

from qcodes import VisaInstrument
from qcodes import ChannelList, InstrumentChannel
from qcodes.utils import validators as vals
from qcodes.instrument.parameter import (
    ParameterWithSetpoints,
)

class ZNBChannel(InstrumentChannel):
    def __init__(
        self,
        parent,
        name: str,
        tracename,
    ) -> None:
        
        super().__init__(parent, name)
        self._tracename = tracename

        self.label_mapping = {
            "MLOG\n": "Magnitude",
            "MLIN\n": "Magnitude",
            "PHAS\n": "Phase",
            "UPH\n": "Unwrapped phase",
            "POL\n": "Complex Magnitude",
            "SMIT\n": "Complex Magnitude",
            "ISM\n": "Complex Magnitude",
            "SWR\n": "Standing Wave Ratio",
            "REAL\n": "Real Magnitude",
            "IMAG\n": "Imaginary Magnitude",
            "GDEL\n": "Delay",
            "COMP\n": "Complex Magnitude",
        }

        self.add_parameter(
            name="format",
            get_cmd=partial(self._get_format, tracename=self._tracename),
            val_mapping={
                "dB": "MLOG\n",
                "Linear Magnitude": "MLIN\n",
                "Phase": "PHAS\n",
                "Unwr Phase": "UPH\n",
                "Polar": "POL\n",
                "Smith": "SMIT\n",
                "Inverse Smith": "ISM\n",
                "SWR": "SWR\n",
                "Real": "REAL\n",
                "Imaginary": "IMAG\n",
                "Delay": "GDEL\n",
                "Complex": "COMP\n",
            },
        )

        self.add_parameter(
            'freq_axis',
            label='Frequency',
            unit='Hz',
            get_cmd=self.get_freq_axis,
            vals=vals.Arrays(shape=(self.parent.npts(),))
        )

        self.add_parameter(
            'trace',
            vals=vals.Arrays(shape=(self.parent.npts(),)),
            setpoints = (self.freq_axis,),
            unit=self.format(),
            label = f'{self.label_mapping[self.format.val_mapping[self.format()]]} {self.name_parts[-1]}',
            get_cmd=self._get_sweep_data,
            parameter_class=ParameterWithSetpoints,
        )

    def get_freq_axis(self):
        start = self.parent.start()
        stop = self.parent.stop()
        npts = self.parent.npts()
        return np.linspace(start=start, stop=stop, num=npts)

    def _get_format(self, tracename: str) -> str:
        n = self._parent._instrument_channel
        self.write(f"CALC{n}:PAR:SEL '{tracename}'")
        return self.ask(f"CALC{n}:FORM?")

    def _get_sweep_data(self, force_polar: bool = False) -> np.ndarray:
        # if force polar is set, the SDAT data format will be used.
        # Here the data will be transferred as a complex number
        # independent of the set format in the instrument.
        if force_polar:
            data_format_command = "SDAT"
        else:
            data_format_command = "FDAT"

            self.write(
                f"CALC{self._parent._instrument_channel}:PAR:SEL "
                f"'{self._tracename}'"
            )
            data_str = self.ask(
                f"CALC{self._parent._instrument_channel}:DATA?"
                f" {data_format_command}"
            )
        data = np.array(data_str.rstrip().split(",")).astype("float64")
        return data


class ZNB(VisaInstrument):
    """
    Args:
        name: instrument name
        address: Address of instrument probably in format
            'TCPIP0::192.168.15.100::inst0::INSTR'
        init_s_params: Automatically setup channels for all S parameters on the
            VNA.
        reset_channels: If True any channels defined on the VNA at the time
            of initialization are reset and removed.
        **kwargs: passed to base class
    """

    def __init__(
        self,
        name: str,
        address: str,
        **kwargs: Any,
    ) -> None:

        super().__init__(name=name, address=address, **kwargs)

        self._instrument_channel = n = 1
        self._min_freq = 10e3
        self._max_freq = 50e9

        self.add_parameter(
            name="power",
            label="Power",
            unit="dBm",
            get_cmd=f"SOUR{n}:POW?",
            set_cmd=f"SOUR{n}:POW {{:.4f}}",
            get_parser=float,
            vals=vals.Numbers(-80, 25),  #  used to be vals=vals.Numbers(-40, 25), before 2022/02/09. Changed by Chris.
        )
        self.add_parameter(
            name="bandwidth",
            label="Bandwidth",
            unit="Hz",
            get_cmd=f"SENS{n}:BAND?",
            set_cmd=f"SENS{n}:BAND {{:.0f}}",  # set cmd added by chris 2022/02/15. Haven't checked what happens when you give it a wrong value.
            get_parser=int,
            vals=vals.Enum(
                *np.append(10 ** 6,
                           np.kron([1, 1.5, 2, 3, 5, 7], 10 ** np.arange(6)))
            ),
            docstring="Measurement bandwidth of the IF filter. "
            "The inverse of this sets the integration "
            "time per point. "
            "There is an 'increased bandwidth option' "
            "(p. 4 of manual) that does not get taken "
            "into account here.",
        )
        self.add_parameter(
            name="avg",
            label="Averages",
            unit="",
            get_cmd=f"SENS{n}:AVER:COUN?",
            set_cmd=f"SENS{n}:AVER:COUN {{:.4f}}",
            get_parser=int,
            vals=vals.Ints(1, 5000),
        )
        self.add_parameter(
            name="start",
            get_cmd=f"SENS{n}:FREQ:START?",
            get_parser=float,
            vals=vals.Numbers(self._min_freq,
                              self._max_freq - 10),
        )
        self.add_parameter(
            name="stop",
            get_cmd=f"SENS{n}:FREQ:STOP?",
            get_parser=float,
            vals=vals.Numbers(self._min_freq + 1,
                              self._max_freq),
        )
        self.add_parameter(
            name="center",
            unit="Hz",
            get_cmd=f"SENS{n}:FREQ:CENT?",
            set_cmd=f"SENS{n}:FREQ:CENT {{:.4f}}",
            get_parser=float,
            vals=vals.Numbers(
                self._min_freq + 0.5, self._max_freq - 10
            ),
        )
        self.add_parameter(
            name="span",
            unit="Hz",
            get_cmd=f"SENS{n}:FREQ:SPAN?",
            set_cmd=f"SENS{n}:FREQ:SPAN {{:.4f}}",
            get_parser=float,
            vals=vals.Numbers(1,
                              self._max_freq - self._min_freq),
        )
        self.add_parameter(
            name="npts",
            get_cmd=f"SENS{n}:SWE:POIN?",
            get_parser=int,
        )

        self.add_parameter(
            name="rf_power",
            get_cmd="OUTP1?",
            set_cmd="OUTP1 {}",
            val_mapping={True: "1\n", False: "0\n"},
        )
        self.add_function("reset", call_cmd="*RST")
        self.add_function("restart_sweep", call_cmd="INITiate:IMMediate") # probably not very general. Works if VNA is set to single sweep mode. 
        self.add_function("tooltip_on", call_cmd="SYST:ERR:DISP ON")
        self.add_function("tooltip_off", call_cmd="SYST:ERR:DISP OFF")
        self.add_function("cont_meas_on", call_cmd="INIT:CONT:ALL ON")
        self.add_function("cont_meas_off", call_cmd="INIT:CONT:ALL OFF")
        self.add_function("update_display_once", call_cmd="SYST:DISP:UPD ONCE")
        self.add_function("update_display_on", call_cmd="SYST:DISP:UPD ON")
        self.add_function("update_display_off", call_cmd="SYST:DISP:UPD OFF")
        self.add_function("rf_off", call_cmd="OUTP1 OFF")
        self.add_function("rf_on", call_cmd="OUTP1 ON")

        self.add_traces()
        self.connect_message()

    def add_traces(self):
        traces = self.ask('CONFigure:TRACe:CATalog?').rstrip().strip("'").split(',')[1::2]


        channels = ChannelList(self, 'matrix', ZNBChannel)

        self.add_submodule('matrix', channels)
        for trace in traces:
            parameter = self.ask(f"CALC{self._instrument_channel}:PAR:MEAS? "
                    f"'{trace}'").rstrip().strip("'")

            channel = ZNBChannel(self, parameter, trace)
            channels.append(channel)
            self.add_submodule(parameter, channel)


