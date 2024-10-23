import logging
import numpy as np
from functools import partial
from typing import Any, Union
from si_prefix import si_format
from qcodes import VisaInstrument
from qcodes import ChannelList, InstrumentChannel
from qcodes.utils import validators as vals
from qcodes.instrument.parameter import (
    ParameterWithSetpoints,
)

log = logging.getLogger(__name__)

class ZNBTrace(InstrumentChannel):
    def __init__(
        self,
        parent: "ZNB",
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

        self.unit_mapping = {
            "MLOG\n": "dB",
            "MLIN\n": "",
            "PHAS\n": "rad",
            "UPH\n": "rad",
            "POL\n": "",
            "SMIT\n": "",
            "ISM\n": "",
            "SWR\n": "U",
            "REAL\n": "U",
            "IMAG\n": "U",
            "GDEL\n": "S",
            "COMP\n": "",
        }

        self.val_mapping = {
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
        }

        self.add_parameter(
            name="format",
            get_cmd=partial(self._get_format, tracename=self._tracename),
            set_cmd=self._set_format,
            val_mapping=self.val_mapping,
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
            unit= self.unit_mapping[self.val_mapping[self.format()]],
            label = f'{self.label_mapping[self.format.val_mapping[self.format()]]} {self.name_parts[-1]}',
            get_cmd=self._get_sweep_data,
            parameter_class=ParameterWithSetpoints,
        )

    def _set_sweep_type(self, tracename: str) -> None:
        channel = self._instrument_channel
        self.write(f"SENS{channel}:SWE:TYPE {tracename}")

    def _set_cw_frequency(self, val: float) -> None:
        channel = self._instrument_channel
        self.write(f"SENS{channel}:FREQ:CW {val:.7f}")

    def _enable_averaging(self, val: str) -> None:
        channel = self._instrument_channel
        self.write(f"SENS{channel}:AVER:STAT {val}")

    def _enable_auto_sweep_time(self, val: str) -> None:
        channel = self._instrument_channel
        self.write(f"SENS{channel}:SWE:TIME:AUTO {val}")

    def get_freq_axis(self):
        start = self.parent.start()
        stop = self.parent.stop()
        npts = self.parent.npts()
        return np.linspace(start=start, stop=stop, num=npts)

    def _get_format(self, tracename: str) -> str:
        n = self.parent._instrument_channel
        self.write(f"CALC{n}:PAR:SEL '{tracename}'")
        val = self.ask(f"CALC{n}:FORM?")

        return val

    def _set_format(self, val: str) -> None:
        n = self.parent._instrument_channel
        self.write(f"CALC{n}:PAR:SEL '{self._tracename}'")
        self.write(f"CALC{n}:FORM {val}")

        self.trace.unit = self.unit_mapping[val]
        self.trace.label = f"{self.short_name} {self.label_mapping[val]}"

    def _get_sweep_data(self, force_polar: bool = False) -> np.ndarray:
        # if force polar is set, the SDAT data format will be used.
        # Here the data will be transferred as a complex number
        # independent of the set format in the instrument.
        if force_polar:
            data_format_command = "SDAT"
        else:
            data_format_command = "FDAT"

        self.write(
            f"CALC{self.parent._instrument_channel}:PAR:SEL "
            f"'{self._tracename}'")

        data_str = self.ask(
            f"CALC{self.parent._instrument_channel}:DATA?"
            f" {data_format_command}")

        if force_polar:
            data_complex = np.array(data_str.rstrip().split(",")).astype("float64") 
            data = []
            for idx in np.arange(0,len(data_complex),2):
                data.append(data_complex[idx]+data_complex[idx+1]*1j)
            data = np.asarray(data)    

        else:
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
        self._min_freq = 10e6
        self._max_freq = 14e9
        self._min_power = -300
        self._max_power = 0

        self.generate_traces()

        self.add_parameter(
            name="power",
            label="Power",
            unit="dBm",
            get_cmd=f"SOUR{n}:POW?",
            set_cmd=f"SOUR{n}:POW {{:.4f}}",
            get_parser=float,
            vals=vals.Numbers(self._min_power, self._max_power),  #  used to be vals=vals.Numbers(-40, 25), before 2022/03/23. Changed by Chris.
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
            set_cmd=f"SENS{n}:FREQ:START {{:.7f}}",
            vals=vals.Numbers(self._min_freq,
                              self._max_freq),
        )
        self.add_parameter(
            name="stop",
            get_cmd=f"SENS{n}:FREQ:STOP?",
            get_parser=float,
            set_cmd=f"SENS{n}:FREQ:STOP {{:.7f}}",
            vals=vals.Numbers(self._min_freq,
                              self._max_freq),
        )
        self.add_parameter(
            name="center",
            unit="Hz",
            get_cmd=f"SENS{n}:FREQ:CENT?",
            set_cmd=f"SENS{n}:FREQ:CENT {{:.4f}}",
            get_parser=float,
            vals=vals.Numbers(
                self._min_freq, self._max_freq
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
            set_cmd=f"SENS{n}:SWE:POIN {{:.7f}}",
            get_parser=int,
        )

        self.add_parameter(
            name="rf_power",
            get_cmd="OUTP1?",
            set_cmd="OUTP1 {}",
            val_mapping={True: "1\n", False: "0\n"},
        )

        self.add_parameter(
            name="S21_complex",
            get_cmd=self._get_S21_complex,
            get_parser=complex,
        )

        self.add_parameter(
            name="S11_complex",
            get_cmd=self._get_S11_complex,
            get_parser=complex,
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
        trace_num = len(traces)
        trace_list = ChannelList(self, 'trace_list', ZNBTrace)
        self.add_submodule('trace_list', trace_list)

        for val, trace in enumerate(traces):
            self.write(f"CALC{self._instrument_channel}:PAR:SEL '{trace}'") # necessary DO NOT DELETE IT!!!!

            if trace_num > 4: # 4 corresponds to S11,S12,S21 and S22 (only one format)
                if val <= trace_num/2-1:
                    self.write(f"CALC{self._instrument_channel}:FORM MLOG")
                if val > trace_num/2-1:
                    self.write(f"CALC{self._instrument_channel}:FORM UPH")

            format = self.ask(f"CALC{self._instrument_channel}:FORM?").rstrip().lower()
            s_param = self.ask(f"CALC{self._instrument_channel}:PAR:MEAS? "
                    f"'{trace}'").rstrip().strip("'")

            trace_str = f"{s_param}_{format}"

            trace_obj = ZNBTrace(self, trace_str, trace)
            trace_list.append(trace_obj)
            self.add_submodule(trace_str, trace_obj)
            self.display_traces(val+1, trace)

    def clear_channels(self) -> None:
        """
        Remove all channels from the instrument and channel list and
        unlock the channel list.
        """
        self.write("CALCulate:PARameter:DELete:ALL")
        for submodule in self.submodules.values():
            if isinstance(submodule, ChannelList):
                submodule._channels = []
                submodule._channel_mapping = {}
                submodule._locked = False

    def generate_traces(self):
        """
        It generates 8 traces corresponding to magnitude and
        phase of: S11,S12,S21 and S22
        """
        self.clear_channels()
        n = self._instrument_channel
        for format_str in ['Mag','PhaseUW']:
            for i in range(1, n + 2):
                for j in range(1, n + 2):
                    param = 'S' + str(i) + str(j) 
                    ch_name = param + '_' + format_str
                    self.write(f"CALC{self._instrument_channel}:PAR:SDEF '{ch_name}','{param}'")

        # if we want 4 traces (S11, S12, S21 and S22) on could use the command: "CALC1:PAR:DEF:SGR 1,2" 

    def display_traces(self, window, trace):
        self.write(f"DISP:WIND{window}:STAT ON")
        self.write(f"DISP:WIND{window}:TRAC{window}:FEED '{trace}'")

    def sweep_init(self,
                    mode: str ='StartStop',
                    start: Union[int,float] = 10e6,
                    stop: Union[int,float] = 14e9,
                    center: Union[int,float] = 3e9,
                    span: Union[int,float] = 10e6,
                    power: Union[int,float] = -80,
                    bandwidth: Union[int,float] = 100,
                    npts: int = 1e3,
                    output_att: Union[int, float] = 0):
        """
        Accepted values fo "mode": "StartStop" or "CenterSpan"
        """
        precision = 4

        if mode == 'StartStop':
            self.start(start)
            self.stop(stop)
            freq1 = f"fstart:{si_format(start, precision=precision).replace(' ','')}Hz"
            freq2 = f"fstop:{si_format(stop, precision=precision).replace(' ','')}Hz"

        if mode == 'CenterSpan':
            self.center(center)
            self.span(span)
            freq1 = f"fcenter:{si_format(center, precision=precision).replace(' ','')}Hz"
            freq2 = f"fspan:{si_format(span, precision=precision).replace(' ','')}Hz"

        self.power(power)
        self.bandwidth(bandwidth)
        self.npts(npts)
        self.rf_on()

        sweep_summary = f"VNAmeasurement_{freq1}_{freq2}_Power:{power}dBm_BW:{si_format(bandwidth,precision=0).replace(' ','')}Hz_pts:{si_format(npts,precision=0).replace(' ','')}_att:{output_att}dB_T=50mK"
        console_print = "\n============ Sweep summary =============\n"+sweep_summary.replace('_','\n')+"\n========================================"
        print(console_print)

        return sweep_summary
    
    def _get_S21_complex(self):
        R = 10**(self.S21_mlog.trace()/20)
        phase = self.S21_uph.trace()*np.pi/180

        S21_complex = R*np.exp(1j*phase)
        
        return S21_complex

    def _get_S11_complex(self):
        R = 10**(self.S11_mlog.trace()/20)
        phase = self.S11_uph.trace()*np.pi/180

        S11_complex = R*np.exp(1j*phase)
        
        return S11_complex

    def sweep_time(self):
        t = self.ask('SWEep:TIME?')
        return t
