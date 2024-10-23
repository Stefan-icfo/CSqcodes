# This Python file uses the following encoding: utf-8
# Loick Le Guevel, 2019
# Etienne Dumur <etienne.dumur@gmail.com>, 2021

from typing import Union, Tuple, Any
from math import ceil
from time import sleep
from tqdm import tqdm
from qcodes.instrument.base import Instrument
from qcodes.instrument.channel import InstrumentChannel, ChannelList
from qcodes.instrument.channel import MultiChannelInstrumentParameter
from qcodes.instrument.visa import VisaInstrument
from qcodes.utils import validators as vals


class iTestChannel(InstrumentChannel):
    """
    A single channel of iTest.
    """

    def __init__(self, parent: Instrument,
                 name: str,
                 chan_num: int) -> None:
        """
        Args:
            parent: The instrument to which the channel is attached.
            name: The name of the channel.
            chan_num: The number of the channel in question.
        """
        super().__init__(parent, name)

        self.chan_num = chan_num
        # Get channel id
        i = ceil(chan_num/4.)
        c = chan_num-(i-1)*4
        self.chan_id = 'i{};c{};'.format(i, c)

        self.add_parameter('v',
                           label='Channel {} voltage'.format(chan_num),
                           unit='V',
                           docstring='Voltage of the channel in volt.',
                           get_cmd=f"{self.chan_id}MEAS:VOLT?",
                           get_parser=float,
                           set_cmd=self.chan_id +
                           "VOLT {:.8f};TRIG:INPUT:INIT",
                           vals=vals.Numbers(-12, 12)
                           )

        self.add_parameter('i',
                           label='Channel {} current'.format(chan_num),
                           unit='A',
                           docstring='Current of the channel in ampere.',
                           get_cmd=f"{self.chan_id}MEAS:CURR?",
                           get_parser=float
                           )

        self.add_parameter('ramp_slope',
                           label='Channel {} ramp slope'.format(chan_num),
                           unit='V/ms',
                           docstring='Slope of the ramp in V/ms.',
                           get_cmd=f"{self.chan_id}VOLT:SLOP?",
                           get_parser=float,
                           set_cmd=self.chan_id + "VOLT:SLOP {:.8f}",
                           vals=vals.Numbers(min_value=1.2e-6, max_value=1.2e-3)  # min value instr is 1.2e-6. max is 0.1.
                           )

        self.add_parameter('output_mode',
                           label='Channel {} output mode'.format(chan_num),
                           unit='',
                           docstring='Mode of the output {exp, ramp, stair, step}.',
                           get_cmd=self.chan_id + 'trig:input?',
                           get_parser=str,
                           set_cmd=self._set_safe_output_mode, #self.chan_id + 'trig:input {}',
                           set_parser=str,
                        #    val_mapping={'exp': 0, 'ramp': 1,
                        #                 'stair': 2, 'step': 3},
                           val_mapping={'exp': 0, 'ramp': 1},
                        #    vals=vals.Enum('ramp', 'exp', 'stair', 'step')
                           vals=vals.Enum('ramp', 'exp')
                           )

        self.add_parameter('v_range',
                           label='Channel {} voltage range'.format(chan_num),
                           unit='V',
                           docstring='Range of the channel in volt.',
                           set_cmd=self.chan_id + "VOLT:RANGE {}",
                           set_parser=float,
                           get_cmd=f"{self.chan_id}VOLT:RANGE?",
                           get_parser=lambda val: val[:-2],
                           vals=vals.Enum(1.2, 12)
                           )

        self.add_parameter('state',
                           docstring='State of the channel {on, off}.',
                           unit='',
                           get_cmd=f"{self.chan_id}OUTP?",
                           get_parser=str,
                           set_cmd=self.chan_id + "OUTP {}",
                           set_parser=str,
                           val_mapping={'on': 1, 'off': 0},
                           vals=vals.Enum('on', 'off')
                           )

        self.add_parameter('pos_sat',
                           get_cmd=f"{self.chan_id}'VOLT:SAT:POS?",
                           get_parser=str,
                           set_cmd=self.chan_id + "VOLT:SAT:POS {}",
                           set_parser=lambda val: f"{val:.8f}" if isinstance(
                               val, (int, float)) else "MAX",
                           )

        self.add_parameter('neg_sat',
                           get_cmd=f"{self.chan_id}'VOLT:SAT:NEG?",
                           get_parser=str,
                           set_cmd=self.chan_id + "VOLT:SAT:NEG {}",
                           set_parser=lambda val: f"{val:.8f}" if isinstance(
                               val, (int, float)) else "MIN",
                           )

        self.add_parameter('bilt_name',
                           docstring='The name of the channel',
                           unit='',
                           set_cmd=self.chan_id + "chan:name '{}'",
                           set_parser=str,
                           initial_value=f'Chan{chan_num:02d}')

        self.add_parameter('trig_delay',
                           label='Channel {} trig dalay'.format(chan_num),
                           unit='ms',
                           docstring='Trigger delay in miliseconds',
                           get_cmd=f"{self.chan_id}'trig:input:delay?",
                           get_parser=int,
                           set_cmd=self.chan_id + "trig:input:delay {}",
                           set_parser=int,
                           vals=vals.Numbers(min_value=0)
                           )

        self.add_parameter('trig_ready_amplitude',
                           label='Channel {} trigger ready amplitude'.format(
                               chan_num),
                           unit='V',
                           docstring='Trigger ready amplitude of the channel in volts.',
                           get_cmd=f"{self.chan_id}'TRIGger:READY:AMPLitude?",
                           get_parser=float,
                           set_cmd=self.chan_id +
                           "TRIGger:READY:AMPLitude {:.8f}",
                           set_parser=float,
                           vals=vals.Numbers(min_value=1.2e-6, max_value=1),
                           )

        self.add_parameter('voltage_status',
                           unit="",
                           label='Channel {} voltage status'.format(chan_num),
                           docstring='Voltage status of the channel.',
                           get_cmd=f"{self.chan_id}'VOLTage:STATus?",
                           get_parser=float,
                           )

        self.add_parameter('ready_status',
                           unit="",
                           label='Channel {} ready status'.format(chan_num),
                           docstring='Ready status of the channel.',
                           get_cmd=f"{self.chan_id}'TRIGger:READY?",
                           get_parser=int,
                           )

        self.add_parameter('step_amplitude',
                           label='Channel {} step amplitude'.format(chan_num),
                           unit='V',
                           docstring='Step amplitude of the channel in volts.',
                           get_cmd=f"{self.chan_id}'VOLTage:STEP:AMPLitude?",
                           get_parser=float,
                           set_cmd=self.chan_id +
                           "VOLTage:STEP:AMPLitude {:.8f}",
                           set_parser=float,
                           vals=vals.Numbers(min_value=1.2e-6, max_value=12),
                           )

        self.add_parameter('step_width',
                           label='Channel {} step width'.format(chan_num),
                           unit='ms',
                           docstring='Step width of the channel in miliseconds.',
                           get_cmd=f"{self.chan_id}'VOLTage:STEP:WIDTH?",
                           get_parser=int,
                           set_cmd=self.chan_id + "VOLTage:STEP:WIDTH {}",
                           set_parser=int,
                           vals=vals.Ints(min_value=5),
                           )
        

    def _set_safe_output_mode(self, mode):
        '''
        This function will set step and delay if we are in the exp mode
        and remove them if we are on other modes
        '''

        if mode=='0': #exp mode
            # delay = 50e-3
            # step = self.ramp_slope() * delay * 1000
            # step = min(step, 1)  # UUUUGE BUG FIX ME. (changed 21/04/21  from min(step, .5e-3))

            self.v.step = 500e-3 # this is a MAXIMUM STEP not necesarily the step size.
            self.v.inter_delay = 50e-3
        else: #all others
            self.v.step = 0
            self.v.inter_delay = 0

        self.write(f"{self.chan_id}trig:input {mode}")



    def block_until_set(self):
        while True:
            status = self.get('voltage_status')
            tqdm.write(f"{self.chan_id}  {status*100:.4f}% done", end='\r')
            if  status == 1.0:
                break

    def block_until_ready(self):
        while True:
            if self.get('ready_status') == 1:
                break

    def set_voltage_and_block(self, value):
            self.set('v', value)
            self.block_until_set()


class iTestMultiChannelParameter(MultiChannelInstrumentParameter):
    """
    """

    def __init__(self, channels, param_name, *args, **kwargs):
        super().__init__(channels, param_name, *args, **kwargs)


class ITest(VisaInstrument):
    """
    This is the QCoDeS python driver for the iTest device from Bilt.
    """

    def __init__(self, name: str,
                 address: str,
                 num_chans: int = 8,
                 **kwargs: Any) -> None:
        """
        Instantiate the instrument.

        Args:
            name: The instrument name used by qcodes
            address: The VISA name of the resource
            num_chans: Number of channels to assign. Default: 16
            init_start: If true set all channels to 0V, 1.2V range and switch
                then on.

        Returns:
            ITest object
        """
        super().__init__(name, address=address,
                         terminator='\n',
                         device_clear=False,
                         **kwargs)

        self.idn = self.get_idn()
        self.num_chans = num_chans
        self.chan_range = range(1, self.num_chans+1)

        # Create the channels
        channels = ChannelList(parent=self,
                               name='Channels',
                               chan_type=iTestChannel,
                               multichan_paramclass=iTestMultiChannelParameter)

        for i in self.chan_range:

            channel = iTestChannel(self, name='chan{:02}'.format(i),
                                   chan_num=i)
            channels.append(channel)
            self.add_submodule('ch{:02}'.format(i), channel)

        channels.lock()
        self.add_submodule('channels', channels)

        self.init_safe_ramp()

        self.connect_message()



    def init_safe_ramp(self, safe_ramp=1.2e-6):
        [chan.ramp_slope(safe_ramp) for chan in self.channels]

        #TODO unsafe because it only works if we were in the exp or ramp mode previously
        [chan.output_mode('ramp') for chan in self.channels]


    def go_to_begging_sweep(self, sweeps_tuple, fixed_tuple):
        
        if not isinstance(sweeps_tuple, tuple):
            raise ValueError

        if not isinstance(fixed_tuple, tuple):
            raise ValueError


        #first set end values...----------------
        for params_sweep in sweeps_tuple:
            params_sweep.set(params_sweep[0])

        for params_sweep in fixed_tuple:
            params_sweep.v(params_sweep.v.fixed_value)
        #----------------------------------------

        try:
            #... then wait do get there-----
            for params_sweep in sweeps_tuple:
                params_sweep.parameter._instrument.block_until_set()

            for params_sweep in fixed_tuple:
                params_sweep.block_until_set()
            #--------------------------------
        except KeyboardInterrupt:
            print('beggining sweep canceled canceled')
        finally:
            for chan in self.channels:
                chan.v(chan.v())


    def ramp_down(self, channels, rate=1.2e-6):
        for chan in channels:
            chan.output_mode('ramp')
            if chan.output_mode() != 'ramp':
                print('problem')
                return
            chan.ramp_slope(rate)
            chan.v(0)

        try:
            for chan in channels:
                chan.block_until_set()
        except KeyboardInterrupt:
            print('ramping down canceled')
        finally:
            for chan in channels:
                chan.v(chan.v())