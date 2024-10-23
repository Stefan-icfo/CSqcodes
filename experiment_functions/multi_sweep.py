from qcodes.instrument.parameter import DelegateParameter
from qcodes.utils.helpers import make_sweep
from si_prefix import si_format


class SweepMultiParam:
    '''
    Very similary to SweepFixedValues, this object will behave like a list of tuples,
    and it will have folowing methods
        set method: this will set all params using the list of tuples 
    '''

    def __init__(self, channel_list):
        self._snapshot = {}
        self._value_snapshot = []

        self._ch_list = []  # list of channels to be swept
        self._ch_sweep_list = []  # list of list for each channel from above

        for chan in channel_list:

            start = chan['sweep_values']['start']
            stop = chan['sweep_values']['stop']

            if start is None:
                raise ValueError('You have to provide start point')
            if stop is None:
                raise ValueError('You have to provide stop point')

            # remember which channel is used
            self._ch_list.append(chan['channel'])

            # remember all the values
            keys = make_sweep(**chan['sweep_values'])
            self._ch_sweep_list.append(keys)

            # validation the values for particular channel
            if hasattr(chan['channel'], 'validate'):
                for value in keys:
                    chan['channel'].validate(value)

        # check if all sweeps have the same lenght
        for value_list in self._ch_sweep_list:
            if len(value_list) != len(self._ch_sweep_list[0]):
                raise ValueError('All channels must have the same size')

        self._zipped_values = list(zip(*self._ch_sweep_list))


        ch_label = ""
        for i, chan in enumerate(self._ch_list):
            ch_label+=f"{chan.instrument.name_parts[-1]}({si_format(self._ch_sweep_list[i][0])},{si_format(self._ch_sweep_list[i][-1])}),"
        ch_label = ch_label[:-1]

        self.parameter = DelegateParameter(name=self._ch_list[0].full_name, source=self._ch_list[0], label=ch_label)

    def reverse(self) -> None:
        """ Reverse SweepFixedValues in place. """
        self._zipped_values.reverse()

    def __iter__(self):
        return iter(self._zipped_values)

    def __getitem__(self, key: slice):
        return self._zipped_values[key]

    def __len__(self):
        return len(self._zipped_values)

    def set(self, values):
        if len(values) != len(self._ch_list):
            raise ValueError('lenght of values has to match channel list')

        for chan, value in zip(self._ch_list, values):
            chan(value)