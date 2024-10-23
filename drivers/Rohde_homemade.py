from typing import Any
from qcodes import VisaInstrument, validators as vals
from qcodes.utils.helpers import create_on_off_val_mapping


class RohdeSchwarz_SMB100A(VisaInstrument):
    """
    This is the QCoDeS driver for the Rohde & Schwarz SGS100A signal generator.

    Status: beta-version.

    .. todo::

        - Add all parameters that are in the manual
        - Add test suite
        - See if there can be a common driver for RS mw sources from which
          different models inherit

    This driver will most likely work for multiple Rohde & Schwarz sources.
    it would be a good idea to group all similar RS drivers together in one
    module.

    Tested working with

    - RS_SMB100A

    This driver does not contain all commands available for the RS_SGS100A but
    only the ones most commonly used.
    """

    def __init__(self, name: str, address: str, **kwargs: Any) -> None:
        super().__init__(name, address, terminator='\n', **kwargs)
        
# --------------------------- RF Frequency Output ----------------------------

        self.add_parameter(name='frequency',
                           label='Frequency',
                           unit='Hz',
                           get_cmd='SOUR:FREQ?',
                           set_cmd='SOUR:FREQ {:.2f}',
                           get_parser=float,
                           vals=vals.Numbers(9e3, 3.2e9))

        self.add_parameter(name='phase',
                           label='Phase',
                           unit='deg',
                           get_cmd='SOUR:PHAS?',
                           set_cmd='SOUR:PHAS {:.2f}',
                           get_parser=float,
                           vals=vals.Numbers(0, 360))

        self.add_parameter(name='power',
                           label='Power',
                           unit='dBm',
                           get_cmd='SOUR:POW?',
                           set_cmd='SOUR:POW {:.2f}',   
                           get_parser=float,
                           vals=vals.Numbers(-120, 25))

        self.add_parameter(name='status',
                           label='RF Output',
                           get_cmd=':OUTP:STAT?',
                           set_cmd=':OUTP:STAT {}',
                           val_mapping=create_on_off_val_mapping(on_val='1',
                                                                 off_val='0'))
# --------------------------- Low Frequency Output ---------------------------

        self.add_parameter(name='lfo_frequency',
                           label='Low frequency output frequency',
                           unit='Hz',
                           get_cmd=':LFO:FREQ?',
                           set_cmd=':LFO:FREQ {:.2f}',
                           get_parser=float,
                           vals=vals.Numbers(1, 1e6))

        self.add_parameter(name='lfo_mode',
                           label='Low frequency output mode',
                           get_cmd=':LFO:FREQ:MODE?',
                           set_cmd=':LFO:FREQ:MODE {}',
                           set_parser=str,
                           vals=vals.Enum('FIXED', 'AUTO', 'SING'))       
        
        self.add_parameter(name='lfo_status',
                           label='Low frequency output status',
                           get_cmd=':LFO:STAT?',
                           set_cmd=':LFO:STAT {}',
                           set_parser=str,
                           val_mapping=create_on_off_val_mapping(on_val='1',
                                                                 off_val='0')) 

# --------------------------- Modulation -------------------------------------

        self.add_parameter(name='FM_deviation',
                           label='Frequency modulated deviation',
                           unit='Hz',
                           get_cmd=':FM?',
                           set_cmd=':FM {:.2f}',
                           get_parser=float,
                           vals=vals.Numbers(1, 16e6))

        self.add_parameter(name='FM_source',
                           label='Frequency modulated source',
                           get_cmd=':FM:SOUR?',
                           set_cmd=':FM:SOUR {}',
                           set_parser=str,
                           vals=vals.Enum('INT', 'EXT'))       

        self.add_parameter(name='FM_status',
                           label='Status of the modulation output',
                           get_cmd=':MOD:STAT?',
                           set_cmd=':MOD:STAT {}',
                           set_parser=str,
                           val_mapping=create_on_off_val_mapping(on_val='1',
                                                                 off_val='0'))   

# # --------------------------- Oscillator -----------------------------------  

        self.add_parameter(name='osc_source',
                           label='Reference oscillator being used',
                           get_cmd='ROSCillator:SOUR?',
                           set_cmd='ROSCillator:SOUR {}',
                           set_parser=str,
                           vals=vals.Enum('INT','EXT'))
                        
        self.add_parameter(name='osc_freq',
                           label='Reference oscillator frequency',
                           get_cmd='ROSCillator:EXT:FREQ?',
                           set_cmd='ROSCillator:EXT:FREQ {:.2}',
                           vals=vals.Numbers(9e3, 3.2e9))
        
# ------------------------------ Other --------------------------------------- 
        
        self.add_parameter(name='wait',
                           label='Waiting time',
                           unit='ms',
                           set_cmd='SYST:WAIT {}',
                           set_parser=int,         # the set value will be an integer
                           vals=vals.Ints(1,1000)) # it make sure you provide integers

                           # Why can I get it??? I can get the last value that I gave but it is not running another SCPI command.
                           # Don't use f-string for get comands!!!

        self.connect_message()

    def on(self) -> None: # it only indicates what do you expect as a return (it helps the compiler if you get errors)
        self.status('on') 

    def off(self) -> None:
        self.status('off')

    def reset(self) -> None:
        self.write('*RST')
    
    def clear_error(self) -> None:
        self.write('*CLS')
    
    def error(self):
        err = self.ask(':SYST:ERR?')
        return err
    
    def shut_down(self) -> None:
        self.write(':SYST:SHUT')
    
    def run_self_tests(self):
        test = self.ask('*TST?')
        return test

    def frequency_up(self) -> None:
        self.wirte('FREQ:UP')
    
    def frequency_down(self) -> None:
        self.wirte('FREQ:DOWN')

    def reset_sweep(self) -> None:
        self.write('SWE:RES')

# --------------------- RF sweep modes ---------------------------------------             
#     def level_sweep(self, start, stop, step=None, dwell=None,mode=None, shape=None):
#         '''Info page 182-188 of the Operating Manual in MANUAL operation and 
#             425-428 for SCPI commands in REMOTE operation. There are different
#             ways of setting the paramaters, here we chose START, STOP, STEP,
#             DWELL TIME; SHAPE and MODE.If we don't specify the mode this will
#             be SINGLE. The mode accepts 'AUTO' for continuos sweep. Functions
#             for *retrace* and *points* also cool.'''
            
#         self.sig_gen.write(':POW:STAR '+str(start))
#         self.sig_gen.write(':POW:STOP '+str(stop))
        
#         if not step:
#             pass
#         if step:    
#             self.sig_gen.write(':SWE:POW:STEP '+str(step))
        
#         if not dwell:
#             pass
#         if dwell:
#             self.sig_gen.write(':SWE:POW:DWEL '+str(dwell))     
            
#         if not shape:
#             self.sig_gen.write(':SWE:SHAP SAWTooth')     
#         if shape=='triangle':
#             self.sig_gen.write(':SWE:SHAP TRIangle')     
#         if shape=='sawtooth':
#             self.sig_gen.write(':SWE:SHAP SAWTooth')    
            
#         if not mode:
#             self.sig_gen.write('TRIG:PSW:SOUR SING')          
#             self.sig_gen.write(':POW:MODE SWEep')
#             self.sig_gen.write('SWE:POW:EXEC')
#         if mode:
#             self.sig_gen.write(':TRIG:PSW:SOUR '+str(mode))  
#             self.sig_gen.write(':POW:MODE SWEep')

    # def freq_sweep(self, start, stop, step=None, dwell=None, mode=None, shape=None):
    #     '''Info page 175-182 of the Operating Manual in MANUAL operation and 
    #         417-425 for SCPI commands in REMOTE operation. There are different
    #         ways of setting the paramaters for the frequency sweep: span, 
    #         center frequency, number of points, etc. Here we chose START,
    #         STOP, STEP (max 1dbm), DWELL TIME, SHAPE and MODE. If we don't 
    #         specify the mode this will be SINGLE. The mode accepts 'AUTO' for
    #         continuos sweep. Functions for *retrace* and *points* also cool.'''
        
    #     self.sig_gen.write(':FREQ:MODE SWEep')
    #     self.sig_gen.write('SWE:FREQ:MODE AUTO') 
    #     self.sig_gen.write(':FREQ:STAR '+str(start))
    #     self.sig_gen.write(':FREQ:STOP '+str(stop))

    #     if not step:
    #         pass
    #     if step:    
    #         self.sig_gen.write(':SWE:STEP:LIN '+str(step))
        
    #     if not dwell:
    #         pass
    #     if dwell:
    #         self.sig_gen.write(':SWE:DWEL '+str(dwell))     
         
    #     if not shape:
    #         self.sig_gen.write(':SWE:SHAP SAWTooth')     
    #     if shape=='triangle':
    #         self.sig_gen.write(':SWE:SHAP TRIangle')     
    #     if shape=='sawtooth':
    #         self.sig_gen.write(':SWE:SHAP SAWTooth')
        
    #     if not mode:
    #         self.sig_gen.write(':TRIG:FSW:SOUR SING')
    #         self.sig_gen.write('SWE:EXEC')
            
    #     if mode:
    #         self.sig_gen.write(':TRIG:FSW:SOUR '+str(mode)) 
        
