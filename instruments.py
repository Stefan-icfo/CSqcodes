import database

import qcodes as qc

from qcodes.tests.instrument_mocks import DummyInstrument, DummyInstrumentWithMeasurement

#from qcodes.instrument_drivers.stanford_research.SR830 import SR830

from qcodes.instrument_drivers.tektronix.Keithley_2000 import Keithley_2000

from qcodes.instrument_drivers.tektronix.Keithley_2400 import Keithley_2400

from qcodes.instrument_drivers.oxford.triton import Triton

from qcodes.instrument_drivers.tektronix.Keithley_2450 import Keithley2450

#from drivers.zurich import MyZurich

from drivers.zurich_CS import MyZurich # from Param's code, if original Triton 2 Zujrich is used again

from drivers.bilt import ITest

from drivers.manual_instr import Manual_Inst

#from drivers.vna3 import ZNB

#import qcodes.instrument_drivers.rohde_schwarz.ZNB as ZNB

#from qcodes.instrument_drivers.rohde_schwarz.ZNB import ZNBChannel

import drivers.QDAC2 

#from zhinst.qcodes import HF2

from drivers.Rohde_homemade import RohdeSchwarz_SMB100A
from Qdac_test_driver import QDAC_test
station = qc.Station()


# #define instuments


#qdac = drivers.QDAC2.QDac2('QDAC', visalib='@py', address="ASRL6:INSTR")
#qdac = drivers.QDAC2.QDac2(name = 'QDAC', address="192.168.1.253" , port=5025)
qdac = drivers.QDAC2_CS.QDac2_CS(name = 'QDAC', address = "TCPIP::192.168.1.253::5025::SOCKET")
#qdac = drivers.QDAC2.QDac2(name = 'QDAC', address = "TCPIP::192.168.1.253::5025::SOCKET")
station.add_component(qdac)

#def jls_extract_def():exit

#    #, tmpfile='./drivers/tritonReg.reg'
#    return 

Triton = Triton(name='Triton', address='192.168.1.20', port=33576) # , tmpfile='./drivers/tritonReg.reg' = jls_extract_def()
station.add_component(Triton)

#vna = ZNB('vna', 'GPIB0::20::INSTR')
#station.add_component(vna)

#define instuments

 #dac = DummyInstrument('dac', gates=['ch1', 'ch2'])
 #station.add_component(dac)

# dmm = DummyInstrumentWithMeasurement('dmm', setter_instr=dac)
#station.add_component(dmm)


manual = Manual_Inst(name='manual')
station.add_component(manual)

zurich = MyZurich(manual,"DEV20039", "localhost", name="zurich")#was Dev 2187 for original Triton 2 zurich, DEV2102 for Triton 1 Zurich
#zurich = MyZurich('zurich', visalib='@py', address="ASRL4:INSTR")#
station.add_component(zurich)
#zurich = HF2("DEV1039", "localhost")
#station.add_component(zurich)

#bilt = ITest(name='bilt', address='GPIB0::5::INSTR', num_chans=4)
#station.add_component(bilt)

#lia1 = SR830(name='lia1', address='GPIB0::8::INSTR')
#station.add_component(lia1)

#k2000 = Keithley_2000(name='k2000', address='GPIB0::11::INSTR')
#station.add_component(k2000)

#k2450 = Keithley2450(name='k2450', address='GPIB0::18::INSTR')
#station.add_component(k2450)

keithley2400 = Keithley_2400(name='keithley2400', address='GPIB0::24::INSTR')

station.add_component(keithley2400)

#rohde = RohdeSchwarz_SMB100A('rohde', address = 'GPIB0::28::INSTR')
#station.add_component(rohde)


###########now import functions
from experiments.Do_GVg_and_adjust_sitpos import do_GVg_and_adjust_sitpos
from experiments.GVg_qdac_zurich_general import GVG_fun
from experiments.cs_mechanics.cs_mechanics_simple_setpoint_adjust_fun import cs_mechanics_simple_setpoint
from experiments.cs_experiment import CSExperiment
exp=CSExperiment()
