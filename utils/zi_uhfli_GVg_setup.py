from instruments import zurich


def zi_uhfli_GVg_setup(source_power, freq,tc):
    '''
    source_power : in Vpk
    freq : frequency of the source and measured frequency at. 
    '''
    # setup the zurich -----------------------------------
    # turn off everything
    zurich.mods.enable(0)
    for i in range(8):
        param1 = getattr(zurich.sigouts.sigouts0.enables, f'enables{i}')
        param2 = getattr(zurich.sigouts.sigouts1.enables, f'enables{i}')
        param1.value(0)
        param2.value(0)
       

    # set the amplitude of 0,0 and 1,1
    zurich.sigouts.sigouts0.amplitudes.amplitudes0.value(source_power)

    # turn on 0,0 and 1,1
    zurich.sigouts.sigouts0.enables.enables0.value(1)

    # turn on outputs
    zurich.sigouts.sigouts0.on(1)


   # turn on 0,0 and 1,1
    zurich.sigouts.sigouts1.enables.enables0.value(3)

       # turn on outputs
    zurich.sigouts.sigouts1.on(3)

    # set demod 0
    zurich.demods.demods3.oscselect(0)
    zurich.demods.demods3.adcselect(0)
    zurich.demods.demods3.timeconstant(tc)
    zurich.demods.demods3.enable(1)
    zurich.oscs.oscs0.freq(freq)

    # set demod 1, 2 (the source and gate drives)
    # zurich.demods.demods1.oscselect(1)
    # zurich.demods.demods2.oscselect(2)

    zurich.mods.mods1.enable(0)
    zurich.mods.mods1.carrier.oscselect(0)
    zurich.mods.mods1.sidebands.sidebands0.mode(0)
    zurich.mods.mods1.sidebands.sidebands0.enable(0)
    zurich.mods.mods1.sidebands.sidebands0.oscselect(0)
    zurich.mods.mods1.sidebands.sidebands1.mode(0)


