from instruments import zurich


def ZI_UFHLI_two_source1D(source_power, gate_power, beat_note_freq):
    '''
    source_power : in Vpk
    gate_power   : in Vpk
    beat_note_freq : The beatnote frequency of the Source and Gate that the ZI will demodulate at. 
    '''
    # setup the zurich -----------------------------------
    # turn off everything
    # zurich.mods.enable(0)
    for i in range(8):
        param1 = getattr(zurich.sigouts.sigouts0.enables, f'enables{i}')
        param2 = getattr(zurich.sigouts.sigouts1.enables, f'enables{i}')
        param1.value(0)
        param2.value(0)
       

    # set the amplitude of 0,0 and 1,1
    zurich.sigouts.sigouts0.amplitudes.amplitudes4.value(source_power)
    zurich.sigouts.sigouts1.amplitudes.amplitudes5.value(gate_power)

    # turn on 0,0 and 1,1
    zurich.sigouts.sigouts0.enables.enables4.value(1)
    zurich.sigouts.sigouts1.enables.enables5.value(1)

    # turn on outputs
    zurich.sigouts.sigouts0.on(1)
    zurich.sigouts.sigouts1.on(1)

    # set demod 0
    zurich.demods.demods3.oscselect(2)
    zurich.demods.demods3.adcselect(0)
    zurich.demods.demods3.enable(1)
    zurich.oscs.oscs0.freq(beat_note_freq)

    # set demod 1, 2 (the source and gate drives)
    # zurich.demods.demods1.oscselect(1)
    # zurich.demods.demods2.oscselect(2)

    zurich.mods.mods0.enable(1)
    # zurich.mods.mods0.carrier.oscselect(0)
    # zurich.mods.mods1.sidebands.sidebands0.mode(1)
    # zurich.mods.mods1.sidebands.sidebands0.enable(0)
    # zurich.mods.mods1.sidebands.sidebands0.oscselect(0)
    # zurich.mods.mods1.sidebands.sidebands1.mode(0)


