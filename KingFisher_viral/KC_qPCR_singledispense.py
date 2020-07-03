import math
from opentrons.types import Point
from opentrons import protocol_api
import time
import os
import numpy as np
from timeit import default_timer as timer
import json
from datetime import datetime
import csv

# metadata
metadata = {
    'protocolName': 'Station C Kingfisher Pathogen qPCR setup Version 2',
    'author': 'Malen Aguirregabiria, Aitor Gastaminza & José Luis Villanueva (jlvillanueva@clinic.cat)',
    'source': 'Hospital Clínic Barcelona, Hospital Universitario Cruces Bilbao',
    'apiLevel': '2.2',
    'description': 'Protocol for Kingfisher sample setup (C) - Viral Kit (ref ...) using CORE script'

}
'''
'technician': '$technician',
'date': '$date'
'''
#Defined variables
##################
NUM_SAMPLES = 16
NUM_SAMPLES = NUM_SAMPLES + 2 # MMIX for PC and NC

air_gap_vol = 20
air_gap_mmix = 5
air_gap_sample = 2
run_id = '43002'

# Tune variables
size_transfer = 4  # Number of wells the distribute function will fill. Deprecated by calculation function
volume_sample = 5  # Volume of the sample
volume_pc = 20
volume_nc = 20
#volume_mmix_available = 50*20 #(NUM_SAMPLES * 1.5 * volume_mmix)  # Total volume of first screwcap
extra_dispensal = 10  # Extra volume for master mix in each distribute transfer
diameter_screwcap = 8.5  # Diameter of the screwcap
temperature = 8  # Temperature of temp module
volume_cone = 60  # Volume in ul that fit in the screwcap cone
tube_type='screwcap_2ml' #'eppendorf_1.5ml'
x_offset = [0,0]
pipette_allowed_capacity=180
#size_transfer = math.floor(pipette_allowed_capacity / volume_mmix)


#############################################################
#Available master mastermixes
#############################################################
MMIX_available={1: 'Seegene', 2: 'Universal', 3: 'Universal_IDT',4: 'Clinic', 5: 'Multiplex'}

mmix_selection = 5 # select the mastermix to be used

MMIX_vol={1: [17,1], 2: [20,1], 3: [20,1], 4: [40,2], 5:[20,1]} # volume of mastermixes per sample and number of wells in which is distributed
MMIX_recipe={1: [5, 5, 5, 2], 2: [8, 5, 1, 2, 2, 1, 1], 3: [12, 5, 1, 1, 1], 4: [1], 5:[6.25,1.25,12.5]} # Reactive volumes for the mmix

size_transfer = math.floor(pipette_allowed_capacity / MMIX_vol[mmix_selection][0]) # Number of wells the distribute function will fill

MMIX_make_location = 9 # Cell C1 in which the first tube for the MMIX will be placed

volume_mmix = MMIX_vol[mmix_selection][0]  # Volume of transfered master mix per well

MMIX_make = {}
for mmix_type in MMIX_recipe.keys():
    MMIX_make[mmix_type] = []
    for needed_vol in MMIX_recipe[mmix_type]:
        MMIX_make[mmix_type].append(needed_vol * NUM_SAMPLES * 1.1)

volume_mmix_available = ((NUM_SAMPLES + 5) * MMIX_vol[mmix_selection][0])  # Total volume of mastermix that will be prepared

#############################################
# Calculated variables
area_section_screwcap = (np.pi * diameter_screwcap**2) / 4
h_cone = (volume_cone * 3 / area_section_screwcap)
num_cols = math.ceil((NUM_SAMPLES-2) / 8)  # Columns we are working on

def run(ctx: protocol_api.ProtocolContext):
    import os
    #from opentrons.drivers.rpi_drivers import gpio
    #gpio.set_rail_lights(False) #Turn off lights (termosensible reagents)
    ctx.comment('Actual used columns: ' + str(num_cols))

    # Define the STEPS of the protocol
    STEP = 0
    STEPS = {  # Dictionary with STEP activation, description, and times
        1: {'Execute': False, 'description': 'Make MMIX'},
        2: {'Execute': False, 'description': 'Transfer MMIX with P300'},
        3: {'Execute': True, 'description': 'Transfer MMIX with P20'},
        4: {'Execute': True, 'description': 'Clean up NC and PC wells'},
        5: {'Execute': True, 'description': 'Transfer elution'},
        6: {'Execute': False, 'description': 'Transfer PC'},
        7: {'Execute': False, 'description': 'Transfer NC'}
    }

    if STEPS[2]['Execute']==True:
        STEPS[3]['Execute']=False # just to make sure if P300 is being used for the MMIX, do not load P20 single

    for s in STEPS:  # Create an empty wait_time
        if 'wait_time' not in STEPS[s]:
            STEPS[s]['wait_time'] = 0

    #Folder and file_path for log time
    folder_path = '/var/lib/jupyter/notebooks/'+str(run_id)
    if not ctx.is_simulating():
        if not os.path.isdir(folder_path):
            os.mkdir(folder_path)
        file_path = folder_path + '/KC_qPCR_time_log.txt'

    # Define Reagents as objects with their properties
    class Reagent:
        def __init__(self, name, flow_rate_aspirate, flow_rate_dispense, rinse,
                     reagent_reservoir_volume, delay, num_wells, h_cono, v_fondo,
                      tip_recycling = 'none'):
            self.name = name
            self.flow_rate_aspirate = flow_rate_aspirate
            self.flow_rate_dispense = flow_rate_dispense
            self.rinse = bool(rinse)
            self.reagent_reservoir_volume = reagent_reservoir_volume
            self.delay = delay
            self.num_wells = num_wells
            self.col = 0
            self.vol_well = 0
            self.h_cono = h_cono
            self.v_cono = v_fondo
            self.unused=[]
            self.tip_recycling = tip_recycling
            self.vol_well_original = reagent_reservoir_volume / num_wells


    # Reagents and their characteristics
    MMIX = Reagent(name = 'Master Mix',
                      rinse = False,
                      flow_rate_aspirate = 2,
                      flow_rate_dispense = 4,
                      reagent_reservoir_volume = volume_mmix_available,
                      num_wells = 1, #change with num samples
                      delay = 0,
                      h_cono = h_cone,
                      v_fondo = volume_cone  # V cono
                      )

    NC = Reagent(name = 'Negative control',
                      rinse = False,
                      flow_rate_aspirate = 1,
                      flow_rate_dispense = 1,
                      reagent_reservoir_volume = 1000, # volume_mmix_available,
                      num_wells = 1, #change with num samples
                      delay = 0,
                      h_cono = h_cone,
                      v_fondo = volume_cone  # V cono
                      )

    PC = Reagent(name = 'Positive control',
                      rinse = False,
                      flow_rate_aspirate = 1,
                      flow_rate_dispense = 1,
                      reagent_reservoir_volume = 1000, # volume_mmix_available,
                      num_wells = 1, #change with num samples
                      delay = 0,
                      h_cono = h_cone,
                      v_fondo = volume_cone  # V cono
                      )

    Elution = Reagent(name='Elution',
                      rinse=False,
                      flow_rate_aspirate = 1,
                      flow_rate_dispense = 2,
                      reagent_reservoir_volume=50,
                      delay=1,
                      num_wells=num_cols,  # num_cols comes from available columns
                      h_cono=0,
                      v_fondo=0
                      )
    tackpath = Reagent(name = 'MMIX_multiplex_tackpath',
                      rinse = False,
                      flow_rate_aspirate = 1,
                      flow_rate_dispense = 1,
                      reagent_reservoir_volume = 1000,
                      num_wells = 1, #change with num samples
                      delay = 0,
                      h_cono = h_cone,
                      v_fondo = volume_cone  # V cono
                      )

    covid_assay = Reagent(name = 'Covid19_assay',
                      rinse = False,
                      flow_rate_aspirate = 1,
                      flow_rate_dispense = 1,
                      reagent_reservoir_volume = 1000,
                      num_wells = 1, #change with num samples
                      delay = 0,
                      h_cono = h_cone,
                      v_fondo = volume_cone  # V cono
                      )

    mmix_water = Reagent(name = 'nuclease_free_water',
                      rinse = False,
                      flow_rate_aspirate = 1,
                      flow_rate_dispense = 1,
                      reagent_reservoir_volume = 1000,
                      num_wells = 1, #change with num samples
                      delay = 0,
                      h_cono = h_cone,
                      v_fondo = volume_cone  # V cono
                      )

    component4 = Reagent(name = 'component4',
                      rinse = False,
                      flow_rate_aspirate = 1,
                      flow_rate_dispense = 1,
                      reagent_reservoir_volume = 1000,
                      num_wells = 1, #change with num samples
                      delay = 0,
                      h_cono = h_cone,
                      v_fondo = volume_cone  # V cono
                      )

    component5 = Reagent(name = 'component5',
                      rinse = False,
                      flow_rate_aspirate = 1,
                      flow_rate_dispense = 1,
                      reagent_reservoir_volume = 1000,
                      num_wells = 1, #change with num samples
                      delay = 0,
                      h_cono = h_cone,
                      v_fondo = volume_cone  # V cono
                      )

    component6 = Reagent(name = 'component6',
                      rinse = False,
                      flow_rate_aspirate = 1,
                      flow_rate_dispense = 1,
                      reagent_reservoir_volume = 1000,
                      num_wells = 1, #change with num samples
                      delay = 0,
                      h_cono = h_cone,
                      v_fondo = volume_cone  # V cono
                      )

    component7 = Reagent(name = 'component7',
                      rinse = False,
                      flow_rate_aspirate = 1,
                      flow_rate_dispense = 1,
                      reagent_reservoir_volume = 1000,
                      num_wells = 1, #change with num samples
                      delay = 0,
                      h_cono = h_cone,
                      v_fondo = volume_cone  # V cono
                      )

    component8 = Reagent(name = 'component8',
                      rinse = False,
                      flow_rate_aspirate = 1,
                      flow_rate_dispense = 1,
                      reagent_reservoir_volume = 1000,
                      num_wells = 1, #change with num samples
                      delay = 0,
                      h_cono = h_cone,
                      v_fondo = volume_cone  # V cono
                      )

    MMIX.vol_well = MMIX.vol_well_original
    PC.vol_well = PC.vol_well_original
    NC.vol_well = NC.vol_well_original
    Elution.vol_well = Elution.vol_well_original
    tackpath.vol_well = tackpath.vol_well_original
    covid_assay.vol_well = covid_assay.vol_well_original
    mmix_water.vol_well = mmix_water.vol_well_original
    component4.vol_well = component4.vol_well_original
    component5.vol_well = component5.vol_well_original
    component6.vol_well = component6.vol_well_original
    component7.vol_well = component7.vol_well_original
    component8.vol_well = component8.vol_well_original

################## Assign class type reactives to a summary
    MMIX_components=[tackpath, covid_assay, mmix_water, component4, component5, component6, component7, component8]
    MMIX_components=MMIX_components[:len(MMIX_make[mmix_selection])]
##################

    ##################
    # Custom functions
    def divide_volume(volume,max_vol):
        num_transfers=math.ceil(volume/max_vol)
        vol_roundup=math.ceil(volume/num_transfers)
        last_vol = volume - vol_roundup*(num_transfers-1)
        vol_list = [vol_roundup for v in range(1,num_transfers)]
        vol_list.append(last_vol)
        return vol_list


    def divide_destinations(l, n): # Divide the list of destinations in size n lists.
        a=[]
        for i in range(0, len(l), n):
            a.append( l[i:i + n])
        return a

    def distribute_custom(pipette, reagent, volume, src, dest, waste_pool, pickup_height, extra_dispensal, disp_height=0):
        # Custom distribute function that allows for blow_out in different location and adjustement of touch_tip
        air_gap=10
        pipette.aspirate((len(dest) * volume) +
                         extra_dispensal, src.bottom(pickup_height), rate = reagent.flow_rate_aspirate)
        pipette.touch_tip(speed=20, v_offset=-5)
        pipette.move_to(src.top(z=5))
        pipette.aspirate(air_gap, rate = reagent.flow_rate_aspirate)  # air gap

        for d in dest:
            pipette.dispense(air_gap, d.top(), rate = reagent.flow_rate_dispense)
            drop = d.top(z = disp_height)
            pipette.dispense(volume, drop, rate = reagent.flow_rate_dispense)
            ctx.delay(seconds = reagent.delay) # pause for x seconds depending on reagent
            pipette.move_to(d.top(z=5))
            pipette.aspirate(air_gap,d.top(z=5), rate = reagent.flow_rate_aspirate)  # air gap

        try:
            pipette.blow_out(waste_pool.wells()[0].bottom(pickup_height + 3))
        except:
            pipette.blow_out(waste_pool.bottom(pickup_height + 3))
        return (len(dest) * volume)

    def move_vol_multichannel(pipet, reagent, source, dest, vol, air_gap_vol, x_offset,
                       pickup_height, rinse, disp_height, blow_out, touch_tip,
                       post_dispense=False, post_dispense_vol=20,
                       post_airgap=False, post_airgap_vol=10):
        '''
        x_offset: list with two values. x_offset in source and x_offset in destination i.e. [-1,1]
        pickup_height: height from bottom where volume
        rinse: if True it will do 2 rounds of aspirate and dispense before the tranfer
        disp_height: dispense height; by default it's close to the top (z=-2), but in case it is needed it can be lowered
        blow_out, touch_tip: if True they will be done after dispensing
        '''
        # Rinse before aspirating
        if rinse == True:
            custom_mix(pipet, reagent, location = source, vol = vol,
                       rounds = 2, blow_out = True, mix_height = 0,
                       x_offset = x_offset)
        # SOURCE
        s = source.bottom(pickup_height).move(Point(x = x_offset[0]))
        pipet.aspirate(vol, s, rate = reagent.flow_rate_aspirate)  # aspirate liquid
        if air_gap_vol != 0:  # If there is air_gap_vol, switch pipette to slow speed
            pipet.aspirate(air_gap_vol, source.top(z = -2),
                           rate = reagent.flow_rate_aspirate)  # air gap
        # GO TO DESTINATION
        drop = dest.top(z = disp_height).move(Point(x = x_offset[1]))
        pipet.dispense(vol + air_gap_vol, drop,
                       rate = reagent.flow_rate_dispense)  # dispense all
        ctx.delay(seconds = reagent.delay) # pause for x seconds depending on reagent
        if blow_out == True:
            pipet.blow_out(dest.top(z = -2))
        if post_dispense == True:
            pipet.dispense(post_dispense_vol, dest.top(z = -2))
        if touch_tip == True:
            pipet.touch_tip(speed = 20, v_offset = -5, radius = 0.9)
        if post_airgap == True:
            pipet.aspirate(post_airgap_vol, dest.top(z = 5), rate = 2)



    def custom_mix(pipet, reagent, location, vol, rounds, blow_out, mix_height,
    x_offset, source_height = 3, post_airgap=False, post_airgap_vol=10,
    post_dispense=False, post_dispense_vol=20,):
        '''
        Function for mixing a given [vol] in the same [location] a x number of [rounds].
        blow_out: Blow out optional [True,False]
        x_offset = [source, destination]
        source_height: height from bottom to aspirate
        mix_height: height from bottom to dispense
        '''
        if mix_height == 0:
            mix_height = 3
        pipet.aspirate(1, location=location.bottom(
            z=source_height).move(Point(x=x_offset[0])), rate=reagent.flow_rate_aspirate)
        for _ in range(rounds):
            pipet.aspirate(vol, location=location.bottom(
                z=source_height).move(Point(x=x_offset[0])), rate=reagent.flow_rate_aspirate)
            pipet.dispense(vol, location=location.bottom(
                z=mix_height).move(Point(x=x_offset[1])), rate=reagent.flow_rate_dispense)
        pipet.dispense(1, location=location.bottom(
            z=mix_height).move(Point(x=x_offset[1])), rate=reagent.flow_rate_dispense)
        if blow_out == True:
            pipet.blow_out(location.top(z=-2))  # Blow out
        if post_dispense == True:
            pipet.dispense(post_dispense_vol, location.top(z = -2))
        if post_airgap == True:
            pipet.dispense(post_airgap_vol, location.top(z = 5))

    def calc_height(reagent, cross_section_area, aspirate_volume, min_height = 0.3, extra_volume = 30):
        nonlocal ctx
        ctx.comment('Remaining volume ' + str(reagent.vol_well) +
                    '< needed volume ' + str(aspirate_volume) + '?')
        if reagent.vol_well < aspirate_volume + extra_volume:
            reagent.unused.append(reagent.vol_well)
            ctx.comment('Next column should be picked')
            ctx.comment('Previous to change: ' + str(reagent.col))
            # column selector position; intialize to required number
            reagent.col = reagent.col + 1
            ctx.comment(str('After change: ' + str(reagent.col)))
            reagent.vol_well = reagent.vol_well_original
            ctx.comment('New volume:' + str(reagent.vol_well))
            height = (reagent.vol_well - aspirate_volume - reagent.v_cono) / cross_section_area
                    #- reagent.h_cono
            reagent.vol_well = reagent.vol_well - aspirate_volume
            ctx.comment('Remaining volume:' + str(reagent.vol_well))
            if height < min_height:
                height = min_height
            col_change = True
        else:
            height = (reagent.vol_well - aspirate_volume - reagent.v_cono) / cross_section_area #- reagent.h_cono
            reagent.vol_well = reagent.vol_well - aspirate_volume
            ctx.comment('Calculated height is ' + str(height))
            if height < min_height:
                height = min_height
            ctx.comment('Used height is ' + str(height))
            col_change = False
        return height, col_change

    ####################################
    # load labware and modules
    # 24 well rack
    tuberack = ctx.load_labware(
        'opentrons_24_aluminumblock_generic_2ml_screwcap', '1',
        'Bloque Aluminio opentrons 24 screwcaps 2000 µL ')

    trash = ctx.load_labware('agilent_1_reservoir_290ml','12', 'trash').wells()[0]

    ############################################
    # tempdeck
    tempdeck = ctx.load_module('tempdeck', '4')
    tempdeck.set_temperature(temperature)

    ##################################
    # qPCR plate - final plate, goes to PCR
    qpcr_plate = tempdeck.load_labware(
        'abi_fast_qpcr_96_alum_opentrons_100ul',
        'chilled qPCR final plate')

    ##################################
    # Sample plate - comes from B
    source_plate = ctx.load_labware(
        "kingfisher_std_96_wellplate_550ul", '2',
        'KF plate with elutions')
    samples = source_plate.wells()[:NUM_SAMPLES]

    ##################################
    # Load Tipracks
    tips20 = [
        ctx.load_labware('opentrons_96_filtertiprack_20ul', slot)
        for slot in ['5']
    ]

    tips200 = [
        ctx.load_labware('opentrons_96_filtertiprack_200ul', slot)
        for slot in ['6']
    ]

    # pipettes
    m20 = ctx.load_instrument(
        'p20_multi_gen2', mount='left', tip_racks=tips20)

    if STEPS[2]['Execute']==True:
        p300 = ctx.load_instrument(
        'p300_single_gen2', mount='right', tip_racks=tips200)
        # used tips counter
        tip_track = {
            'counts': {p300: 0, m20: 0},
            'maxes': {p300: len(tips200) * 96, m20: len(tips20)*96}
        }
    else:
        tips20mmix = [ ctx.load_labware('opentrons_96_filtertiprack_20ul', slot)
                        for slot in ['8'] ]
        p20 = ctx.load_instrument(
        'p20_single_gen2', mount='right', tip_racks=tips20mmix)
        # used tips counter
        tip_track = {
            'counts': {p20: 0, m20: 0},
            'maxes': {p20: len(tips200) * 96, m20: len(tips20)*96}
        }



    ################################################################################
    # Declare which reagents are in each reservoir as well as deepwell and elution plate
    MMIX.reagent_reservoir = tuberack.rows()[0][:MMIX.num_wells] # 1 row, 2 columns (first ones)
    MMIX_components_location=tuberack.wells()[MMIX_make_location:(MMIX_make_location + len(MMIX_make[mmix_selection]))]
    ctx.comment('Wells in: '+ str(tuberack.rows()[0][:MMIX.num_wells]) + ' element: '+str(MMIX.reagent_reservoir[MMIX.col]))
    PC.reagent_reservoir = tuberack.rows()[1][:1]
    NC.reagent_reservoir = tuberack.rows()[2][:1]
    # setup up sample sources and destinations
    samples = source_plate.wells()[:NUM_SAMPLES]
    samples_multi = source_plate.rows()[0][:num_cols]
    pcr_wells = qpcr_plate.wells()[:(NUM_SAMPLES)]
    pcr_wells_multi = qpcr_plate.rows()[0][:num_cols]
    pc_well = qpcr_plate.wells()[(NUM_SAMPLES-2):(NUM_SAMPLES-1)][0]
    nc_well = qpcr_plate.wells()[(NUM_SAMPLES-1):(NUM_SAMPLES)][0]
    pc_well_clean = source_plate.wells()[(NUM_SAMPLES-2):(NUM_SAMPLES-1)][0]
    nc_well_clean = source_plate.wells()[(NUM_SAMPLES-1):(NUM_SAMPLES)][0]
    # Divide destination wells in small groups for P300 pipette
    dests = list(divide_destinations(pcr_wells, size_transfer))

    ##########
    # pick up tip and if there is none left, prompt user for a new rack
    def pick_up(pip):
        nonlocal tip_track
        if not ctx.is_simulating():
            if tip_track['counts'][pip] == tip_track['maxes'][pip]:
                ctx.pause('Replace ' + str(pip.max_volume) + 'µl tipracks before \
                resuming.')
                pip.reset_tipracks()
                tip_track['counts'][pip] = 0

        if not pip.hw_pipette['has_tip']:
            pip.pick_up_tip()
    ##########
    ############################################################################
    # STEP 1: Make Master MIX
    ############################################################################
    ctx._hw_manager.hardware.set_lights(rails=False) # set lights off when using MMIX
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()
        # Check if among the pipettes, p300_single is installed
        used_vol=[]
        ctx.comment('#######################################################')
        ctx.comment('Selected MMIX: '+MMIX_available[mmix_selection])
        ctx.comment('#######################################################')

        for i,[source, vol] in enumerate(zip(MMIX_components_location, MMIX_make[mmix_selection])):
            pick_up(p300)
            ctx.comment('#######################################################')
            ctx.comment('Add component: '+MMIX_components[i].name)
            ctx.comment('#######################################################')
            if (vol + air_gap_vol) > pipette_allowed_capacity: # because 200ul is the maximum volume of the tip we will choose 180
            # calculate what volume should be transferred in each step
                vol_list=divide_volume(vol, pipette_allowed_capacity)
                for vol in vol_list:
                    move_vol_multichannel(p300, reagent=MMIX_components[i], source=source, dest=MMIX.reagent_reservoir[0],
                    vol=vol, air_gap_vol=air_gap_vol, x_offset = x_offset,pickup_height=1, # should be changed with picku_up_height calculation
                    rinse=False, disp_height=-10,blow_out=True, touch_tip=True)
            else:
                move_vol_multichannel(p300, reagent=MMIX_components[i], source=source, dest=MMIX.reagent_reservoir[0],
                vol=vol, air_gap_vol=air_gap_vol, x_offset=x_offset, pickup_height=1, # should be changed with picku_up_height calculation
                rinse=False, disp_height=-10,blow_out=True, touch_tip=True)

            if i+1<len(MMIX_components):
                p300.drop_tip()
            else:
                ctx.comment('#######################################################')
                ctx.comment('Final mix')
                ctx.comment('#######################################################')
                custom_mix(p300, reagent = MMIX, location = MMIX.reagent_reservoir[0], vol = 180, rounds = 5,
                blow_out = True, mix_height = 2, x_offset = x_offset)
                p300.drop_tip()

            tip_track['counts'][p300]+=1

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('#######################################################')
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        ctx.comment('#######################################################')
        STEPS[STEP]['Time:'] = str(time_taken)


    ############################################################################
    # STEP 2: Transfer Master MIX with P300
    ############################################################################
    ctx._hw_manager.hardware.set_lights(rails=False) # set lights off when using MMIX
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()
        p300.pick_up_tip()

        for dest in pcr_wells:
            [pickup_height, col_change] = calc_height(MMIX, area_section_screwcap, volume_mmix)

            move_vol_multichannel(p300, reagent = MMIX, source = MMIX.reagent_reservoir[MMIX.col],
            dest = dest, vol = volume_mmix, air_gap_vol = air_gap_mmix, x_offset = x_offset,
                   pickup_height = pickup_height, disp_height = -10, rinse = False,
                   blow_out=True, touch_tip=False)

            #used_vol.append(used_vol_temp)

        p300.drop_tip()
        tip_track['counts'][p300]+=1
        #MMIX.unused_two = MMIX.vol_well

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('#######################################################')
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        ctx.comment('#######################################################')
        STEPS[STEP]['Time:'] = str(time_taken)

    ############################################################################
    # STEP 3: Transfer Master MIX with P20
    ############################################################################
    ctx._hw_manager.hardware.set_lights(rails=False) # set lights off when using MMIX
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()
        p20.pick_up_tip()

        for dest in pcr_wells:
            [pickup_height, col_change] = calc_height(MMIX, area_section_screwcap, volume_mmix)

            move_vol_multichannel(p20, reagent = MMIX, source = MMIX.reagent_reservoir[MMIX.col],
            dest = dest, vol = volume_mmix, air_gap_vol = 0, x_offset = x_offset,
                   pickup_height = pickup_height, disp_height = -10, rinse = False,
                   blow_out=True, touch_tip=True)

            #used_vol.append(used_vol_temp)

        p20.drop_tip()
        tip_track['counts'][p20]+=1
        #MMIX.unused_two = MMIX.vol_well

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('#######################################################')
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        ctx.comment('#######################################################')
        STEPS[STEP]['Time:'] = str(time_taken)
        ctx.pause('Put samples please')

    ############################################################################
    # STEP 4: Clean up PC and NC well
    ############################################################################
    ctx._hw_manager.hardware.set_lights(rails=False) # set lights off when using MMIX
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()
        clean_up_wells=[pc_well_clean,nc_well_clean]
        p20.pick_up_tip()
        for src in clean_up_wells:
            for vol in range(3):
                p20.aspirate(20,src)
                p20.dispense(20,trash)
        p20.drop_tip()
        tip_track['counts'][p20]+=1

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('#######################################################')
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        ctx.comment('#######################################################')
        STEPS[STEP]['Time:'] = str(time_taken)

    ############################################################################
    # STEP 5: TRANSFER Samples
    ############################################################################
    ctx._hw_manager.hardware.set_lights(rails=False) # set lights off when using MMIX
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()
        ctx.comment('pcr_wells')
        #Loop over defined wells
        for s, d in zip(samples_multi, pcr_wells_multi):
            m20.pick_up_tip()
            #Source samples
            move_vol_multichannel(m20, reagent = Elution, source = s, dest = d,
            vol = volume_sample, air_gap_vol = air_gap_sample, x_offset = x_offset,
                   pickup_height = 0.3, disp_height = -10, rinse = False,
                   blow_out=True, touch_tip=True, post_airgap=False)
            ## ADD Custom mix
            m20.drop_tip()
            tip_track['counts'][m20]+=8

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('#######################################################')
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        ctx.comment('#######################################################')
        STEPS[STEP]['Time:'] = str(time_taken)
        tempdeck.deactivate()

    ############################################################################
    # STEP 6: Transfer PC with P20
    ############################################################################
    ctx._hw_manager.hardware.set_lights(rails=False) # set lights off when using MMIX
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()
        p20.pick_up_tip()

        [pickup_height, col_change] = calc_height(PC, area_section_screwcap, volume_mmix)

        move_vol_multichannel(p20, reagent = PC, source = PC.reagent_reservoir[PC.col],
        dest = pc_well, vol = volume_pc, air_gap_vol = 0, x_offset = x_offset,
               pickup_height = pickup_height, disp_height = -10, rinse = False,
               blow_out=True, touch_tip=True)

        p20.drop_tip()
        tip_track['counts'][p20]+=1
        #MMIX.unused_two = MMIX.vol_well

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('#######################################################')
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        ctx.comment('#######################################################')
        STEPS[STEP]['Time:'] = str(time_taken)

    ############################################################################
    # STEP 7: Transfer NC with P20
    ############################################################################
    ctx._hw_manager.hardware.set_lights(rails=False) # set lights off when using MMIX
    STEP += 1
    if STEPS[STEP]['Execute'] == True:
        start = datetime.now()
        p20.pick_up_tip()

        [pickup_height, col_change] = calc_height(NC, area_section_screwcap, volume_mmix)

        move_vol_multichannel(p20, reagent = NC, source = NC.reagent_reservoir[NC.col],
        dest = nc_well, vol = volume_nc, air_gap_vol = 0, x_offset = x_offset,
               pickup_height = pickup_height, disp_height = -10, rinse = False,
               blow_out=True, touch_tip=True)

        p20.drop_tip()
        tip_track['counts'][p20]+=1
        #MMIX.unused_two = MMIX.vol_well

        end = datetime.now()
        time_taken = (end - start)
        ctx.comment('#######################################################')
        ctx.comment('Step ' + str(STEP) + ': ' +
                    STEPS[STEP]['description'] + ' took ' + str(time_taken))
        ctx.comment('#######################################################')
        STEPS[STEP]['Time:'] = str(time_taken)

    ############################################################################
    # Export the time log to a tsv file
    if not ctx.is_simulating():
        with open(file_path, 'w') as f:
            f.write('STEP\texecution\tdescription\twait_time\texecution_time\n')
            for key in STEPS.keys():
                row = str(key)
                for key2 in STEPS[key].keys():
                    row += '\t' + format(STEPS[key][key2])
                f.write(row + '\n')
        f.close()

    ############################################################################
    # Light flash end of program

    time.sleep(2)
    import os
    #os.system('mpg123 -f -8000 /etc/audio/speaker-test.mp3 &')

    '''if STEPS[1]['Execute'] == True:
        total_used_vol = np.sum(used_vol)
        total_needed_volume = total_used_vol
        ctx.comment('Total Master Mix used volume is: ' + str(total_used_vol) + '\u03BCl.')
        ctx.comment('Needed Master Mix volume is ' +
                    str(total_needed_volume + extra_dispensal*len(dests)) +'\u03BCl')
        ctx.comment('Used Master Mix volumes per run are: ' + str(used_vol) + '\u03BCl.')
        ctx.comment('Master Mix Volume remaining in tubes is: ' +
                    format(np.sum(MMIX.unused)+extra_dispensal*len(dests)+MMIX.vol_well) + '\u03BCl.')
        ctx.comment('200 ul Used tips in total: ' + str(tip_track['counts'][p300]))
        ctx.comment('200 ul Used racks in total: ' + str(tip_track['counts'][p300] / 96))'''

    if (STEPS[3]['Execute'] == True & STEPS[4]['Execute'] == True):
        ctx.comment('20 ul Used tips in total: ' + str(tip_track['counts'][m20]+tip_track['counts'][p20]))
        ctx.comment('20 ul Used racks in total: ' + str((tip_track['counts'][m20]+tip_track['counts'][p20]) / 96))

    if (STEPS[2]['Execute'] == True & STEPS[4]['Execute'] == True):
        ctx.comment('200 ul Used tips in total: ' + str(tip_track['counts'][p300]))
        ctx.comment('200 ul Used racks in total: ' + str(tip_track['counts'][p300] / 96))
        ctx.comment('20 ul Used tips in total: ' + str(tip_track['counts'][m20]))
        ctx.comment('20 ul Used racks in total: ' + str((tip_track['counts'][m20] / 96)))

    if (STEPS[2]['Execute'] == True & STEPS[4]['Execute'] == False):
        ctx.comment('200 ul Used tips in total: ' + str(tip_track['counts'][p300]))
        ctx.comment('200 ul Used racks in total: ' + str(tip_track['counts'][p300] / 96))

    if (STEPS[3]['Execute'] == False & STEPS[4]['Execute'] == True):
        ctx.comment('20 ul Used tips in total: ' + str(tip_track['counts'][m20]))
        ctx.comment('20 ul Used racks in total: ' + str((tip_track['counts'][m20] / 96)))


    for i in range(3):
        ctx._hw_manager.hardware.set_lights(rails=False)
        #ctx._hw_manager.hardware.set_lights(button=(1,0,0))
        time.sleep(0.3)
        ctx._hw_manager.hardware.set_lights(rails=True)
        #ctx._hw_manager.hardware.set_button_light(0,0,1)
        time.sleep(0.3)
        ctx._hw_manager.hardware.set_lights(rails=False)
    #ctx._hw_manager.hardware.set_button_light(0,1,0)

    ctx.comment('Finished! \nMove plate to PCR')
