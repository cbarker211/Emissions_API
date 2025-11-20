def update_mass_info(temp_dict,name,variant):
    ''' For details on rocket approximations, see Launch and Rocket Checks.xlsx, sheet: Missing Propellant Info.
        All rockets with the same number of stages are compared using prop type, LEO payload and rocket length.
        
        All rocket mass information (propellant, stage, and fairing mass) is initially compiled using DISCOSweb.
        However, DISCOSweb can contain incomplete, out of date, or incorrect information.
        This file updates the propellant, stage, and fairing mass for each rocket from 2020-2022. 
        Priority for information is as follows:
            -   Rocket manuals / Mission guides 
            -   Space Launch Report / Spaceflight101
            -   DISCOSweb / DISCOSweb (wet-dry)
            -   Norbert Brugge / Jonathan's Space Report GCAT (dry mass only).
    '''
    
    # We could rearrange all of this to get the rocket information directly from GCAT.
    #   Pros: GCAT labels stage number and we could quickly and easily query mass.
    #   Cons: GCAT also has missing information, and would need to be cross-referenced with other sources anyway.
     
    if name in ["Angara A5","Angara A5 Persei","Angara A5 Orion"]:
        # https://web.archive.org/web/20220406013831/http://www.spacelaunchreport.com/angara.html
        # https://discosweb.esoc.esa.int/launch-vehicles/107449
        # https://discosweb.esoc.esa.int/launch-vehicles/130
        # https://spaceflight101.com/spacerockets/angara-a5/
        # https://www.mach5lowdown.com/wp-content/uploads/PUG/Angara-Mission-Planners-Guide-Rev-0-2002-12.pdf
        
        # Prop masses:
        # DISCOSweb: 530400, 132600, 36000, 19800 (Briz-M), 15000 (Persei)
        # SLR      : 530400, 132600, 35800, 19800 (Briz-M), 18700 (Persei) 
        # Sp101    : 512000, 128000, 35800, 19800 (Briz-M), N/A   (Persei)
        # Manual   : 530400, 132600, 35800, 19800 (Briz-M), N/A   (Persei)
        
        # Stage masses:
        # DISCOSweb: 39200, 9800, 2000, 2370 (Briz-M), 2140 (Persei), N/A  (Fairing Briz-M), N/A (Fairing Persei)
        # SLR      : N/A,   N/A,  N/A,  2370 (Briz-M), 2900 (Persei), N/A  (Fairing Briz-M), N/A (Fairing Persei)
        # Sp101    : 39200, 9800, 4000, 2370 (Briz-M), N/A  (Persei), 1600 (Fairing Briz-M), N/A (Fairing Persei)
        # Manual   : N/A,   N/A,  N/A,  2370 (Briz-M), N/A  (Persei), 2600?? (Fairing Briz-M), N/A (Fairing Persei)
        
        temp_dict[f"Booster Number"] = 4
        temp_dict[f"Stage2 Propellant Mass"]  = 35800 # Manual
        temp_dict[f"Stage2 Stage Mass"]       = 4000  # Sp101

        # Not much information on Orion, but it seems similar to Persei, based on 11S861-03 Phase I Version 2 upgrades to Blok DM-03.
        if name in ["Angara A5 Persei","Angara A5 Orion"]:
            temp_dict[f"Stage3 Propellant Mass"]  = 18700 # SLR
            temp_dict[f"Stage3 Stage Mass"]       = 2900  # SLR
            temp_dict[f"Fairing Mass"] = 1600 # No info, approximating using Briz-M.
        elif name == "Angara A5":
            temp_dict[f"Fairing Mass"] = 1600 # Sp101
            
    if name == "Angara 1.2":
        # Stage 1 is same as Angara A5, but no boosters.
        # Different stage 2.
        # Manual: https://www.mach5lowdown.com/wp-content/uploads/PUG/Angara-Mission-Planners-Guide-Rev-0-2002-12.pdf
        temp_dict[f"Stage1 Propellant Mass"]  = 128800 # Manual
        temp_dict[f"Stage1 Stage Mass"]       = (9800+10500) / 2  # SLR/Sp101
        temp_dict[f"Stage2 Propellant Mass"]  = 25700  # Manual
        temp_dict[f"Stage2 Stage Mass"]       = 2355   # Sp101
        # https://www.russianspaceweb.com/angara1.html There is an 'Aggregate Module' upper stage according to this.
        # Setting using DW info.
        temp_dict[f"Stage3 Propellant Mass"]  = 2000
        temp_dict[f"Stage3 Stage Mass"]       = 1000
        temp_dict[f"Stage3 Propellant Name"]  = 'UDMH (Unsymmetrical Dimethyl Hydrazine)/N2O4'
        temp_dict[f"Stage3 Fuel Type"]        = 'Hypergolic'
        temp_dict[f"Fairing Mass"] = (500+810)/2 # Two options on Sp101.
        
    if name == "Ariane 5ECA":
        
        # Manaul is for ESC-A upper stage, whereas launches since 2019 have used ESC-D, and improved version.
        # However due to lack of information, we will use ESC-A details.
        # https://www.arianespace.com/wp-content/uploads/2016/10/Ariane5-users-manual-Jun2020.pdf
        # Prop masses:
        # DISCOSweb: 480450 (two), 173440, 14540
        # Manual   : 240000 (one), 170000, 14900   
             
        # Stage masses:
        # DISCOSweb: 76400, 14700, 5000
        # Manual   : N/A,   14700, 4540
        
        temp_dict[f"Booster Number"]            = 2 # Manual
        temp_dict[f"Stage0 Propellant Mass"]    = 240000 * int(temp_dict[f"Booster Number"]) # Manual
        temp_dict[f"Stage0 Stage Mass"]         = 31000     # SLR/Sp101
        temp_dict[f"Stage1 Propellant Mass"]    = 170000    # Manual 
        temp_dict[f"Stage1 Stage Mass"]         = 14700     # Manual
        temp_dict[f"Stage2 Propellant Mass"]    = 14900     # Manual
        temp_dict[f"Stage2 Stage Mass"]         = 4540      # Manual
        temp_dict[f"Fairing Mass"]              = 2675      # Manual        
    
    if name[:7] == "Atlas V":
        # Booster propellant and stage mass changed on November 13th 2020 when updated from AJ-60A to GEM63.  
        # https://web.archive.org/web/20180918143456/http://www.northropgrumman.com/Capabilities/GEM/Documents/GEM_63_GEM_63XL.pdf  
        
        # The Centaur upper stage has used the RL-10C-1 engine since 2014.
        # Two launches have used the RL-10C-1-1A engine since (https://en.wikipedia.org/wiki/RL10 and NB), on 18/05/2021 and 04/08/2022.
        # https://web.archive.org/web/20190629113839/https://www.rocket.com/sites/default/files/documents/Capabilities/PDFs/RL10_data_sheet.pdf
        # https://ntrs.nasa.gov/api/citations/20110015783/downloads/20110015783.pdf
        # Dry mass reduced by 2kg.
        
        # Fairing mass:
        #   - https://www.ulalaunch.com/docs/default-source/rockets/atlasvusersguide2010.pdf            2305 kg (EPF) / 3524 kg (5m short)
        #   - https://spaceflight101.com/spacerockets/atlas-v-531/                                      2305 kg (EPF) / 3524 kg (5m short) 
        #   - https://web.archive.org/web/20220406013818/http://www.spacelaunchreport.com/atlas5.html   2260 kg (EPF)                                    
                
        temp_dict[f"Booster Number"] = name[9]
        if "G" in variant: 
            temp_dict[f"Stage0 Propellant Mass"]   = 44200 * int(temp_dict[f"Booster Number"]) # Manual for GEM63
            temp_dict[f"Stage0 Stage Mass"]        = 5100 * int(temp_dict[f"Booster Number"])  # Manual for GEM63
        else:
            temp_dict[f"Stage0 Propellant Mass"]   = 42630 * int(temp_dict[f"Booster Number"]) # SLR for AJ-60A
            temp_dict[f"Stage0 Stage Mass"]        = 3630 * int(temp_dict[f"Booster Number"])  # SLR for AJ-60A
        
        # The original RL-10A-4-2 engine has the same masses as the RL-10C-1 engine, so we can apply this to everything.
        temp_dict[f"Stage2 Propellant Mass"]    = 20800 # SLR for RL-10C-1

        # But the newer RL-10C-1-1A engine has a slightly lower dry mass.
        if "1" in variant: 
            temp_dict[f"Stage2 Stage Mass"] = 2028  # Manual for RL-10C-1-1A
        else:
            temp_dict[f"Stage2 Stage Mass"] = 2030  # SLR for RL-10C-1
        
        # The 400 series.
        if name[8] == "4":    
            temp_dict[f"Stage1 Stage Mass"] = 21054  # https://www.ulalaunch.com/docs/default-source/rockets/atlasvusersguide2010a.pdf?sfvrsn=f84bb59e_2
            temp_dict[f"Fairing Mass"]      = 2305   # https://www.ulalaunch.com/docs/default-source/rockets/atlasvusersguide2010a.pdf?sfvrsn=f84bb59e_2
        # The 500 series.
        elif name[8] == "5":  
            temp_dict[f"Stage1 Stage Mass"] = 21351  # https://www.ulalaunch.com/docs/default-source/rockets/atlasvusersguide2010a.pdf?sfvrsn=f84bb59e_2  
            temp_dict[f"Fairing Mass"]      = 3524   # https://www.ulalaunch.com/docs/default-source/rockets/atlasvusersguide2010a.pdf?sfvrsn=f84bb59e_2 
        # The N Series, no fairing and uses DEC not SEC. No info so using DISCOSweb for stage masses.
        else:
            temp_dict[f"Stage2 Propellant Mass"] = 20800 # DW wet-dry, same as SEC.
            
    if name in ["Antares 230","Antares 230+"]:
        
        # Prop masses:
        # DISCOSweb: 242000, 12837
        # SLR      : 242400, 23100 
        # Sp101    : 242000, N/A
        
        # Stage masses:
        # DISCOSweb: 20600, 1185
        # SLR      : 18800, 2100
        # Sp101    : 20600, 2100
        
        # Antares 230, second stage solid propellant has a typo in DISCOSweb.
        # Looking at user manuals, this value is also wrong, it seems to be for the Castor 30B second stage.
        # This vehicle uses the Castor 30XL, which has double the propellant mass.
        # https://www.northropgrumman.com/wp-content/uploads/CASTOR-Motor-Series.pdf Had to convert LBM to kg.
        temp_dict[f"Stage1 Propellant Mass"] = (242400+242000) / 2
        temp_dict[f"Stage1 Stage Mass"]      = (18800+20600)   / 2
        temp_dict[f"Stage2 Propellant Mass"] = 24924 # CASTOR documentation
        temp_dict[f"Stage2 Stage Mass"]      = 1392  # CASTOR documentation
        # Fairing mass:
        #   - https://web.archive.org/web/20240614065100/https://spaceflight101.com/spacerockets/antares-200-series/ 970kg
        #   - https://web.archive.org/web/20220406013828/http://www.spacelaunchreport.com/taurus2.html               972kg
        temp_dict[f"Fairing Mass"] = 971
    
    if name == "Astra Rocket 3":
        # Astra Rocket (USA) - 3 failed launches in 2020.
        # 2-stage rocket (Kerosene/LOX both stages)
        # First stage powered by 5 Delphin engines, second stage powered by 1 Aether engine.
        # No other 2-stage rockets with similar LEO payload or length.
        # 25 kg LEO payload capacity according to SLR. Looking for other rockets with <1000 kg LEO payload and similar size.
        # Electron, Pegasus XL, Qased, Simorgh all <1000kg. However, only Electron (300kg) has the same propellant.
        # Electron (17m) has three stages, but also has a similar height to Astra (11.6m).
        temp_dict[f"Stage1 Propellant Mass"] = 9250             # Using value for Electron from SLR/DISCOSweb(wet/dry) 
        temp_dict[f"Stage2 Propellant Mass"] = (2050+2150) / 2  # Using value for Electron from SLR/Sp101
        temp_dict[f"Stage1 Stage Mass"]      = 950              # Using value for Electron
        temp_dict[f"Stage2 Stage Mass"]      = 250              # Using value for Electron
        temp_dict[f"Fairing Mass"]           = 44               # Using value for Electron

    if name == "Ceres-1":
        # Ceres-1 
        # Will assume hypergolic, as wikipedia and following link suggest hydrazine.
        # https://www.nasaspaceflight.com/2021/12/chinas-galactic-energy-launches-second-ceres-1-rocket-successfully/
        # https://www.galactic-energy.cn/index.php/En/List/cid/14
        # Website lists fourth stage as 'advanced liquid upper stage', will assume hypergolic.
        # See excel file for comparison details. Best represented by Shavit.
        # Was originally represented by Kuaizhou-1, but comparing based on LEO payload and rocket length has shown Shavit is more similar.

        # For fairing mass, approximating using Kuaizhou-1 (which approximates using Shavit). All similar LEO payload mass and size.
        
        temp_dict[f"Stage1 Propellant Mass"] = 12750 # SLR Shavit
        temp_dict[f"Stage2 Propellant Mass"] = 12750 # SLR Shavit
        temp_dict[f"Stage3 Propellant Mass"] = 1890  # SLR Shavit
        temp_dict[f"Stage4 Propellant Mass"] = 166   # SLR Shavit
        temp_dict[f"Stage1 Stage Mass"]      = 1240  # SLR Shavit
        temp_dict[f"Stage2 Stage Mass"]      = 1376  # SLR Shavit
        temp_dict[f"Stage3 Stage Mass"]      = 684   # SLR Shavit
        temp_dict[f"Stage4 Stage Mass"]      = 71    # SLR Shavit
        temp_dict[f"Fairing Mass"]           = 57    # Shavit http://www.b14643.de/Spacerockets_1/Rest_World/Shavit/Description/Frame.htm
        
    if name == "Delta 4H":
        
        # https://web.archive.org/web/20240422092922/https://spaceflight101.com/spacerockets/delta-iv-heavy-rs-68a/
        # Prop masses:
        # DISCOSweb: 399300 (two), 199650, 27200
        # SLR      : 204000 (one), 204000, 27200
        # Sp101    : 202000 (one), 202000, 27220
        # Manual   : N/A,          N/A,    27200
                
        # Stage masses:
        # DISCOSweb: 28000 (???), 3510,  3510
        # SLR      : 28000 (one), 28000, 3350  
        # Sp101    : 26400 (one), 26400, 3490

        # Fairing mass:
        #   - https://web.archive.org/web/20220406013823/http://www.spacelaunchreport.com/delta4.html   3550
        temp_dict[f"Booster Number"] = 2
        temp_dict[f"Stage0 Propellant Mass"]  = 203000 * int(temp_dict[f"Booster Number"]) # SLR/Sp101
        temp_dict[f"Stage1 Propellant Mass"]  = 203000              # SLR/Sp101  
        temp_dict[f"Stage0 Stage Mass"]       = (28000+26400) / 2 * int(temp_dict[f"Booster Number"]) # SLR/Sp101
        temp_dict[f"Stage1 Stage Mass"]       = (28000+26400) / 2   # SLR/Sp101
        temp_dict[f"Stage2 Stage Mass"]       = (3350+3490) / 2     # SLR/Sp101
        temp_dict[f"Fairing Mass"]            = 3550                # SLR

    if name.startswith("Diamant"):
        # https://www.russianspaceweb.com/diamant.html
        if name == "Diamant A":
            temp_dict[f"Stage1 Propellant Mass"]  = 12800
            temp_dict[f"Stage2 Propellant Mass"]  = 2230
            temp_dict[f"Stage3 Propellant Mass"]  = 640
            temp_dict[f"Stage1 Stage Mass"]       = 1900
            temp_dict[f"Stage2 Stage Mass"]       = 670
            temp_dict[f"Stage3 Stage Mass"]       = 70
        elif name == "Diamant B":
            temp_dict[f"Stage1 Propellant Mass"]  = 18031
            temp_dict[f"Stage2 Propellant Mass"]  = 2230
            temp_dict[f"Stage3 Propellant Mass"]  = 700
            temp_dict[f"Stage1 Stage Mass"]       = 1969
            temp_dict[f"Stage2 Stage Mass"]       = 670
            #temp_dict[f"Stage3 Stage Mass"] only available in JSR
        elif name == "Diamant BP.4":
            temp_dict[f"Stage1 Propellant Mass"]  = 17000
            temp_dict[f"Stage2 Propellant Mass"]  = 4000
            temp_dict[f"Stage3 Propellant Mass"]  = 680
            temp_dict[f"Stage1 Stage Mass"]       = 1969 # assuming same 1st stage mass as predecssor
            temp_dict[f"Stage2 Stage Mass"]       = 780
            temp_dict[f"Stage3 Stage Mass"]       = 100

    if name.startswith("Electron"):
        
        # Prop masses:
        # DISCOSweb: 9250, 2050, 245
        # SLR      : 9250, 2050, N/A
        # Sp101    : 9250, 2150, 
        # Manual   : N/A,  N/A,  N/A
                
        # Stage masses:
        # DISCOSweb: 950, 250, 55
        # SLR      : 950, 250, N/A  
        # Sp101    : 950, 250, N/A
        # Manual   : N/A, N/A, 40
        
        # https://www.rocketlabusa.com/assets/Uploads/Rocket-Lab-Launch-Payload-Users-Guide-6.5.pdf
        # https://www.rocketlabusa.com/assets/Uploads/Electron-Payload-User-Guide-7.0.pdf

        # 3rd stage is 'kick-stage', and has a 'liquid bi-propellant' which wikipedia says is hypergolic.
        # Kick-stage is placed into low elliptical orbit, and burns some propellant to move into final orbit.
        # Seems like it is specific to mission, so would have to take an educated guess to the average prop.

        temp_dict[f"Stage1 Propellant Mass"] = 9250             # Using value from SLR/DISCOSweb(wet/dry)/Sp101
        temp_dict[f"Stage2 Propellant Mass"] = (2050+2150) / 2  # SLR/Sp101
        temp_dict[f"Stage3 Propellant Mass"] = 245              # Using value from DISCOSweb(wet/dry)
        temp_dict[f"Stage1 Stage Mass"]      = 950
        temp_dict[f"Stage2 Stage Mass"]      = 250
        temp_dict[f"Stage3 Stage Mass"]      = 40               # https://rocketlabcorp.com/assets/Electron-Payload-User-Guide-7.0-v6.pdf 
        temp_dict[f"Fairing Mass"]           = 44               # https://rocketlabcorp.com/assets/Electron-Payload-User-Guide-7.0-v6.pdf
      
    if name == "Epsilon-2 CLPS":
        # Manual: https://global.jaxa.jp/projects/rockets/epsilon/pdf/EpsilonUsersManual_e.pdf
        # Manual and Sp101 have details for PBS not CLPS upper stage, filling in with SLR.
        temp_dict[f"Stage1 Propellant Mass"] = 66300 # Manual
        temp_dict[f"Stage2 Propellant Mass"] = 15000 # Manual 
        temp_dict[f"Stage3 Propellant Mass"] = 2500  # Manual 
        temp_dict[f"Stage4 Propellant Mass"] = 145   # SLR
        temp_dict[f"Stage1 Stage Mass"]      = 8700  # Manual
        temp_dict[f"Stage2 Stage Mass"]      = 2000  # Manual 
        temp_dict[f"Stage3 Stage Mass"]      = 800   # Manual 
        temp_dict[f"Stage4 Stage Mass"]      = 155   # SLR
        temp_dict[f"Fairing Mass"]           = 1000  # Manual 
        
        temp_dict[f"Stage2 Propellant Name"] = "HTPB"
        temp_dict[f"Stage3 Propellant Name"] = "HTPB"
        temp_dict[f"Stage4 Propellant Name"] = "Hydrazine"
        temp_dict[f"Stage2 Fuel Type"]       = "Solid"
        temp_dict[f"Stage3 Fuel Type"]       = "Solid"
        temp_dict[f"Stage4 Fuel Type"]       = "Hypergolic"
        
    if name == "Falcon Heavy":
    
        # Prop masses:
        # DISCOSweb: 294300 (???), 395700, 108670
        # SLR      : 407600 (one), 407600, 107200
        # Sp101    : 411000 (one), 411000, 107500
                
        # Stage masses:
        # DISCOSweb: 33100 (???), 25600, 3230
        # SLR      : 17000 (one), 17000, 4500
        # Sp101    : 22500 (one), 25600, 4000
        
        # Manual: https://www.mach5lowdown.com/wp-content/uploads/PUG/falcon-users-guide-2021-09.pdf
        # "Falcon Heavy’s first-stage comprises three Falcon 9 first stages with enhancements provided to strengthen the cores. 
        # Furthermore, Falcon Heavy utilizes the same second stage and same payload fairing as flown on Falcon 9..."
        
        # Fairing mass:
        #   - https://spaceflight101.com/spacerockets/falcon-9-ft/                                          1750 kg
        #   - https://web.archive.org/web/20220411235211/http://www.spacelaunchreport.com/falcon9ft.html    2000 kg
        
        temp_dict[f"Booster Number"] = 2 
        temp_dict[f"Stage0 Propellant Mass"]  = (407600+411000) / 2 * 2
        temp_dict[f"Stage0 Stage Mass"]       = (17000+22500) / 2 * 2
        temp_dict[f"Stage1 Propellant Mass"]  = (407600+411000) / 2
        temp_dict[f"Stage2 Propellant Mass"]  = (107200+107500) / 2
        temp_dict[f"Stage1 Stage Mass"]       = (17000+25600) / 2
        temp_dict[f"Stage2 Stage Mass"]       = (4500+4000) / 2
        temp_dict[f"Fairing Mass"]            = (2000+1750) / 2 # Average of SLR/Sp101
        
    if name == "Falcon 9 v1.2" or (name == "Falcon 9" and variant == "FT5"):
        # This is listed incorrectly on DISCOSweb, as it should actually be Falcon 9 Block 5, not v1.2. 
        
        # Prop masses:
        # DISCOSweb: 294300, 108700
        # SLR      : 418700, 111500
        # Sp101    : 411000, 107500
                
        # Stage masses:
        # DISCOSweb: 33100, 4300
        # SLR      : 27200, 4500
        # Sp101    : 22200, 4000

        # Fairing mass:
        #   - https://spaceflight101.com/spacerockets/falcon-9-ft/                                          1750 kg
        #   - https://web.archive.org/web/20220411235211/http://www.spacelaunchreport.com/falcon9ft.html    2000 kg
        
        temp_dict[f"Stage1 Propellant Mass"] = (418700+411000) / 2
        temp_dict[f"Stage2 Propellant Mass"] = (111500+107500) / 2
        temp_dict[f"Stage1 Stage Mass"]      = (27200+22200) / 2
        temp_dict[f"Stage2 Stage Mass"]      = (4500+4000) / 2
        temp_dict[f"Fairing Mass"]           = (2000+1750) / 2 # Average of SLR/Sp101
    
    if name == "Firefly Alpha":
        # Dry masses on DISCOSweb match those in user manual
        # https://www.mach5lowdown.com/wp-content/uploads/PUG/Firefly_Alpha_PUG_20190830.pdf
        temp_dict[f"Stage1 Propellant Mass"] = 44800-2895 # Using wet-dry value from DISCOSweb.
        temp_dict[f"Stage2 Propellant Mass"] = 8000-909   # Using wet-dry value from DISCOSweb.
        temp_dict[f"Fairing Mass"] = 0
    
    if name == "GSLV Mk II":
        # Its Mk 2+ for 2021 onward launches according to NB and SLR. 
        # DISCOSweb and Sp101 stage 3 details are for CUS-12 older 3rd stage.
        # The Mk 2+ uses the newer CUS-15 3rd stage, which SLR lists.
        # Booster values from source are for each booster.
        
        # Prop masses:
        # DISCOSweb: N/A,   N/A,    N/A,   12800 (Wet-dry)
        # SLR      : 42600, 138200, 39500, 14996 (Mk2+ CUS-15)
        # Sp101    : 42000, 138000, 39400, 12800
        
        # Stage masses:
        # DISCOSweb: N/A,  N/A,   N/A,  2500
        # SLR      : 5600, 28300, 5500, 2601 
        # Sp101    : 5600, 28300, 5500, 2500
        
        temp_dict[f"Booster Number"] = 4
        temp_dict[f"Stage0 Propellant Mass"]  = (42600+42000) / 2 * int(temp_dict[f"Booster Number"])   # SLR/Sp101
        temp_dict[f"Stage1 Propellant Mass"]  = (138200+138000) / 2                                     # SLR/Sp101
        temp_dict[f"Stage2 Propellant Mass"]  = (39500+39400) / 2                                       # SLR/Sp101
        temp_dict[f"Stage3 Propellant Mass"]  = 14996                                                   # SLR
        temp_dict[f"Stage0 Stage Mass"]       = 5600 * int(temp_dict[f"Booster Number"])                # SLR/Sp101
        temp_dict[f"Stage1 Stage Mass"]       = 28300                                                   # SLR/Sp101
        temp_dict[f"Stage2 Stage Mass"]       = 5500                                                    # SLR/Sp101
        temp_dict[f"Stage3 Stage Mass"]       = 2601                                                    # SLR
        temp_dict[f"Fairing Mass"] = 0
        
        temp_dict[f"Stage0 Propellant Name"]  = "UDMH (Unsymmetrical Dimethyl Hydrazine)/N2O4"
        temp_dict[f"Stage1 Propellant Name"]  = "HTPB"
        temp_dict[f"Stage2 Propellant Name"]  = "UDMH (Unsymmetrical Dimethyl Hydrazine)/N2O4"
        temp_dict[f"Stage3 Propellant Name"]  = "LH2 (Liquid Hydrogen)/LOX"
        
    if name == "GSLV Mk III":        
        # Prop masses:
        # DISCOSweb: 205000 (one), 115400 (W-D), 27700 (W-D)
        # SLR      : 410000 (two), 116000,       28600
        # Sp101    : 206690 (one), 115000,       28000
        
        # Stage masses:
        # DISCOSweb: 31000 (one), 9600,  5300
        # SLR      : 62000 (two), 9600,  4700
        # Sp101    : 31300 (one), 10600, 5000
        
        temp_dict[f"Booster Number"] = 2
        temp_dict[f"Stage0 Propellant Mass"] = (206690+205000)
        temp_dict[f"Stage1 Propellant Mass"]  = 115500
        temp_dict[f"Stage2 Propellant Mass"]  = (28600+28000) / 2
        temp_dict[f"Stage0 Stage Mass"]      = (31300+31000)
        temp_dict[f"Stage1 Stage Mass"]       = (9600+10600) / 2
        temp_dict[f"Stage2 Stage Mass"]       = (4700+5000) / 2
        temp_dict[f"Fairing Mass"] = 0
    
    if "H-IIA" in name:
        # Details from https://global.jaxa.jp/projects/rockets/h2a/ 
        temp_dict[f"Stage1 Propellant Mass"] = 101100
        temp_dict[f"Stage2 Propellant Mass"] = 16900 
        temp_dict[f"Stage1 Stage Mass"]      = 12900 
        temp_dict[f"Stage2 Stage Mass"]      = 3100  
        temp_dict[f"Fairing Mass"]           = 1400 # website/Sp101/SLR
        
        if name == "H-IIA 202":
            temp_dict[f"Booster Number"] = 2
        elif name == "H-IIA 204":
            temp_dict[f"Booster Number"] = 4
        
        temp_dict[f"Stage0 Propellant Mass"] = 130000 / 2 * int(temp_dict[f"Booster Number"]) 
        temp_dict[f"Stage0 Stage Mass"]      = 21000 / 2 * int(temp_dict[f"Booster Number"])
        
    if name == "H-IIB":
        # Details from https://global.jaxa.jp/projects/rockets/h2b/
        temp_dict[f"Booster Number"] = 4
        temp_dict[f"Stage0 Propellant Mass"] = 263800
        temp_dict[f"Stage1 Propellant Mass"]  = 177800
        temp_dict[f"Stage2 Propellant Mass"]  = 16600 
        temp_dict[f"Stage0 Stage Mass"]      = 306000-263800 
        temp_dict[f"Stage1 Stage Mass"]       = 202000-177800 
        temp_dict[f"Stage2 Stage Mass"]       = 20000-16600  
        temp_dict[f"Fairing Mass"]            = 3200 # website/Sp101/SLR
     
    if name == "Jielong-3":
        # Apparently this uses the same motors as the Zhongke-1 rocket (AKA PR-1, Kinetica-1, Lijian-1) 
        # Details of the Zhongke-1 rocket are at http://www.cas-space.com/product?t=1
        # See excel file for comparison details. Best represented by Epsilon-2 CPLS.
        temp_dict[f"Stage1 Propellant Mass"] = 66300 # Epsilon-2 CPLS.
        temp_dict[f"Stage2 Propellant Mass"] = 15000 # Epsilon-2 CPLS.
        temp_dict[f"Stage3 Propellant Mass"] = 2500  # Epsilon-2 CPLS.
        temp_dict[f"Stage4 Propellant Mass"] = 145   # Epsilon-2 CPLS.
        temp_dict[f"Stage1 Stage Mass"]      = 8700  # Epsilon-2 CPLS.
        temp_dict[f"Stage2 Stage Mass"]      = 2000  # Epsilon-2 CPLS.
        temp_dict[f"Stage3 Stage Mass"]      = 800   # Epsilon-2 CPLS.
        temp_dict[f"Stage4 Stage Mass"]      = 155   # Epsilon-2 CPLS.
        temp_dict[f"Fairing Mass"]           = 1000  # Epsilon-2 CPLS.
    
    if name == "Juno II":
        # http://www.astronautix.com/j/junoii.html
        # https://en.wikipedia.org/wiki/Juno_II#Specifications 
        # Wiki assumes the masses by using the fact the 4th stage is one baby seargent engine, and stages 2-3 are multiples of this.
        # http://www.astronautix.com/c/castorengine.html This is one baby seargent engine.
        # http://www.astronautix.com/j/jupiterstage.html This is the first stage.
        temp_dict[f"Stage1 Propellant Mass"] = 48988
        temp_dict[f"Stage2 Propellant Mass"] = 231
        temp_dict[f"Stage3 Propellant Mass"] = 63
        temp_dict[f"Stage1 Stage Mass"]      = 5443
        temp_dict[f"Stage2 Stage Mass"]      = 231
        temp_dict[f"Stage3 Stage Mass"]      = 63

        if variant == "(3)":
            temp_dict[f"Stage4 Propellant Mass"] = 0
            temp_dict[f"Stage4 Stage Mass"]      = 0
        else:
            temp_dict[f"Stage4 Propellant Mass"] = 21
            temp_dict[f"Stage4 Stage Mass"]      = 21
    
    if name == "Jupiter C":
        # This is another name for Juno-I.
        # Same as Juno II but with a different first stage http://www.astronautix.com/j/jupitercstage.html
        temp_dict[f"Stage1 Propellant Mass"] = 24540
        temp_dict[f"Stage2 Propellant Mass"] = 231
        temp_dict[f"Stage3 Propellant Mass"] = 63
        temp_dict[f"Stage4 Propellant Mass"] = 21

        temp_dict[f"Stage1 Stage Mass"]      = 3890
        temp_dict[f"Stage2 Stage Mass"]      = 231
        temp_dict[f"Stage3 Stage Mass"]      = 63
        temp_dict[f"Stage4 Stage Mass"]      = 21

        if variant == "JunoI":
            temp_dict[f"Stage4 Propellant Mass"] = 21
            temp_dict[f"Stage4 Stage Mass"]      = 21
        else:
            temp_dict[f"Stage4 Propellant Mass"] = 0
            temp_dict[f"Stage4 Stage Mass"]      = 0

    if name in ["Kuaizhou","Kuaizhou-1","Kuaizhou-1A"]:
        # No data on SLR or Sp101 or CSR. Have to use wet-dry from DISCOSweb. Dry mass already loaded from DW.
        
        temp_dict[f"Stage1 Propellant Mass"] = 16000 # DW
        temp_dict[f"Stage2 Propellant Mass"] = 8000 # DW
        temp_dict[f"Stage3 Propellant Mass"] = 3000 # DW
        temp_dict[f"Stage4 Propellant Mass"] = 1000 # DW (wet-dry)

        temp_dict[f"Stage1 Stage Mass"]      = 621 # DW
        temp_dict[f"Stage2 Stage Mass"]      = 686 # DW
        temp_dict[f"Stage3 Stage Mass"]      = 183 # DW
        temp_dict[f"Stage4 Stage Mass"]      = 200 # DW

        # Kuaizhou-1 - 2 successes and 1 failure.
        # https://web.archive.org/web/20220406013827/http://www.spacelaunchreport.com/kz.html
        # https://chinaspacereport.wordpress.com/launch-vehicles/kuaizhou/
        # https://www.nasaspaceflight.com/2014/11/china-launches-kuaizhou-2-second-launch-24-hours/ 'liquid upper stage'
        # Ryan22 approximated the the final mass using the mean of all four-stage solid rockets.
        # Norbert Brügge gives stage4 propellant mass as ~1.1t, so DISCOSweb figure is reasonable.
        temp_dict[f"Fairing Mass"] = 0
    
    if name == "Kuaizhou-11":
        # Kuaizhou-11 (Stage 1-3, no prop mass or type)
        # Might be a liquid fourth stage for propulsion control, looks like its part of the satellite.
        # Comparing to other three-stage rockets. See excel file for comparison details. 
        # Best represented by Long March 6. 
        temp_dict[f"Stage1 Propellant Name"] = 'Solid'
        temp_dict[f"Stage2 Propellant Name"] = 'Solid'
        temp_dict[f"Stage3 Propellant Name"] = 'Solid'
        temp_dict[f"Stage1 Fuel Type"] = "Solid"
        temp_dict[f"Stage2 Fuel Type"] = "Solid"
        temp_dict[f"Stage3 Fuel Type"] = "Solid"
        temp_dict[f"Stage1 Propellant Mass"] = 76000             # Using CZ-6
        temp_dict[f"Stage2 Propellant Mass"] = (15000+15150) / 2 # Using CZ-6
        temp_dict[f"Stage3 Propellant Mass"] = 480               # Using CZ-6
        temp_dict[f"Stage1 Stage Mass"]      = 7530              # Using CZ-6
        temp_dict[f"Stage2 Stage Mass"]      = 1490              # Using CZ-6
        temp_dict[f"Stage3 Stage Mass"]      = 520               # Using CZ-6
        temp_dict[f"Fairing Mass"]           = 1500              # Using CZ-6
    
    if name == "LauncherOne":
        temp_dict[f"Stage1 Propellant Mass"] = 20000 # NB
        temp_dict[f"Stage2 Propellant Mass"] = 2750  # NB/DISCOSweb(wet/dry)
        temp_dict[f"Stage1 Stage Mass"]      = 3000  # NB
        temp_dict[f"Stage2 Stage Mass"]      = 250   # NB/DISCOSweb(wet/dry)
        temp_dict[f"Fairing Mass"]           = 100   # NB

    if name == "Lambda 4S":
        # https://www.isas.jaxa.jp/e/enterp/rockets/vehicles/l-4s/index.shtml
        temp_dict["Booster Number"]          = 2
        temp_dict[f"Stage0 Propellant Mass"] = 310.0
        temp_dict[f"Stage1 Propellant Mass"] = 3887.0
        temp_dict[f"Stage2 Propellant Mass"] = 1845.0
        temp_dict[f"Stage3 Propellant Mass"] = 547.5
        temp_dict[f"Stage4 Propellant Mass"] = 87.95
        temp_dict[f"Stage0 Stage Mass"]      = 190.0
        temp_dict[f"Stage1 Stage Mass"]      = 2145.55 # stage 1 dry mass - dry mass of other stages (see website)
        temp_dict[f"Stage2 Stage Mass"]      = 1562.5
        temp_dict[f"Stage3 Stage Mass"]      = 394.9
        temp_dict[f"Stage4 Stage Mass"]      = 23.05

    if "Long March (CZ) 2C" in name:
        # Long March (CZ) 2C 
        # https://discosweb.esoc.esa.int/launch-vehicles/157
        # https://en.wikipedia.org/wiki/Long_March_2C
        # http://cgwic.com/Launchservice/LM2C.html
        temp_dict[f"Stage1 Propellant Mass"] = 162706  # http://cgwic.com/Launchservice/LM2C.html
        temp_dict[f"Stage2 Propellant Mass"] = 54667   # http://cgwic.com/Launchservice/LM2C.html
        temp_dict[f"Stage1 Stage Mass"]      = 10000   # SLR
        temp_dict[f"Stage2 Stage Mass"]      = 4000    # SLR
        temp_dict[f"Fairing Mass"]           = 800 # http://www.b14643.de/Spacerockets_1/China/CZ-2C/Description/Frame.htm

        if name == "Long March (CZ) 2C/YZ-1S":
            # No clear info available. DISCOSweb dry mass matches JSR reentry data, so keeping DW values.
            temp_dict[f"Stage3 Propellant Mass"] = 1350 # DW
            temp_dict[f"Stage3 Stage Mass"]      = 2500 # DW
            temp_dict[f"Stage3 Propellant Name"]  = "UDMH (Unsymmetrical Dimethyl Hydrazine)/N2O4"
            temp_dict[f"Stage3 Fuel Type"]        = "Hypergolic"
    
    if "Long March (CZ) 2D" in name:
        # Long March (CZ) 2D 
        # DISCOSweb dry match SLR
        # https://discosweb.esoc.esa.int/launch-vehicles/158
        # https://en.wikipedia.org/wiki/Long_March_2D
        # http://cgwic.com/Launchservice/LM2D.html
        temp_dict[f"Stage1 Propellant Mass"] = 182000  # http://cgwic.com/Launchservice/LM2D.html
        temp_dict[f"Stage2 Propellant Mass"] = 52700   # http://cgwic.com/Launchservice/LM2D.html
        temp_dict[f"Fairing Mass"]           = 750 # http://www.b14643.de/Spacerockets_1/China/CZ-2D/Description/Frame.htm

        if name == "Long March (CZ) 2D/YZ-3":
            # Basing this off YZ-1S
            temp_dict[f"Stage3 Propellant Mass"]   = 1350
            temp_dict[f"Stage3 Stage Mass"]        = 2500
    
    if name == "Long March (CZ) 2F":
        
        # Prop masses:
        # DISCOSweb : 41501 (one), 186310, 86150
        # SLR       : 37800 (one), 187000, 86000
        # Sp101     : 40800 (one), 186500, 91500
                
        # Stage masses:
        # DISCOSweb : 3200 (one), 1100,  5500
        # SLR       : 3200 (one), 9500,  5500
        # Sp101     : 3000 (one), 12500, 5500
        
        temp_dict[f"Booster Number"] = 4
        temp_dict[f"Stage0 Propellant Mass"] = (37800+40800) / 2 * 4   # SLR/Sp101
        temp_dict[f"Stage0 Stage Mass"]      = (3000+3200) / 2 * 4     # SLR/Sp101
        temp_dict[f"Stage1 Stage Mass"]       = (9500+12500) / 2        # SLR/Sp101
        temp_dict[f"Stage1 Propellant Mass"]  = (187000+186500) / 2     # SLR/Sp101
        temp_dict[f"Stage2 Propellant Mass"]  = (86000+91500) / 2       # SLR/Sp101 
        temp_dict[f"Fairing Mass"]            = 6000                    # http://www.b14643.de/Spacerockets_1/China/CZ-2EF/Description/Frame.htm
        
    if "Long March (CZ) 3" in name:
        #These should really be the 3BE/3CE versions, so the 1998/1999 manuals are now outdated.
        
        # Manual: https://www.mach5lowdown.com/wp-content/uploads/PUG/LM-3B-User-Manual-v1999.pdf
        # Manual: https://www.mach5lowdown.com/wp-content/uploads/PUG/LM-3C-User-Manual-v1998.pdf
    
        if name.startswith("Long March (CZ) 3B"):
            temp_dict[f"Booster Number"] = 4
        elif name.startswith("Long March (CZ) 3C"):
            temp_dict[f"Booster Number"] = 2
            
        # Prop masses:
        # DISCOSweb : 37800 (one), 180775, 50000, 18199
        # SLR       : 41100 (one), 186200, 32600, 18200
        # Sp101     : 41100 (one), 186200, 49400, 18193         
    
        # Stage masses:
        # DISCOSweb : N/A,        9000,  5000, 2800
        # SLR       : 3900 (one), 9800,  4000, 2800
        # Sp101     : 3200 (one), 12600, 3848, 2740, 1500 (PLF)
            
        temp_dict[f"Stage0 Propellant Mass"] = 41100 * int(temp_dict[f"Booster Number"])   # SLR/Sp101
        temp_dict[f"Stage1 Propellant Mass"]  = 186200                                      # SLR/Sp101
        temp_dict[f"Stage2 Propellant Mass"]  = (32600+49400) / 2                           # SLR/Sp101
        temp_dict[f"Stage3 Propellant Mass"]  = (18193+18200) / 2                           # SLR/Sp101
        
        temp_dict[f"Stage0 Stage Mass"]     = (3900+3200) / 2 * int(temp_dict[f"Booster Number"]) # SLR/Sp101
        temp_dict[f"Stage1 Stage Mass"]      = (9800+12600)/ 2                                     # SLR/Sp101
        temp_dict[f"Stage2 Stage Mass"]      = (4000+3848) / 2                                     # SLR/Sp101
        temp_dict[f"Stage3 Stage Mass"]      = (2800+2740) / 2                                     # SLR/Sp101
        temp_dict[f"Fairing Mass"]           = 1500 # Sp101

        if name.endswith("YZ-1"):
            temp_dict[f"Stage4 Stage Mass"]       = 5000
            temp_dict[f"Stage4 Propellant Mass"]  = 50000
            temp_dict[f"Stage4 Propellant Name"]  = 'UDMH (Unsymmetrical Dimethyl Hydrazine)/N2O4'
            temp_dict[f"Stage4 Fuel Type"]        = 'Hypergolic'
        
    if name in ["Long March (CZ) 4B","Long March (CZ) 4C"]: 
        
        # Prop masses:
        # DISCOSweb : 183200, 35600, 16198
        # SLR       : 183200, 35550, 14000
              
        # Stage masses:
        # DISCOSweb : 9500, 4000, 1000
        # SLR       : 9500, 4000, 2000
        
        temp_dict[f"Stage2 Propellant Mass"]  = 35550   # SLR
        temp_dict[f"Stage3 Propellant Mass"]  = 14000   # SLR
        temp_dict[f"Stage3 Stage Mass"]       = 2000    # SLR
        temp_dict[f"Fairing Mass"]            = 1800    # http://www.b14643.de/Spacerockets_1/China/CZ-4/Description/Frame.htm
     
    if "Long March (CZ) 5" in name:  
        
        # Prop masses:
        # DISCOSweb : 174000, 158500, 22900
        # SLR       : 152000, 158000, 22900
        # Sp101     : 135000, 158500, 26500
                
        # Stage masses:
        # DISCOSweb : 12000, 18000, 3400
        # SLR       : 12500, 17000, 3520
        # Sp101     : 12000, 18000, 3400
        
        temp_dict[f"Booster Number"] = 4 
        temp_dict[f"Stage0 Stage Mass"]      = (12500+12000) / 2 * 4  
        temp_dict[f"Stage1 Stage Mass"]       = (17000+18000) / 2 
        
        temp_dict[f"Stage0 Propellant Mass"] = (152000+135000) / 2 * 4 
        temp_dict[f"Stage1 Propellant Mass"]  = (158000+158500) / 2    
        
        if name == "Long March (CZ) 5B":
            temp_dict[f"Fairing Mass"]  = 2400       # http://www.b14643.de/Spacerockets_1/China/CZ-5/Description/Frame.htm
        if name in ["Long March (CZ) 5","Long March (CZ) 5/YZ-2"]:
            temp_dict[f"Stage2 Propellant Mass"]  = (22900+26500) / 2
            temp_dict[f"Stage2 Stage Mass"]       = (3520+3400) / 2  
            if "YZ-2" in name: 
                temp_dict[f"Stage3 Propellant Mass"]  = 7000 # DW
                temp_dict[f"Stage3 Stage Mass"]       = 5000 # DW   
            temp_dict[f"Fairing Mass"]  = 3000 # http://www.b14643.de/Spacerockets_1/China/CZ-5/Description/Frame.htm
        
    if name == "Long March (CZ) 6":
        # 3rd stage prop:
        # DISCOSweb:    'PBV' with LH2/LOX propellant.
        # Wikipedia:     N2O4/UDMH propellant.
        # GSP/CSR/CLR:   H2O2/Kerosene -  4 × YF-85.
        temp_dict[f"Stage3 Propellant Name"] = 'Kerosene/H2O2'
        temp_dict[f"Stage3 Fuel Type"]       = "Kerosene"
        
        # Prop masses:
        # DISCOSweb : 75800, 15150, N/A
        # SLR       : 76000, 15000, N/A
        # Sp101     : 76000, 15150, N/A
                
        # Stage masses:
        # DISCOSweb : 8200, 1850, 1000
        # SLR       : 7530, 1490, N/A
        # Sp101     : N/A,  N/A,  N/A
        
        # NB gives propellant weights for gross (93.7t), stage 1 (77.02), and stage 2 (16.2).
        # Subtracting stages 1+2 from gross gives 93.7-77.02-16.2 = 0.48t
        # With stage 3 liftoff weight of 1000, that means empty is 1000-480 = 520
        temp_dict[f"Stage1 Propellant Mass"] = 76000              # SLR/Sp101
        temp_dict[f"Stage2 Propellant Mass"] = (15000+15150) / 2  # SLR/Sp101
        temp_dict[f"Stage3 Propellant Mass"] = 480                # NB
        temp_dict[f"Stage1 Stage Mass"]      = 7530               # SLR
        temp_dict[f"Stage2 Stage Mass"]      = 1490               # SLR
        temp_dict[f"Stage3 Stage Mass"]      = 520                # NB
        temp_dict[f"Fairing Mass"]           = 1500 # http://www.b14643.de/Spacerockets_1/China/CZ-6/Description/Frame.htm
    
    if name == "Long March (CZ) 6A":
        # No info for prop mass or stage mass anywhere.
        # NB: The CZ-6A is completely different to the the CZ-6. The CZ-6A is basically a CZ-7A without the third stage.
        # NB: Four new solid-fuel boosters have been added.
        # The LM 7 lacks the upper stage of the 7A, so we can estimate the fairing from the CZ-7.
        temp_dict[f"Fairing Mass"]            = 2800 # http://www.b14643.de/Spacerockets_1/China/CZ-7/Description/Frame.htm
        
        # From CZ-7A details.
        temp_dict[f"Stage1 Propellant Mass"]  = 174000
        temp_dict[f"Stage2 Propellant Mass"]  = 65000
        temp_dict[f"Stage1 Stage Mass"]       = (12000+12500) / 2
        temp_dict[f"Stage2 Stage Mass"]       = (5500+6000) / 2  
        
        # The 6A uses 4 solid rocket boosters, but no info is available and no other LM use solid boosters.
        # Going to approximate using CZ-7A again.
        temp_dict[f"Booster Number"]           = 4
        temp_dict[f"Stage0 Propellant Mass"]  = (75500+75000) / 2 * int(temp_dict[f"Booster Number"])
        temp_dict[f"Stage0 Stage Mass"]       = 6000 * int(temp_dict[f"Booster Number"])
        
        temp_dict[f"Stage0 Propellant Name"]  = "Solid"
        temp_dict[f"Stage0 Fuel Type"]        = "Solid"
             
    if "Long March (CZ) 7" in name:
        
        # Prop masses:
        # DISCOSweb: 75500 (each), 174000, 65000 (7), 65500 (7A) , 18199 (H18)
        # SLR      : 75500 (each), 174000, 65000                 , 18200 (H18)
        # Sp101    : 75000 (each), 174000, 65000                 , N/A   (H18)
        
        # Stage masses:
        # DISCOSweb: 6000 (each), 12000, 6000 (7), 5900 (7A) , 2800 (H18)   
        # SLR      : 6000 (each), 12500, 5500                , 2800 (H18)   
        # Sp101    : 6000 (each), 12000, 6000                , N/A  (H18)   
        
        temp_dict[f"Booster Number"]          = 4
        temp_dict[f"Stage0 Propellant Mass"]  = (75500+75000) / 2 * int(temp_dict[f"Booster Number"])
        temp_dict[f"Stage0 Stage Mass"]       = 6000 * int(temp_dict[f"Booster Number"])
        temp_dict[f"Stage1 Stage Mass"]       = (12000+12500) / 2
        temp_dict[f"Stage2 Stage Mass"]       = (5500+6000) / 2 
        
        if name == "Long March (CZ) 7A":
            temp_dict[f"Stage2 Propellant Mass"]  = 65000      # SLR/Sp101
            temp_dict[f"Stage3 Propellant Mass"]  = 18200      # SLR
            temp_dict[f"Fairing Mass"]            = 2500 # http://www.b14643.de/Spacerockets_1/China/CZ-7A/Description/Frame.htm
        elif name == "Long March (CZ) 7":
            temp_dict[f"Fairing Mass"]            = 2800 # http://www.b14643.de/Spacerockets_1/China/CZ-7/Description/Frame.htm
        
        if name.endswith("YZ-1A"):
            temp_dict["Stage3 Stage Mass"]       = 3100 # DW for HO stage
            temp_dict["Stage3 Propellant Mass"]  = 26000 - 3100 # DW for HO stage
            temp_dict[f"Fairing Mass"]            = 2800 # http://www.b14643.de/Spacerockets_1/China/CZ-7/Description/Frame.htm
        
    if "Long March (CZ) 8" in name:
        # Long March (CZ) 8
        # https://discosweb.esoc.esa.int/launch-vehicles/102354
        # https://web.archive.org/web/20220411235209/http://www.spacelaunchreport.com/cz5.html
        # https://en.wikipedia.org/wiki/Long_March_8
        temp_dict[f"Booster Number"] = 2
        temp_dict[f"Stage0 Propellant Mass"] = 75500 * 2  # SLR
        temp_dict[f"Stage0 Stage Mass"]      = 6000   * 2 # SLR
        temp_dict[f"Stage1 Stage Mass"]      = 12500      # SLR   
        temp_dict[f"Fairing Mass"] = 2500 # http://www.b14643.de/Spacerockets_1/China/CZ-8/Description/Frame.htm        

        # https://spacenews.com/first-launch-of-long-march-8a-sends-second-group-of-guowang-megaconstellation-satellites-into-orbit/
        # The Long March 8A is an upgraded variant of the standard Long March 8, which debuted in December 2020. 
        # It features the same first stage and side boosters as the original but includes a newly designed 
        # 3.35-meter-diameter hydrogen-oxygen second stage, allowing a wider, 5.2-meter-diameter payload fairing.      

    if name == "Long March (CZ) 11":
        # Long March (CZ) 11.
        # Ryan22 approximated the propellant masses using the masses from the Vega rocket.
        # Fourth Stage is Solid (Wiki, GSP, CSR/NB) or Liquid (SLR)?
        # See excel file for comparison details. Best represented by Minotaur-1.
        # Was originally represented by Minotaur-4, but since adding 2021/2022 rockets, Minotaur-1 is a much closer match.
        temp_dict[f"Stage1 Propellant Mass"] = 35000              # SLR
        temp_dict[f"Stage2 Propellant Mass"] = 6237               # Minotaur-1
        temp_dict[f"Stage3 Propellant Mass"] = 3915               # Minotaur-1
        temp_dict[f"Stage4 Propellant Mass"] = (782+770.2) / 2    # Minotaur-1
        temp_dict[f"Stage1 Stage Mass"]      = 4770               # SLR
        temp_dict[f"Stage2 Stage Mass"]      = 795                # Minotaur-1
        temp_dict[f"Stage3 Stage Mass"]      = (416+391) / 2      # Minotaur-1
        temp_dict[f"Stage4 Stage Mass"]      = (203+102.1) / 2    # Minotaur-1
        temp_dict[f"Fairing Mass"]           = 300                # Minotaur-1
    
    if name == "Minotaur 1":
        # Prop masses:
        # DISCOSweb: 20785 (W-D), 6237(W-D), 3924, 782
        # SLR      : 20785,       6237,      3915, 782
        # Sp101    : 20785 (W-D), 6237(W-D), 3915, 770.2
        
        # Stage masses:
        # DISCOSweb: 2292, 795, 202, 202
        # SLR      : 2292, 795, 416, 203 
        # Sp101    : 2292, 795, 391, 102.1
        
        temp_dict[f"Stage1 Propellant Name"] = 'Solid'
        temp_dict[f"Stage2 Propellant Name"] = 'Solid'
        temp_dict[f"Stage1 Fuel Type"]       = "Solid"
        temp_dict[f"Stage2 Fuel Type"]       = "Solid"
        temp_dict[f"Stage1 Propellant Mass"] = 20785
        temp_dict[f"Stage2 Propellant Mass"] = 6237
        temp_dict[f"Stage3 Propellant Mass"] = 3915
        temp_dict[f"Stage4 Propellant Mass"] = (782+770.2) / 2
        temp_dict[f"Stage3 Stage Mass"]      = (416+391) / 2
        temp_dict[f"Stage4 Stage Mass"]      = (203+102.1) / 2
        temp_dict[f"Fairing Mass"]           = 300
        
    if name == "Minotaur-4":
        # Prop masses:
        # DISCOSweb : N/A,   24900, 6800, 782
        # SLR       : 45370, 24490, 7070, 770
        # Sp101     : 45400, 24500, 7080, 770.2
                
        # Stage masses:
        # DISCOSweb : N/A,  2900, 1400, 202
        # SLR       : 3620, 3180, 640,  410 
        # Sp101     : 3600, 3200, 650,  102.1

        # Fairing mass:
        #   - https://spaceflight101.com/spacerockets/minotaur-iv/                                          450 kg
        #   - https://web.archive.org/web/20220406013816/http://www.spacelaunchreport.com/mintaur4.html     400 kg
        temp_dict[f"Stage1 Propellant Name"] = 'HTPB'
        temp_dict[f"Stage2 Propellant Name"] = 'HTPB'
        temp_dict[f"Stage3 Propellant Name"] = 'NEPE'
        temp_dict[f"Stage4 Propellant Name"] = 'HTPB'
        temp_dict[f"Stage1 Fuel Type"] = "Solid"
        temp_dict[f"Stage2 Fuel Type"] = "Solid"
        temp_dict[f"Stage3 Fuel Type"] = "Solid"
        temp_dict[f"Stage4 Fuel Type"] = "Solid"
        temp_dict[f"Stage1 Propellant Mass"] = (45370+45400) / 2 
        temp_dict[f"Stage2 Propellant Mass"] = (24490+24500) / 2 
        temp_dict[f"Stage3 Propellant Mass"] = (7070+7080) / 2 
        temp_dict[f"Stage4 Propellant Mass"] = 770 
        temp_dict[f"Stage1 Stage Mass"]      = (3620+3600) / 2 
        temp_dict[f"Stage2 Stage Mass"]      = (3180+3200) / 2  
        temp_dict[f"Stage3 Stage Mass"]      = (640+650) / 2  
        temp_dict[f"Stage4 Stage Mass"]      = (410+102) / 2  
        temp_dict[f"Fairing Mass"] = (400+450) / 2
         
    if name == "Nuri":  
        # No info anywhere else, so using DISCOSweb wet mass - dry mass.
        temp_dict[f"Stage1 Propellant Mass"] = 147100-16000   # DW (W-D)
        temp_dict[f"Stage2 Propellant Mass"] = 41200-4400     # DW (W-D)
        temp_dict[f"Stage3 Propellant Mass"] = 12100-1290     # DW (W-D)
        temp_dict[f"Fairing Mass"] = 0
        
        temp_dict[f"Stage1 Propellant Name"] = 'Kerosene/LOX'
        temp_dict[f"Stage2 Propellant Name"] = 'Kerosene/LOX'
        temp_dict[f"Stage3 Propellant Name"] = 'Kerosene/LOX'
        temp_dict[f"Stage1 Fuel Type"] = "Kerosene"
        temp_dict[f"Stage2 Fuel Type"] = "Kerosene"
        temp_dict[f"Stage3 Fuel Type"] = "Kerosene"
     
    if "PSLV-" in name or (name == "PSLV" and variant == "CA"):
        # https://www.isro.gov.in/media_isro/pdf/Missions/PSLVC41/PSLVC41.pdf
        # https://www.isro.gov.in/media_isro/pdf/Missions/PSLV_C44_LaunchKit.pdf
        # https://www.isro.gov.in/media_isro/pdf/Missions/PSLVC60/PSLVC60-mission-brochure-english.pdf

        # CA has less prop in 4th stage
        
        # Prop masses:
        # DISCOSweb    : 72000 (???),    138000, 40500, 7300, 2100
        # SLR          : 12000 (one),    138000, 41500, 7600, 2500 (XL/DL), 2100 (CA)
        # Sp101        : 12000 (one-XL), 138000, 40700, 6700, 2000 (XL/DL), 1600 (CA)
        # Manual (DL)  : 12200 (one),    139000, 4100,  7650, 1600
        # Manual (XL)  : 12200 (one),    138200, 42000, 7600, 2500
                
        # Stage masses:
        # DISCOSweb    : 12000 (???), 30000, 5300, 920,  960
        # SLR          : 2700 (one),  30000, 5400, 700,  420
        # Sp101        : 2010 (one),  30200, 5300, 1100, 920 
        
        # Fairing mass:
        #   - https://spaceflight101.com/spacerockets/pslv/                                         1150 kg
        #   - https://web.archive.org/web/20220406013825/http://www.spacelaunchreport.com/pslv.html 1100 kg    
        
        if name == "PSLV-XL":
            temp_dict[f"Booster Number"] = 6
            temp_dict[f"Stage4 Propellant Mass"]    = (2500+2000) / 2 # SLR/Sp101
        if name == "PSLV-QL":
            temp_dict[f"Booster Number"] = 4
            temp_dict[f"Stage4 Propellant Mass"]    = (2500+2000) / 2 # SLR/Sp101
        elif name == "PSLV-DL":
            temp_dict[f"Booster Number"] = 2
            temp_dict[f"Stage4 Propellant Mass"]    = (2500+2000) / 2 # SLR/Sp101
        elif name == "PSLV-CA" or (variant == "CA"):
            temp_dict[f"Booster Number"] = 0
            temp_dict[f"Stage4 Propellant Mass"]    = (2100+1600) / 2 # SLR/Sp101
        
        temp_dict[f"Stage0 Propellant Mass"]    = 12200 * int(temp_dict[f"Booster Number"]) # Manual
        temp_dict[f"Stage1 Propellant Mass"]    = (139000+138200) / 2                       # Manual
        temp_dict[f"Stage2 Propellant Mass"]    = 41500                                     # Manual
        temp_dict[f"Stage3 Propellant Mass"]    = (7650+7600) / 2                           # Manual

        temp_dict[f"Stage0 Stage Mass"]         = (2700+2010) / 2 * int(temp_dict[f"Booster Number"])   # SLR/Sp101
        temp_dict[f"Stage1 Stage Mass"]         = (30000+30200) / 2 # SLR/Sp101
        temp_dict[f"Stage2 Stage Mass"]         = (5400+5300) / 2   # SLR/Sp101
        temp_dict[f"Stage3 Stage Mass"]         = (700+1100) / 2    # SLR/Sp101
        temp_dict[f"Stage4 Stage Mass"]         = (420+920) / 2     # SLR/Sp101
        temp_dict[f"Fairing Mass"]              = (1150+1100) / 2   # SLR/Sp101

    if "Proton-K" in name:
        # https://web.archive.org/web/20220406013821/http://www.spacelaunchreport.com/proton.html

        temp_dict[f"Stage1 Propellant Mass"] = 419410 # SLR
        temp_dict[f"Stage2 Propellant Mass"] = 156113 # SLR
        temp_dict[f"Stage3 Propellant Mass"] = 46562  # SLR
        
        temp_dict[f"Stage1 Stage Mass"]      = 31000  # SLR
        temp_dict[f"Stage2 Stage Mass"]      = 11750  # SLR
        temp_dict[f"Stage3 Stage Mass"]      = 4185   # SLR
        
        if name == "Proton-K/DM-2":
            temp_dict[f"Stage4 Propellant Mass"] = 15050  # SLR
            temp_dict[f"Stage4 Stage Mass"]      = 2440   # SLR

        temp_dict[f"Fairing Mass"] = 0

    if "Proton-M" in name:
        # DISCOSweb hasn't updated to include the 'Enhanced' variant, which has been standard since 2007.
        # Unsure whether Sp101 has either, so using SLR.
        
        # Prop masses:
        # DISCOSweb: 419400, 157300, 43062, 19800 (Briz-M), 15000 (DM-03)
        # SLR      : 428300, 157300, 46562, 19800 (Briz-M), N/A   (DM-03)
        # Sp101    : 419400, 156113, 46562, 19800 (Briz-M), 18600 (DM-03)
        
        # Stage masses:
        # DISCOSweb: 31000, 11000, 3500, 2370 (Briz-M), 2140 (DM-03)
        # SLR      : 30600, 11000, 3500, 2370 (Briz-M), N/A  (DM-03)
        # Sp101    : 31000, 11715, 4185, 2370 (Briz-M), 3500 (DM-03)
        
        # Fairing masses:
        # Sp101    :  2000 (Briz-M Proton), 1600 (Briz-M Angara A5), 2000 (DM-03)

        temp_dict[f"Fairing Mass"]           = 2000   # Sp101 has this value for Proton Medium and Briz-M / DM-03. Using here.

        if variant == "Enhanced":
            temp_dict[f"Stage1 Propellant Mass"] = 428300 # SLR
            temp_dict[f"Stage2 Propellant Mass"] = 157300 # SLR
            temp_dict[f"Stage3 Propellant Mass"] = 46562  # SLR
            temp_dict[f"Stage1 Stage Mass"]      = 30600  # SLR
            temp_dict[f"Stage2 Stage Mass"]      = 11000  # SLR
            temp_dict[f"Stage3 Stage Mass"]      = 3500   # SLR
            
            if name.endswith("Briz-M"):
                temp_dict[f"Stage4 Propellant Mass"] = 19800  # SLR
                temp_dict[f"Stage4 Stage Mass"]      = 2370   # SLR
        else:
            temp_dict[f"Stage1 Propellant Mass"] = 419410 # SLR
            temp_dict[f"Stage2 Propellant Mass"] = 156113 # SLR
            temp_dict[f"Stage3 Propellant Mass"] = 46562  # SLR
            temp_dict[f"Stage1 Stage Mass"]      = 30600  # SLR
            temp_dict[f"Stage2 Stage Mass"]      = 11400  # SLR
            temp_dict[f"Stage3 Stage Mass"]      = 3700   # SLR

            if name.endswith("Briz-M"):
                temp_dict[f"Stage4 Propellant Mass"] = 19800 # SLR
                temp_dict[f"Stage4 Stage Mass"]      = 2500  # SLR
            
        if name.endswith("DM-3"):
            temp_dict[f"Stage4 Propellant Mass"] = 18600 # Sp101
            temp_dict[f"Stage4 Stage Mass"]      = 3500  # Sp101
       
    if "Pegasus XL" in name:
        # Manual: https://web.archive.org/web/20221230022508/https://www.mach5lowdown.com/wp-content/uploads/PUG/pegasus-user-guide-2007-1.pdf
        # Prop masses:
        # DISCOSweb: 15024, 3924, 782
        # SLR      : 15010, 3930, 770
        # Sp101    : 15014, 3915, 770.2
        # Manual   : 15032, 3923, 770
        
        # Stage masses:
        # DISCOSweb: 2886, 202, 202
        # SLR      : 1370, 420, 110
        # Sp101    : 1369, 391, 102.1, 
        # Manual   : 1386, 416, 108  
        
        temp_dict[f"Stage1 Propellant Mass"] = 15032
        temp_dict[f"Stage2 Propellant Mass"] = 3923
        temp_dict[f"Stage3 Propellant Mass"] = 770
        temp_dict[f"Stage1 Stage Mass"]      = 1386 
        temp_dict[f"Stage2 Stage Mass"]      = 416
        temp_dict[f"Stage3 Stage Mass"]      = 108 
        temp_dict[f"Fairing Mass"] = 0
    
    if name == "Qased":

        # First stage uses UDMH/N2O4, second/third are solid fueled.
        # First stage says oxidiser is AK27 on DISCOSweb/SLR/NB, N2O4 on Wiki.
        
        # Upper stage is misassigned in all sources (wiki)
        # Approximating third stage by averaging values from each source.
        
        # Prop masses:
        # DISCOSweb : 12330, 1500, 770
        # SLR       : 14400, N/A,  240
        # NB        : 12330, 1500, 185
                
        # Stage masses:
        # DISCOSweb : 2790, 400, 30
        # SLR       : 2700, N/A, 55
        # NB        : 2370, 200, 315

        temp_dict[f"Stage1 Propellant Name"] = 'UDMH (Unsymmetrical Dimethyl Hydrazine)/AK-27I'
        temp_dict[f"Stage2 Propellant Name"] = 'Solid'
        temp_dict[f"Stage3 Propellant Name"] = 'Solid'
        temp_dict[f"Stage1 Fuel Type"]       = "Solid"
        temp_dict[f"Stage2 Fuel Type"]       = "Solid"
        temp_dict[f"Stage3 Fuel Type"]       = "Solid"
        temp_dict[f"Stage1 Propellant Mass"] = 14400 # SLR
        temp_dict[f"Stage2 Propellant Mass"] = 1500  # NB 
        temp_dict[f"Stage3 Propellant Mass"] = 398   # Misassigned in all sources, using average of SLR and NB.
        temp_dict[f"Stage1 Stage Mass"]      = 2700  # SLR
        temp_dict[f"Stage2 Stage Mass"]      = 200   # NB
        temp_dict[f"Stage3 Stage Mass"]      = (30+55+315)/ 3 # Misassigned in all sources, using average of SLR and NB.
        temp_dict[f"Fairing Mass"]           = 100 # SLR
    
    if name == "Safir":
        # https://web.archive.org/web/20220411235210/http://www.spacelaunchreport.com/safir.html
        temp_dict[f"Stage1 Propellant Mass"] = 18200 # SLR
        temp_dict[f"Stage2 Propellant Mass"] = 3050  # SLR
        temp_dict[f"Stage1 Stage Mass"]      = 3500  # SLR
        temp_dict[f"Stage2 Stage Mass"]      = 350   # SLR
        temp_dict[f"Fairing Mass"]           = 100   # SLR

    if name == "Scout F-1":
        # http://www.astronautix.com/s/scoutf-1.html
        temp_dict[f"Stage3 Propellant Mass"] = 1100
        temp_dict[f"Stage3 Stage Mass"]      = 300

    #if name in ["Scout X-2M","Scout X-3M"]: 
    #    temp_dict[f"Stage4 Propellant Mass"] = 
    #    temp_dict[f"Stage4 Stage Mass"]      = 

    if "Shavit" in name:
        # Shavit (Stage 1-3, no prop mass or type)
        # https://discosweb.esoc.esa.int/launch-vehicles/76
        # https://space.skyrocket.de/doc_lau/shavit.htm
        # https://web.archive.org/web/20220406013832/http://www.spacelaunchreport.com/shavit.html
        # https://www.rafael.co.il/wp-content/uploads/2021/07/RAFAEL-SPACE-PROPULSION-2021-CATALOGUE-2.pdf
        # https://www.iai.co.il/sites/default/files/2020-07/Shavit%20Brochure.pdf
        # The stages are listed on DISCOSweb as the RSA-3 vehicle from South Africa, however the launch in 2020 was for Shavit-2.
        # Fourth stage isn't included on DISCOSweb or NB, but is included on SLR and mentioned on wiki.
        # DISCOSweb(wet/dry) are similar to values for Shavit/Shavit-1. Will use SLR values as these are correct for Shavit-2.

        temp_dict[f"Stage1 Propellant Name"] = 'HTPB'
        temp_dict[f"Stage2 Propellant Name"] = 'HTPB'
        temp_dict[f"Stage3 Propellant Name"] = 'HTPB'
        temp_dict[f"Stage1 Fuel Type"] = "Solid"
        temp_dict[f"Stage2 Fuel Type"] = "Solid"
        temp_dict[f"Stage3 Fuel Type"] = "Solid"

        if name == "Shavit":
            temp_dict[f"Stage1 Propellant Mass"] = 9100 # SLR
            temp_dict[f"Stage1 Stage Mass"]      = 1115  # SLR
            temp_dict[f"Stage2 Propellant Mass"] = 9100 # SLR
            temp_dict[f"Stage2 Stage Mass"]      = 1288  # SLR
        elif name == "Shavit 1":
            temp_dict[f"Stage1 Propellant Mass"] = 12750 # SLR
            temp_dict[f"Stage1 Stage Mass"]      = 1240  # SLR
            temp_dict[f"Stage2 Propellant Mass"] = 9100 # SLR
            temp_dict[f"Stage2 Stage Mass"]      = 1288  # SLR
        elif name == "Shavit 2":
            temp_dict[f"Stage1 Propellant Mass"] = 12750 # SLR
            temp_dict[f"Stage1 Stage Mass"]      = 1240  # SLR
            temp_dict[f"Stage2 Propellant Mass"] = 12750 # SLR
            temp_dict[f"Stage2 Stage Mass"]      = 1376  # SLR

        temp_dict[f"Stage3 Propellant Mass"] = 1890  # SLR
        temp_dict[f"Stage3 Stage Mass"]      = 684   # SLR
        temp_dict[f"Fairing Mass"]           = 57 # http://www.b14643.de/Spacerockets_1/Rest_World/Shavit/Description/Frame.htm
      
    if name in ["Shuang Quxian-1","Shuang Quxian 1"]: 
        # Info only available through DISCOSweb. However, there is a small typo in the stage 1 prop mass.
        # Adding DW values here so it works when getting info from JSR.

        temp_dict[f"Stage1 Propellant Mass"] = 14356 # SLR
        temp_dict[f"Stage2 Propellant Mass"] = 8776  # SLR
        temp_dict[f"Stage3 Propellant Mass"] = 2618  # SLR
        temp_dict[f"Stage4 Propellant Mass"] = 1050  # SLR
        temp_dict[f"Stage1 Stage Mass"]      = 1744  # SLR
        temp_dict[f"Stage2 Stage Mass"]      = 974   # SLR
        temp_dict[f"Stage3 Stage Mass"]      = 282   # SLR
        temp_dict[f"Stage4 Stage Mass"]      = 300   # SLR
        temp_dict[f"Fairing Mass"] = 0
        # https://webcache.googleusercontent.com/search?q=cache:rNwod1W3xxcJ:www.i-space.com.cn/statics/ispace/doc/Hyperbola-1%2520User%2520Manual.pdf&hl=en&gl=uk
        # According to this cached user manual and NB, all stages are solid. 
        
    if name == "Simorgh":
        # Simorgh
        # https://discosweb.esoc.esa.int/launch-vehicles/88324
        # https://web.archive.org/web/20220406013833/http://www.spacelaunchreport.com/simorgh.html
        # https://en.wikipedia.org/wiki/Saman-1_(rocket_stage)
        # Saman-1 is an upper-stage, used as a space tug to transfer satellites to higher orbits.
        temp_dict[f"Stage3 Propellant Name"] = 'Solid'
        temp_dict[f"Stage1 Fuel Type"] = "Hypergolic"
        temp_dict[f"Stage2 Fuel Type"] = "Hypergolic"
        temp_dict[f"Stage3 Fuel Type"] = "Solid"
        temp_dict[f"Stage1 Propellant Mass"] = 63000 # SLR (v. similar to DISCOSwb(wet/dry))
        temp_dict[f"Stage2 Propellant Mass"] = 8000  # SLR (v. similar to DISCOSwb(wet/dry))
        temp_dict[f"Stage3 Propellant Mass"] = 185   # SLR (v. similar to DISCOSwb(wet/dry)) 
        temp_dict[f"Stage1 Stage Mass"]      = 12700 # SLR (v. similar to DISCOSwb(wet/dry))
        temp_dict[f"Stage2 Stage Mass"]      = 1600  # SLR (v. similar to DISCOSwb(wet/dry))
        temp_dict[f"Stage3 Stage Mass"]      = 230   # SLR (v. similar to DISCOSwb(wet/dry)) 
        temp_dict[f"Fairing Mass"]           = 120   # SLR (same value in NB)
    
    if name == "Space Launch System - Block 1 Crew":
        # https://www.nasa.gov/wp-content/uploads/2022/03/sls-reference-guide-2022-v2-print-0.pdf
        
        # Prop masses:
        # DISCOSweb : 1295300 (???), 978900,  27200
        # SLR (2012): 631495  (one), 979452,  26853, 8165 (PLF)
        # Manual    : 626700  (one), 1002700, 29000
                
        # Stage masses:
        # DISCOSweb: 192000 (???), 89400,  7300
        # SLR      : 100390 (one), 112000, 4355 
        # Manual   : 99300  (one), 99300,  3700
        
        # User Manual values.
        temp_dict[f"Booster Number"]            = 2
        temp_dict[f"Stage0 Propellant Mass"]   = 626700 * int(temp_dict[f"Booster Number"])
        temp_dict[f"Stage1 Propellant Mass"]    = 1002700
        temp_dict[f"Stage2 Propellant Mass"]    = 29000
        temp_dict[f"Stage0 Stage Mass"]        = 99300 * int(temp_dict[f"Booster Number"])
        temp_dict[f"Stage1 Stage Mass"]         = 99300
        temp_dict[f"Stage2 Stage Mass"]         = 3700
        temp_dict[f"Fairing Mass"]              = 8165 # SLR
    
    if name == "Strela":
        # https://web.archive.org/web/20230614172003/https://spaceflight101.com/spacerockets/strela/
        temp_dict[f"Stage1 Stage Mass"]         = 5700
        temp_dict[f"Stage2 Stage Mass"]         = 1500
        temp_dict[f"Stage3 Stage Mass"]         = 725
        temp_dict[f"Stage1 Propellant Mass"]    = 71500
        temp_dict[f"Stage2 Propellant Mass"]    = 10700
        temp_dict[f"Stage3 Propellant Mass"]    = 375
    
    if name == "SSLV":
        # https://www.isro.gov.in/mission_SSLV_D2.html
        # https://www.isro.gov.in/media_isro/pdf/Missions/SSLV/SSLV_D2_EOS_07_DigitalBrochure.pdf
        
        # Stage masses:
        # DISCOSweb: 12000, 300,  1500, 450
        # JSR      : 8000,  1000, 400,  200
        # NB       : 13000, 4300, 2500, 100
        # DISCOSweb stage 2 looks wrong, lets average JSR/NB. 
        temp_dict[f"Stage1 Stage Mass"]         = (8000+13000) / 2
        temp_dict[f"Stage2 Stage Mass"]         = (1000+4300) / 2
        temp_dict[f"Stage3 Stage Mass"]         = (400+2500) / 2
        temp_dict[f"Stage4 Stage Mass"]         = (200+100) / 2
        temp_dict[f"Stage1 Propellant Mass"]    = 87000
        temp_dict[f"Stage4 Propellant Mass"]    = 50
        temp_dict[f"Fairing Mass"]              = 90 # NB
    
    if "Soyuz" in name:
        
        # Soyuz-2-1A Fregat
        # Soyuz-2-1A
        # Soyuz-2-1A Fregat-M
        # Soyuz-ST-A Fregat-M	
        
        # Soyuz-2-1B
        # Soyuz-ST-B Fregat-MT
        # Soyuz-2-1B Fregat
        # Soyuz-2-1B Fregat-M
        
        # Soyuz-2-1V	
        # Soyuz-2-1V Volga
        	
        # https://web.archive.org/web/20230424155119/https://www.mach5lowdown.com/wp-content/uploads/PUG/Soyuz-Users-Manual-March-2012.pdf Details for 2.1b ST
        # https://web.archive.org/web/20221230022919/https://www.mach5lowdown.com/wp-content/uploads/PUG/soyuz_users_manual_190401-1.pdf Details for 2.1a    
        # A to B version only changes Stage 2 engine from RD-0110 to RD-0124.
        # V Version keeps B stages but updates the 1st stage engine.
        
        # Sort out boosters
        if name in ["Soyuz-2-1V","Soyuz-2-1V Volga"]:
            temp_dict[f"Booster Number"] = 0
        else:
            temp_dict[f"Booster Number"] = 4
            temp_dict[f"Stage0 Propellant Mass"] = (27900 + 11260) * 4     # Manual
            temp_dict[f"Stage0 Stage Mass"]      = 3784 * 4                # Manual
        
        # Sort out first stage
        if name in ["Soyuz-2-1V","Soyuz-2-1V Volga"]:
            temp_dict[f"Stage1 Propellant Mass"]  = 119000 # SLR
            temp_dict[f"Stage1 Stage Mass"]       = 9300   # Sp101
        else:
            temp_dict[f"Stage1 Propellant Mass"]  = 63800 + 26300       # Manual
            temp_dict[f"Stage1 Stage Mass"]       = 6545                # Manual
        
        # Sort out second stage
        if name in ["Soyuz-2-1B Fregat","Soyuz-2-1B Fregat-M","Soyuz-2-1V","Soyuz-2-1V Volga"]:
            temp_dict[f"Stage2 Stage Mass"]       = 2355            # Manual
            temp_dict[f"Stage2 Propellant Mass"]  = 17800 + 7600    # Manual
        elif name in ["Soyuz-2-1A","Soyuz-2-1A Fregat-M","Soyuz-2-1A Fregat"]: 
            temp_dict[f"Stage2 Stage Mass"]       = 2410            # Manual
            temp_dict[f"Stage2 Propellant Mass"]  = 22790           # Manual
        elif name in ["Soyuz-ST-A Fregat-M"]: 
            temp_dict[f"Stage2 Stage Mass"]       = 2470            # Manual
            temp_dict[f"Stage2 Propellant Mass"]  = 22830           # Manual
            
        # Sort out third stage
        if name in ["Soyuz-2-1A Fregat","Soyuz-2-1B Fregat"]:
            temp_dict[f"Stage3 Stage Mass"]       = 1000            # Manual and SLR
            temp_dict[f"Stage3 Propellant Mass"]  = 5350            # Manual and SLR
        elif name in ["Soyuz-2-1A Fregat-M","Soyuz-ST-A Fregat-M","Soyuz-2-1B Fregat-M"]:
            temp_dict[f"Stage3 Stage Mass"]       = 960             # Sp101
            temp_dict[f"Stage3 Propellant Mass"]  = 5750            # Sp101
        elif name in ["Soyuz-ST-B Fregat-MT"]:
            temp_dict[f"Stage3 Stage Mass"]       = 1050            # Sp101
            temp_dict[f"Stage3 Propellant Mass"]  = 7100            # Sp101
        elif name in ["Soyuz-2-1V Volga"]:
            temp_dict[f"Stage3 Stage Mass"]       = 890 # Sp101
            # Variable from 300 - 900 kg (SLR and Sp101)
            temp_dict[f"Stage3 Propellant Mass"]  = 900
            
        # Sort out fairing
        # Fairing mass:
        #   - http://www.b14643.de/Spacerockets_1/East_Europe_1/Soyuz-2/Description/Frame.htm           1000 kg (Soyuz-2 (1a)	Soyuz-2K (1b))
        #   - http://www.b14643.de/Spacerockets_1/East_Europe_1/Soyuz-2/Description/Frame.htm           1800 kg (Soyuz-2/Fregat	Soyuz-2K/Fregat)
        #   - http://www.b14643.de/Spacerockets_1/East_Europe_1/Soyuz-ST/Description/Frame.htm          1800 kg (Soyuz-2K/Fregat)
        #   - http://www.b14643.de/Spacerockets_1/West_Europe/Soyuz-ESA/Description/Frame.htm           1800 kg (Soyuz STK-A, Soyuz STK-B - Fregat and Fregat-M/MT)
        #   - https://www.mach5lowdown.com/wp-content/uploads/PUG/Soyuz-Users-Manual-March-2012.pdf     1700 kg (ST Fairing for Arianespace with Fregat)
        if name in ["Soyuz-ST-A Fregat-M","Soyuz-ST-B Fregat-MT"]:
            temp_dict[f"Fairing Mass"] = 1700
        elif name in ["Soyuz-2-1A","Soyuz-2-1B","Soyuz-2-1V","Soyuz-2-1V Volga"]:
            temp_dict[f"Fairing Mass"] = 1000 # NB
        elif name in ["Soyuz-2-1A Fregat","Soyuz-2-1A Fregat-M","Soyuz-2-1B Fregat","Soyuz-2-1B Fregat-M"]:
            temp_dict[f"Fairing Mass"] = 1800 # NB
    
    if name in ["Taurus 3110","Taurus 3210","Minotaur-C 3210"]:
        # Stage 2 - Orion 50SXLG
        # https://cdn.northropgrumman.com/-/media/wp-content/uploads/Orion-Motor-Series.pdf?v=1.0.0
        temp_dict[f"Stage2 Propellant Mass"]    = 33145 * 0.453592 # Manual (weights in lbm)
        temp_dict[f"Stage2 Stage Mass"]         = 2237 * 0.453592 # Manual (weights in lbm)

    if name == "Vega":            
        # Manual - https://www.arianespace.com/wp-content/uploads/2018/05/Vega-Users-Manual_Issue-04_April-2014.pdf
        temp_dict[f"Stage1 Propellant Mass"]    = 87710 # Manual
        temp_dict[f"Stage2 Propellant Mass"]    = 23814 # Manual
        temp_dict[f"Stage3 Propellant Mass"]    = 10567 # Manual
        temp_dict[f"Stage4 Propellant Mass"]    = 577   # Manual
        temp_dict[f"Stage1 Stage Mass"]         = 8533  # Manual
        temp_dict[f"Stage2 Stage Mass"]         = 2486  # Manual
        temp_dict[f"Stage3 Stage Mass"]         = 1433  # Manual
        temp_dict[f"Stage4 Stage Mass"]         = 688   # Manual
        temp_dict[f"Fairing Mass"]              = 540   # Manual
        
    if name == "Vega C":            
        # Manual - https://www.mach5lowdown.com/wp-content/uploads/PUG/vega-c-user-manual-issue-0-revision-0_20180705.pdf
        temp_dict[f"Stage1 Propellant Mass"]    = 141634 # Manual
        temp_dict[f"Stage2 Propellant Mass"]    = 36239  # Manual
        temp_dict[f"Stage3 Propellant Mass"]    = 10567  # Manual
        temp_dict[f"Stage4 Propellant Mass"]    = 740    # Manual
        temp_dict[f"Stage1 Stage Mass"]         = 13393  # Manual
        temp_dict[f"Stage2 Stage Mass"]         = 4238   # Manual
        temp_dict[f"Stage3 Stage Mass"]         = 1433   # Manual
        temp_dict[f"Stage4 Stage Mass"]         = 698    # Manual
        temp_dict[f"Fairing Mass"]              = 860    # Manual
    
    if name == "Zhongke 1A":      
        # See excel file for comparison details. Best represented by Vega C.
        temp_dict[f"Stage1 Fuel Type"] = "Solid"
        temp_dict[f"Stage2 Fuel Type"] = "Solid"
        temp_dict[f"Stage3 Fuel Type"] = "Solid"
        temp_dict[f"Stage4 Fuel Type"] = 'Solid'
        temp_dict[f"Stage1 Propellant Mass"]    = 141634 # Vega-C
        temp_dict[f"Stage2 Propellant Mass"]    = 36239  # Vega-C
        temp_dict[f"Stage3 Propellant Mass"]    = 10567  # Vega-C
        temp_dict[f"Stage4 Propellant Mass"]    = 740    # Vega-C
        temp_dict[f"Stage1 Stage Mass"]         = 13393  # Vega-C
        temp_dict[f"Stage2 Stage Mass"]         = 4238   # Vega-C
        temp_dict[f"Stage3 Stage Mass"]         = 1433   # Vega-C
        temp_dict[f"Stage4 Stage Mass"]         = 698    # Vega-C
        temp_dict[f"Fairing Mass"]              = 860    # Vega-C
    
    if "Zhuque-2" in name:    
        # See excel file for comparison details. Best represented by Antares 230.
        temp_dict[f"Stage1 Propellant Mass"] = (242400+242000) / 2 # Antares 230
        temp_dict[f"Stage1 Stage Mass"]      = (18800+20600)   / 2 # Antares 230
        temp_dict[f"Stage2 Propellant Mass"] = 24924               # Antares 230
        temp_dict[f"Stage2 Stage Mass"]      = 1392                # Antares 230
        temp_dict[f"Fairing Mass"]           = 971                 # Antares 230
        
        # https://spacenews.com/landspace-puts-2-satellites-in-orbit-with-enhanced-zhuque-2-rocket/
        # "vernier thrusters replaced by a vector control system, saving 400 kilograms in mass."
        # So we will just use the same masses as Antares 230 but with a smaller stage 2 mass.
        if name == "Zhuque-2E":
            temp_dict[f"Stage2 Stage Mass"] = temp_dict[f"Stage2 Stage Mass"] - 400

    ############################
    # New Vehicles in 2023/2024
    ############################

    if name == "Goche Yeollyo Uju Balsache (GYUB) - TV2":
        # No Information available for this vehicle.
        #Used Information Similar to Qased
        temp_dict[f"Stage1 Propellant Mass"] = 14400 # SLR
        temp_dict[f"Stage2 Propellant Mass"] = 1500  # NB 
        temp_dict[f"Stage3 Propellant Mass"] = 398   # Misassigned in all sources, using average of SLR and NB.
        temp_dict[f"Stage1 Stage Mass"]      = 2700  # SLR
        temp_dict[f"Stage2 Stage Mass"]      = 200   # NB
        temp_dict[f"Stage3 Stage Mass"]      = (30+55+315)/ 3 # Misassigned in all sources, using average of SLR and NB.
        temp_dict[f"Fairing Mass"]           = 100 # SLR

    if name == "H-III 22":
        # Manual - https://www.mhi.com/products/space/launch_srv_lineup.html#pdh3
        # Approximating using H2B, as SLR says that "In many ways, H3 appears to be an improved H-2B."
        temp_dict[f"Booster Number"] = 2
        temp_dict[f"Stage0 Propellant Mass"] = 263800 / 4 * temp_dict[f"Booster Number"]
        temp_dict[f"Stage1 Propellant Mass"]  = 177800
        temp_dict[f"Stage2 Propellant Mass"]  = 16600 
        temp_dict[f"Stage0 Stage Mass"]      = (306000-263800) / 4 * temp_dict[f"Booster Number"]
        temp_dict[f"Stage1 Stage Mass"]       = 202000-177800 
        temp_dict[f"Stage2 Stage Mass"]       = 20000-16600  
        temp_dict[f"Fairing Mass"]            = 3200 # website/Sp101/SLR
    
    if name == "RS1":
        # Manual - https://ablspacesystems.com/wp-content/uploads/2022/06/ABL-Payload-Users-Guide-2022-V1.pdf
        # Used Information Similar to Firefly Alpha
        temp_dict[f"Stage1 Propellant Name"] = 'Liquid Oxygen (LOX)/RP-1' # Manual
        temp_dict[f"Stage2 Propellant Name"] = 'Liquid Oxygen (LOX)/RP-1' # Manual 
        temp_dict[f"Stage1 Fuel Type"] = 'Kerosene' 
        temp_dict[f"Stage2 Fuel Type"] = 'Kerosene'
        temp_dict[f"Stage1 Propellant Mass"]    = 44800-289
        temp_dict[f"Stage2 Propellant Mass"]    = 8000-909
        temp_dict[f"Stage1 Stage Mass"]         = 2895
        temp_dict[f"Stage2 Stage Mass"]         = 909

    if name == "Starship":
        pass
        # No information available for this vehicle.

    if name == "Terran-1":
        # https://www.relativityspace.com/terran-r
        # Used Information Similar to Zhuque-2, which itself is based on Antares 230
        temp_dict[f"Stage1 Fuel Type"] = 'Methane'
        temp_dict[f"Stage2 Fuel Type"] = 'Methane'
        temp_dict[f"Stage1 Propellant Name"]    = 'Liquid Oxygen (LOX)/Subcooled Methane' # Manual
        temp_dict[f"Stage2 Propellant Name"]    = 'Liquid Oxygen (LOX)/Subcooled Methane' # Manual
        temp_dict[f"Stage1 Propellant Mass"]    = (242400+242000) / 2
        temp_dict[f"Stage2 Propellant Mass"]    = 24924 # CASTOR documentation
        temp_dict[f"Stage1 Stage Mass"]         = (18800+20600)   / 2
        temp_dict[f"Stage2 Stage Mass"]         = 1392  # CASTOR documentation

    if name == "Tianlong 2":
        # No information available for this vehicle.
        #https://space.skyrocket.de/doc_lau/tianlong-2.htm
        # Used Information Similar to Electron (Curie)
        temp_dict[f"Stage3 Fuel Type"] = 'Hypergolic' #SLR
        temp_dict[f"Stage1 Propellant Mass"]    = 9250
        temp_dict[f"Stage2 Propellant Mass"]    = (2050+2150) / 2
        temp_dict[f"Stage3 Propellant Mass"]    = 245
        temp_dict[f"Stage1 Stage Mass"]         = 950
        temp_dict[f"Stage2 Stage Mass"]         = 250
        temp_dict[f"Stage3 Stage Mass"]         = 40 

    if name == "Ariane 62":
        # No other information available for this vehicle.
        # Boosters are P120C from Vega C
        # Manual - https://www.mach5lowdown.com/wp-content/uploads/PUG/vega-c-user-manual-issue-0-revision-0_20180705.pdf
        # Stage Mass information taken from H-IIA 202
        temp_dict[f"Booster Number"] = 2
        temp_dict[f"Stage0 Stage Mass"]        = 13393  * int(temp_dict[f"Booster Number"]) # Manual
        temp_dict[f"Stage0 Propellant Mass"]   = 141634 * int(temp_dict[f"Booster Number"]) # Manual
        temp_dict[f"Stage1 Stage Mass"]         = 12900
        temp_dict[f"Stage2 Stage Mass"]         = 3100

    if name in ["Chollima-1","NK Kerolox LV"]:
        # Using information from Long March (CZ) 2C/YZ-1S
        # NK Kerolox LV is based on Chollima-1.
        temp_dict[f"Stage1 Propellant Mass"]    = 162706
        temp_dict[f"Stage2 Propellant Mass"]    = 54667
        temp_dict[f"Stage3 Propellant Mass"]    = 1350
        temp_dict[f"Stage1 Stage Mass"]         = 10000
        temp_dict[f"Stage2 Stage Mass"]         = 4000
        temp_dict[f"Stage3 Stage Mass"]         = 2500

    if name == "Gravity-1":
        #No other information available for this vehicle.
        #Used Information Similar to Kuaizhou-11, which uses CZ-6.
        temp_dict[f"Stage1 Propellant Mass"]    = 76000
        temp_dict[f"Stage2 Propellant Mass"]    = (15000+15150) / 2
        temp_dict[f"Stage3 Propellant Mass"]    = 480 
        temp_dict[f"Stage1 Stage Mass"]         = 7530
        temp_dict[f"Stage2 Stage Mass"]         = 1490
        temp_dict[f"Stage3 Stage Mass"]         = 520
        temp_dict[f"Fairing Mass"]              = 1500

    if name == "Long March (CZ) 12":
        # No other information available for this vehicle.
        #Used Information Similar to Falcon 9 v1.2
        temp_dict[f"Stage1 Propellant Mass"]    = (418700+411000) / 2
        temp_dict[f"Stage2 Propellant Mass"]    = (111500+107500) / 2
        temp_dict[f"Stage1 Stage Mass"]         = (27200+22200) / 2
        temp_dict[f"Stage2 Stage Mass"]         = (4500+4000) / 2
        temp_dict[f"Fairing Mass"]              = (2000+1750) / 2 # Average of SLR/Sp101

    if name == "Long March (CZ) 6C":
        # No information available for this vehicle.
        #Used Information Similar to Long March (CZ) 6
        temp_dict[f"Stage1 Propellant Mass"]    = 76000
        temp_dict[f"Stage2 Propellant Mass"]    = (15000+15150) / 2 
        temp_dict[f"Stage1 Stage Mass"]         = 7530
        temp_dict[f"Stage2 Stage Mass"]         = 1490
        
    if name == "Vulcan Centaur VC2S":
        # Using stage 2 details from DISCOSweb.
        # https://web.archive.org/web/20180918143456/http://www.northropgrumman.com/Capabilities/GEM/Documents/GEM_63_GEM_63XL.pdf
        temp_dict[f"Booster Number"] = 2 
        temp_dict[f"Stage0 Propellant Mass"] = 48000 * int(temp_dict[f"Booster Number"]) # Manual
        temp_dict[f"Stage0 Stage Mass"]      = 5400 * int(temp_dict[f"Booster Number"]) # Manual
        temp_dict[f"Stage0 Fuel Type"] = 'Solid'
        temp_dict[f"Stage1 Propellant Mass"]    = 368000 # SLR
        temp_dict[f"Stage1 Stage Mass"]         = 400000-368000 # SLR
    
    if name == "Qaem-100 ":
        #Information Similar to Qased
        temp_dict[f"Stage1 Propellant Mass"] = 14400 # SLR
        temp_dict[f"Stage2 Propellant Mass"] = 1500  # NB 
        temp_dict[f"Stage3 Propellant Mass"] = 398   # .
        temp_dict[f"Stage1 Stage Mass"]      = 2700  # SLR
        temp_dict[f"Stage2 Stage Mass"]      = 200   # NB
        temp_dict[f"Stage3 Stage Mass"]      = (30+55+315)/ 3 # 
        temp_dict[f"Fairing Mass"]           = 100 # SLR

    if name == "Spectrum":
        pass 
        # This is 2025-F04. Rocket failed early into flight and crashed into the sea, so ignoring.
    
    if temp_dict["Fairing Mass"] == 0 and name != "Atlas V N22 v2020":
        #print(f"Setting fairing mass to average (1756) for {name}")
        temp_dict["Fairing Mass"] = 1756
       
    temp_dict["Fairing Mass"] = round(temp_dict["Fairing Mass"],1)
    
    return temp_dict