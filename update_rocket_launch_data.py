import sys
from pprint import pprint

def update_mass_info(temp_dict,name,variant):
    '''   
        All launch and dry stage mass information is initially compiled using GCAT.
        This file updates the propellant, stage, and fairing mass for each rocket. 
        Priority for information is as follows:
            - Rocket manuals / Mission guides 
            - Space Launch Report / Spaceflight101 / Jonathan's Space Report GCAT.
            - DISCOSweb
    '''
    
    # We could rearrange all of this to get the rocket information directly from GCAT.
    #   Pros: GCAT labels stage number and we could quickly and easily query mass.
    #   Cons: GCAT also has missing information, and would need to be cross-referenced with other sources anyway.

    # TODO: Perhaps this should be a YAML/TOML field, where there is a notes section.
     
    if name in ["Angara A5","Angara A5/Persei","Angara A5/Orion"]:
        # https://web.archive.org/web/20220406013831/http://www.spacelaunchreport.com/angara.html
        # https://discosweb.esoc.esa.int/launch-vehicles/107449
        # https://discosweb.esoc.esa.int/launch-vehicles/130
        # https://web.archive.org/web/20240730160859/https://spaceflight101.com/spacerockets/angara-a5/
        # https://web.archive.org/web/20221230022503/https://www.mach5lowdown.com/wp-content/uploads/PUG/Angara-Mission-Planners-Guide-Rev-0-2002-12.pdf
        # GCAT masses are identical to the manual for S0-S2.
        
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

        # Not much information on Orion, but it seems similar to Persei, based on 11S861-03 Phase I Version 2 upgrades to Blok DM-03.
        if name in ["Angara A5/Persei","Angara A5/Orion"]:
            temp_dict[f"Stage3 Propellant Mass"]  = 18700 # SLR
            temp_dict[f"Stage3 Stage Mass"]       = 2900  # SLR
            temp_dict[f"Fairing Mass"]            = 1600  # No info, approximating using Briz-M SLR info.
        elif name == "Angara A5": # Briz-M upper stage
            temp_dict["Stage3 Propellant Mass"]   = 19800 # SLR
            temp_dict[f"Stage3 Stage Mass"]       = 2370  # SLR
            temp_dict[f"Fairing Mass"]            = 1600  # Sp101
            
    if name in ["Angara-1.2","Angara 1.2PP"]:
        # Stage 1 is same as Angara A5, but no boosters.
        # Different stage 2.
        # Manual: https://web.archive.org/web/20230424155119/https://www.mach5lowdown.com/wp-content/uploads/PUG/Angara-Mission-Planners-Guide-Rev-0-2002-12.pdf
        temp_dict[f"Stage1 Propellant Mass"]  = 128800 # Manual
        temp_dict[f"Stage1 Stage Mass"]       = (9800+10500) / 2  # SLR/Sp101 (GCAT is 9800)
        temp_dict[f"Stage2 Propellant Mass"]  = 25700  # Manual
        temp_dict[f"Stage2 Stage Mass"]       = 2355   # Sp101 (GCAT and SLR are 4000)

        # https://www.russianspaceweb.com/angara1.html There is an 'Aggregate Module' upper stage according to this.
        # No info on GCAT, SLR or Sp101. Setting using DW info.
        temp_dict[f"Stage3 Propellant Mass"]  = 2000
        temp_dict[f"Stage3 Stage Mass"]       = 1000
        temp_dict[f"Fairing Mass"]            = (500+810)/2 # Two options on Sp101.
        
    if name.startswith("Ariane 5ECA"):
        
        # Manual is for ESC-A upper stage, whereas launches since 2019 have used ESC-D, an improved version.
        # https://web.archive.org/web/20240924193307/https://www.arianespace.com/wp-content/uploads/2016/10/Ariane5-users-manual-Jun2020.pdf
        # Using SLR for ESC-D propellant mass.
        temp_dict[f"Stage0 Propellant Mass"]     = 240000 * int(temp_dict[f"Booster Number"]) # Manual
        temp_dict[f"Stage0 Stage Mass"]          = 31000     # SLR/Sp101 (GCAT is 40000)
        temp_dict[f"Stage1 Propellant Mass"]     = 170000    # Manual 
        temp_dict[f"Stage1 Stage Mass"]          = 14700     # Manual
        temp_dict[f"Stage2 Stage Mass"]          = 4540      # Manual / SLR

        if name == "Ariane 5ECA":
            temp_dict[f"Stage2 Propellant Mass"] = 14900     # Manual
        elif name == "Ariane 5ECA+":
            temp_dict[f"Stage2 Propellant Mass"] = 15260     # SLR for ESC-D
        temp_dict[f"Fairing Mass"]               = 2675      # Manual        
    
    if name.startswith("Ariane 6") and variant != "1":
        # The '-' variants use the original P120C boosters from Vega C.
        # Only the latest '1' variant uses the updated P160C boosters.
        # Updated information from ESA
        # https://www.esa.int/Enabling_Support/Space_Transportation/Europe_s_Spaceport/Europe_s_next_solid_propellant_rocket_motor_passes_review

        if variant == '-':
            temp_dict[f"Stage0 Stage Mass"]        = 11000  * int(temp_dict[f"Booster Number"]) # ESA
            temp_dict[f"Stage0 Propellant Mass"]   = 142000 * int(temp_dict[f"Booster Number"]) # ESA
        elif variant == '1':
            temp_dict[f"Stage0 Stage Mass"]        = 11000  * int(temp_dict[f"Booster Number"]) # ESA
            temp_dict[f"Stage0 Propellant Mass"]   = 156000 * int(temp_dict[f"Booster Number"]) # ESA
    
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
        
        temp_dict[f"Stage1 Propellant Mass"]    = 284089 # https://www.ulalaunch.com/docs/default-source/rockets/atlasvusersguide2010a.pdf?sfvrsn=f84bb59e_2
        # The 400 series.
        if name[8] == "4":    
            temp_dict[f"Stage1 Stage Mass"] = 21054  # https://www.ulalaunch.com/docs/default-source/rockets/atlasvusersguide2010a.pdf?sfvrsn=f84bb59e_2
            temp_dict[f"Fairing Mass"]      = 2305   # https://www.ulalaunch.com/docs/default-source/rockets/atlasvusersguide2010a.pdf?sfvrsn=f84bb59e_2
        # The 500 series.
        elif name[8] == "5":  
            temp_dict[f"Stage1 Stage Mass"] = 21351  # https://www.ulalaunch.com/docs/default-source/rockets/atlasvusersguide2010a.pdf?sfvrsn=f84bb59e_2  
            temp_dict[f"Fairing Mass"]      = 3524   # https://www.ulalaunch.com/docs/default-source/rockets/atlasvusersguide2010a.pdf?sfvrsn=f84bb59e_2 
        # The N Series, no fairing and uses DEC not SEC. No info in manual so using GCAT.
        else:
            pass
            
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

        temp_dict[f"Stage0 Propellant Mass"]  = 203000 * int(temp_dict[f"Booster Number"]) # SLR/Sp101
        temp_dict[f"Stage1 Propellant Mass"]  = 203000              # SLR/Sp101
        temp_dict[f"Stage2 Propellant Mass"]  = 27200               # Manual  
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
            temp_dict[f"Stage1 Stage Mass"]       = 1969 # assuming same 1st stage mass as predecessor
            temp_dict[f"Stage2 Stage Mass"]       = 780
            temp_dict[f"Stage3 Stage Mass"]       = 100

    if name == "Electron":
        
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

        # 3rd stage is 'kick-stage', and has a 'liquid bi-propellant' which wikipedia says is hypergolic.
        # Kick-stage is placed into low elliptical orbit, and burns some propellant to move into final orbit.
        # Seems like it is specific to mission, so would have to take an educated guess to the average prop.

        temp_dict[f"Stage1 Propellant Mass"]     = 9250             # SLR/Sp101
        temp_dict[f"Stage2 Propellant Mass"]     = (2050+2150) / 2  # SLR/Sp101
        temp_dict[f"Stage1 Stage Mass"]          = 950              # SLR/Sp101
        temp_dict[f"Stage2 Stage Mass"]          = 250              # SLR/Sp101
        temp_dict[f"Fairing Mass"]               = 44               # https://rocketlabcorp.com/assets/Rocket-Lab-Electron-Payload-User-Guide-8.0.pdf

        if variant != "2":
            temp_dict[f"Stage3 Propellant Mass"] = 245              # Using value from DISCOSweb(wet/dry)
            temp_dict[f"Stage3 Stage Mass"]      = 40               # https://rocketlabcorp.com/assets/Rocket-Lab-Electron-Payload-User-Guide-8.0.pdf
      
    if name == "Epsilon":
        # Manual: https://global.jaxa.jp/projects/rockets/epsilon/pdf/EpsilonUsersManual_e.pdf
        # Manual for newer Epsilon S: https://www.ihi.co.jp/en/all_news/pdf/ip01.pdf
        # Manual has details for upper stage, but says total and prop mass are 1.0 tons (so zero dry mass).
        # Filling in with SLR.

        # Same first stage for all variants.
        temp_dict[f"Stage1 Propellant Mass"] = 66300 # Manual
        temp_dict[f"Stage1 Stage Mass"]      = 8700  # Manual

        # 2nd Stage M-34c for original Epsilon, M35 for Epsilon 2 (Enhanced Epsilon).
        # 3rd Stage KM-V2b for original Epsilon, KM-V2c for Epsilon 2 (Enhanced Epsilon).
        if variant in ["2", "2PBS"]:
            temp_dict[f"Stage2 Propellant Mass"] = 15000 # Manual 
            temp_dict[f"Stage2 Stage Mass"]      = 2000  # Manual 
            temp_dict[f"Stage3 Propellant Mass"] = 2500  # Manual
            temp_dict[f"Stage3 Stage Mass"]      = 800   # Manual
        else:
            temp_dict[f"Stage2 Propellant Mass"] = 10800 # SLR
            temp_dict[f"Stage2 Stage Mass"]      = 800   # SLR
            temp_dict[f"Stage3 Propellant Mass"] = 2500  # SLR
            temp_dict[f"Stage3 Stage Mass"]      = 500   # SLR

        if variant in ["2PBS",""]:
            temp_dict[f"Stage4 Propellant Mass"] = 145   # SLR         
            temp_dict[f"Stage4 Stage Mass"]      = 155   # SLR
        
        temp_dict[f"Fairing Mass"]               = 1000  # Manual
        
    if name == "Falcon Heavy":
    
        # Prop masses:
        # DISCOSweb: 294300 (???), 395700, 108670
        # SLR      : 407600 (one), 407600, 107200
        # Sp101    : 411000 (one), 411000, 107500
                
        # Stage masses:
        # DISCOSweb: 33100 (???), 25600, 3230
        # SLR      : 17000 (one), 17000, 4500
        # Sp101    : 22500 (one), 25600, 4000
        
        # Manual: https://www.spacex.com/assets/media/falcon-users-guide-2025-05-09.pdf
        # "Falcon Heavy’s first-stage comprises three Falcon 9 first stages with enhancements provided to strengthen the cores. 
        # Furthermore, Falcon Heavy utilizes the same second stage and same payload fairing as flown on Falcon 9..."
        
        # Fairing mass:
        #   - https://spaceflight101.com/spacerockets/falcon-9-ft/                                          1750 kg
        #   - https://web.archive.org/web/20220411235211/http://www.spacelaunchreport.com/falcon9ft.html    2000 kg

        temp_dict[f"Stage0 Propellant Mass"]  = (407600+411000) / 2 * 2 # SLR/Sp101
        temp_dict[f"Stage0 Stage Mass"]       = (17000+22500) / 2 * 2   # SLR/Sp101
        temp_dict[f"Stage1 Propellant Mass"]  = (407600+411000) / 2     # SLR/Sp101
        temp_dict[f"Stage2 Propellant Mass"]  = (107200+107500) / 2     # SLR/Sp101
        temp_dict[f"Stage1 Stage Mass"]       = (17000+25600) / 2       # SLR/Sp101
        temp_dict[f"Stage2 Stage Mass"]       = (4500+4000) / 2         # SLR/Sp101
        temp_dict[f"Fairing Mass"]            = (2000+1750) / 2         # SLR/Sp101
        
    if name == "Falcon 9" and variant == "FT5":
        
        # Prop masses:
        # SLR      : 418700, 111500
        # Sp101    : 411000, 107500
                
        # Stage masses:
        # SLR      : 27200, 4500
        # Sp101    : 22200, 4000

        # Fairing mass:
        #   - https://spaceflight101.com/spacerockets/falcon-9-ft/                                          1750 kg
        #   - https://web.archive.org/web/20220411235211/http://www.spacelaunchreport.com/falcon9ft.html    2000 kg
        
        temp_dict[f"Stage1 Propellant Mass"] = (418700+411000) / 2 # SLR/Sp101
        temp_dict[f"Stage2 Propellant Mass"] = (111500+107500) / 2 # SLR/Sp101
        temp_dict[f"Stage1 Stage Mass"]      = (27200+22200) / 2   # SLR/Sp101
        temp_dict[f"Stage2 Stage Mass"]      = (4500+4000) / 2     # SLR/Sp101
        temp_dict[f"Fairing Mass"]           = (2000+1750) / 2     # SLR/Sp101
    
    if name == "GSLV Mk II":

        # Prop masses:
        # SLR      : 42600, 138200, 39500, 14996 (Mk2+ CUS-15)
        # Sp101    : 42000, 138000, 39400, 12800 (CUS15)
        
        # Stage masses:
        # SLR      : 5600, 28300, 5500, 2601  
        # Sp101    : 5600, 28300, 5500, 2500 (CUS15)

        temp_dict[f"Stage0 Propellant Mass"]  = (42600+42000) / 2 * int(temp_dict[f"Booster Number"])   # SLR/Sp101
        temp_dict[f"Stage0 Stage Mass"]       = 5600 * int(temp_dict[f"Booster Number"])                # SLR/Sp101
        
        temp_dict[f"Stage1 Propellant Mass"]  = (138200+138000) / 2                                     # SLR/Sp101
        temp_dict[f"Stage1 Stage Mass"]       = 28300                                                   # SLR/Sp101

        if variant == "-": # GS2 second stage
            temp_dict[f"Stage2 Propellant Mass"]  = (39500+39400) / 2                                   # SLR/Sp101
            temp_dict[f"Stage2 Stage Mass"]       = 5500                                                # SLR/Sp101
        else: # GS2-HTVE second stage
            # No info, using GCAT.
            pass

        if variant in ["-","a"]: # CUS-12 upper stage
            temp_dict[f"Stage3 Propellant Mass"]  = 12800  # SLR/Sp101
            temp_dict[f"Stage3 Stage Mass"]       = 2500   # SLR/Sp101
        else: # CUS-15 upper stage
            temp_dict[f"Stage3 Propellant Mass"]  = 14996  # SLR
            temp_dict[f"Stage3 Stage Mass"]       = 2601   # SLR
        
    if name in ["GSLV Mk III","LVM3"]:        
        # Prop masses:
        # DISCOSweb: 205000 (one), 115400 (W-D), 27700 (W-D)
        # SLR      : 205000 (one), 116000,       28600
        # Sp101    : 206690 (one), 115000,       28000
        
        # Stage masses:
        # DISCOSweb: 31000 (one), 9600,  5300
        # SLR      : 31000 (one), 9600,  4700
        # Sp101    : 31300 (one), 10600, 5000
        
        temp_dict[f"Stage0 Propellant Mass"] = (206690+205000)
        temp_dict[f"Stage1 Propellant Mass"] = 115500
        temp_dict[f"Stage2 Propellant Mass"] = (28600+28000) / 2
        temp_dict[f"Stage0 Stage Mass"]      = (31300+31000)
        temp_dict[f"Stage1 Stage Mass"]      = (9600+10600) / 2
        temp_dict[f"Stage2 Stage Mass"]      = (4700+5000) / 2
    
    if name.startswith("H-IIA"):
        # Details from https://global.jaxa.jp/projects/rockets/h2a/ 
        temp_dict[f"Stage1 Propellant Mass"] = 101100
        temp_dict[f"Stage2 Propellant Mass"] = 16900 
        temp_dict[f"Stage1 Stage Mass"]      = 12900 
        temp_dict[f"Stage2 Stage Mass"]      = 3100  
        temp_dict[f"Fairing Mass"]           = 1400 # website/Sp101/SLR
        
        temp_dict[f"Booster Number"] = name[8]
        
        temp_dict[f"Stage0 Propellant Mass"] = 130000 / 2 * int(temp_dict[f"Booster Number"]) 
        temp_dict[f"Stage0 Stage Mass"]      = 21000 / 2 * int(temp_dict[f"Booster Number"])
        if len(name) == 10:
            temp_dict[f"Stage-1 Propellant Mass"] = 52400 / 4 * int(name[9])
            temp_dict[f"Stage-1 Stage Mass"]      = (31/2-52.4/4) * 1000 * int(name[9])
        
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
    
    if name in ["Juno II", "Jupiter C"]:

        # Jupiter C is another name for Juno I, with a different first stage. http://www.astronautix.com/j/jupitercstage.html
        # http://www.astronautix.com/j/junoii.html
        # https://en.wikipedia.org/wiki/Juno_II#Specifications 
        # Wiki assumes the masses by using the fact the 4th stage is one baby seargent engine, and stages 2-3 are multiples of this.
        # http://www.astronautix.com/c/castorengine.html This is one baby seargent engine.

        if name == "Juno II": 
            temp_dict[f"Stage1 Propellant Mass"] = 48988 # http://www.astronautix.com/j/jupiterstage.html
            temp_dict[f"Stage1 Stage Mass"]      = 5443  # http://www.astronautix.com/j/jupiterstage.html
        elif name == "Jupiter C": # Redstone first stage
            temp_dict[f"Stage1 Propellant Mass"] = 24540 # http://www.astronautix.com/j/jupitercstage.html
            temp_dict[f"Stage1 Stage Mass"]      = 3890  # http://www.astronautix.com/j/jupitercstage.html

        temp_dict[f"Stage2 Propellant Mass"] = 231 # 11 Baby Sergeant
        temp_dict[f"Stage2 Stage Mass"]      = 231 # 11 Baby Sergeant
        temp_dict[f"Stage3 Propellant Mass"] = 63  # 3 Baby Sergeant
        temp_dict[f"Stage3 Stage Mass"]      = 63  # 3 Baby Sergeant

        if (name == "Juno II" and variant == "-") or (name == "Jupiter C" and variant == "JunoI"):
            temp_dict[f"Stage4 Propellant Mass"] = 21 # 1 Baby Sergeant
            temp_dict[f"Stage4 Stage Mass"]      = 21 # 1 Baby Sergeant
    
    if name == "LauncherOne":
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

    if name.startswith("Chang Zheng 2C"):

        # Chang Zheng 2C 
        temp_dict[f"Stage1 Propellant Mass"] = 162706  # http://cgwic.com/Launchservice/LM2C.html
        temp_dict[f"Stage2 Propellant Mass"] = 54667   # http://cgwic.com/Launchservice/LM2C.html
        temp_dict[f"Stage1 Stage Mass"]      = 10000   # SLR
        temp_dict[f"Stage2 Stage Mass"]      = 4000    # SLR
        temp_dict[f"Fairing Mass"]           = 800     # http://www.b14643.de/Spacerockets_1/China/CZ-2C/Description/Frame.htm

        if name.endswith("YZ-1S"):
            # No info available from GCAT, SLR or Sp101. DISCOSweb dry mass matches JSR reentry data, so keeping DW values.
            temp_dict[f"Stage3 Propellant Mass"] = 1350 # DW
            temp_dict[f"Stage3 Stage Mass"]      = 2500 # DW
        
        # GCAT has data for the SM/SMA stages now. 
    
    if name.startswith("Chang Zheng 2D"):
        # Chang Zheng 2D 
        temp_dict[f"Stage1 Stage Mass"] = 9500    # SLR
        temp_dict[f"Stage2 Stage Mass"] = 4000    # SLR
        temp_dict[f"Stage1 Propellant Mass"] = 182000  # http://cgwic.com/Launchservice/LM2D.html
        temp_dict[f"Stage2 Propellant Mass"] = 52700   # http://cgwic.com/Launchservice/LM2D.html
        temp_dict[f"Fairing Mass"]           = 750     # http://www.b14643.de/Spacerockets_1/China/CZ-2D/Description/Frame.htm

        if name.endswith("YZ-3"):
            # No info from GCAT, SLR or Sp101. Basing this off YZ-1S
            temp_dict[f"Stage3 Propellant Mass"]   = 1350
            temp_dict[f"Stage3 Stage Mass"]        = 2500
    
    if name == "Chang Zheng 2F":
        
        # Prop masses:
        # DISCOSweb : 41501 (one), 186310, 86150
        # SLR       : 37800 (one), 187000, 86000
        # Sp101     : 40800 (one), 186500, 91500
                
        # Stage masses:
        # DISCOSweb : 3200 (one), 1100,  5500
        # SLR       : 3200 (one), 9500,  5500
        # Sp101     : 3000 (one), 12500, 5500
        temp_dict[f"Stage0 Propellant Mass"]  = (37800+40800) / 2 * 4 # SLR/Sp101
        temp_dict[f"Stage0 Stage Mass"]       = (3000+3200) / 2 * 4   # SLR/Sp101
        temp_dict[f"Stage1 Stage Mass"]       = (9500+12500) / 2      # SLR/Sp101
        temp_dict[f"Stage1 Propellant Mass"]  = (187000+186500) / 2   # SLR/Sp101
        temp_dict[f"Stage2 Propellant Mass"]  = (86000+91500) / 2     # SLR/Sp101 
        temp_dict[f"Fairing Mass"]            = 6000                  # http://www.b14643.de/Spacerockets_1/China/CZ-2EF/Description/Frame.htm
        
    if name.startswith("Chang Zheng 3A") or name.startswith("Chang Zheng 3B") or name.startswith("Chang Zheng 3C"):
        
        # These are for the original version, the upgraded version has updated boosters and first stage.
        # Manual: https://web.archive.org/web/20221230022505/https://www.mach5lowdown.com/wp-content/uploads/PUG/LM-3B-User-Manual-v1999.pdf
        # Manual: https://web.archive.org/web/20221230022531/https://www.mach5lowdown.com/wp-content/uploads/PUG/LM-3C-User-Manual-v1998.pdf

        if variant.startswith("G"): # Updated boosters and first stage for the upgraded version.
            temp_dict[f"Stage0 Propellant Mass"]  = 41100 * int(temp_dict[f"Booster Number"])           # SLR/Sp101
            temp_dict[f"Stage1 Propellant Mass"]  = 186200                                              # SLR/Sp101
            temp_dict[f"Stage0 Stage Mass"]       = (3900+3200) / 2 * int(temp_dict[f"Booster Number"]) # SLR/Sp101
            temp_dict[f"Stage1 Stage Mass"]       = (9800+12600)/ 2                                     # SLR/Sp101
        else: # Original version.
            temp_dict[f"Stage0 Propellant Mass"]  = 37746 * int(temp_dict[f"Booster Number"])           # Manual
            temp_dict[f"Stage1 Propellant Mass"]  = 171775                                              # Manual
            temp_dict[f"Stage0 Stage Mass"]       = (3000+3330) / 2 * int(temp_dict[f"Booster Number"]) # SLR/Sp101
            temp_dict[f"Stage1 Stage Mass"]       = (12125+9200)/ 2                                     # SLR/Sp101

        if name == "Chang Zheng 3A":
            temp_dict[f"Stage0 Propellant Mass"]  = 0
            temp_dict[f"Stage0 Stage Mass"]       = 0
        else: # 3B and 3C have updated second stage, but 3A, 3B, and 3C share the same first and third stage.
            temp_dict[f"Stage2 Propellant Mass"]  = 49605                                               # Manual
            temp_dict[f"Stage2 Stage Mass"]       = (3848+4000) / 2                                     # SLR/Sp101
        
        temp_dict[f"Stage3 Propellant Mass"]  = 18193                                                   # Manual
        temp_dict[f"Stage3 Stage Mass"]       = (2740+2800) / 2                                         # SLR/Sp101
        temp_dict[f"Fairing Mass"]            = 1500                                                    # Sp101

        if name.endswith("YZ-1"):
            # No info available from GCAT, SLR or Sp101. Using DW.
            temp_dict[f"Stage4 Stage Mass"]       = 5000  # DW
            temp_dict[f"Stage4 Propellant Mass"]  = 50000 # DW
            temp_dict[f"Stage4 Fuel Type"]        = 'Hypergolic'
        
    if name in ["Chang Zheng 4B","Chang Zheng 4C"]: 
        
        temp_dict[f"Stage2 Propellant Mass"]  = 35550   # SLR
        temp_dict[f"Stage3 Propellant Mass"]  = 14000   # SLR
        temp_dict[f"Stage3 Stage Mass"]       = 2000    # SLR
        temp_dict[f"Fairing Mass"]            = 1800    # http://www.b14643.de/Spacerockets_1/China/CZ-4/Description/Frame.htm
     
    if name.startswith("Chang Zheng 5"):  
        
        # Prop masses:
        # DISCOSweb : 174000, 158500, 22900
        # SLR       : 152000, 158000, 22900
        # Sp101     : 135000, 158500, 26500
                
        # Stage masses:
        # DISCOSweb : 12000, 18000, 3400
        # SLR       : 12500, 17000, 3520
        # Sp101     : 12000, 18000, 3400

        temp_dict[f"Stage0 Stage Mass"]       = (12500+12000) / 2 * 4  
        temp_dict[f"Stage1 Stage Mass"]       = (17000+18000) / 2 
        temp_dict[f"Stage0 Propellant Mass"]  = (152000+135000) / 2 * 4 
        temp_dict[f"Stage1 Propellant Mass"]  = (158000+158500) / 2  

        # No info for YZ-2 in GCAT, SLR or Sp101. Using DW values.  
        if "5B" in name:
            temp_dict[f"Fairing Mass"]  = 2400       # http://www.b14643.de/Spacerockets_1/China/CZ-5/Description/Frame.htm
            if name.endswith("YZ-2"): 
                # No info available from GCAT, SLR or Sp101. Using DW.
                temp_dict[f"Stage2 Propellant Mass"]  = 7000 # DW
                temp_dict[f"Stage2 Stage Mass"]       = 5000 # DW  
        else:
            temp_dict[f"Stage2 Propellant Mass"]  = (22900+26500) / 2
            temp_dict[f"Stage2 Stage Mass"]       = (3520+3400) / 2  
            if name.endswith("YZ-2"): 
                temp_dict[f"Stage3 Propellant Mass"]  = 7000 # DW
                temp_dict[f"Stage3 Stage Mass"]       = 5000 # DW   
            temp_dict[f"Fairing Mass"]  = 3000 # http://www.b14643.de/Spacerockets_1/China/CZ-5/Description/Frame.htm
        
    if name == "Chang Zheng 6":     
        # Prop masses:
        # SLR       : 76000, 15000, N/A
        # Sp101     : 76000, 15150, N/A
                
        # Stage masses:
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
             
    if name.startswith("Chang Zheng 7"):
        
        # Prop masses:
        # SLR      : 75500 (each), 174000, 65000                 , 18200 (H18)
        # Sp101    : 75000 (each), 174000, 65000                 , N/A   (H18)
        
        # Stage masses:
        # SLR      : 6000 (each), 12500, 5500                , 2800 (H18)   
        # Sp101    : 6000 (each), 12000, 6000                , N/A  (H18)
        
        temp_dict[f"Stage0 Propellant Mass"]  = (75500+75000) / 2 * int(temp_dict[f"Booster Number"]) # SLR/Sp101
        temp_dict[f"Stage1 Propellant Mass"]  = 174000 
        temp_dict[f"Stage2 Propellant Mass"]  = 65000  # SLR/Sp101

        temp_dict[f"Stage0 Stage Mass"]       = 6000 * int(temp_dict[f"Booster Number"]) # SLR/Sp101
        temp_dict[f"Stage1 Stage Mass"]       = (12000+12500) / 2 # SLR/Sp101
        temp_dict[f"Stage2 Stage Mass"]       = (5500+6000) / 2   # SLR/Sp101
        
        if name.endswith("YZ-1A"):
            # No info available from GCAT, SLR or Sp101. Using DW.
            temp_dict["Stage3 Stage Mass"]        = 3100         # DW for HO stage
            temp_dict["Stage3 Propellant Mass"]   = 26000 - 3100 # DW for HO stage
            temp_dict[f"Fairing Mass"]            = 2800 # http://www.b14643.de/Spacerockets_1/China/CZ-7/Description/Frame.htm
        elif name == "Chang Zheng 7":
            temp_dict[f"Fairing Mass"]            = 2800 # http://www.b14643.de/Spacerockets_1/China/CZ-7/Description/Frame.htm
        elif name == "Chang Zheng 7A":
            temp_dict[f"Stage3 Propellant Mass"]  = (18200+18193) / 2 # SLR/Sp101
            temp_dict[f"Stage3 Stage Mass"]       = (2740+2800)   / 2 # SLR/Sp101
            temp_dict[f"Fairing Mass"]            = 2500 # http://www.b14643.de/Spacerockets_1/China/CZ-7A/Description/Frame.htm
        
    if "Chang Zheng 8" in name:
        # https://web.archive.org/web/20220411235209/http://www.spacelaunchreport.com/cz5.html

        if variant == "-":
            temp_dict[f"Booster Number"] = 2
            temp_dict[f"Stage0 Propellant Mass"] = 75500 * 2  # SLR
            temp_dict[f"Stage0 Stage Mass"]      = 6000  * 2  # SLR

        temp_dict[f"Stage1 Propellant Mass"] = 174000         # SLR
        temp_dict[f"Stage2 Propellant Mass"] = 18200          # SLR
        temp_dict[f"Stage1 Stage Mass"]      = 12500          # SLR   
        temp_dict[f"Stage2 Stage Mass"]      = 2800           # SLR 
        temp_dict[f"Fairing Mass"] = 2500 # http://www.b14643.de/Spacerockets_1/China/CZ-8/Description/Frame.htm        

        # https://spacenews.com/first-launch-of-long-march-8a-sends-second-group-of-guowang-megaconstellation-satellites-into-orbit/
        # The Chang Zheng 8A is an upgraded variant of the standard Chang Zheng 8, which debuted in December 2020. 
        # It features the same first stage and side boosters as the original but includes a newly designed 
        # 3.35-meter-diameter hydrogen-oxygen second stage, allowing a wider, 5.2-meter-diameter payload fairing.      
    
    if name == "Chang Zheng 12":
        
        # No other information available for 12, 12A or 12B.
        # Approximating using Falcon 9 FT5 (see find_proxy.ipynb).
        temp_dict[f"Stage1 Propellant Mass"] = (418700+411000) / 2 # SLR/Sp101
        temp_dict[f"Stage2 Propellant Mass"] = (111500+107500) / 2 # SLR/Sp101
        temp_dict[f"Stage1 Stage Mass"]      = (27200+22200) / 2   # SLR/Sp101
        temp_dict[f"Stage2 Stage Mass"]      = (4500+4000) / 2     # SLR/Sp101
        temp_dict[f"Fairing Mass"]           = (2000+1750) / 2     # SLR/Sp101
        temp_dict["proxy"]                   = "Falcon 9 FT5"

    if name == "Chang Zheng 12A":
        
        # No other information available for 12, 12A or 12B.
        # Approximating using Zhuque-2E (see find_proxy.ipynb).
        temp_dict[f"Stage1 Propellant Mass"] = 135000.0 # GCAT
        temp_dict[f"Stage2 Propellant Mass"] = 45000.0  # GCAT
        temp_dict[f"Stage1 Stage Mass"]      = 15000.0  # GCAT
        temp_dict[f"Stage2 Stage Mass"]      = 5000.0   # GCAT
        temp_dict[f"Fairing Mass"]           = 1756     # GCAT
        temp_dict["proxy"]                   = "Zhuque-2E"
        
    if name in ["Minotaur I"]:
        # Prop masses:
        # SLR      : 20785,       6237,      3915, 782
        # Sp101    : 20785 (W-D), 6237(W-D), 3915, 770.2
        
        # Stage masses:
        # SLR      : 2292, 795, 416, 203 
        # Sp101    : 2292, 795, 391, 102.1
        
        temp_dict[f"Stage1 Propellant Mass"] = 20785
        temp_dict[f"Stage2 Propellant Mass"] = 6237
        temp_dict[f"Stage3 Propellant Mass"] = 3915
        temp_dict[f"Stage1 Stage Mass"]      = 2292
        temp_dict[f"Stage2 Stage Mass"]      = 795
        temp_dict[f"Stage3 Stage Mass"]      = (416+391) / 2
        temp_dict[f"Stage4 Propellant Mass"] = (782+770.2) / 2
        temp_dict[f"Stage4 Stage Mass"]      = (203+102.1) / 2
        temp_dict[f"Fairing Mass"]           = 300
        
    if name.startswith("Minotaur IV"):
        # Prop masses:
        # SLR       : 45370, 24490, 7070, 770   (Orion), 2009 (Star 48BV)
        # Sp101     : 45400, 24500, 7080, 770.2 (Orion), 2010 (Star 48BV)
                
        # Stage masses:
        # SLR       : 3620, 3180, 640,  410   (Orion), 155   (Star 48BV)
        # Sp101     : 3600, 3200, 650,  102.1 (Orion), 154.5 (Star 48BV)

        # Fairing mass:
        #   - https://spaceflight101.com/spacerockets/minotaur-iv/                                          450 kg
        #   - https://web.archive.org/web/20220406013816/http://www.spacelaunchreport.com/mintaur4.html     400 kg
        
        temp_dict[f"Stage1 Propellant Mass"] = (45370+45400) / 2 
        temp_dict[f"Stage2 Propellant Mass"] = (24490+24500) / 2 
        temp_dict[f"Stage3 Propellant Mass"] = (7070+7080) / 2 
        temp_dict[f"Stage1 Stage Mass"]      = (3620+3600) / 2 
        temp_dict[f"Stage2 Stage Mass"]      = (3180+3200) / 2  
        temp_dict[f"Stage3 Stage Mass"]      = (640+650) / 2  
        
        temp_dict[f"Fairing Mass"] = (400+450) / 2

        if name == "Minotaur IV": # GCA/Orion 38 stage 4, also GCA/Orion 38 stage 5 for 2 variant.
            temp_dict[f"Stage4 Propellant Mass"] = 770 # SLR/Sp101
            temp_dict[f"Stage4 Stage Mass"]      = (410+102) / 2 # SLR/Sp101 
        elif name == "Minotaur IV+": # Star 48BV stage 4.
            temp_dict[f"Stage4 Propellant Mass"] = 2010 # SLR/Sp101
            temp_dict[f"Stage4 Stage Mass"]      = 154.5  # SLR/Sp101
        elif name == "Minotaur IV Lite": # HAPS stage 4.
            temp_dict[f"Stage4 Propellant Mass"] = 72   # Sp101
            temp_dict[f"Stage4 Stage Mass"]      = 0.74 # Sp101
     
    if name == "Nuri":  
        # No info anywhere else, so using DISCOSweb wet mass - dry mass.
        temp_dict[f"Stage1 Stage Mass"] = 16000   # DW
        temp_dict[f"Stage2 Stage Mass"] = 4400    # DW
        temp_dict[f"Stage3 Stage Mass"] = 1290    # DW
        temp_dict[f"Stage1 Propellant Mass"] = 147100-16000   # DW (W-D)
        temp_dict[f"Stage2 Propellant Mass"] = 41200-4400     # DW (W-D)
        temp_dict[f"Stage3 Propellant Mass"] = 12100-1290     # DW (W-D)
     
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
        
        if name in ["PSLV-XL","PSLV-QL","PSLV-DL"]:
            temp_dict[f"Stage4 Propellant Mass"]    = (2500+2000) / 2 # SLR/Sp101
        elif variant == "CA":
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
       
    if "Pegasus XL" in name:
        # Manual: https://web.archive.org/web/20221230022508/https://www.mach5lowdown.com/wp-content/uploads/PUG/pegasus-user-guide-2007-1.pdf
        # Prop masses:
        # SLR      : 15010, 3930, 770
        # Sp101    : 15014, 3915, 770.2
        # Manual   : 15032, 3923, 770
        
        # Stage masses:
        # SLR      : 1370, 420, 110
        # Sp101    : 1369, 391, 102.1, 
        # Manual   : 1386, 416, 108  
        
        temp_dict[f"Stage1 Propellant Mass"] = 15032
        temp_dict[f"Stage2 Propellant Mass"] = 3923
        temp_dict[f"Stage3 Propellant Mass"] = 770
        temp_dict[f"Stage1 Stage Mass"]      = 1386 
        temp_dict[f"Stage2 Stage Mass"]      = 416
        temp_dict[f"Stage3 Stage Mass"]      = 108 
    
    if "Proton-M" in name:
        # Fairing masses:
        # Sp101    :  2000 (Briz-M Proton), 1600 (Briz-M Angara A5), 2000 (DM-03)
        temp_dict[f"Fairing Mass"]           = 2000 

    if name == "Safir":
        # https://web.archive.org/web/20220411235210/http://www.spacelaunchreport.com/safir.html
        temp_dict[f"Stage1 Propellant Mass"] = 18200 # SLR
        temp_dict[f"Stage2 Propellant Mass"] = 3050  # SLR
        temp_dict[f"Stage1 Stage Mass"]      = 3500  # SLR
        temp_dict[f"Stage2 Stage Mass"]      = 350   # SLR
        temp_dict[f"Fairing Mass"]           = 100   # SLR

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

        if name == "Shavit":
            temp_dict[f"Stage1 Propellant Mass"] = 9100  # SLR
            temp_dict[f"Stage1 Stage Mass"]      = 1115  # SLR
            temp_dict[f"Stage2 Propellant Mass"] = 9100  # SLR
            temp_dict[f"Stage2 Stage Mass"]      = 1288  # SLR
        elif name == "Shavit 1":
            temp_dict[f"Stage1 Propellant Mass"] = 12750 # SLR
            temp_dict[f"Stage1 Stage Mass"]      = 1240  # SLR
            temp_dict[f"Stage2 Propellant Mass"] = 9100  # SLR
            temp_dict[f"Stage2 Stage Mass"]      = 1288  # SLR
        elif name == "Shavit 2":
            temp_dict[f"Stage1 Propellant Mass"] = 12750 # SLR
            temp_dict[f"Stage1 Stage Mass"]      = 1240  # SLR
            temp_dict[f"Stage2 Propellant Mass"] = 12750 # SLR
            temp_dict[f"Stage2 Stage Mass"]      = 1376  # SLR

        temp_dict[f"Stage3 Propellant Mass"] = 1890      # SLR
        temp_dict[f"Stage3 Stage Mass"]      = 684       # SLR
        temp_dict[f"Fairing Mass"]           = 57        # http://www.b14643.de/Spacerockets_1/Rest_World/Shavit/Description/Frame.htm
        
    if name == "Simorgh":
        # https://web.archive.org/web/20220406013833/http://www.spacelaunchreport.com/simorgh.html
        # Saman-1 is an upper-stage, used as a space tug to transfer satellites to higher orbits.
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
    
    if name == "SLS Block 1":
        # Manual: https://www.nasa.gov/wp-content/uploads/2022/03/sls-reference-guide-2022-v2-print-0.pdf

        temp_dict[f"Stage0 Propellant Mass"]   = 626700 * int(temp_dict[f"Booster Number"]) # Manual
        temp_dict[f"Stage1 Propellant Mass"]   = 1002700 # Manual
        temp_dict[f"Stage2 Propellant Mass"]   = 29000   # Manual
        temp_dict[f"Stage0 Stage Mass"]        = 99300 * int(temp_dict[f"Booster Number"]) # Manual
        temp_dict[f"Stage1 Stage Mass"]        = 99300 # Manual
        temp_dict[f"Stage2 Stage Mass"]        = 3700 # Manual
        temp_dict[f"Fairing Mass"]             = 8165 # SLR

    if name == "Strela":
        # https://web.archive.org/web/20230614172003/https://spaceflight101.com/spacerockets/strela/
        temp_dict[f"Stage1 Stage Mass"]      = 5700  # Sp101
        temp_dict[f"Stage2 Stage Mass"]      = 1500  # Sp101
        temp_dict[f"Stage3 Stage Mass"]      = 725   # Sp101
        temp_dict[f"Stage1 Propellant Mass"] = 71500 # Sp101
        temp_dict[f"Stage2 Propellant Mass"] = 10700 # Sp101
        temp_dict[f"Stage3 Propellant Mass"] = 375   # Sp101
    
    if name == "SSLV":
        # https://www.isro.gov.in/mission_SSLV_D2.html
        # https://www.isro.gov.in/media_isro/pdf/Missions/SSLV/SSLV_D2_EOS_07_DigitalBrochure.pdf
        
        temp_dict[f"Stage1 Propellant Mass"]    = 87000 # Mission guides
        temp_dict[f"Stage2 Propellant Mass"]    = 7700  # Mission guides
        temp_dict[f"Stage3 Propellant Mass"]    = 4500  # Mission guides
        temp_dict[f"Stage4 Propellant Mass"]    = 50    # Mission guides
        temp_dict[f"Fairing Mass"]              = 90 # NB
    
    if "Soyuz-2" in name or "Soyuz-ST" in name:
        
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
        if "2-1V" in name:
            temp_dict[f"Booster Number"] = 0
        else:
            temp_dict[f"Booster Number"] = 4
            temp_dict[f"Stage0 Propellant Mass"] = (27900 + 11260) * 4     # Soyuz-Users-Manual-March-2012.pdf
            temp_dict[f"Stage0 Stage Mass"]      = 3784 * 4                # Soyuz-Users-Manual-March-2012.pdf
        
        # Sort out first stage
        if name == "Soyuz-2-1V":
            temp_dict[f"Stage1 Propellant Mass"]  = 119000 # SLR
            temp_dict[f"Stage1 Stage Mass"]       = 9300   # Sp101
        else:
            temp_dict[f"Stage1 Propellant Mass"]  = 63800 + 26300       # Soyuz-Users-Manual-March-2012.pdf
            temp_dict[f"Stage1 Stage Mass"]       = 6545                # Soyuz-Users-Manual-March-2012.pdf
        
        # Sort out second stage
        if "2-1B" in name or "2-1V" in name:
            temp_dict[f"Stage2 Stage Mass"]       = 2355            # Soyuz-Users-Manual-March-2012.pdf
            temp_dict[f"Stage2 Propellant Mass"]  = 17800 + 7600    # Soyuz-Users-Manual-March-2012.pdf
        elif "2-1A" in name: 
            temp_dict[f"Stage2 Stage Mass"]       = 2410            # soyuz_users_manual_190401-1.pdf
            temp_dict[f"Stage2 Propellant Mass"]  = 22790           # soyuz_users_manual_190401-1.pdf
        elif "ST-A" in name: 
            temp_dict[f"Stage2 Stage Mass"]       = 2470            # soyuz_users_manual_190401-1.pdf
            temp_dict[f"Stage2 Propellant Mass"]  = 22830           # soyuz_users_manual_190401-1.pdf
            
        # Sort out third stage
        if variant == "Fregat":
            temp_dict[f"Stage3 Stage Mass"]       = 1000            # Manual and SLR
            temp_dict[f"Stage3 Propellant Mass"]  = 5350            # Manual and SLR
        elif variant == "Fregat-M":
            temp_dict[f"Stage3 Stage Mass"]       = 960             # Sp101
            temp_dict[f"Stage3 Propellant Mass"]  = 5750            # Sp101
        elif variant == "Fregat-MT":
            temp_dict[f"Stage3 Stage Mass"]       = 1050            # Sp101
            temp_dict[f"Stage3 Propellant Mass"]  = 7100            # Sp101
        elif variant == "Volga":
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
        if "ST" in name:
            temp_dict[f"Fairing Mass"] = 1700
        elif "Fregat" in variant:
            temp_dict[f"Fairing Mass"] = 1800 # NB
        else:
            temp_dict[f"Fairing Mass"] = 1000 # NB
        
    if name in ["Taurus 3110","Taurus 3210","Minotaur-C 3210"]:
        # Stage 2 - Orion 50SXLG
        # https://cdn.northropgrumman.com/-/media/wp-content/uploads/Orion-Motor-Series.pdf?v=1.0.0
        temp_dict[f"Stage2 Propellant Mass"]    = 33145 * 0.453592 # Manual (weights in lbm)
        temp_dict[f"Stage2 Stage Mass"]         = 2237  * 0.453592 # Manual (weights in lbm)

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
    
    if name.startswith("Vulcan Centaur VC"):
        
        # Boosters are GEM-63XL
        # https://web.archive.org/web/20180918143456/http://www.northropgrumman.com/Capabilities/GEM/Documents/GEM_63_GEM_63XL.pdf
        temp_dict[f"Stage0 Propellant Mass"] = 48000 * int(temp_dict[f"Booster Number"]) # Manual
        temp_dict[f"Stage0 Stage Mass"]      = 5400 * int(temp_dict[f"Booster Number"])  # Manual
        temp_dict[f"Stage1 Propellant Mass"] = 368000        # SLR
        temp_dict[f"Stage1 Stage Mass"]      = 400000-368000 # SLR
    
    if temp_dict["Fairing Mass"] == 0 and name != "Atlas V N22 v2020":
        # This is just an average fairing mass.
        temp_dict["Fairing Mass"] = 1756
       
    temp_dict["Fairing Mass"] = round(temp_dict["Fairing Mass"],1)
    
    return temp_dict