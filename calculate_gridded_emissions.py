#!/usr/bin/python
""" Save out rocket emissions as JSON for output to the API.

    Input Files:
        launchfile     (str): The path of the launch database, a netcdf file containing information for each rocket launch.
        reentryfile    (str): The path of the reentry database, a netcdf file containing information for each object reentry.
        rocketinfofile (str): The path of the rocket database, a netcdf file containing information for each rocket.
        peifile        (str): The path of the pei database, a csv file containing primary emission indices of major emission species for each propellant type.

    Returns:
        JSON output file (json): A database of launch and reentry information and emissions.
"""

# Import modules:
import argparse
from netCDF4 import Dataset # type: ignore
import numpy as np
import pandas as pd
import geopandas as gpd
import fiona
import json
import calendar
import sys
import requests
import xarray as xr
import cProfile
import pstats

from python_modules.distribute_emis_func import make_grid_LL, read_gc_box_height, get_ross_profiles, interp_prop_mass
from python_modules.alt_emis_func import calculate_bc_ei, calculate_nox_ei, calculate_co_ei, calculate_cl_ei
from python_modules.web_scrape_func import scrape_jsr
class InputData:
    '''Read rocket launch and re-entry activity data.'''
    def __init__(self, launchfile: str,reentryfile: str, rocketinfofile: str, peifile: str, start_year):

        def open_file(filepath):
            ds = xr.open_dataset(filepath).to_dataframe().reset_index(drop=True)
            ds = ds.rename(columns={"Time(UTC)": "Time_UTC", "Time (UTC)": "Time_UTC"})
            ds["Date"] = pd.to_datetime(ds["Date"], format="%Y%m%d")
            return ds

        self.dsl  = open_file(launchfile)
        self.handle_falcon_landings()
        if start_year >= 2020:
            self.dsre = open_file(reentryfile)
        self.dsr = xr.open_dataset(rocketinfofile).to_dataframe().reset_index(drop=True)
        self.define_pei(peifile)
    
    def handle_falcon_landings(self):    

        session = requests.Session() 
        df_rcat = scrape_jsr("https://planet4589.org/space/gcat/tsv/cat/rcat.tsv",session)
        self.falcon_landing_dict = {}
        falcon_df = self.dsl[self.dsl["Rocket_Name"].isin(["Falcon 9", "Falcon Heavy"])].copy()

        for row in falcon_df.itertuples(index=False):
            cospar_id = str(row.COSPAR_ID)
                
            df_launch = df_rcat[df_rcat["Launch_Tag"].str.contains(cospar_id, na=False)]
            df_launch = df_launch[df_launch["Type"].str[0].isin(["R"])]

            def falcon_dest(stage):
                stage_row = df_launch[df_launch["PLName"].str.contains(stage, case=False, na=False)]
                if "LZ" in stage_row["Dest"].values[0]:
                    return "ground"
                elif stage_row["Status"].values[0] in ["L","LF"] or stage_row['Dest'].values[0] == "OCISLY":
                    return "ocean"
                elif stage_row["Status"].values[0] == "S":
                    return "expended"
                else:
                    raise ValueError(f"Unknown reentry method for {cospar_id} - {stage_row['Status'].values} - {stage_row['Dest'].values}")

            if len(df_launch) == 1: # Falcon 9
                self.falcon_landing_dict[cospar_id] = falcon_dest("St1")

            elif len(df_launch) == 3: # Falcon Heavy. Sometimes different stages re-enter differently.
                
                self.falcon_landing_dict[cospar_id] = {}
                self.falcon_landing_dict[cospar_id]["center"] = falcon_dest("center")
                self.falcon_landing_dict[cospar_id]["left"]   = falcon_dest("left")
                self.falcon_landing_dict[cospar_id]["right"]  = falcon_dest("right")

            else:
                raise ValueError(f"Unexpected number of Falcon lower stages for {cospar_id} - {len(df_launch)}")
    
    def define_pei(self, peifile):
        """ Define the primary emission indices.
            Each array should contain 5 values, corresponding to the emission index for :
                Hypergolic, Kerosene, Solid, Hydrogen, Methane

        Args:
            peifile (str): The path of the pei database.
        """        
        
        # Define propellant types.
        self.pei_fuel_type     = np.array(['Hypergolic', 'Kerosene', 'Solid', 'Hydrogen', 'Methane'])
        
        # Load pei values.
        pei_file = pd.read_csv(peifile)
        self.h2o_pei   = pei_file["H2O"].to_numpy() # Each variable must be a numpy array in order to make the np.where routine work later.
        self.h2_pei    = pei_file["H2"].to_numpy()          
        self.co_pei    = pei_file["CO"].to_numpy()
        self.co2_pei   = pei_file["CO2"].to_numpy()
        self.bc_pei    = pei_file["BC"].to_numpy()
        self.nox_pei   = pei_file["NOx"].to_numpy()
        self.al2o3_pei = pei_file["Al2O3"].to_numpy()
        self.cly_pei   = pei_file["Cly"].to_numpy()
              
class OutputEmis:
    def __init__(self, event_alts, months, dataset, year, stage_alt_dict):

        ####################################
        # Initial Setup
        ####################################

        # Setup input variables.
        self.event_alts = event_alts
        self.year = year
        self.stage_alt_dict = stage_alt_dict
        self.total_landing_prop = 0
        self.included_prop = 0
        self.missing_prop_total = 0
        self.model_alt = MODEL_ALT
        events_data = {}
        self.csv_count, self.csv_count_2 = 0, 0
        
        # Build list of rocket names and variants
        rocket_pairs = zip(input_data.dsr["Rocket_Name"].values, input_data.dsr["Rocket_Variant"].values)
        
        self.rocket_index_map = {
            (name, variant): i
            for i, (name, variant) in enumerate(rocket_pairs)
        }
        
        # Create variables for totals for a future sanity check:
        # BC launch, CO launch, CO2 launch, NOx launch, H2O launch, Al2O3 launch, Cl launch, HCl launch, Cl2 launch 0-8
        # NOx reentry, Al2O3 reentry, BC reentry, Cl reentry, HCl reentry 9-13
        self.emission_totals = np.zeros(14)
        self.missing_emis = np.zeros(14)
        self.prop_consumed = np.zeros((3,len(input_data.dsl)))
        self.total_prop_consumed = np.zeros((len(input_data.h2o_pei),3))
        self.launch_count = 0
        self.mass_survive_total = 0

        #########################
        # Set up grid
        #########################

        # First get Ross vertical profiles of propellant burned on a fine grid.
        # The smallest GEOS-Chem grid box has a height of 0.129 km, so setting fine grid res to 0.1 km.
        # Setting maximum to 100 km, this is the top of the highest bin in the Ross distribution. 
        self.ross_alt_edge, self.ross_prop_mass, self.ross_cumulative_mass = get_ross_profiles()
        fine_grid_res, fine_grid_top = 0.1,100.0
        self.fine_grid_bot_alt = np.arange(0,100,fine_grid_res)*1e3
        self.fine_grid_mid_alt = np.arange(fine_grid_res/2,fine_grid_top,fine_grid_res)*1e3
        self.fine_grid_top_alt = np.arange(fine_grid_res,fine_grid_top+fine_grid_res,fine_grid_res)*1e3
        self.prop_in_fine_grid, self.fine_grid_mass = interp_prop_mass(self.fine_grid_bot_alt, 
                                                                       self.fine_grid_mid_alt, 
                                                                       self.fine_grid_top_alt, 
                                                                       self.ross_alt_edge, 
                                                                       self.ross_cumulative_mass)
        
        # TODO: Will need to put this in the loop later if we want to use meteorological box heights.
        layer_data = np.loadtxt("./input_files/" + f"geoschem_vertical_grid_{LEVELS}.csv",delimiter=",")    
        bot_alt_base = layer_data[::-1,0]*1000  # shape (LEVELS,)
        mid_alt_base = layer_data[::-1,1]*1000
        top_alt_base = layer_data[::-1,2]*1000

        # Reshape for broadcasting: (LEVELS,1,1)
        bot_alt_3d = bot_alt_base[:, None, None]
        mid_alt_3d = mid_alt_base[:, None, None]
        top_alt_3d = top_alt_base[:, None, None]
        
        # Place re-entry emissions between 60-80 km.
        self.bot_reenter = np.searchsorted(top_alt_base, 60000, side='right')
        self.top_reenter = np.searchsorted(top_alt_base, 80000, side='right')
        self.n_reenter_levs = self.top_reenter - self.bot_reenter + 1

        # Make grid at input resolution (no bounds supplied so full global).
        # Can use this with different in and out edges.
        if GRID_RES == "2x25":
            grid = make_grid_LL("2x2.5")  
        else:
            grid = make_grid_LL(GRID_RES)  
        #Extract useful variables from the grid.    
        self.nlon = len(grid['lon'])
        self.nlat = len(grid['lat'])
        self.lon =  grid['lon']
        self.lat =  grid['lat']

        # Area can also be precomputed
        area_2d = np.ones((self.nlat, self.nlon))

        ########################################
        # Filter based on the year and dataset
        ########################################

        launch_mask = (input_data.dsl["Date"].dt.year == year)
        if dataset == 1:
            launch_mask &= ~input_data.dsl["Megaconstellation_Flag"].astype(bool)
        elif dataset == 2:
            launch_mask &= input_data.dsl["Megaconstellation_Flag"].astype(bool)
        launch_length  = np.sum(launch_mask)

        total_length = launch_length

        if year >= 2020:
            reentry_mask = (input_data.dsre["Date"].dt.year == year)
            if dataset == 1:
                reentry_mask &= ~input_data.dsre["Megaconstellation_Flag"].astype(bool)
            elif dataset == 2:
                reentry_mask &= input_data.dsre["Megaconstellation_Flag"].astype(bool)
            reentry_length = np.sum(reentry_mask)
            total_length += reentry_length

        if launch_length > 0:
            self.output_csv_launch_prop  = np.zeros((launch_length*3,LEVELS))
        if total_length > 0:
            self.output_csv_emis  = np.zeros(((total_length)*10,LEVELS))
        
        #############################################
        # Build a list of ocean landing coordinates.
        #############################################

        # Build the ocean landing list for the year in question.
        raul_data = gpd.read_file('./databases/reentry/General_SpaceX_Map_Raul.kml', driver='KML', layer = 2) # Falcon landing data.  
        ocean_landings = raul_data[
            (raul_data["Name"].str.contains("ASDS")) &
            (raul_data["Description"].str.contains("Landing -"))].copy()
        ocean_landings = ocean_landings[~ocean_landings["Name"].str.contains("planned")].reset_index(drop=True)
        ocean_landings = ocean_landings[~ocean_landings["Name"].str.contains("planed")].reset_index(drop=True)
        date_col = []
        for landing in range(len(ocean_landings["Name"])):
            year_ind = ocean_landings["Description"][landing].index("Landing -")
            date_col.append(ocean_landings["Description"][landing][year_ind:].replace(" ","")[8:17].replace("-",""))
        ocean_landings['Date'] = pd.to_datetime(date_col,format='%d%b%y')
        self.ocean_landings = ocean_landings[
            (ocean_landings["Date"] < f'{year+1}-01-01') &
            (ocean_landings["Date"] >= f'{year}-01-01')
        ].reset_index(drop=True)

        #################################
        # Loop over each day and month.
        #################################
        for m in months:
            
            # Process month data
            ndays = calendar.monthrange(self.year, m)[1]
            self.strmon=str(m).zfill(2)
            print('MONTH = ', self.strmon)   

            month_launch_mask = (launch_mask & (input_data.dsl["Date"].dt.month == m))
            if year >= 2020:
                month_reentry_mask = (reentry_mask & (input_data.dsre["Date"].dt.month == m))           
            
            #Loop over all days:
            for d in range(ndays):

                # Find launches and reentries.
                day_launch_mask = ( month_launch_mask &
                    (np.array(input_data.dsl["Date"].dt.day)  == d+1) &
                    (~np.isnan(np.array(input_data.dsl['Time_UTC'])))
                )
                daily_launches_df = input_data.dsl[day_launch_mask].reset_index(drop=True).copy()

                if year >= 2020:

                    day_reentry_mask = ( month_reentry_mask &
                        (np.array(input_data.dsre["Date"].dt.day) == d+1) &
                        (~np.isnan(np.array(input_data.dsre['Time_UTC'])))
                    )                
                    daily_reentries_df = input_data.dsre[day_reentry_mask].reset_index(drop=True).copy()
                else:
                    daily_reentries_df = pd.DataFrame()

                if len(daily_launches_df) + len(daily_reentries_df) == 0:
                    continue

                self.strday = str(d+1).zfill(2) # Process day data
                #if not np.any(day_launch_mask) and not np.any(day_reentry_mask):
                #    continue  # skip to next day
                
                self.bot_alt = bot_alt_3d + np.zeros((LEVELS, self.nlat, self.nlon))
                self.mid_alt = mid_alt_3d + np.zeros((LEVELS, self.nlat, self.nlon))
                self.top_alt = top_alt_3d + np.zeros((LEVELS, self.nlat, self.nlon))
                self.area = area_2d.copy()  # ensure a fresh array per day

                # Initialize:
                self.pmin,self.pmax=np.nan,np.nan
                self.qmin,self.qmax=np.nan,np.nan
                
                ################################################################################
                # Call grid_emis function to calculate distribution from launches and re-entries
                ################################################################################
                daily_launches, daily_reentries = [], []
                if len(daily_launches_df)>0:
                    daily_launches = self.grid_emis(daily_launches_df,'launch')
       
                if len(daily_reentries_df)>0:
                    daily_reentries = self.grid_emis(daily_reentries_df,'reentry')
                
                daily_events_data = {"launches": daily_launches,
                                     "reentries": daily_reentries}
                events_data[f"{self.year}-{self.strmon}-{self.strday}"] = daily_events_data 
                
        ####################################
        # Conservation of mass check
        ####################################
        
        # At the end of the year, check that the emissions going out to the netcdf files are the same as those from the fine grid.
        # This is done by comparing the emis_distribution output (equivalent to netcdf without the conversion to kgm-2s-1) and the totals.

        # BC launch, CO launch, CO2 launch, NOx launch, H2O launch, Al2O3 launch, Cl launch, HCl launch, Cl2 launch 0-8
        # NOx reentry, Al2O3 reentry, BC reentry, Cl reentry, HCl reentry 9-13
        final_emis = np.asarray([self.emission_totals[0]+self.emission_totals[11],
                                 self.emission_totals[1],
                                 self.emission_totals[2],
                                 self.emission_totals[3]+self.emission_totals[9],
                                 self.emission_totals[4], 
                                 self.emission_totals[5]+self.emission_totals[10],
                                 self.emission_totals[6]+self.emission_totals[12],
                                 self.emission_totals[7]+self.emission_totals[13],
                                 self.emission_totals[8]])
        
        spec_names = ["BC", "CO", "CO2", "NOx", "H2O", "Al2O3", "Cl", "HCl", "Cl2"]
        diff_out = np.zeros(9)
        for event in range(int(np.shape(self.output_csv_emis)[0]/10)):
            for spec in range(9):
                diff_out[spec] += np.sum(self.output_csv_emis[event*10+spec,:]) 
        
        for i_emis in range(9):
            with np.errstate(divide='ignore', invalid='ignore'):
                if np.abs(diff_out[i_emis]-final_emis[i_emis])/final_emis[i_emis]*100.0 > 0.1:
                    raise RuntimeError(f"Error with {spec_names[i_emis]} emissions - {np.abs(diff_out[i_emis]-final_emis[i_emis])/final_emis[i_emis]*100.0:.2f}%.") 

        if dataset == 3:
            filename = f'./out_files/{(year // 10) * 10}/data_{self.year}.json'    
        else:
            filename = f'./out_files/{(year // 10) * 10}/data_{self.year}_{dataset}.json'                            
        with open(filename, 'w') as json_file:
            json.dump(events_data, json_file, indent=4)
             
    def process_launch_event_altitudes(self, roc_ind, launch_rocket, launch_variant, launch_id):
        
        ################################################################
        # Define the propellant saved for landing for reusable rockets.
        ################################################################

        REENTRY = 5.6
        LANDING_BURN = 1.2
        BOOSTBACK = 5.6

        landing_percent_dict = {
            "ocean": 100.0 - REENTRY - LANDING_BURN,
            "ground": 100.0 - REENTRY - BOOSTBACK - LANDING_BURN,
            "expended": 100.0
        }

        if launch_rocket == "Falcon 9":
            booster_percent = 0.0
            landing_info = input_data.falcon_landing_dict[launch_id]
            s1_percent = landing_percent_dict.get(landing_info, None)
            if s1_percent is None: 
                raise ValueError(f"Unexpected landing site for {launch_id}")
        elif launch_rocket == "Falcon Heavy":
            landing_info_center = input_data.falcon_landing_dict[launch_id]["center"]
            landing_info_booster_left = input_data.falcon_landing_dict[launch_id]["left"]
            landing_info_booster_right = input_data.falcon_landing_dict[launch_id]["right"]
            if landing_info_booster_left != landing_info_booster_right:
                raise ValueError(f"Different landing sites for boosters of {launch_id}")
            s1_percent = landing_percent_dict.get(landing_info_center, None)
            booster_percent = landing_percent_dict.get(landing_info_booster_left, None)
            if s1_percent is None or booster_percent is None:
                raise ValueError(f"Unexpected landing site for {launch_id} - {s1_percent} {booster_percent}")
        else:
            booster_percent = 100.0
            s1_percent = 100.0

        ############################################
        # Determine event altitudes for each stage.
        ############################################
                    
        stage_keys = ['BECO', 'MECO', 'SEI1', 'SECO']
        stage_alts = {key: self.stage_alt_dict.get(f"{launch_rocket} {launch_variant} {key}", np.nan) for key in stage_keys}

        stage_alt_beco = stage_alts['BECO']
        stage_alt_meco = stage_alts['MECO']
        stage_alt_sei  = stage_alts['SEI1']
        stage_alt_seco = stage_alts['SECO']

        row = input_data.dsr.iloc[roc_ind]
        stages = tuple(bool(row[f"Stage{i}_Fuel_Type"]) for i in range(6))
        
        if np.isnan(stage_alt_meco) and np.isnan(stage_alt_sei) and np.isnan(stage_alt_seco):

            config_map = {
                # booster, s1, s2, s3, s4, s5
                (True,  True,  False, False, False, False): 0,
                (True,  True,  True,  False, False, False): 0,
                (True,  True,  True,  True,  False, False): 1,
                (True,  True,  True,  True,  True , False): 2,
                (True,  True,  True,  True,  True , True):  2,
                (False, True,  True,  False, False, False): 3,
                (False, True,  True,  True,  False, False): 4,
                (False, True,  True,  True,  True , False): 5,
                (False, True,  True,  True,  True , True):  5
            }
            # Look up configuration safely
            rocket_config_type = config_map.get(stages)

            if rocket_config_type is None:
                raise IndexError(f"Incorrect rocket configuration for {launch_rocket} (stages={stages})")
                    
            stage_alt_beco = self.event_alts[0][rocket_config_type]
            stage_alt_meco = self.event_alts[1][rocket_config_type]
            stage_alt_sei  = self.event_alts[2][rocket_config_type]    
            stage_alt_seco = self.event_alts[3][rocket_config_type]

        ################################################################
        # Distribute propellant mass over fine grid based on altitudes.
        ################################################################

        def get_alt_index(alt_km, fine_top, fine_bot):
            """Return index of layer containing altitude, or None if above grid."""
            if alt_km*1e3 > fine_top[-1]:
                return None
            else:
                idx = np.searchsorted(fine_top, alt_km*1e3, side='right') 
                if not (np.round(fine_bot[idx],1) <= alt_km*1e3 <= np.round(fine_top[idx],1)):
                    raise IndexError(f"Altitude indexing error {fine_bot[idx]} {alt_km*1e3} {fine_top[idx]}")
                return idx + 1

        def normalize_mass(mass_array, percent):
            mass_array = np.asarray(mass_array)
            total = np.sum(mass_array)
            if total == 0:
                print(mass_array)
                raise ValueError(f"Zero fine-grid mass {percent}")
            return mass_array / total * percent
        
        fei_alt_dict = {
            "GSLV Mk III": 46.534,
            "Pegasus": 11.9, # All Pegasus values from SLR.
            "Pegasus H": 11.9,
            "Pegasus XL": 11.9,
            "Pegasus XL/HAPS": 11.9,
            "Pegasus/HAPS": 11.9,
            "LauncherOne": 12.37,
        }  
        # Get altitude, default to 0 if rocket not in dictionary
        fei_alt = fei_alt_dict.get(launch_rocket, 0)

        self.fine_grid_mass_stages = [np.array([]) for _ in range(8)]
        
        ###########
        # Boosters
        ###########

        if stages[0]:
            
            self.booster_alt_index = get_alt_index(stage_alt_beco, self.fine_grid_top_alt, self.fine_grid_bot_alt)

            if stage_alt_beco * 1e3 > self.fine_grid_top_alt[-1]:
                self.fine_grid_mass_stages[0] = np.asarray(self.fine_grid_mass)
            else:
                self.fine_grid_mass_stages[0] = normalize_mass(self.fine_grid_mass[:self.booster_alt_index].copy(), booster_percent)
                         
            if launch_rocket == "Falcon Heavy":
                self.total_landing_prop += (100.0 - booster_percent) * input_data.dsr["Stage0_PropMass"][roc_ind] / 100.0

        ###########
        # Stage 1
        ###########  
        # TODO: Needs reworking if running for different model ceilings above 80km.            
        self.fei_alt_index = np.searchsorted(self.fine_grid_top_alt, fei_alt*1e3)
        self.MECO_alt_index = get_alt_index(stage_alt_meco, self.fine_grid_top_alt, self.fine_grid_bot_alt)
        if stage_alt_meco * 1e3 > self.fine_grid_top_alt[-1]:
            self.fine_grid_mass_stages[1] = self.fine_grid_mass[self.fei_alt_index:].copy()
            if launch_rocket == "GSLV Mk III":
                fei_mass = np.interp(fei_alt,self.ross_alt_edge, self.ross_cumulative_mass)
                fei_percent = (self.prop_in_fine_grid - fei_mass) / (100.0 - fei_mass) * 100.0
                self.fine_grid_mass_stages[1] = normalize_mass(self.fine_grid_mass_stages[1], fei_percent)
        else:
            fine_grid_mass_stage1 = self.fine_grid_mass[self.fei_alt_index:self.MECO_alt_index].copy()
            self.fine_grid_mass_stages[1] = normalize_mass(fine_grid_mass_stage1, s1_percent)
                        
        if launch_rocket == "Falcon Heavy":
            self.total_landing_prop += (100.0 - s1_percent) * input_data.dsr["Stage1_PropMass"][roc_ind] / 100.0

        ###########
        # Stage 2
        ###########
                               
        # If SEI occurs above the fine grid, then we can ignore the second stage.
        if stage_alt_sei * 1e3 >= self.fine_grid_top_alt[-1]:
            self.fine_grid_mass_stages[2] = np.asarray([0])
            self.sei_alt_index, self.seco_alt_index = None, None   
        
        # If SEI is within fine grid, we need to differentiate based on SECO.
        else:
            self.sei_alt_index = np.searchsorted(self.fine_grid_top_alt, stage_alt_sei*1e3, side='right')
            self.seco_alt_index = get_alt_index(stage_alt_seco, self.fine_grid_top_alt, self.fine_grid_bot_alt)
            
            # If SECO occurs above the fine grid, then we only use some of the stage 2 emissions.
            if stage_alt_seco * 1e3 > self.fine_grid_top_alt[-1]:
                # Interpolate the cumulative mass to optimise what percentage of propellant is within GEOS-Chem.
                sei_mass = np.interp(stage_alt_sei,self.ross_alt_edge, self.ross_cumulative_mass)
                seco_percent = (self.prop_in_fine_grid - sei_mass) / (100.0 - sei_mass) * 100.0
            else:
                seco_percent = 100.0
            self.fine_grid_mass_stages[2] = normalize_mass(self.fine_grid_mass[self.sei_alt_index:self.seco_alt_index].copy(), seco_percent)
            
        ###########
        # Stage 3
        ###########
                
        # TODO: Check if any other rocket third stages start below 100 km.        
        # Calculate third stage emissions for Minotaur 1 only.    
        if launch_rocket == "Minotaur 1":
            # Minotaur 1 third stage ignition occurs below the model boundary, so we can need to add some of the stage 3 emissions.
            self.TEI_alt_index = np.searchsorted(self.fine_grid_top_alt, 73.47*1e3, side='right')
            # Now interpolate the cumulative mass to optimise what percentage of propellant is within GEOS-Chem.
            tei_mass = np.interp(73.47,self.ross_alt_edge, self.ross_cumulative_mass)
            tei_percent = (self.prop_in_fine_grid - tei_mass) / (100.0 - tei_mass) * 100.0
            self.fine_grid_mass_stages[3] = normalize_mass(self.fine_grid_mass[self.TEI_alt_index:].copy(), tei_percent)
         
        ##################
        # Falcon Reusable
        ##################
           
        # Finally create an extra stage just for the landing emissions of Falcon 9 v1.2. 
        # TODO: Needs reworking if running for different model ceilings above 80km.
        # Would need to add boostback emissions if top is above 80.
        if launch_rocket in ["Falcon 9","Falcon Heavy"]:
            if input_data.falcon_landing_dict[launch_id] != "expended":

                # 5.6% are used in the entry burn, over 70-54.7 km.
                self.entry_top = np.searchsorted(self.fine_grid_top_alt, 70000, side='right')
                self.entry_bot = np.searchsorted(self.fine_grid_top_alt, 54700, side='right')
                self.fine_grid_mass_stages[6] = normalize_mass(self.fine_grid_mass[self.entry_bot:self.entry_top+1].copy(), 5.6)
                
                # 1.2% are used in the landing burn, over 3.3-0 km.
                self.landing_top = np.searchsorted(self.fine_grid_top_alt, 3300, side='right')
                self.fine_grid_mass_stages[7] = normalize_mass(self.fine_grid_mass[:self.landing_top+1].copy(), 1.2)
        
        return stage_alt_beco, stage_alt_meco, stages
    
    def calc_emis(self, start_ind, stop_ind, total_vertical_propellant, stage, launch_tuple, launch_details):
        '''
        Calculate the emissions over the range of the stag within the fine grid (0-100 km).
        '''

        time_index, q, p, pei_indices, prop_masses, falcon_flag = launch_tuple
        pei_index = pei_indices[stage]
        prop_mass = prop_masses[stage]
        species_keys = ["launch_bc","co","co2","launch_nox","fuel_nox","h2o","launch_al","launch_cl","launch_hcl","cl2"]

        if falcon_flag:
            if start_ind != None:
                vertical_profile = self.fine_grid_mass_stages[6]
            else:
                vertical_profile = self.fine_grid_mass_stages[7]
        else:
            vertical_profile = self.fine_grid_mass_stages[stage]

        if np.sum(1e-2 * vertical_profile) > 1.01:
            raise ValueError("Error with Stage {stage}. Propellant distribution exceeds unity.")

        seg_slice = slice(start_ind, stop_ind)
        seg_mid_alts_km = self.fine_grid_mid_alt[seg_slice] * 1e-3  # km

        # Calculate the emission indices for each species.
        ei_bc                 = calculate_bc_ei (seg_mid_alts_km, input_data.bc_pei[pei_index])
        ei_co, ei_co2         = calculate_co_ei (seg_mid_alts_km, input_data.co_pei[pei_index], input_data.co2_pei[pei_index])
        ei_sec_nox            = calculate_nox_ei(seg_mid_alts_km)
        ei_cl, ei_hcl, ei_cl2 = calculate_cl_ei (seg_mid_alts_km, input_data.cly_pei[pei_index])
        ei_h2o                = input_data.h2o_pei[pei_index] + input_data.h2_pei[pei_index] * (1.008 * 2 + 16) / (1.008 * 2)

        # Calculate the emissions for each species in g.
        emis_full = np.zeros((len(self.fine_grid_mid_alt),11))
        factor = prop_mass * 1e-2 * vertical_profile
        for index, emission_index in enumerate([ei_bc, ei_co, ei_co2, ei_sec_nox, input_data.nox_pei[pei_index], ei_h2o, input_data.al2o3_pei[pei_index], ei_cl, ei_hcl, ei_cl2]):
            emis_full[seg_slice,index] = emission_index * factor
        emis_full[seg_slice,10] = vertical_profile

        # Place the emissions into a larger array covering the whole fine grid (0-100km).
        selected_alts = []
        bound_arr = np.append(self.bot_alt[:, q, p], self.top_alt[-1, q, p])
        bound_idx = np.searchsorted(self.fine_grid_bot_alt, bound_arr)

        # Regrid the vertical_profile to the desired model profile.
        for i in range(len(bound_arr)-1):

            bot_ind = bound_idx[i]
            top_ind = bound_idx[i+1]
            selected_alts.extend(np.arange(bot_ind, top_ind))

            # Slice and sum once to save time.
            emis_slice = emis_full[bot_ind:top_ind,:]
            emis_sum   = np.sum(emis_slice, axis=0)
            
            # Sum the emissions in this range for each species and place in an array.
            for i, key in enumerate(species_keys):
                self.rocket_data_arrays[key][time_index,i,q,p]  += (emis_sum[i] * 1e-6)
                launch_details["emissions"][key] += (emis_sum[i] * 1e-6)
            total_vertical_propellant[i,0] += emis_sum[10] * prop_mass
            total_vertical_propellant[i,1] += emis_sum[0] 
            total_vertical_propellant[i,2] += emis_sum[1] 
            total_vertical_propellant[i,3] += emis_sum[2]
            total_vertical_propellant[i,4] += np.sum(emis_sum[3:5]) 
            total_vertical_propellant[i,5] += emis_sum[5]
            total_vertical_propellant[i,6] += emis_sum[6]
            total_vertical_propellant[i,7] += emis_sum[7]
            total_vertical_propellant[i,8] += emis_sum[8]
            total_vertical_propellant[i,9] += emis_sum[9]
             
        if len(list(set(selected_alts))) != len(selected_alts):
            raise IndexError("Error in fine grid indexing.")
        
        # Slice and sum once to save time.
        total_slice = emis_full[selected_alts[0]:selected_alts[-1]+1,:]
        total_sum   = np.sum(total_slice, axis=0)
        
        # BC launch, CO launch, CO2 launch, NOx launch, H2O launch, Al2O3 launch, Cl launch, HCl launch, Cl2 launch
        # NOx reentry, Al2O3 reentry, BC reentry, HCl reentry, Cl reentry
        self.emission_totals[0] += total_sum[0]       
        self.emission_totals[1] += total_sum[1]         
        self.emission_totals[2] += total_sum[2]   
        self.emission_totals[3] += np.sum(total_sum[3:5])     
        self.emission_totals[4] += total_sum[5]   
        self.emission_totals[5] += total_sum[6]   
        self.emission_totals[6] += total_sum[7]   
        self.emission_totals[7] += total_sum[8]   
        self.emission_totals[8] += total_sum[9]   
        self.included_prop      += (total_sum[10] * prop_mass * 1e-2)

        if falcon_flag:
            self.missing_prop[stage] += np.round((np.sum(emis_full[selected_alts[0]:selected_alts[-1]+1,10])),2)
        else:
            self.missing_prop[stage] += np.round((100-np.sum(emis_full[selected_alts[0]:selected_alts[-1]+1,10])),2)
        
        return total_vertical_propellant        

    def calc_missing_emis(self, pei_index, prop_mass, percent_included):

        ei_bc = calculate_bc_ei([self.model_alt], input_data.bc_pei[pei_index])
        ei_co, ei_co2 = calculate_co_ei([self.model_alt], input_data.co_pei[pei_index], input_data.co2_pei[pei_index])
        ei_sec_nox = calculate_nox_ei([self.model_alt])
        ei_h2o =  input_data.h2o_pei[pei_index] + input_data.h2_pei[pei_index] * (1.008 * 2 + 16) / (1.008 * 2) 
        ei_cl, ei_hcl, ei_cl2 = calculate_cl_ei([self.model_alt], input_data.cly_pei[pei_index]) 
        
        t_bc_emis = ei_bc * prop_mass * 1e-2 * percent_included
        t_co_emis = ei_co * prop_mass * 1e-2 * percent_included
        t_co2_emis = ei_co2 * prop_mass * 1e-2 * percent_included
        t_launch_nox_emis = ei_sec_nox * prop_mass * 1e-2 * percent_included
        t_fuel_nox_emis = input_data.nox_pei[pei_index] * prop_mass * 1e-2 * percent_included
        t_h2o_emis = ei_h2o * prop_mass * 1e-2 * percent_included
        t_al2o3_emis = input_data.al2o3_pei[pei_index] * prop_mass * 1e-2 * percent_included
        t_cl_emis = ei_cl * prop_mass * 1e-2 * percent_included
        t_hcl_emis = ei_hcl * prop_mass * 1e-2 * percent_included
        t_cl2_emis = ei_cl2 * prop_mass * 1e-2 * percent_included
        
        self.missing_prop_total += prop_mass * 1e-2 * percent_included
        self.missing_emis[0] += np.sum(t_bc_emis) 
        self.missing_emis[1] += np.sum(t_co_emis)
        self.missing_emis[2] += np.sum(t_launch_nox_emis)
        self.missing_emis[2] += np.sum(t_fuel_nox_emis)
        self.missing_emis[3] += np.sum(t_h2o_emis)
        self.missing_emis[4] += np.sum(t_al2o3_emis)
        self.missing_emis[5] += np.sum(t_cl_emis)
        self.missing_emis[6] += np.sum(t_hcl_emis)
        self.missing_emis[7] += np.sum(t_cl2_emis)
        self.missing_emis[8] += np.sum(t_co2_emis)    

    def launch_emis(self,row,q,p,time_index):
        '''Calculate the launch emissions'''

        # Interpolate propellant mass distribution.
        # Linearly interpolate Ross profile proportion of propellant mass to GEOS-Chem vertical grid.
        # This is just to double check that the propellant consumption profile is calculated correctly.
        self.propellant_in_model, self.gc_relative_mass = interp_prop_mass(self.bot_alt[:,q,p], self.mid_alt[:,q,p], self.top_alt[:,q,p],
                                                                            self.ross_alt_edge, self.ross_cumulative_mass)

        ############################################
        # Check the rocket type and get info.
        ############################################

        key = (row.Rocket_Name, row.Rocket_Variant)
        if key in self.rocket_index_map:
            roc_ind = self.rocket_index_map[key]
        else:
            raise IndexError(f"No propellant mass found for {key}")
            
        launch_details = {
            "date":      f"{self.year}-{self.strmon}-{self.strday}",
            "id":        row.COSPAR_ID,
            "time":      row.Time_UTC,
            "site":      "",
            "rocket":    row.Rocket_Name,
            "variant":   row.Rocket_Variant,
            "lat":       row.Latitude,
            "lon":       row.Longitude,
            "smc":       bool(row.Megaconstellation_Flag),
            "location":  row.Site,
            "emissions": {sp: 0 for sp in ["launch_bc","co","co2","launch_nox","fuel_nox","h2o","launch_al","launch_cl","launch_hcl","cl2"]}
        }
        
        ############################################
        # Process the launch event altitudes.
        ############################################
        
        stage_alt_beco, stage_alt_meco, stages = self.process_launch_event_altitudes(roc_ind, row.Rocket_Name, row.Rocket_Variant, row.COSPAR_ID)
        prop_masses, pei_indices = np.zeros(6), np.full(6, -1, dtype=int)
        for i, active in enumerate(stages):
            if active:
                prop_masses[i] = input_data.dsr[f"Stage{i}_PropMass"][roc_ind]

                idx = np.where(input_data.pei_fuel_type == input_data.dsr[f"Stage{i}_Fuel_Type"][roc_ind])[0]
                if len(idx) != 1:
                    raise RuntimeError(f"Expected one PEI fuel match for Stage {i}, got {idx}") 
                pei_indices[i] = idx[0]
        # TODO: Needs reworking if running for different model ceilings above 80km.

        ############################################
        # Deal with failed launches.
        ############################################

        # Most failures are for upper stages, and so can be treated as normal here.
        # Full information is provided in source_info/failed_launch_info.txt.
                    
        # Skip launches where the launch failed close to the launch pad, and all failed launches before 2020.
        if row.COSPAR_ID in ['2020-F04','2021-F04','2023-F02','2024-F01'] or (row.COSPAR_ID[5] == "F" and self.year < 2020):
            self.csv_count += 3
            self.csv_count_2 += 10
            return
        
        # This is where stage 1 failed during ascent.
        failed_alt_events = {
            '2020-F07': 900,    # Astra Rocket 3 launch, rocket shut off at 0.9km.
            '2021-F01': 10700,  # Shuang Quxian-1 launch, rocket disintegrated at Max-Q. Approximating altitude using Proton-M and Minotaur-IV max-q alts.
            '2021-F07': 31000   # Astra Rocket 3 launch, rocket shut off at 31km.
        }
        if row.COSPAR_ID in failed_alt_events:
            failed_alt = failed_alt_events[row.COSPAR_ID]
            cutoff_ind = np.searchsorted(self.fine_grid_mid_alt, failed_alt, side='right')
            self.fine_grid_mass_stages[1][cutoff_ind:] = 0.0 # Zero above the cutoff_ind.

        # This is where stage 2 never ignited.
        zero_stage2_events = {
            '2020-F02','2020-F05',
            '2021-F02','2021-F07','2021-F08','2022-F01',
            '2022-F02','2022-F03','2023-F01','2023-F04',
            '2023-F05','2023-F06','2023-F07','2023-F09',
            '2024-F02','2024-F04'
        }
        if row.COSPAR_ID in zero_stage2_events:
            self.fine_grid_mass_stages[2] = np.asarray([0])
            self.sei_alt_index = None
        
        ##################################################
        # Sanity Checks for Propellant Mass Distributions.
        ##################################################

        # The propellant consumed should never be bigger than the total propellant.
        error_lim = 0.001 # Maximum error from rounding / floating point errors.
        def check_error(stage,cospar,error_lim):
            with np.errstate(divide='ignore', invalid='ignore'):
                per_error = ((self.prop_consumed[stage,self.launch_count] - prop_masses[stage]) / prop_masses[stage] * 100.0)
                if per_error > error_lim:
                    raise ValueError(f'Error with emissions for stage {stage} of {cospar}. Prop consumed: {self.prop_consumed[stage,self.launch_count]}, Prop mass: {prop_masses[stage]}, Per Error: {per_error}.') 
            
        if stages[0]: # If there are boosters.
            self.prop_consumed[0,self.launch_count] = np.sum(self.fine_grid_mass_stages[0] * prop_masses[0]  * 1e-2)
            if (stage_alt_beco < self.model_alt) and (row.Rocket_Name != "Falcon Heavy"): # Suppressed for Falcon Heavy, as this has reusable boosters.
                check_error(0,row.COSPAR_ID,error_lim)
                
        self.prop_consumed[1,self.launch_count]  = np.sum(self.fine_grid_mass_stages[1] * prop_masses[1] * 1e-2)
        self.prop_consumed[2,self.launch_count]  = np.sum(self.fine_grid_mass_stages[2] * prop_masses[2] * 1e-2)

        # For all rockets, the propellant consumed for each stage should never be bigger than the total propellant in each stage.
        check_error(1,row.COSPAR_ID,error_lim)
        if row.Rocket_Name != "Long March (CZ) 5B":
            check_error(2,row.COSPAR_ID,error_lim)

        # When MECO occurs in the model, the consumed propellant should be within 1% of the total propellant mass of stage 1.  
        # The error is suppressed for rockets that failed early into the launch.
        # Also suppressed for Falcon 9, as this has a reusable first stage.
        # A check that the Falcon landing distribution is no greater than 7% of the stage 1 emissions is undertaken later in the grid_emis function.
        if stage_alt_meco < self.model_alt:
            cospar = row.COSPAR_ID

            if (cospar not in failed_alt_events): 
                rocket = input_data.dsr["Rocket_Name"][roc_ind]

                if rocket != "Falcon 9":
                    check_error(1,cospar,1)
                elif input_data.falcon_landing_dict[cospar] != "expended":
                    check_error(1,cospar,1)
                
        self.total_prop_consumed[pei_indices[0],0]  += self.prop_consumed[0,self.launch_count]
        self.total_prop_consumed[pei_indices[1],1]  += self.prop_consumed[1,self.launch_count]
        self.total_prop_consumed[pei_indices[2],2]  += self.prop_consumed[2,self.launch_count]                  
        
        ##############################################
        # Calculate the emissions for each species.
        ##############################################     
        
        # Creating a 2d array for the prop output.
        # Total prop, bc, co, co2, nox, h2o, al2o3, cl, hcl, cl2
        total_vertical_propellant = np.zeros((len(self.mid_alt[:,q,p]),10)) 
        self.missing_prop = np.zeros(6) # An array for propellant consumed >80 km.

        launch_tuple = (time_index, q, p, pei_indices, prop_masses, False)  
                        
        # Check whether there is a booster:
        if stages[0]:
            total_vertical_propellant = self.calc_emis(None,
                        self.booster_alt_index,
                        total_vertical_propellant,
                        0,
                        launch_tuple,
                        launch_details)
                
        # Every rocket has a first stage.
        total_vertical_propellant = self.calc_emis(self.fei_alt_index,
                                                    self.MECO_alt_index,
                                                    total_vertical_propellant,
                                                    1,
                                                    launch_tuple,
                                                    launch_details
                                                    )
        
        # Check whether there is a second stage:
        if stages[2] and self.sei_alt_index != None:
            total_vertical_propellant = self.calc_emis(self.sei_alt_index,
                        self.seco_alt_index,
                        total_vertical_propellant,
                        2,
                        launch_tuple,
                        launch_details)
            
        # NOTE: This section needs to be tweaked if wanting to run for different vertical heights above 80km.
        # Check for more rockets that have third stage emissions within model.    
        # Add third stage emissions for Minotaur 1.
        if row.Rocket_Name == "Minotaur 1":
            total_vertical_propellant = self.calc_emis(self.TEI_alt_index,
                        None,
                        total_vertical_propellant,
                        3,
                        launch_tuple,
                        launch_details)
            
        # If the rocket is a Falcon 9, then add the landing emissions. Kerosene, so no Al2O3 or Cly.     
        
        if row.Rocket_Name in ["Falcon 9","Falcon Heavy"]:
            if input_data.falcon_landing_dict[row.COSPAR_ID] == "expended":
                pass
            else: 
                falcon_p, falcon_q = None, None

                if input_data.falcon_landing_dict[row.COSPAR_ID] == "ground":
                    falcon_p = p
                    falcon_q = q           
                else:
                    matching = self.ocean_landings[self.ocean_landings["Date"].isin([f"{self.year}-{self.strmon}-{self.strday}"])].reset_index(drop=True)
                    # There is a typo in the database, a 2023 launch is listed as 2022.      
                    if (matching.shape[0] == 1) or (f"{self.year}-{self.strmon}-{self.strday}" == "2022-04-27"):

                        falcon_lat = np.array(matching["geometry"].y)[0]
                        falcon_lon = np.array(matching["geometry"].x)[0]
                        
                    # There are two launches on the same day.
                    elif row.COSPAR_ID in ["2022-124","2023-037"]:
                        falcon_lat = np.array(matching["geometry"].y)[0]
                        falcon_lon = np.array(matching["geometry"].x)[0]
                        
                    elif row.COSPAR_ID in ["2022-125","2023-038"]:
                        falcon_lat = np.array(matching["geometry"].y)[1]
                        falcon_lon = np.array(matching["geometry"].x)[1]
                        
                    elif matching.shape[0] == 0:
                        # The database hasn't been well updated for 2022, so lets just fill in based on most common geolocation for all other 2020-2022 launches.
                        if np.isclose(row.Longitude, -81.0, atol=1.0) and np.isclose(row.Latitude, 29.0, atol=1.0):
                            falcon_lon = -75
                            falcon_lat = 32
                        elif np.isclose(row.Longitude, -120.6, atol=1.0) and np.isclose(row.Latitude, 34.7, atol=1.0):
                            falcon_lon = -122.5
                            falcon_lat = 30
                        elif np.isclose(row.Longitude, 100.3, atol=1.0) and np.isclose(row.Latitude, 41.3, atol=1.0):
                            falcon_lon = 100.3
                            falcon_lat = 41.3
                        else:
                            raise RuntimeError("Launch not from assigned site.",row.Longitude,row.Latitude)
                    
                    elif matching.shape[0] > 1: 
                        raise RuntimeError(f"Multiple ocean entries for Falcon Stage 1 landing for {row.COSPAR_ID}.")
                    else:
                        raise RuntimeError(f"Problem geolocating Falcon Stage 1 landing for {row.COSPAR_ID}, {self.lon[falcon_p],self.lat[falcon_q]}") 

                    falcon_p = np.argmin(abs(falcon_lon-self.lon))
                    falcon_q = np.argmin(abs(falcon_lat-self.lat))
                        
                    if falcon_p > self.pmax:
                        self.pmax = falcon_p
                    if falcon_q > self.qmax:
                        self.qmax = falcon_q
                
                if np.sum(self.fine_grid_mass_stages[6]) + np.sum(self.fine_grid_mass_stages[7]) > 7:
                    raise RuntimeError("Error with Stage 1 Ocean Landing. Propellant distribution exceeds what is expected.")
                
                falcon_map = {
                    "Falcon 9": 1,
                    "Falcon Heavy": 0,
                }

                try:
                    falcon_stage = falcon_map[row.Rocket_Name]
                except KeyError:
                    raise RuntimeError(f"Error identifying Falcon rocket type: {row.Rocket_Name}")
                
                # NOTE: This section needs to be tweaked if wanting to run for different vertical heights above 80km.
                # Should implement the boostback burn if the vertical height is increased to 100 km.
                # First the entry emissions.
                falcon_tuple = (time_index, falcon_q, falcon_p, pei_indices, prop_masses, True) 
                total_vertical_propellant = self.calc_emis(self.entry_bot,
                            self.entry_top+1,
                            total_vertical_propellant,
                            falcon_stage,
                            falcon_tuple,
                            launch_details)
                # Now the landing emissions.
                total_vertical_propellant = self.calc_emis(None,
                            self.landing_top+1,
                            total_vertical_propellant,
                            falcon_stage,
                            falcon_tuple,
                            launch_details)

        # NOTE: This section may need to be tweaked if wanting to run for different vertical heights above 80km. 
        
        # Calculate the missing emissions from above model.
        # Check if there are missing emissions from the boosters stage.
        if stages[0] and self.missing_prop[0] > 0:
            if 100 - np.sum(self.fine_grid_mass_stages[0]) > 100:
                sys.exit("Error with Booster emissions above model.") 
            self.calc_missing_emis(pei_indices[0],prop_masses[0],self.missing_prop[0])                  
    
        # Check if there are missing emissions from the first stage.
        if stages[1] and self.missing_prop[1] > 0 and row.COSPAR_ID not in ['2020-F04','2020-F07','2021-F01','2021-F07']:
            if 100 - np.sum(self.fine_grid_mass_stages[1]) > 100:
                sys.exit("Error with Stage 1 emissions above model.") 

            self.calc_missing_emis(pei_indices[1],prop_masses[1],self.missing_prop[1] - self.missing_prop[5])
            
        # Check if there are missing emissions from the second stage.
        if stages[2] and self.missing_prop[2] > 0 and row.COSPAR_ID not in ['2020-F02','2020-F04','2020-F05','2020-F07','2021-F01',
                                                                            '2021-F02','2021-F07','2021-F08','2022-F01','2022-F02']:
            if 100 - np.sum(self.fine_grid_mass_stages[2]) > 100:
                sys.exit("Error with Stage 2 emissions above model.") 
            
            if row.COSPAR_ID == "2022-F03":
                included_emis = 240/315*100 - (100-self.missing_prop[2])
            else:
                included_emis = self.missing_prop[2]
                
            self.calc_missing_emis(pei_indices[2],prop_masses[2],included_emis)
            
        # Check if there are missing emissions from the third stage.
        if stages[3] and row.COSPAR_ID not in ['2020-F02','2020-F03','2020-F05','2020-F06','2021-F01',
                                                '2021-F06','2021-F10','2022-F02','2022-F05','2022-F07']:
            if row.COSPAR_ID == "2021-F09":
                included_emis = 475/521*100
            elif row.Rocket_Name == "Minotaur 1":
                included_emis = self.missing_prop[3]
            else:
                included_emis = 100
                
            self.calc_missing_emis(pei_indices[3],prop_masses[3],included_emis)
        
        # Check if there are missing emissions from the fourth stage.
        if stages[4] and row.COSPAR_ID[-3:] not in ['2020-F08','2020-F09','2021-F01','2022-F02',
                                                    '2022-F04','2022-F05','2022-F07']:
            self.calc_missing_emis(pei_indices[4],prop_masses[4],100)

        if stages[5]:
            self.calc_missing_emis(pei_indices[5],prop_masses[5],100)
        new_dict = {
            "BC":    f'{launch_details["emissions"]["launch_bc"]:.6f}',
            "CO":    f'{launch_details["emissions"]["co"]:.6f}',
            "CO2":   f'{launch_details["emissions"]["co2"]:.6f}',
            "NOx":   f'{(launch_details["emissions"]["launch_nox"] + launch_details["emissions"]["fuel_nox"]):.6f}',
            "H2O":   f'{launch_details["emissions"]["h2o"]:.6f}',
            "Cly":   f'{(launch_details["emissions"]["launch_cl"] + launch_details["emissions"]["cl2"] + launch_details["emissions"]["launch_hcl"]):.6f}',
            "Al2O3": f'{launch_details["emissions"]["launch_al"]:.6f}'
        }
        launch_details["emissions"] = new_dict                 
        
        ##############################################
        # Output the emissions to a file for viewing.
        ##############################################                  

        self.output_csv_launch_prop[self.csv_count,:] = total_vertical_propellant[:,0]
        self.csv_count += 1
        self.output_csv_launch_prop[self.csv_count,:] = self.mid_alt[:,q,p]*1e-3
        self.csv_count += 1  
        self.output_csv_launch_prop[self.csv_count,:] = self.gc_relative_mass
        self.csv_count += 1 
        
        if np.sum(self.gc_relative_mass) == 0 and np.sum(total_vertical_propellant[:,0]) == 0:
            print(f'??? {input_data.dsr["Rocket_Name"][roc_ind]}')
        
        for i in range(1,10):
            self.output_csv_emis[self.csv_count_2,:] = total_vertical_propellant[:,i]
            self.csv_count_2 += 1 
        self.output_csv_emis[self.csv_count_2,:] = self.mid_alt[:,q,p]*1e-3
        self.csv_count_2 += 1

        self.launch_count += 1

        return launch_details

    def reentry_emis(self,row,q,p,time_index):
        '''Calculate the reentry emissions'''

        reentry_details = {
            "date": f"{self.year}-{self.strmon}-{self.strday}",
            "id": row.COSPAR_ID,
            "time": row.Time_UTC,
            "reusability": "",
            "name": row.Object_Name,
            "category": row.Category, 
            "lat": row.Latitude,
            "lon": row.Latitude,
            "smc": bool(row.Megaconstellation_Flag),
            "location": int(row.Location_Constraint),
            "burnup": row.Burnup, 
            "emissions": {sp: 0 for sp in ["reentry_nox","reentry_al","reentry_bc","reentry_hcl","reentry_cl"]}
        }
        
        if row.COSPAR_ID[:8] in ["2021-F09","2022-065","2023-72"] and row.Category == "S1":
            reentry_details["lat"] = 34.43194444
            reentry_details["lon"] = 127.535
        elif row.COSPAR_ID[:8] in ["2020-065","2022-167","2022-046","2022-126","2023-135","2024-102","2024-153","2024-173","2024-245","2025-007","2025-105"] and row.Category == "S1":
            reentry_details["lat"] = 34.9
            reentry_details["lon"] = 121.2
        else:
            reentry_details["lat"] = row.Latitude
            reentry_details["lon"] = row.Longitude

        total_vertical_propellant = np.zeros((len(self.mid_alt[:,q,p]),10))
        reentry_ei = np.zeros(5) # Al2O3, NOx, BC, Cl, HCl
        if np.ma.is_masked(row.Ablatable_Mass) or np.ma.is_masked(row.Other_Mass):
            pass
        else:
        
            reentry_ei[0]   = row.Ablation_Degree * row.Percent_Aluminium
            if row.Ablation_Degree == 0:
                reentry_ei[1] = 0.175
            else:
                reentry_ei[1] = 0.4
                
            # Calculate the total mass surviving re-entry in tonnes
            # We are adding all mass that comes back, so this includes reusable stages and fairings.
            self.mass_survive = ((row.Ablatable_Mass* (1-row.Ablation_Degree)) + row.Other_Mass) * 1e-3
                
            # Add chlorine and bc reentry emissions from ATISPADE Report. Don't add for lower stages, only upper (>S2).
            if row.Category in ["P","C"]:
                reentry_ei[2:] = [0.041,0.015,0.008] # Worst Case
            elif row.Category in ["S2","S3","S4"]:
                reentry_ei[2:] = [0.029,0.011,0.005] # Worst Case
                
            # For consistency with launch emissions, the totals are kept in g units. 
            # NOx reentry, Al2O3 reentry, BC reentry, Cl reentry, HCl reentry 9-13           
            t_nox_reentry = (row.Ablatable_Mass + row.Other_Mass) * reentry_ei[1] * 1000
            t_al2o3_reentry = row.Ablatable_Mass * reentry_ei[0] * 1000

            def output_reentry_emis(key,temp,i):
                self.emission_totals[i] += temp
                reentry_details["emissions"][key] += temp * 1e-6   
                self.rocket_data_arrays[key][time_index,self.bot_reenter:self.top_reenter+1,q,p] += temp / self.n_reenter_levs * 1e-6

            output_reentry_emis("reentry_nox",t_nox_reentry,9)
            output_reentry_emis("reentry_al",t_al2o3_reentry,10)

            total_vertical_propellant[self.bot_reenter:self.top_reenter+1,4]   += np.full((self.n_reenter_levs),t_nox_reentry/self.n_reenter_levs)
            total_vertical_propellant[self.bot_reenter:self.top_reenter+1,6]   += np.full((self.n_reenter_levs),t_al2o3_reentry/self.n_reenter_levs)
            
            if row.Category in ["P","C","S2","S3","S4"]:
                
                prop_dict = {"reentry_bc": 1,"reentry_cl":7,"reentry_hcl": 8}
                for i, key in enumerate(["reentry_bc","reentry_cl","reentry_hcl"]):
                    t_reentry = (row.Ablatable_Mass + row.Other_Mass) * reentry_ei[i+2] * 1000
                    self.emission_totals[i+11] += t_reentry
                    reentry_details["emissions"][key] += t_reentry * 1e-6
                    self.rocket_data_arrays[key][time_index,self.bot_reenter:self.top_reenter+1,q,p] += t_reentry / self.n_reenter_levs * 1e-6
                    total_vertical_propellant[self.bot_reenter:self.top_reenter+1,prop_dict[key]]   += np.full((self.n_reenter_levs),t_reentry/self.n_reenter_levs)

        new_dict = {
            "NOx":   f'{reentry_details["emissions"]["reentry_nox"]:.6f}',
            "Al2O3": f'{reentry_details["emissions"]["reentry_al"]:.6f}',
            "BC":    f'{reentry_details["emissions"]["reentry_bc"]:.6f}',
            "HCl":   f'{reentry_details["emissions"]["reentry_hcl"]:.6f}',
            "Cl":    f'{reentry_details["emissions"]["reentry_cl"]:.6f}',
            "Unablated_Mass": self.mass_survive,
        }
        reentry_details["emissions"] = new_dict                
        self.mass_survive_total += self.mass_survive

        for i in range(1,10):
            self.output_csv_emis[self.csv_count_2,:] = total_vertical_propellant[:,i]
            self.csv_count_2 += 1
        self.output_csv_emis[self.csv_count_2,:] = self.mid_alt[:,q,p]*1e-3
        self.csv_count_2 += 1

        return reentry_details
    
    def grid_emis(self, df, emis_type):
        """Grid the data onto the GEOS-Chem horizontal and vertical grid"""
        
        # Pre allocate the emission array.
        daily_info = []
        species = ['launch_nox', 'fuel_nox',  'h2o',       'launch_bc',  'co',
                   'co2',        'launch_al', 'launch_hcl','launch_cl',  'cl2',
                   'reentry_nox','reentry_al','reentry_bc','reentry_hcl','reentry_cl']

        # Set up the emission arrays for this launch/reentry.
        self.rocket_data_arrays = {}
        for sp in species:
            self.rocket_data_arrays[sp] = np.zeros((HOURS, LEVELS, self.nlat, self.nlon))

        #Loop over each launch/reentry.
        for row in df.itertuples():
            
            # TODO: Add this back in when re-entries are readded.
            #if pd.isna(row.COSPAR_ID) or pd.isna(row.Longitude) or pd.isna(row.Latitude):
            #    raise ValueError(f"{row.COSPAR_ID} {row.Longitude} {row.Latitude}")
            
            ################
            # Setup grid.
            ################
            # Get grid horizontal indices (p is lon index; q is lat index).
            # This works out the nearest latitude and longitude on the grid.
            p,q = np.argmin(np.abs(row.Longitude-self.lon)), np.argmin(np.abs(row.Latitude-self.lat))
            
            # Usually, any reentries from launches are placed at the same lat/lon as the launch itself.
            # This code is only used where a non launch related reentry occurs on the same day as a launch.
            # It updates the grid range (lat and lon edges) to account for the new location and to prepare for saving to file:
            if np.isnan(self.pmin):
                self.pmin,self.pmax=p,p
                self.qmin,self.qmax=q,q
            if ~np.isnan(self.pmin) and self.pmin>p:
                self.pmin=p
            if ~np.isnan(self.pmax) and self.pmax<p:
                self.pmax=p
            if ~np.isnan(self.qmin) and self.qmin>q:
                self.qmin=q
            if ~np.isnan(self.qmax) and self.qmax<q:
                self.qmax=q
                                
            # Get factor to convert from kg to kg/m2/s if needed for GEOS-Chem outputs:
            #self.kg_to_kgm2s        = self.area[q,p]*TIMESTEP
            #self.falcon_kg_to_kgm2s = self.area[falcon_q,falcon_p]*TIMESTEP # TODO: Where do we calculate Falcon q and Falcon p
            
            # Format the time of the event.
            time_index = int(row.Time_UTC*60*60)//TIMESTEP
            
            # Work out vertical distribution for launch or reentry emissions.
            if emis_type=='launch':
                launch_details = self.launch_emis(row,q,p,time_index)
                if launch_details is not None:
                    daily_info.append(launch_details)
            elif emis_type=='reentry':
                reentry_details = self.reentry_emis(row,q,p,time_index)
                daily_info.append(reentry_details)
            else:
                raise ValueError("Unexpected event type")

        return daily_info

def check_total_emissions(year,dataset,res,levels,emis_data):
    """ Print total emissions of each species.
        The total emissions including all afterburning are compared to a scenario where only primary emission indices are used.
    """

    # BC launch, CO launch, CO2 launch, NOx launch, H2O launch, Al2O3 launch, Cl launch, HCl launch, Cl2 launch 0-8
    # NOx reentry, Al2O3 reentry, BC reentry, Cl reentry, HCl reentry 9-13

    total_inc_emis = np.sum(emis_data.emission_totals[:]*1e-9)
    total_inc_prop = emis_data.included_prop*1e-6
     
    data = {'Species':                ['BC (Launch)',   'CO (Launch)',     'CO2 (Launch)', 'NOx (Launch)', 'H2O (Launch)', 'Al2O3 (Launch)',
                                       'Cl (Launch)',   'HCl (Launch)',    'Cl2 (Launch)', 
                                       'NOx (Reentry)', 'Al2O3 (Reentry)', 'BC (Reentry)', 'Cl (Reentry)', 'HCl (Reentry)', 
                                       'Cly (Total)',
                                       'NOx (Total)',
                                       'Total Emis',
                                       'Total Prop',
                                       'Surviving Mass'],
            'Emissions 0-80 km [Gg]': [*(emis_data.emission_totals * 1e-9).tolist(),
                                       np.sum(emis_data.emission_totals[6:9]*1e-9) + np.sum(emis_data.emission_totals[12:]*1e-9),
                                       emis_data.emission_totals[3]*1e-9 + emis_data.emission_totals[9]*1e-9,
                                       total_inc_emis,
                                       total_inc_prop,
                                       emis_data.mass_survive_total * 1e-3]}
    df = pd.DataFrame(data)
    df_rounded = df.round(9)
    print(df_rounded)  
    df_rounded.to_csv(f"./out_files/{(year // 10) * 10}/emis_stats_{year}_{dataset}_{res}_{levels}.csv",sep=',',index=False)    

# Main section of the program
if __name__ == "__main__":

    # Configure the script arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument('-sm', "--start_month", default = "1", choices=str(np.arange(1,13)), help='Start Month (will override final month if greater than final month).')
    parser.add_argument('-fm', "--final_month", default = "12", choices=str(np.arange(1,13)), help='Final Month.')
    parser.add_argument('-sd', "--start_dataset", default = "1", choices=str(np.arange(1,4)), help='Dataset. 1=Non-SMC, 2=SMC, 3=All')
    parser.add_argument('-fd', "--final_dataset", default = "3", choices=str(np.arange(1,4)), help='Dataset. 1=Non-SMC, 2=SMC, 3=All')
    parser.add_argument('-sy', "--start_year", default = "2023", choices=str(np.arange(1957,2025)), help='Start Year.')
    parser.add_argument('-fy', "--final_year", default = "2024", choices=str(np.arange(1957,2025)), help='Final Year.')
    args = parser.parse_args()

    ######################################   
    # Define resolutions and set timings.
    ###################################### 
    
    TIMESTEP  = 3600
    LEVELS    = 72
    MODEL_ALT = 80
    GRID_RES  = "2x25"    
    HOURS     = 24
    print(f"Timestep: {TIMESTEP}s. Resolution: {GRID_RES}x{LEVELS}. Vertical Range: 0-{MODEL_ALT} km.")

    def main():
        # Define the ranges of years, months and datasets to process.
        def make_range(start, end):
            start, end = int(start), int(end)
            if start > end:
                end = start + 1
            return np.arange(start, end + 1), start, end
        
        years,    start_year,    final_year    = make_range(args.start_year, args.final_year)
        months,   start_month,   final_month   = make_range(args.start_month, args.final_month)
        datasets, start_dataset, final_dataset = make_range(args.start_dataset, args.final_dataset)
        
        print(f"Years: {start_year}-{final_year}. Months: {start_month}-{final_month}.")

        #################################  
        # Define launch event altitudes.
        #################################

        # BECO, MECO and SEI values from literature sources are used wherever possible. 
        # When not available, the average of other rockets with the same configuration is used.

        stage_alt_dict = {}
        stage_alt_rockets = np.genfromtxt("./input_files/launch_event_altitudes.csv",dtype=str,skip_header=1,usecols=[0,1],delimiter=",")
        stage_alt_data = np.genfromtxt("./input_files/launch_event_altitudes.csv",dtype=np.float64,skip_header=1,usecols=[2,3,4,5],delimiter=",")

        stages = ["BECO", "MECO", "SEI1", "SECO"]
        for (name, variant), row in zip(stage_alt_rockets, stage_alt_data):
            for stage, value in zip(stages, row):
                stage_alt_dict[f"{name} {variant} {stage}"] = None if value == "" else np.float64(value)
        
        # BECO, MECO, SEI1, SECO
        # B+1/2S, B+3S, B+4S, 2S, 3S, 4S 
        event_alts = [[66,55,29,0,0,0],
                    [220,120,64,90,56,52],
                    [229,120,64,103,61,59],
                    [356,232,216,312,176,149]]
            
        ################
        # Import files.
        ################
        
        fiona.drvsupport.supported_drivers['KML'] = 'rw' # type: ignore
        if start_year < 2020:
            launch_path       = f'./databases/launch_activity_data_1957-2019.nc'
            rocket_info_path  = f'./databases/rocket_attributes_1957-2019.nc'
        elif start_year >= 2020 and final_year <= 2022:
            launch_path       = f'./databases/launch_activity_data_2020-2022.nc'
            rocket_info_path  = f'./databases/rocket_attributes_2020-2022.nc'
        elif start_year >= 2023 and final_year <= 2024:
            launch_path       = f'./databases/launch_activity_data_2023-2024.nc'
            rocket_info_path  = f'./databases/rocket_attributes_2023-2024.nc'
        else: 
            raise ImportError(f"Error: Unsupported time range for {start_year}-{final_year}")
        
        if start_year < 2020:
            reentry_path  = f'./databases/reentry_activity_data_1957-2019.nc'
        elif start_year >= 2020 and final_year <= 2022:
            reentry_path  = f'./databases/reentry_activity_data_2020-2022_moredatacorrectlocations.nc'
        elif start_year >= 2023 and final_year <= 2024:
            reentry_path  = f'./databases/reentry_activity_data_2023-2024.nc'
        else:
            raise ImportError(f"Error: Unsupported time range for {start_year}-{final_year}")

        pei_path          = './input_files/primary_emission_indices.csv'  
        global input_data
        input_data       = InputData(launch_path, reentry_path, rocket_info_path, pei_path, start_year)
        print("Successfully loaded input databases.")
        
        #Loop over all years and run functions depending on input arguments.
        for year in years:
            for dataset in datasets:
                print(f"Year: {year} Dataset: {dataset}")    
                emis_data = OutputEmis(event_alts, months, dataset, year, stage_alt_dict)
                check_total_emissions(year,dataset,GRID_RES,LEVELS,emis_data)

    # Profile the main function
    with cProfile.Profile() as pr:
        main()

    # Save and print profiling stats
    #stats = pstats.Stats(pr)
    #stats.strip_dirs()
    #stats.sort_stats(pstats.SortKey.TIME)
    #stats.print_stats(50)  # Print the top 20 time-consuming functions