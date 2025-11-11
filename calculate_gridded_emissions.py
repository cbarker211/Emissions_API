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
pd.set_option('display.max_colwidth', None)
import geopandas as gpd
import fiona
import json
import calendar
import sys

from python_modules.distribute_emis_func import make_grid_LL, read_gc_box_height, get_ross_profiles, interp_prop_mass
from python_modules.alt_emis_func import calculate_bc_ei, calculate_nox_ei, calculate_co_ei, calculate_cl_ei

class RocketData:
    '''Read rocket launch and re-entry activity data.'''
    def __init__(self, launchfile: str,reentryfile: str, rocketinfofile: str, peifile: str):
        
        self.read_launch_activity(launchfile)
        #self.read_reentry_activity(reentryfile) # TODO: Re-enable when reentry emissions are added.
        self.read_rocket_info(rocketinfofile)
        self.define_pei(peifile)
    
    def read_launch_activity(self, launchfile):
        """Read launch activity from the database.

        Args:
            launchfile (str): The path of the launch database.
        """        

        #Open Dataset.
        with Dataset(launchfile, mode='r') as launchdata:
        
            #Extract variables.
            self.launch_lons    = launchdata.variables['Longitude'][:]
            self.launch_lats    = launchdata.variables['Latitude'][:]
            self.launch_time    = launchdata.variables['Time(UTC)'][:]  # UTC hours
            self.launch_datestr = launchdata.variables['Date'][:]
            self.launch_id      = launchdata.variables['COSPAR_ID'][:].astype(str)
            self.launch_rocket  = launchdata.variables['Rocket_Name'][:]
            try:
                self.launch_variant = launchdata.variables['Rocket_Variant'][:]
            except:
                self.launch_variant = np.array(['-']*len(self.launch_id))
            self.launch_smc     = launchdata.variables['Megaconstellation_Flag'][:]
            self.launch_site    = launchdata.variables['Site'][:]
            
            self.launch_year, self.launch_month, self.launch_day = [], [], []
            for datestring in self.launch_datestr:
                self.launch_year.append(np.int64(datestring[:4]))
                self.launch_month.append(np.int64(datestring[4:6]))
                self.launch_day.append(np.int64(datestring[6:]))
        
    def read_reentry_activity(self, reentryfile):
        """Read re-entry activity from the database.

        Args:
            reentryfile (str): The path of the reentry database.
        """        

        # Open Dataset
        with Dataset(reentryfile, mode='r') as reentrydata:
        
            # Extract variables:
            self.reentry_id           = reentrydata.variables['COSPAR_ID'][:]
            self.reentry_name         = reentrydata.variables['Object_Name'][:]
            self.reentry_category     = reentrydata.variables['Category'][:]
            self.reentry_lats         = reentrydata.variables['Latitude'][:]
            self.reentry_lons         = reentrydata.variables['Longitude'][:]
            self.reentry_datestr      = reentrydata.variables['Date'][:]
            self.reentry_time         = reentrydata.variables['Time (UTC)'][:]
            self.reentry_abl_mass     = reentrydata.variables['Ablatable_Mass'][:]
            self.reentry_abl_deg      = reentrydata.variables['Ablation_Degree'][:]
            self.reentry_per_alu      = reentrydata.variables['Percent_Aluminium'][:]
            self.reentry_other_mass   = reentrydata.variables['Other_Mass'][:]
            self.reentry_smc          = reentrydata.variables['Megaconstellation_Flag'][:]
            self.reentry_location     = reentrydata.variables['Location_Constraint'][:]
            self.reentry_burnup       = reentrydata.variables['Burnup'][:]
            
            self.reentry_year, self.reentry_month, self.reentry_day = [], [], []
            for datestring in self.reentry_datestr:
                self.reentry_year.append(np.int64(datestring[:4]))
                self.reentry_month.append(np.int64(datestring[4:6]))
                self.reentry_day.append(np.int64(datestring[6:]))
              
    def read_rocket_info(self, rocketinfofile):
        """Read information about the propellant mass and type for each rocket stage.
        The stage mass is also included in the database for reentry purposes, but not required here.

        Args:
            rocketinfofile (str): The path of the rocket database.
        """        

        # Open Dataset.
        with Dataset(rocketinfofile, mode='r') as vehicledata:
        
            # Extract variables:
            self.rocket_name         = vehicledata.variables['Rocket_Name'][:]
            try:
                self.rocket_variant = vehicledata.variables['Rocket_Variant'][:]
            except:
                self.rocket_variant = np.array(['-']*len(self.rocket_name))
            try:
                self.booster_prop_mass   = vehicledata.variables['Stage0_PropMass'][:]
                self.booster_prop_type   = vehicledata.variables['Stage0_Fuel_Type'][:]
            except:
                self.booster_prop_mass   = vehicledata.variables['Booster_PropMass'][:]
                self.booster_prop_type   = vehicledata.variables['Booster_Fuel_Type'][:]
            self.stage1_prop_mass    = vehicledata.variables['Stage1_PropMass'][:]
            self.stage1_prop_type    = vehicledata.variables['Stage1_Fuel_Type'][:]
            self.stage2_prop_mass    = vehicledata.variables['Stage2_PropMass'][:]
            self.stage2_prop_type    = vehicledata.variables['Stage2_Fuel_Type'][:]
            self.stage3_prop_mass    = vehicledata.variables['Stage3_PropMass'][:]
            self.stage3_prop_type    = vehicledata.variables['Stage3_Fuel_Type'][:]
            self.stage4_prop_mass    = vehicledata.variables['Stage4_PropMass'][:]
            self.stage4_prop_type    = vehicledata.variables['Stage4_Fuel_Type'][:]
    
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
    def __init__(self, event_alts, reentry_levs, months, dataset, year, ground_landing_list,stage_alt_dict):

        ####################################
        # Initial Setup
        ####################################
        
        # Setup input variables.
        self.event_alts = event_alts
        self.bot_reenter = min(reentry_levs)
        self.top_reenter = max(reentry_levs)
        self.n_reenter_levs = self.top_reenter - self.bot_reenter + 1
        self.include_landing = True
        self.dataset = dataset
        self.year = year
        self.ground_landing_list = ground_landing_list
        self.stage_alt_dict = stage_alt_dict
        self.total_landing_prop = 0
        self.included_prop = 0
        self.model_alt = MODEL_ALT
        self.rocket_index_map = {r: i for i, r in enumerate(rocket_data.rocket_name)}
        vert_filepath = f"geoschem_vertical_grid_{LEVELS}.csv"

        # Define Molecular Weights:
        self.mw_h   = 1.008
        self.mw_h2  = self.mw_h * 2
        self.mw_h2o = self.mw_h2 + 16.00
        events_data = {}

        # Get Ross vertical profiles of propellant burned on a fine grid.
        # The smallest GEOS-Chem grid box has a height of 0.129 km, so setting find grid res to 0.1 km.
        # Setting maximum to 100 km, this is the top of the highest bin in the Ross distribution. 
        self.ross_alt_edge, self.ross_prop_mass, self.ross_cumulative_mass = get_ross_profiles()
        fine_grid_res = 0.1
        self.fine_grid_bot_alt = np.arange(0,100,fine_grid_res)*1e3
        self.fine_grid_mid_alt = np.arange(fine_grid_res/2,100,fine_grid_res)*1e3
        self.fine_grid_top_alt = np.arange(fine_grid_res,100.1,fine_grid_res)*1e3
        self.prop_in_fine_grid, self.fine_grid_mass = interp_prop_mass(self.fine_grid_bot_alt, self.fine_grid_mid_alt, self.fine_grid_top_alt, 
                                                                         self.ross_alt_edge, self.ross_cumulative_mass)

        # Create variables for totals for a future sanity check:
        # BC launch, CO launch, CO2 launch, NOx launch, H2O launch, Al2O3 launch, Cl launch, HCl launch, Cl2 launch 0-8
        # NOx reentry, Al2O3 reentry, BC reentry, Cl reentry, HCl reentry 9-13
        self.emission_totals = np.zeros(14)
        self.prop_consumed = np.zeros((3,len(rocket_data.launch_id)))
        self.total_prop_consumed = np.zeros((len(rocket_data.h2o_pei),3))
        self.launch_count = 0
        self.mass_survive_total = 0
        
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
        
        # Create the arrays for the output file to check prop/species distribution.
        self.csv_count, self.csv_count_2 = 0, 0

        launch_mask = np.zeros(len(rocket_data.launch_year), dtype=bool)
        #reentry_mask = np.zeros(len(rocket_data.reentry_year), dtype=bool) 
        # TODO: Re-enable when reentry emissions are added.
        if self.dataset == 1:
            launch_mask  = (rocket_data.launch_year  == year) & (~(rocket_data.launch_smc.astype(bool)))
            #reentry_mask = (rocket_data.reentry_year == year) & (~(rocket_data.reentry_smc.astype(bool)))
        elif self.dataset == 2:
            launch_mask  = (rocket_data.launch_year  == year) & (rocket_data.launch_smc.astype(bool))
            #reentry_mask = (rocket_data.reentry_year == year) & (rocket_data.reentry_smc.astype(bool))
        elif self.dataset == 3:
            launch_mask  = (rocket_data.launch_year  == year)
            #reentry_mask = (rocket_data.reentry_year == year)

        launch_length  = np.sum(launch_mask)
        reentry_length = 0#np.sum(reentry_mask)
                                
        if launch_length > 0:
            self.output_csv_launch_prop  = np.zeros((launch_length*3,LEVELS))
        if launch_length+reentry_length > 0:
            self.output_csv_emis  = np.zeros(((launch_length+reentry_length)*10,LEVELS))
        
        # Build the ocean landing list for the year in question.
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

        layer_data = np.loadtxt("./input_files/" + vert_filepath ,delimiter=",") 
        # TODO: Will need to put this in the loop later if we want to use meteorological box heights.
        bot_alt_base = layer_data[::-1,0]*1000  # shape (LEVELS,)
        mid_alt_base = layer_data[::-1,1]*1000
        top_alt_base = layer_data[::-1,2]*1000

        # Reshape for broadcasting: (LEVELS,1,1)
        bot_alt_3d = bot_alt_base[:, None, None]
        mid_alt_3d = mid_alt_base[:, None, None]
        top_alt_3d = top_alt_base[:, None, None]

        # Area can also be precomputed
        area_2d = np.ones((self.nlat, self.nlon))

        #################################
        # Loop over each day and month.
        #################################
        for m in months:
            
            # Process month data
            ndays = calendar.monthrange(self.year, m)[1]
            self.strmon=str(m).zfill(2)
            print('MONTH = ', self.strmon)             
            
            #Loop over all days:
            for d in range(ndays):

                # Find launches and reentries.
                launch_mask = (
                    (np.array(rocket_data.launch_month) == m) &
                    (np.array(rocket_data.launch_day)   == d + 1) &
                    (np.array(rocket_data.launch_year)  == self.year) &
                    (~np.isnan(np.array(rocket_data.launch_time)))
                )

                #reentry_mask = ( # TODO: Re-enable when reentry emissions are added.
                #    (np.array(rocket_data.reentry_month) == m) &
                #    (np.array(rocket_data.reentry_day) == d+1) &
                #    (np.array(rocket_data.reentry_year) == self.year) &
                #    (~np.isnan(np.array(rocket_data.reentry_time)))
                #)

                # Apply dataset-specific SMC filter. No extra filter needed for dataset 3 as all launches/reentries included.
                if self.dataset == 1:
                    launch_mask  &= ~(rocket_data.launch_smc.astype(bool))
                    #reentry_mask &= ~(rocket_data.reentry_smc.astype(bool))
                elif self.dataset == 2:
                    launch_mask  &= (rocket_data.launch_smc.astype(bool))
                    #reentry_mask &= (rocket_data.reentry_smc.astype(bool))

                self.strday = str(d+1).zfill(2) # Process day data
                if not np.any(launch_mask): # and not np.any(reentry_mask):
                    continue  # skip to next day
                
                self.bot_alt = bot_alt_3d + np.zeros((LEVELS, self.nlat, self.nlon))
                self.mid_alt = mid_alt_3d + np.zeros((LEVELS, self.nlat, self.nlon))
                self.top_alt = top_alt_3d + np.zeros((LEVELS, self.nlat, self.nlon))
                self.area = area_2d.copy()  # ensure a fresh array per day

                # Initialize:
                self.pmin,self.pmax=np.nan,np.nan
                self.qmin,self.qmax=np.nan,np.nan

                # Get indices
                l_ind = np.where(launch_mask)[0]
                r_ind = []#np.where(reentry_mask)[0]
                
                ########################################################################
                # Call grid_emis function to calculate distribution from launches
                ########################################################################
                daily_launches = []
                if len(l_ind)>0:
                    daily_launches = self.grid_emis(l_ind,
                                   rocket_data.launch_lons[l_ind],
                                   rocket_data.launch_lats[l_ind],
                                   rocket_data.launch_time[l_ind],
                                   'launch',
                                   rocket_data.launch_id[l_ind],
                                   rocket_data.launch_rocket[l_ind],
                                   rocket_data.launch_smc[l_ind],
                                   '',
                                   rocket_data.launch_time[l_ind],
                                   rocket_data.launch_site[l_ind],
                                   '',
                                   )

                ########################################################################
                # Call grid_emis function to calculate distribution from re-entries
                ########################################################################              
                daily_reentries = []
                if len(r_ind)>0:
                    daily_reentries = self.grid_emis(r_ind,
                                   rocket_data.reentry_lons[r_ind],
                                   rocket_data.reentry_lats[r_ind],
                                   rocket_data.reentry_time[r_ind],
                                   'reentry',
                                   rocket_data.reentry_id[r_ind],
                                   rocket_data.reentry_name[r_ind],
                                   rocket_data.reentry_smc[r_ind],
                                   rocket_data.reentry_category[r_ind],
                                   rocket_data.reentry_time[r_ind],
                                   rocket_data.reentry_location[r_ind],
                                   rocket_data.reentry_burnup[r_ind],
                                   )
                
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
        
        for date in events_data:
            for value in events_data[date]:
                if len(events_data[date][value]) > 0:
                    for event in events_data[date][value]:
                        for key in event:
                            if type(event[key]) not in [str, np.float64, bool, int, np.str_, float]:
                                if type(event[key]) == dict:
                                    for species in event[key]:
                                        if type(event[key][species]) not in [np.float64, float]:
                                            print(type(event[key][species]),species)
                                else:
                                    print(type(event[key]), event, key)

        if self.dataset == 3:
            filename = f'./out_files/{(year // 10) * 10}/data_{self.year}.json'    
        else:
            filename = f'./out_files/{(year // 10) * 10}/data_{self.year}_{self.dataset}.json'                            
        with open(filename, 'w') as json_file:
            json.dump(events_data, json_file, indent=4)
             
    def process_launch_event_altitudes(self, valid_index, launch_rocket,launch_id):
        
        # Define the propellant saved for landing for reusable rockets.
        if launch_id in self.ground_landing_list:
            landing_prop_percent = 100.0-5.6-5.6-1.2
        else:
            landing_prop_percent = 100.0-5.6-1.2
                    
        stage_keys = ['BECO', 'MECO', 'SEI1', 'SECO']
        stage_alts = {key: stage_alt_dict.get(f"{launch_rocket} {key}", np.nan) for key in stage_keys}

        stage_alt_beco = stage_alts['BECO']
        stage_alt_meco = stage_alts['MECO']
        stage_alt_sei  = stage_alts['SEI1']
        stage_alt_seco = stage_alts['SECO']
        
        if np.isnan(stage_alt_meco) and np.isnan(stage_alt_sei) and np.isnan(stage_alt_seco):

            stages = [
                rocket_data.booster_prop_type[valid_index] != '',
                rocket_data.stage1_prop_type[valid_index] != '',
                rocket_data.stage2_prop_type[valid_index] != '',
                rocket_data.stage3_prop_type[valid_index] != '',
                rocket_data.stage4_prop_type[valid_index] != ''
            ]

            config_map = {
                # booster, s1, s2, s3, s4
                (True,  True,  False, False, False): 0,
                (True,  True,  True,  False, False): 0,
                (True,  True,  True,  True,  False): 1,
                (True,  True,  True,  True,  True ): 2,
                (False, True,  True,  False, False): 3,
                (False, True,  True,  True,  False): 4,
                (False, True,  True,  True,  True ): 5
            }

            # Look up configuration safely
            stages_tuple = tuple(stages)
            rocket_config_type = config_map.get(stages_tuple)

            if rocket_config_type is None:
                raise IndexError(f"Incorrect rocket configuration for {launch_rocket} (stages={stages_tuple})")
                    
            stage_alt_beco = self.event_alts[0][rocket_config_type]
            stage_alt_meco = self.event_alts[1][rocket_config_type]
            stage_alt_sei  = self.event_alts[2][rocket_config_type]    
            stage_alt_seco = self.event_alts[3][rocket_config_type]

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
        
        # TODO: Add FEI for all air-launched rockets. 
        fei_alt_dict = {
            "GSLV Mk III": 46.534,
            "Pegasus XL": 11.79,
            "LauncherOne": 12.37
        }   
        # Get altitude, default to 0 if rocket not in dictionary
        fei_alt = fei_alt_dict.get(launch_rocket, 0)
        
        if (launch_rocket in ["Falcon 9 v1.2", "Falcon Heavy"] and self.include_landing):
            percent = landing_prop_percent
        else:
            percent = 100.0
        
        ###########
        # Boosters
        ###########
        
        if rocket_data.booster_prop_type[valid_index] != '':
            
            self.booster_alt_index = get_alt_index(stage_alt_beco, self.fine_grid_top_alt, self.fine_grid_bot_alt)

            if stage_alt_beco * 1e3 > self.fine_grid_top_alt[-1]:
                self.fine_grid_mass_booster = np.asarray(self.fine_grid_mass)
            else:
                self.fine_grid_mass_booster = normalize_mass(self.fine_grid_mass[:self.booster_alt_index].copy(), percent)
                         
            if launch_rocket == "Falcon Heavy" and self.include_landing == True:
                self.total_landing_prop += (100.0 - landing_prop_percent) * rocket_data.booster_prop_mass[valid_index] / 100.0

        ###########
        # Stage 1
        ###########  
        # TODO: Needs reworking if running for different model ceilings above 80km.            
        self.fei_alt_index = np.searchsorted(self.fine_grid_top_alt, fei_alt*1e3)
        self.MECO_alt_index = get_alt_index(stage_alt_meco, self.fine_grid_top_alt, self.fine_grid_bot_alt)
        if stage_alt_meco * 1e3 > self.fine_grid_top_alt[-1]:
            self.fine_grid_mass_stage1 = self.fine_grid_mass[self.fei_alt_index:].copy()
            if launch_rocket == "GSLV Mk III":
                fei_mass = np.interp(fei_alt,self.ross_alt_edge, self.ross_cumulative_mass)
                fei_percent = (self.prop_in_fine_grid - fei_mass) / (100.0 - fei_mass) * 100.0
                self.fine_grid_mass_stage1 = normalize_mass(self.fine_grid_mass_stage1, fei_percent)
        else:
            fine_grid_mass_stage1 = self.fine_grid_mass[self.fei_alt_index:self.MECO_alt_index].copy()
            self.fine_grid_mass_stage1 = normalize_mass(fine_grid_mass_stage1, percent)
                        
        if launch_rocket == "Falcon Heavy" and self.include_landing == True:
            self.total_landing_prop += (100.0 - landing_prop_percent) * rocket_data.stage1_prop_mass[valid_index] / 100.0

        ###########
        # Stage 2
        ###########
                               
        # If SEI occurs above the fine grid, then we can ignore the second stage.
        if stage_alt_sei * 1e3 >= self.fine_grid_top_alt[-1]:
            self.fine_grid_mass_stage2 = np.asarray([0])
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
            self.fine_grid_mass_stage2 = normalize_mass(self.fine_grid_mass[self.sei_alt_index:self.seco_alt_index].copy(), seco_percent)
            
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
            self.fine_grid_mass_stage3 = normalize_mass(self.fine_grid_mass[self.TEI_alt_index:].copy(), tei_percent)
         
        ##################
        # Falcon Reusable
        ##################
           
        # Finally create an extra stage just for the landing emissions of Falcon 9 v1.2. 
        # TODO: Needs reworking if running for different model ceilings above 80km.
        # Would need to add boostback emissions if top is above 80.
        if launch_rocket in ["Falcon 9 v1.2","Falcon Heavy"]:

            # 5.6% are used in the entry burn, over 70-54.7 km.
            self.entry_top = np.searchsorted(self.fine_grid_top_alt, 70000, side='right')
            self.entry_bot = np.searchsorted(self.fine_grid_top_alt, 54700, side='right')
            self.fine_grid_mass_entry = normalize_mass(self.fine_grid_mass[self.entry_bot:self.entry_top+1].copy(), 5.6)
            
            # 1.2% are used in the landing burn, over 3.3-0 km.
            self.landing_top = np.searchsorted(self.fine_grid_top_alt, 3300, side='right')
            self.fine_grid_mass_landing = normalize_mass(self.fine_grid_mass[:self.landing_top+1].copy(), 1.2)

        else:
            self.fine_grid_mass_landing = 0
            self.fine_grid_mass_entry = 0
        
        return stage_alt_beco, stage_alt_meco
    
    def calc_emis(self,start_ind,stop_ind,pei_index,prop_mass,vertical_profile,time_index, q, p, total_vertical_propellant,stage):
        '''
        Calculate the emissions over the range of the stag within the fine grid (0-100 km).
        '''
        
        # Calculate the emission indices for each species.
        ei_bc                 = calculate_bc_ei (self.fine_grid_mid_alt[start_ind:stop_ind]*1e-3, rocket_data.bc_pei[pei_index])
        ei_co, ei_co2         = calculate_co_ei (self.fine_grid_mid_alt[start_ind:stop_ind]*1e-3, rocket_data.co_pei[pei_index], rocket_data.co2_pei[pei_index])
        ei_sec_nox            = calculate_nox_ei(self.fine_grid_mid_alt[start_ind:stop_ind]*1e-3)
        ei_cl, ei_hcl, ei_cl2 = calculate_cl_ei (self.fine_grid_mid_alt[start_ind:stop_ind]*1e-3, rocket_data.cly_pei[pei_index])
        ei_h2o                = rocket_data.h2o_pei[pei_index] + rocket_data.h2_pei[pei_index] * self.mw_h2o / self.mw_h2 

        # Calculate the emissions for each species in g.
        emis_full = np.zeros((len(self.fine_grid_mid_alt),11))
        for index, emission_index in enumerate([ei_bc, ei_co, ei_co2, ei_sec_nox, rocket_data.nox_pei[pei_index], ei_h2o, rocket_data.al2o3_pei[pei_index], ei_cl, ei_hcl, ei_cl2]):
            emis_full[start_ind:stop_ind,index] = emission_index * prop_mass * 1e-2 * vertical_profile
        emis_full[start_ind:stop_ind,10] = vertical_profile

        # Place the emissions into a larger array covering the whole fine grid (0-100km).
        vertical_propellant = np.zeros(len(self.mid_alt[:,q,p]))
        selected_alts = []
        # Regrid the vertical_profile to the desired model profile.
        for i, alt in enumerate(self.mid_alt[:,q,p]):
            if i == 0:
                bot_ind = 0
            else:
                bot_ind = np.argmin(np.abs(self.fine_grid_mid_alt - self.bot_alt[i,q,p])) + 1
            top_ind = np.argmin(np.abs(self.fine_grid_mid_alt - self.top_alt[i,q,p])) + 1
            selected_alts.extend(np.arange(bot_ind,top_ind))
            
            # Sum the emissions in this range for each species and place in an array.
            for i, key in enumerate(["launch_bc","co","co2","launch_nox","fuel_nox","h2o","launch_al","launch_cl","launch_hcl","cl2"]):
                self.rocket_data_arrays[key][time_index,i,q,p]  += (np.sum(emis_full[bot_ind:top_ind,i]) * 1e-6)
            total_vertical_propellant[i,0] += np.sum(emis_full[bot_ind:top_ind,10]) * prop_mass
            total_vertical_propellant[i,1] += np.sum(emis_full[bot_ind:top_ind,0]) 
            total_vertical_propellant[i,2] += np.sum(emis_full[bot_ind:top_ind,1]) 
            total_vertical_propellant[i,3] += np.sum(emis_full[bot_ind:top_ind,2])
            total_vertical_propellant[i,4] += np.sum(emis_full[bot_ind:top_ind,3:5]) 
            total_vertical_propellant[i,5] += np.sum(emis_full[bot_ind:top_ind,5]) 
            total_vertical_propellant[i,6] += np.sum(emis_full[bot_ind:top_ind,6]) 
            total_vertical_propellant[i,7] += np.sum(emis_full[bot_ind:top_ind,7]) 
            total_vertical_propellant[i,8] += np.sum(emis_full[bot_ind:top_ind,8]) 
            total_vertical_propellant[i,9] += np.sum(emis_full[bot_ind:top_ind,9]) 
             
        if len(list(set(selected_alts))) != len(selected_alts):
            raise IndexError("Error in fine grid indexing.")
        
        # BC launch, CO launch, CO2 launch, NOx launch, H2O launch, Al2O3 launch, Cl launch, HCl launch, Cl2 launch
        # NOx reentry, Al2O3 reentry, BC reentry, HCl reentry, Cl reentry
        self.emission_totals[0] += np.sum(emis_full[selected_alts[0]:selected_alts[-1]+1,0])       
        self.emission_totals[1] += np.sum(emis_full[selected_alts[0]:selected_alts[-1]+1,1])         
        self.emission_totals[2] += np.sum(emis_full[selected_alts[0]:selected_alts[-1]+1,2])   
        self.emission_totals[3] += np.sum(emis_full[selected_alts[0]:selected_alts[-1]+1,3:5])     
        self.emission_totals[4] += np.sum(emis_full[selected_alts[0]:selected_alts[-1]+1,5])   
        self.emission_totals[5] += np.sum(emis_full[selected_alts[0]:selected_alts[-1]+1,6])   
        self.emission_totals[6] += np.sum(emis_full[selected_alts[0]:selected_alts[-1]+1,7])   
        self.emission_totals[7] += np.sum(emis_full[selected_alts[0]:selected_alts[-1]+1,8])   
        self.emission_totals[8] += np.sum(emis_full[selected_alts[0]:selected_alts[-1]+1,9])   
        self.included_prop      += (np.sum(emis_full[selected_alts[0]:selected_alts[-1]+1,10]) * prop_mass * 1e-2)
        
        return total_vertical_propellant        
                  
    def grid_emis(self,index,lon,lat,hour,emis_type,event_id,name,smc,category,time,location, burnup):
        """Grid the data onto the GEOS-Chem horizontal and vertical grid"""
        
        daily_info = []                
        #Loop over each launch/reentry.
        for w in range(len(lon)):

            # Set up the emission arrays for this launch/reentry.
            self.rocket_data_arrays = {}
            species = ['launch_nox', 'fuel_nox',  'h2o',       'launch_bc',  'co',
                       'co2',        'launch_al', 'launch_hcl','launch_cl',  'cl2',
                       'reentry_nox','reentry_al','reentry_bc','reentry_hcl','reentry_cl'
            ]

            for sp in species:
                self.rocket_data_arrays[sp] = np.zeros((HOURS, LEVELS, self.nlat, self.nlon))
            
            ################
            # Setup grid.
            ################
            # Get grid horizontal indices (p is lon index; q is lat index).
            # This works out the nearest latitude and longitude on the grid.
            p,q = np.argmin(np.abs(lon[w]-self.lon)), np.argmin(np.abs(lat[w]-self.lat))
            
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
            time_index = int(hour[w]*60*60)//TIMESTEP
            
            #########################################################
            # Work out vertical distribution for launch emissions.
            #########################################################

            if emis_type=='launch':
                
                ############################################
                # Interpolate propellant mass distribution.
                ############################################

                # Linearly interpolate Ross profile proportion of propellant mass to GEOS-Chem vertical grid.
                # This is just to double check that the propellant consumption profile is calculated correctly.
                self.propellant_in_model, self.gc_relative_mass = interp_prop_mass(self.bot_alt[:,q,p], self.mid_alt[:,q,p], self.top_alt[:,q,p],
                                                                                    self.ross_alt_edge, self.ross_cumulative_mass)

                ############################################
                # Check the rocket type 
                ############################################
                
                valid_index = 0
                if name[w] in self.rocket_index_map:
                    valid_index = self.rocket_index_map[name[w]]
                else:
                    raise IndexError(f"No propellant mass found for {name[w]}")
                    
                launch_details = {
                    "date": f"{self.year}-{self.strmon}-{self.strday}",
                    "id": event_id[w],
                    "time": time[w],
                    "site": "",
                    "rocket": name[w],
                    "lat": lat[w],
                    "lon": lon[w],
                    "smc": bool(smc[w]),
                    "location": location[w],
                }
                
                if location[w] == "Naro Space Center":
                    launch_details["lat"] = 34.43194444
                    launch_details["lon"] = 127.535
                elif location[w] == "China Sea Launch":
                    launch_details["lat"] = 34.9
                    launch_details["lon"] = 121.2
                    
                ############################################################################
                # Process the launch event altitudes.
                ############################################################################
                
                if pd.isna(event_id[w]) or pd.isna(lon[w]) or pd.isna(lat[w]):
                    raise ValueError(f"Event ID: {event_id[w]} Latitude: {lat[w]}, Longitude: {lon[w]}, Site: {location[w]}")
                                                    
                # Next find index of fuel type in primary emissions index data:
                pei_booster_index  = np.where( rocket_data.pei_fuel_type == rocket_data.booster_prop_type[valid_index])[0]
                pei_stage1_index   = np.where( rocket_data.pei_fuel_type == rocket_data.stage1_prop_type[valid_index] )[0]
                pei_stage2_index   = np.where( rocket_data.pei_fuel_type == rocket_data.stage2_prop_type[valid_index] )[0]
                pei_stage3_index   = np.where( rocket_data.pei_fuel_type == rocket_data.stage3_prop_type[valid_index] )[0]
                
                stage_alt_beco, stage_alt_meco = self.process_launch_event_altitudes(valid_index, name[w], event_id[w])
                # TODO: Needs reworking if running for different model ceilings above 80km.
                # Most failures are for upper stages, and so can be treated as normal here.
                # Full information is provided in source_info/failed_launch_info.txt.

                ############################################################################
                # Deal with failed launches.
                ############################################################################
                         
                # TODO: Failed launches back to 1957.
                # This is where the launch failed close to the launch pad.
                if event_id[w] in ['2020-F04','2021-F04','2023-F02','2024-F01'] or (event_id[w][5] == "F" and self.year < 2020):
                    self.csv_count += 3
                    self.csv_count_2 += 10
                    continue
                
                # This is where stage 1 failed during ascent.
                failed_alt_events = {
                    '2020-F07': 900,    # Astra Rocket 3 launch, rocket shut off at 0.9km.
                    '2021-F01': 10700,  # Shuang Quxian-1 launch, rocket disintegrated at Max-Q. Approximating altitude using Proton-M and Minotaur-IV max-q alts.
                    '2021-F07': 31000   # Astra Rocket 3 launch, rocket shut off at 31km.
                }

                if event_id[w] in failed_alt_events:
                    failed_alt = failed_alt_events[event_id[w]]
                    cutoff_ind = np.searchsorted(self.fine_grid_mid_alt, failed_alt, side='right')
                    self.fine_grid_mass_stage1[cutoff_ind:] = 0.0 # Zero above the cutoff_ind.

                # This is where stage 2 never ignited.
                zero_stage2_events = {
                    '2020-F02','2020-F05',
                    '2021-F02','2021-F07','2021-F08','2022-F01',
                    '2022-F02','2022-F03','2023-F01','2023-F04',
                    '2023-F05','2023-F06','2023-F07','2023-F09',
                    '2024-F02','2024-F04'
                }

                if event_id[w] in zero_stage2_events:
                    self.fine_grid_mass_stage2 = np.asarray([0])
                    self.sei_alt_index = None
                
                ##################################################
                # Sanity Checks for Propellant Mass Distributions
                ##################################################
                
                # Each time, use 0.001% as a the maximum error from rounding / floating point errors.
                error_lim = 0.001
                
                # For all rockets with boosters, the propellant consumed for the boosters should never be bigger than the total booster propellant.
                # Suppressed for Falcon Heavy, as this has reusable boosters. 
                with np.errstate(divide='ignore', invalid='ignore'):
                    if rocket_data.booster_prop_type[valid_index] != '':
                        self.prop_consumed[0,self.launch_count] = np.sum(self.fine_grid_mass_booster * rocket_data.booster_prop_mass[valid_index] * 1e-2)
                        if stage_alt_beco < self.model_alt:
                            if ((np.abs(self.prop_consumed[0,self.launch_count] - rocket_data.booster_prop_mass[valid_index]) / rocket_data.booster_prop_mass[valid_index] * 100.0) > error_lim) and (rocket_data.rocket_name[valid_index] != "Falcon Heavy"):
                                print(np.abs(self.prop_consumed[0,self.launch_count]),rocket_data.booster_prop_mass[valid_index],np.abs(self.prop_consumed[0,self.launch_count] - rocket_data.booster_prop_mass[valid_index]) / rocket_data.booster_prop_mass[valid_index] * 100.0)
                                raise ValueError(f"Error with booster emissions -2. {event_id[w]}") 
                            
                    self.prop_consumed[1,self.launch_count]  = np.sum(self.fine_grid_mass_stage1 * rocket_data.stage1_prop_mass[valid_index] * 1e-2)
                    self.prop_consumed[2,self.launch_count]  = np.sum(self.fine_grid_mass_stage2 * rocket_data.stage2_prop_mass[valid_index] * 1e-2)

                    # For all rockets, the propellant consumed for each stage should never be bigger than the total propellant in each stage.
                    if ((self.prop_consumed[1,self.launch_count] - rocket_data.stage1_prop_mass[valid_index]) / rocket_data.stage1_prop_mass[valid_index] * 100.0 ) > error_lim:
                        raise ValueError(f"Error with Stage 1 emissions - 1.")
                    if (rocket_data.rocket_name[valid_index] != "Long March (CZ) 5B"):
                        if (((self.prop_consumed[2,self.launch_count] - rocket_data.stage2_prop_mass[valid_index]) / rocket_data.stage2_prop_mass[valid_index] * 100.0 ) > error_lim):
                            raise ValueError("Error with Stage 2 emissions.")

                    # When MECO occurs in the model, the consumed propellant should be within 1% of the total propellant mass of stage 1.  
                    # The error is suppressed for 2020-F07, 2021-F01, and 2021-F07 where rockets failed early. 
                    # Also suppressed for Falcon 9, as this has a reusable first stage.
                    # A check that the Falcon landing distribution is no greater than 7% of the stage 1 emissions is undertaken later in the grid_emis function.
                    if stage_alt_meco < self.model_alt:
                        if (event_id[w] not in ['2020-F07','2021-F01','2021-F07']) and rocket_data.rocket_name[valid_index] != "Falcon 9 v1.2" and ((np.abs(self.prop_consumed[1,self.launch_count] - rocket_data.stage1_prop_mass[valid_index]) / rocket_data.stage1_prop_mass[valid_index] * 100) > 1):
                            print(event_id[w],self.prop_consumed[1,self.launch_count],rocket_data.stage1_prop_mass[valid_index],(np.abs(self.prop_consumed[1,self.launch_count] - rocket_data.stage1_prop_mass[valid_index]) / rocket_data.stage1_prop_mass[valid_index] * 100.0))
                            raise ValueError("Error with Stage 1 emissions - 2.")  
                        
                self.total_prop_consumed[pei_booster_index,0] += self.prop_consumed[0,self.launch_count]
                self.total_prop_consumed[pei_stage1_index,1]  += self.prop_consumed[1,self.launch_count]
                self.total_prop_consumed[pei_stage2_index,2]  += self.prop_consumed[2,self.launch_count]                  

                total_prop_mass = rocket_data.booster_prop_mass[valid_index] + rocket_data.stage1_prop_mass[valid_index] + rocket_data.stage2_prop_mass[valid_index]
                total_prop_mass += rocket_data.stage3_prop_mass[valid_index] + rocket_data.stage4_prop_mass[valid_index]
                
                ##############################################
                # Calculate the emissions for each species.
                ##############################################     
                
                # Creating a 2d array for the prop output.
                # Total prop, bc, co, co2, nox, h2o, al2o3, cl, hcl, cl2
                total_vertical_propellant = np.zeros((len(self.mid_alt[:,q,p]),10))    
                                
                # Check whether there is a booster:
                if (rocket_data.booster_prop_type[valid_index] != ''):
                    if np.sum(1e-2 * self.fine_grid_mass_booster) > 1.01:
                        raise ValueError("Error with Boosters. Propellant distribution exceeds unity.")
                    total_vertical_propellant = self.calc_emis(None,
                                self.booster_alt_index,
                                pei_booster_index,
                                rocket_data.booster_prop_mass[valid_index],
                                self.fine_grid_mass_booster,
                                time_index, 
                                q, 
                                p, 
                                total_vertical_propellant,
                                0)
                        
                # Every rocket has a first stage.
                if np.sum(1e-2 * self.fine_grid_mass_stage1) > 1.01:
                    raise ValueError("Error with Stage 1. Propellant distribution exceeds unity.")  
                total_vertical_propellant = self.calc_emis(self.fei_alt_index,
                                                           self.MECO_alt_index,
                                                           pei_stage1_index,
                                                           rocket_data.stage1_prop_mass[valid_index],
                                                           self.fine_grid_mass_stage1,
                                                           time_index, 
                                                           q, 
                                                           p, 
                                                           total_vertical_propellant,
                                                           1
                                                           )

                # Check whether there is a second stage:
                if rocket_data.stage2_prop_type[valid_index] != '' and self.sei_alt_index != None:
                    if np.sum(1e-2 * self.fine_grid_mass_stage2) > 1.01:
                        raise ValueError("Error with Stage 2. Propellant distribution exceeds unity.")
                    total_vertical_propellant = self.calc_emis(self.sei_alt_index,
                                self.seco_alt_index,
                                pei_stage2_index,
                                rocket_data.stage2_prop_mass[valid_index],
                                self.fine_grid_mass_stage2,
                                time_index, 
                                q, 
                                p, 
                                total_vertical_propellant,
                                2)
                    
                # NOTE: This section needs to be tweaked if wanting to run for different vertical heights above 80km.
                # Check for more rockets that have third stage emissions within model.    
                # Add third stage emissions for Minotaur 1.
                if rocket_data.rocket_name[valid_index] == "Minotaur 1":
                    if np.sum(1e-2 * self.fine_grid_mass_stage3) > 1.01:
                        raise ValueError("Error with Stage 3 for Minotaur 1. Propellant distribution exceeds unity.")
                    total_vertical_propellant = self.calc_emis(self.TEI_alt_index,
                                None,
                                pei_stage3_index,
                                rocket_data.stage3_prop_mass[valid_index],
                                self.fine_grid_mass_stage3,
                                time_index, 
                                q, 
                                p, 
                                total_vertical_propellant,
                                3)
                    
                # If the rocket is a Falcon 9, then add the landing emissions. Kerosene, so no Al2O3 or Cly.     
                
                if rocket_data.rocket_name[valid_index] in ["Falcon 9 v1.2","Falcon Heavy"] and self.include_landing == True:
                    falcon_p, falcon_q = None, None

                    if event_id[w] in self.ground_landing_list:
                        falcon_p = p
                        falcon_q = q
                        
                    else:
                        matching = self.ocean_landings[self.ocean_landings["Date"].isin([f"{self.year}-{self.strmon}-{self.strday}"])].reset_index(drop=True)
                        # There is a typo in the database, a 2023 launch is listed as 2022.      
                        if (matching.shape[0] == 1) or (f"{self.year}-{self.strmon}-{self.strday}" == "2022-04-27"):

                            falcon_lat = np.array(matching["geometry"].y)[0]
                            falcon_lon = np.array(matching["geometry"].x)[0]
                            
                        # There are two launches on the same day.
                        elif event_id[w] in ["2022-124","2023-037"]:
                            falcon_lat = np.array(matching["geometry"].y)[0]
                            falcon_lon = np.array(matching["geometry"].x)[0]
                            
                        elif event_id[w] in ["2022-125","2023-038"]:
                            falcon_lat = np.array(matching["geometry"].y)[1]
                            falcon_lon = np.array(matching["geometry"].x)[1]
                            
                        elif matching.shape[0] == 0:
                            # The database hasn't been well updated for 2022, so lets just fill in based on most common geolocation for all other 2020-2022 launches.
                            if lon[w] == -81.0 and lat[w] == 28.5:
                                falcon_lon = -75
                                falcon_lat = 32
                            elif lon[w] == -120.6 and lat[w] == 34.7:
                                falcon_lon = -122.5
                                falcon_lat = 30
                            elif lon[w] == 100.3 and lat[w] == 41.3:
                                falcon_lon = 100.3
                                falcon_lat = 41.3
                            else:
                                raise RuntimeError("Launch not from assigned site.")
                        
                        elif matching.shape[0] > 1: 
                            raise RuntimeError(f"Multiple ocean entries for Falcon Stage 1 landing for {event_id[w]}.")
                        else:
                            raise RuntimeError(f"Problem geolocating Falcon Stage 1 landing for {event_id[w]}, {self.lon[falcon_p],self.lat[falcon_q]}") 

                        falcon_p = np.argmin(abs(falcon_lon-self.lon))
                        falcon_q = np.argmin(abs(falcon_lat-self.lat))
                            
                        if falcon_p > self.pmax:
                            self.pmax = falcon_p
                            #print(f"Falcon9 stage 1 reentry longitude out of bounds for {event_id[w]}. Updating bounds.")
                        if falcon_q > self.qmax:
                            self.qmax = falcon_q
                            #print(f"Falcon9 stage 1 reentry latitude out of bounds for {event_id[w]}. Updating bounds.")                  
                    
                    if np.sum(self.fine_grid_mass_entry) + np.sum(self.fine_grid_mass_landing) > 7:
                        raise RuntimeError("Error with Stage 1 Ocean Landing. Propellant distribution exceeds what is expected.")
                        
                    if rocket_data.rocket_name[valid_index] == "Falcon 9 v1.2":
                        pei = pei_stage1_index
                        prop_mass = rocket_data.stage1_prop_mass[valid_index]
                    elif rocket_data.rocket_name[valid_index] == "Falcon Heavy":
                        pei = pei_booster_index
                        prop_mass = rocket_data.booster_prop_mass[valid_index]
                    else:
                        pei = None
                        prop_mass = None
                        raise RuntimeError("Error identifying Falcon rocket type.")
                    
                    # NOTE: This section needs to be tweaked if wanting to run for different vertical heights above 80km.
                    # Should implement the boostback burn if the vertical height is increased to 100 km.
                    # First the entry emissions.
                    total_vertical_propellant = self.calc_emis(self.entry_bot,
                                self.entry_top+1,
                                pei,
                                prop_mass,
                                self.fine_grid_mass_entry,
                                time_index, 
                                falcon_q, 
                                falcon_p, 
                                total_vertical_propellant,
                                5)
                    # Now the landing emissions.
                    total_vertical_propellant = self.calc_emis(None,
                                self.landing_top+1,
                                pei,
                                prop_mass,
                                self.fine_grid_mass_landing,
                                time_index, 
                                falcon_q, 
                                falcon_p, 
                                total_vertical_propellant,
                                5)
                    
                launch_details["emissions"] = {
                    "BC":    np.sum(self.rocket_data_arrays["launch_bc"]),
                    "CO":    np.sum(self.rocket_data_arrays["co"]),
                    "CO2":   np.sum(self.rocket_data_arrays["co2"]),
                    "NOx":   np.sum(self.rocket_data_arrays["launch_nox"]) + np.sum(self.rocket_data_arrays["fuel_nox"]),
                    "H2O":   np.sum(self.rocket_data_arrays["h2o"]),
                    "Cly":   np.sum(self.rocket_data_arrays["launch_cl"]) + np.sum(self.rocket_data_arrays["cl2"]) + np.sum(self.rocket_data_arrays["launch_hcl"]),
                    "Al2O3": np.sum(self.rocket_data_arrays["launch_al"]),
                    }   
                
                daily_info.append(launch_details) 
                
                ##############################################
                # Output the emissions to a file for viewing.
                ##############################################                  
                #with np.errstate(divide='ignore', invalid='ignore'):
                #    total_vertical_propellant[:,0] = total_vertical_propellant[:,0] / total_prop_mass
                self.output_csv_launch_prop[self.csv_count,:] = total_vertical_propellant[:,0]
                self.csv_count += 1
                self.output_csv_launch_prop[self.csv_count,:] = self.mid_alt[:,q,p]*1e-3
                self.csv_count += 1  
                self.output_csv_launch_prop[self.csv_count,:] = self.gc_relative_mass
                self.csv_count += 1 
                
                if np.sum(self.gc_relative_mass) == 0 and np.sum(total_vertical_propellant[:,0]) == 0:
                    print(f"??? {rocket_data.rocket_name[valid_index]}")
                
                for i in range(1,10):
                    self.output_csv_emis[self.csv_count_2,:] = total_vertical_propellant[:,i]
                    self.csv_count_2 += 1 
                self.output_csv_emis[self.csv_count_2,:] = self.mid_alt[:,q,p]*1e-3
                self.csv_count_2 += 1

                self.launch_count += 1
                
            #########################################################
            # Work out vertical distribution for reentry emissions.
            #########################################################
            
            if emis_type=='reentry':
                reentry_details = {
                    "date": f"{self.year}-{self.strmon}-{self.strday}",
                    "id": event_id[w],
                    "time": time[w],
                    "reusability": "",
                    "name": name[w],
                    "category": category[w], 
                    "lat": lat[w],
                    "lon": lon[w],
                    "smc": bool(smc[w]),
                    "location": int(location[w]),
                    "burnup": burnup[w], 
                }
                
                if event_id[w][:8] in ["2021-F09","2022-065","2023-72"] and category[w] == "S1":
                    reentry_details["lat"] = 34.43194444
                    reentry_details["lon"] = 127.535
                elif event_id[w][:8] in ["2020-065","2022-167","2022-046","2022-126","2023-135","2024-102","2024-153","2024-173","2024-245","2025-007","2025-105"] and category[w] == "S1":
                    reentry_details["lat"] = 34.9
                    reentry_details["lon"] = 121.2
                else:
                    reentry_details["lat"] = lat[w]
                    reentry_details["lon"] = lon[w]
                
                total_vertical_propellant = np.zeros((len(self.mid_alt[:,q,p]),10))
                reentry_ei = np.zeros(5) # Al2O3, NOx, BC, Cl, HCl
                if np.ma.is_masked(rocket_data.reentry_abl_mass[index[w]]) or np.ma.is_masked(rocket_data.reentry_other_mass[index[w]]):
                    pass
                else:
                
                    reentry_ei[0]   = rocket_data.reentry_abl_deg[index[w]] * rocket_data.reentry_per_alu[index[w]]
                    if rocket_data.reentry_abl_deg[index[w]] == 0:
                        reentry_ei[1] = 0.175
                    else:
                        reentry_ei[1] = 0.4
                        
                    # Calculate the total mass surviving re-entry in tonnes
                    # We are adding all mass that comes back, so this includes reusable stages and fairings.
                    self.mass_survive = ((rocket_data.reentry_abl_mass[index[w]]* (1-rocket_data.reentry_abl_deg[index[w]])) + rocket_data.reentry_other_mass[index[w]]) * 1e-3
                     
                    # Add chlorine and bc reentry emissions from ATISPADE Report. Don't add for lower stages, only upper (>S2).
                    if rocket_data.reentry_category[index[w]] in ["P","C"]:
                        reentry_ei[2:] = [0.041,0.015,0.008] # Worst Case
                    elif rocket_data.reentry_category[index[w]] in ["S2","S3","S4"]:
                        reentry_ei[2:] = [0.029,0.011,0.005] # Worst Case
                     
                    # For consistency with launch emissions, the totals are kept in g units. 
                    # NOx reentry, Al2O3 reentry, BC reentry, Cl reentry, HCl reentry 9-13           
                    t_nox_reentry = (rocket_data.reentry_abl_mass[index[w]] + rocket_data.reentry_other_mass[index[w]]) * reentry_ei[1] * 1000
                    self.emission_totals[9] += t_nox_reentry
                    self.rocket_data_arrays["reentry_nox"][time_index,self.bot_reenter:self.top_reenter+1,q,p] += t_nox_reentry / self.n_reenter_levs * 1e-6 
                    
                    t_al2o3_reentry = rocket_data.reentry_abl_mass[index[w]] * reentry_ei[0] * 1000
                    self.emission_totals[10] += t_al2o3_reentry
                    self.rocket_data_arrays["reentry_al"][time_index,self.bot_reenter:self.top_reenter+1,q,p] += t_al2o3_reentry / self.n_reenter_levs * 1e-6

                    total_vertical_propellant[self.bot_reenter:self.top_reenter+1,4]   += np.full((self.n_reenter_levs),t_nox_reentry/self.n_reenter_levs)
                    total_vertical_propellant[self.bot_reenter:self.top_reenter+1,6]   += np.full((self.n_reenter_levs),t_al2o3_reentry/self.n_reenter_levs)
                    
                    if rocket_data.reentry_category[index[w]] in ["P","C","S2","S3","S4"]:
                        
                        prop_dict = {"reentry_bc": 1,"reentry_cl":7,"reentry_hcl": 8}
                        for i, key in enumerate(["reentry_bc","reentry_cl","reentry_hcl"]):
                            t_reentry = (rocket_data.reentry_abl_mass[index[w]] + rocket_data.reentry_other_mass[index[w]]) * reentry_ei[i+2] * 1000
                            self.emission_totals[i+11] += t_reentry
                            self.rocket_data_arrays[key][time_index,self.bot_reenter:self.top_reenter+1,q,p] += t_reentry / self.n_reenter_levs * 1e-6
                            total_vertical_propellant[self.bot_reenter:self.top_reenter+1,prop_dict[key]]   += np.full((self.n_reenter_levs),t_reentry/self.n_reenter_levs)
                    
                reentry_details["emissions"] = {"NOx":   np.sum(self.rocket_data_arrays["reentry_nox"]),
                                                "Al2O3": np.sum(self.rocket_data_arrays["reentry_al"]),
                                                "BC":    np.sum(self.rocket_data_arrays["reentry_bc"]),
                                                "HCl":   np.sum(self.rocket_data_arrays["reentry_hcl"]),
                                                "Cl":    np.sum(self.rocket_data_arrays["reentry_cl"]),
                                                "Unablated_Mass": self.mass_survive,
                             }   
                
                self.mass_survive_total += self.mass_survive
                daily_info.append(reentry_details) 

                for i in range(1,10):
                    self.output_csv_emis[self.csv_count_2,:] = total_vertical_propellant[:,i]
                    self.csv_count_2 += 1
                self.output_csv_emis[self.csv_count_2,:] = self.mid_alt[:,q,p]*1e-3
                self.csv_count_2 += 1
        
        return daily_info

def check_total_emissions(year,dataset,res,levels):
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
    print(df.round(4))  
    df.to_csv(f"./out_files/{(year // 10) * 10}/emis_stats_{year}_{dataset}_{res}_{levels}.csv",sep=',')    

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
      
    # Define launch event altitudes.

    # BECO, MECO and SEI values from literature sources are used wherever possible. 
    # When not available, the average of other rockets with the same configuration is used.

    stage_alt_dict = {}
    stage_alt_rockets = np.genfromtxt("./input_files/launch_event_altitudes.csv",dtype=str,skip_header=1,usecols=[0],delimiter=",")
    stage_alt_data = np.genfromtxt("./input_files/launch_event_altitudes.csv",dtype=np.float64,skip_header=1,usecols=[1,2,3,4],delimiter=",")

    stages = ["BECO", "MECO", "SEI1", "SECO"]
    for rocket, row in zip(stage_alt_rockets, stage_alt_data):
        for stage, value in zip(stages, row):
            stage_alt_dict[f"{rocket} {stage}"] = None if value == "" else np.float64(value)
    
    # BECO, MECO, SEI1, SECO
    # B+1/2S, B+3S, B+4S, 2S, 3S, 4S 
    event_alts = [[66,55,29,0,0,0],
                  [220,120,64,90,56,52],
                  [229,120,64,103,61,59],
                  [356,232,216,312,176,149]]

    # TODO: Check for Falcon Heavy expendable/reusable lower stages.
    ground_landing_list = ["2020-016","2020-059","2020-086","2020-101","2021-059","2022-002","2022-008","2022-009",
                           "2022-040","2022-057","2022-063","2022-144","2022-166","2022-168","2022-173","2022-179",
                           "2023-001","2023-004","2023-008","2023-029","2023-050","2023-054","2023-070","2023-084",
                           "2023-108","2023-128","2023-133","2023-157","2023-173","2023-174","2023-185","2023-204",
                           "2023-210","2024-003","2024-014","2024-021","2024-025","2024-028","2024-030","2024-042",
                           "2024-043","2024-054","2024-066","2024-070","2024-081","2024-101","2024-119","2024-139",
                           "2024-146","2024-149","2024-163","2024-178","2024-188","2024-200","2024-206","2024-247"]
    
    ######################################   
    # Define resolutions and set timings.
    ###################################### 
    
    TIMESTEP  = 3600
    LEVELS    = 72
    MODEL_ALT = 80
    GRID_RES  = "2x25"    
    HOURS     = 24    
    print(f"Timestep: {TIMESTEP}s. Resolution: {GRID_RES}x{LEVELS}. Vertical Range: 0-{MODEL_ALT} km.")
    
    # Re-entry emissions are placed at ~68 km (reentry burnup location).
    # The emissions are evenly distributed across two vertical layers in the model (equivalent to 60-80 km).
    # TODO: Should make it so this automatically places from 60-80 km.
    # This should probably be specific to the box heights of the day.
    if LEVELS==47:
        reentry_levs = [45,46] # 46-47 for the 47-layer model
    elif LEVELS==72:
        reentry_levs = [64,71] # 65-72 for the 72-layer model
    else:
        raise ValueError("Invalid number of levels.")
        
    ################
    # Import files.
    ################
    
    fiona.drvsupport.supported_drivers['KML'] = 'rw' # type: ignore
    raul_data = gpd.read_file('./databases/reentry/General_SpaceX_Map_Raul.kml', driver='KML', layer =2) # Falcon landing data.  
    if start_year == 1957:
        source = "_jsr"
        launch_path       = f'./databases/launch_activity_data_1957-2019_jsr.nc'
        rocket_info_path  = f'./databases/rocket_attributes_1957-2019_jsr.nc'
    else: 
        launch_path       = f'./databases/launch_activity_data_{start_year}-{final_year}.nc'
        rocket_info_path  = f'./databases/rocket_attributes_{start_year}-{final_year}.nc'
    

    # TODO: Removing re-entries for now, re-enable when data is fixed.
    #if start_year == 2020:
    #    reentry_path  = f'./databases/reentry_activity_data_{start_year}-{final_year}_moredatacorrectlocations.nc'
    #else:
    #    reentry_path  = f'./databases/reentry_activity_data_{start_year}-{final_year}.nc'
    reentry_path = ""

    pei_path          = './input_files/primary_emission_indices.csv'  
    rocket_data       = RocketData(launch_path, reentry_path, rocket_info_path, pei_path)
    print("Successfully loaded input databases.")
    
    #Loop over all years and run functions depending on input arguments.
    for year in years:
        for dataset in datasets:
            print(f"Year: {year} Dataset: {dataset}")    
            # Go through process of gridding and saving rocket emissions:
            emis_data = OutputEmis(event_alts, 
                                   reentry_levs, 
                                   months,
                                   dataset,
                                   year,
                                   ground_landing_list,
                                   stage_alt_dict,
                                   )

            check_total_emissions(year,dataset,GRID_RES,LEVELS)