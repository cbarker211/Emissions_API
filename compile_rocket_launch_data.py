from datetime import datetime, timedelta
from time import sleep
import numpy as np
import xarray as xr
import argparse


import sys
sys.path.append('./python_modules/')
from discosweb_api_func import server_request, wait_function, response_error_handler
from update_rocket_launch_data import update_mass_info

""" Script to request data from the DISCOSweb database. 
NB: Each request is limited to 30 results, and there is a maximum limit of requests in a certain timeframe (20 requests per 60s).
To solve this, the script makes multiple requests over different filters, and automatically calculates the time required to wait."""

class import_launches_from_discos:
    
    def __init__(self,start_year,final_year):
        """This function sets up the class.

        Args:
            year (_type_): The year.
        """        
        
        self.start_year = start_year
        self.final_year = final_year
        
    def get_launch_list(self): 
        """This function is called in all of the following functions.
        It gets the launch information for an entire year, which will be saved as self.full_launch_list

        Args:
            year (int): The year to get launches for.
        """        
        
        #Initialize variables.
        self.full_launch_list = []
        start_epoch = f'{str(self.start_year)}-01-01'
        end_epoch = f'{str(self.final_year+1)}-01-01'
        
        #This loop will run until all launches are found.
        page_number = 1
        while True:
            #Add the filtering here. Currently set up to show successful and unsuccessful launches per year.
            #The filtering params currently get all launches between start_epoch and end_epoch, sorted by epoch.
            params={
                    'filter': f"ge(epoch,epoch:'{start_epoch}')&lt(epoch,epoch:'{end_epoch}')",
                    'sort' : 'epoch',  
                    'page[number]' : page_number
                }
            #Now we need a while loop to deal with the rate limit. 
            #If the response is ok (code <400), then it breaks the while loop, and if not then the while loop continues after a delay.
            response = server_request(params, '/launches')
            
            if response.ok:            
                #This converts the reponse to a json object and exracts the data.
                doc = response.json()
                temp_launch_list = doc['data']
                self.full_launch_list.extend(temp_launch_list)
                
                #Each request gives a 'page' of results, and the maximum is 30 results per page. 
                #The script loops to request the next page.
                if len(temp_launch_list) == 30:
                    page_number += 1
                    continue
                else:
                    break
            elif response.status_code == 429:
                message = ""
                response_error_handler(response,message)
                continue
            else:
                response_error_handler(response,"")
            break
                    
        # Skip the Chang'e 5 launch which was the 'ascender' launching from the moon to dock with the 'orbiter'.
        # It eventually 'crashed' back onto the moon to avoid space debris, so is irrelevant for this inventory.
        change_found = False
        for count, launch in enumerate(self.full_launch_list):
            if launch["attributes"]["cosparLaunchNo"] == "2020-F11":
                change_found = True
                break
        if change_found == True:
            self.full_launch_list= np.delete(self.full_launch_list, count)
        
        #This section simply goes through each launch and counts the numbers of successes and failures.
        #It then prints the totals to the screen. Can be switched off if not needed.    
        launch_success = 0
        launch_failure = 0
        for launch in self.full_launch_list:
            if launch["attributes"]["failure"] == True:
                launch_failure+=1
            elif launch["attributes"]["failure"] == False:
                launch_success+=1
            else:
                print("Error: Launch not recorded as success or failure.")         
        print(f'\nFrom {self.start_year} to {self.final_year}, there were {launch_success} successful launches and {launch_failure} unsuccessful launches.')

        return self.full_launch_list

    def get_launch_info(self):
        """This function gathers all of the required launch information that HEMCO needs.

        Returns:
            launch_info: A list of all of the extracted launch information for the year.
        """        
        self.get_launch_list()
                
        self.Cospar_Id = []
        self.launch_time = []
        self.launch_datestr = []
        self.latitude = []
        self.longitude = []
        self.rocket_name = []
        self.vehicle_id = []
        self.mcs_check = []
        
        for count, launch in enumerate(self.full_launch_list):
            
            self.Cospar_Id.append(launch["attributes"]["cosparLaunchNo"])
            temp_launch_epoch = launch["attributes"]["epoch"]
            self.launch_datestr.append(temp_launch_epoch[:10].replace("-",""))
            
            temp_launch_time = temp_launch_epoch[10:]
            hour = int(temp_launch_time[1:3])
            minute = int(temp_launch_time[4:6])
            sec = int(temp_launch_time[7:9])
            time_utc = np.float64(hour + minute / 60 + sec / 3600)
            self.launch_time.append(time_utc)
            
            # Mark as MCS any launches containing payloads with the name:
            #   - Starlink
            #   - Oneweb
            #   - Yinhe
            #   - Lynk
            #   - E-Space
            params={
                    'filter' : "(eq(objectClass,Payload))&(contains(name,Starlink)|contains(name,OneWeb)|contains(name,Oneweb)|contains(name,Yinhe)|contains(name,Lynk)|contains(name,E-Space))",
                    'sort' : 'id'         
                }
            while True:
                response = server_request(params,f'/launches/{launch["id"]}/objects')
                
                if response.ok:
                    doc = response.json()
                    temp_object_list = doc['data']
                    # The Lynk 04 satellite is misassigned to 2020-011 instead of 2020-016 (confirmed used JSR). Fixed here.
                    if len(temp_object_list) > 0 and self.Cospar_Id[-1] != "2020-011":
                        self.mcs_check.append(True)
                    elif self.Cospar_Id[-1] == "2020-016": 
                        self.mcs_check.append(True)  
                    else:
                        self.mcs_check.append(False)   
                elif response.status_code == 429:
                    message = f"On launch {count} of {len(self.full_launch_list)} in {self.Cospar_Id[-1][:4]}."
                    response_error_handler(response,message)
                    continue
                else:
                    response_error_handler(response,"")
                break
            
            # Get the launch site info.
            while True:
                response = server_request({},f'/launches/{launch["id"]}/site')
                if response.ok:
                    doc = response.json()
                    self.latitude.append(np.float64(doc['data']["attributes"]["latitude"]))
                    self.longitude.append(np.float64(doc['data']["attributes"]["longitude"]))
                    break
                elif response.status_code == 429:
                    message = f"On launch {count} of {len(self.full_launch_list)} in {self.Cospar_Id[-1][:4]}."
                    response_error_handler(response,message)
                    continue
                else:
                    response_error_handler(response,"")
                break
                
            # Get the vehicle info.
            while True:
                response = server_request({},f'/launches/{launch["id"]}/vehicle')
                if response.ok:
                    doc = response.json()
                    if doc['data']["attributes"]["name"][:5] == "Atlas":
                        if datetime(int(self.launch_datestr[-1][:4]),int(self.launch_datestr[-1][4:6]),int(self.launch_datestr[-1][6:8])) < datetime(2020,11,13):
                            self.rocket_name.append(doc['data']["attributes"]["name"])
                        elif self.Cospar_Id[-1] in ["2021-042","2022-092"]:
                            self.rocket_name.append(doc['data']['attributes']["name"] + " v2021")
                        else:
                            self.rocket_name.append(doc['data']['attributes']["name"] + " v2020")
                    else:
                        self.rocket_name.append(doc['data']["attributes"]["name"])
                    self.vehicle_id.append(int(doc["data"]["id"]))
                elif response.status_code == 429:
                    message = f"On launch {count} of {len(self.full_launch_list)} in {self.Cospar_Id[-1][:4]}."
                    response_error_handler(response,message)
                    continue
                else:
                    response_error_handler(response,"")
                break
                    
        return [self.Cospar_Id, self.launch_time, self.launch_datestr, self.latitude, self.longitude, self.rocket_name, self.mcs_check]
            
    def launch_info_to_netcdf(self):
        """This saves the launch information as a NetCDF file for later processing for GEOS-Chem.
        This is only needed if using GEOS-Chem. Recommended to ignore if you don't need this.
        """        
        
        #Set up the dimensions of the netcdf file.
        #launches = np.arange(len(self.Cospar_Id), dtype=np.int64)
        dims = ('launches')
    
        #Create the DataArrays.
        data_da_cospar_id = xr.DataArray(self.Cospar_Id, dims=dims,
            attrs=dict(long_name="COSPAR ID",
                       short_name="COSPAR ID"))
        
        data_da_time = xr.DataArray(self.launch_time, dims=dims,
            attrs=dict(long_name="Time (UTC)",
                       short_name="Time (UTC)",
                       ))
        
        data_da_date = xr.DataArray(self.launch_datestr, dims=dims,
            attrs=dict(long_name="Date in form YYYYMMDD",
                       short_name="Date",
                       ))
        
        data_da_lon = xr.DataArray(self.longitude, dims=dims,
            attrs=dict(long_name="Longitude",
                       short_name="Longitude",
                       units="Degrees"))
        
        data_da_lat = xr.DataArray(self.latitude, dims=dims,
            attrs=dict(long_name="Latitude",
                       short_name="Latitude",
                       units="Degrees"))
        
        data_da_rocket_name = xr.DataArray(self.rocket_name, dims=dims,
            attrs=dict(long_name="Rocket Type",
                       short_name="Rocket Type",
                       ))
        
        data_da_vehicle_id = xr.DataArray(self.vehicle_id, dims=dims,
            attrs=dict(long_name="Vehicle ID",
                       short_name="Vehicle ID",
                       ))
        
        data_da_launch_MCS = xr.DataArray(self.mcs_check, dims=dims,
            attrs=dict(long_name="MCS Launch Check",
                       short_name="MCS Launch Check",
                       ))
        
        # Create an xarray Dataset from the DataArrays.
        ds = xr.Dataset()
        ds['COSPAR_ID'] = data_da_cospar_id
        ds['Time(UTC)'] = data_da_time
        ds['Date'] = data_da_date
        ds['Longitude'] = data_da_lon
        ds['Latitude'] = data_da_lat
        ds['Rocket_Name'] = data_da_rocket_name
        ds['DISCOSweb_Rocket_ID'] = data_da_vehicle_id
        ds['Megaconstellation_Flag'] = data_da_launch_MCS
             
        #Save to file and close the DataSet     
        ds.to_netcdf(f'./databases/launch_activity_data_{self.start_year}-{self.final_year}.nc')
        ds.close()
        
    def handle_stage_info(self, count, stage, i, temp_dict, unique_vehicle_name_list):
        
        ##########################
        # Manually rearrage stages
        ##########################
        
        # Sometimes the stages are not in the right order on DISCOSweb, so we need to manually rearrange.
        if unique_vehicle_name_list[i] in ["Angara A5", "Angara A5 Persei", "Ariane 5ECA", 
                                           "Atlas V N22 v2020", "Atlas V 401 v2020", "Atlas V 411", "Atlas V 421 v2021",
                                           "Atlas V 511 v2020", "Atlas V 531 v2020", 
                                           "Atlas V 541", "Atlas V 541 v2020", "Atlas V 551", "Atlas V 551 v2020",
                                           "Delta 4H", "Epsilon-2 CLPS", "H-IIA 202", "H-IIA 204", "H-IIB", "GSLV Mk II ", "GSLV Mk III",
                                           "Jielong-3", "PSLV-DL", "PSLV-XL", "PSLV-CA", "Space Launch System - Block 1 Crew", "Falcon Heavy",
                                           "Long March (CZ) 11", "Long March (CZ) 2F", "Long March (CZ) 3B", "Long March (CZ) 3C", "Long March (CZ) 5", 
                                           "Long March (CZ) 5B", "Long March (CZ) 6", "Long March (CZ) 6A", "Long March (CZ) 7", "Long March (CZ) 7A", 
                                           "Long March (CZ) 8", "Minotaur 1", "Pegasus XL", "Vega C",
                                           "Soyuz-2-1A", "Soyuz-2-1A Fregat-M", "Soyuz-2-1B Fregat", "Soyuz-2-1B Fregat-M", "Soyuz-ST-A Fregat-M",
                                           "Soyuz-2-1A Fregat","Soyuz-2-1B","Soyuz-2-1V Volga","Soyuz-ST-B Fregat-MT"]:
            
            stage_name = stage["attributes"]["name"]

            if "booster" in stage_name.lower() or stage_name in ["EPS P241", "Soyuz-FG Blok-B,V,G,D","K2-1","GSLV-SOM","Falcon 9 Merlin-1D+","S-200"]:
                stage_number = "Booster"
                
            # For stage ids 773 and 775, this is stages with the same name for the Delta 4H rocket.  
            # Same for 112791 and 112794 for Jielong-3.
            elif (stage_name in ["URM-1","EPC H173", "Atlas V CCB", "H-II LE-7A", "H-II widebody LE-7A", "PS1", "Soyuz-FG Blok-A", 	"Soyuz-2.1V Blok-A",
                                 "Long March (CZ) 11 Stage 1", "L-186 (YF21B)", "L-172 (YF21C)", "CZ-5-500", "L-165 (H5-1)", "L-76 (YF100)",
                                 "GSLV-1", "Minuteman 1", "ORION 50SXL","SLS-B1 First Stage","P120C","Falcon Heavy Centre Core","L-110","JIE-3 first stage"]) or stage["id"] == "775":
                stage_number = "Stage1"
                
            elif (stage_name in ["URM-2","ESC-A","Centaur-5 SEC","Centaur-5 DEC", "H-II LE-5B", "PS2", "Blok-I", "Long March (CZ) 11 Stage 2", 
                                 "L-84 (YF26)", "L-45 (YF-24E)", "CZ-5-HO", "L-15 (YF115)", "L-15 (YF-115)", "K3-2 short",
                                 "M-35","GSLV-2","SR19","SLS  iCPS","Zefiro 40C","Falcon Heavy second stage","JIE-3 second stage","C-25"]) or stage["id"] == "773":
                stage_number = "Stage2"
                
            elif (stage_name in ["Briz-M", "PS3", "Fregat-M", "Fregat", "Fregat-MT", "Long March (CZ) 11 Stage 3", "PBV",
                                "Blok-DM-3","KM-V2c","GS3 (GSLV Mk II)","Zefiro 9","Volga (Soyuz-2-1V Volga)"]) or stage["id"] == "112791":
                stage_number = "Stage3"
                
            elif (stage_name in ["PSLV fourth stage (PS4)", "Long March (CZ) 11 Stage 4","CLPS","AVUM+"]) or stage["id"] == "112794":
                stage_number = "Stage4"
                
            # The Pegasus, Antares, Minotaur and Taurus rockets all use similar stages 
            elif stage_name == "ORION 50XL":
                if unique_vehicle_name_list[i] == "Pegasus XL":
                    stage_number = "Stage2"
                elif unique_vehicle_name_list[i] == "Minotaur 1":
                    stage_number = "Stage3"
            elif stage_name == "ORION 38 (Pegasus XL)":
                if unique_vehicle_name_list[i] == "Pegasus XL":
                    stage_number = "Stage3"
                elif unique_vehicle_name_list[i] == "Minotaur 1":
                    stage_number = "Stage4"
            elif stage_name == "H-II SRB-A":
                if unique_vehicle_name_list[i] == "Epsilon-2 CLPS":
                    stage_number = "Stage1"
                else:
                    stage_number = "Booster"
            
            # Long March often uses the same stage for different rockets, but not always at the same position.    
            elif stage_name == "H-18 (Long March (CZ) YF)":
                if unique_vehicle_name_list[i] in ["Long March (CZ) 3B", "Long March (CZ) 3C", "Long March (CZ) 7A"]:
                    stage_number = "Stage3"
                elif unique_vehicle_name_list[i] == "Long March (CZ) 8":
                    stage_number = "Stage2"

            elif stage_name == "K3-1":
                if unique_vehicle_name_list[i] in ["Long March (CZ) 5", "Long March (CZ) 5B"]:
                    stage_number = "Booster"
                elif unique_vehicle_name_list[i] in ["Long March (CZ) 8", "Long March (CZ) 7", "Long March (CZ) 7A"]:
                    stage_number = "Stage1"
                    
            else:
                print(f"Warning: Error manually assigning stages and boosters. Stage {stage_name} not in check lists.")
        else:
            stage_number = f"Stage{count+1}"
            
        if unique_vehicle_name_list[i] in ["Soyuz-ST-A Fregat-M","Soyuz-2-1A Fregat"]:
            print(unique_vehicle_name_list[i],stage_name, stage_number)
            
        #############################
        # Obtain the propellant mass
        #############################
        
        # Firstly get whatever value DISCOSweb has for the mass.
        # Also convert any none values to 0.
        temp_fuel_mass = stage["attributes"]["fuelMass"]
        if temp_fuel_mass == None:
            temp_fuel_mass = 0
        temp_oxidiser_mass = stage["attributes"]["oxidiserMass"]
        if temp_oxidiser_mass == None:
            temp_oxidiser_mass = 0
        temp_solid_mass = stage["attributes"]["solidPropellantMass"]
        if temp_solid_mass == None:
            temp_solid_mass = 0
        
        # Calculate total prop mass and check if its missing completely.
        temp_prop_mass = temp_fuel_mass+temp_oxidiser_mass+temp_solid_mass
        if temp_prop_mass == 0:
            pass
            #print(f"Warning: No propellant mass found for {stage_number} of {unique_vehicle_name_list[i]}")
        
        # Get the propellant name and assign type.
        while True:
            sleep(0.2) # Sometimes DISCOSweb will refuse the request if you try too many times rapidly, so set a small delay.
            response = server_request({},f'/launch-vehicles/stages/{stage["id"]}/propellant')
            if response.ok:
                doc = response.json()
                if doc['data'] != None:
                    propellant_info = doc['data']["attributes"]
                    
                    if (propellant_info["fuel"] != None) and (propellant_info["oxidiser"] != None) and (propellant_info["solidPropellant"] != None):
                        print("HYBRID ROCKET DETECTED")
                    #First get the name.
                    if (propellant_info["fuel"] != None) and (propellant_info["oxidiser"] != None):
                        temp_dict[f"{stage_number} Propellant Name"] = propellant_info["fuel"] + "/" + propellant_info["oxidiser"]
                    elif propellant_info["solidPropellant"] != None:
                        temp_dict[f"{stage_number} Propellant Name"] = propellant_info["solidPropellant"]
                        #The else case is handled below in the next else block.
                        
                    # Fix minor spelling errors from database.
                    if stage["attributes"]["name"] == "Antares-200 first stage (RD-181)":
                        temp_dict[f"Stage1 Propellant Name"] = "RP1 (Rocket Propellant 1)/LOX"
                    
                    if stage["attributes"]["name"] == "AVUM":
                        temp_dict[f"Stage4 Propellant Name"] = "UDMH (Unsymmetrical Dimethyl Hydrazine)/N2O4"
                    
                    #Now assign to a fuel type.
                    if temp_dict[f"{stage_number} Propellant Name"] in ['HTPB','Solid','HTPB-1912','TP-H8299',"PBAN-Al/NH4ClO4"]:
                        temp_dict[f"{stage_number} Fuel Type"] = "Solid"  
                    elif temp_dict[f"{stage_number} Propellant Name"] in ['RP1 (Rocket Propellant 1)/LOX','Kerosene/LOX']:
                        temp_dict[f"{stage_number} Fuel Type"] = "Kerosene" 
                    elif "Hydrazine" in temp_dict[f"{stage_number} Propellant Name"] or temp_dict[f"{stage_number} Propellant Name"] == "UH25 (75%UDMH+25%N2H4)/N2O4":
                        temp_dict[f"{stage_number} Fuel Type"] = "Hypergolic"  
                    elif temp_dict[f"{stage_number} Propellant Name"] in ['LH2 (Liquid Hydrogen)/LOX']:
                        temp_dict[f"{stage_number} Fuel Type"] = "Hydrogen"  
                    elif temp_dict[f"{stage_number} Propellant Name"] in ['LNG (Liquid Natural Gas)/LOX']:
                        temp_dict[f"{stage_number} Fuel Type"] = "Methane"  
                    else:
                        pass
                        print(f"Warning: Couldn't assign propellant <{temp_dict[f'{stage_number} Propellant Name']}> for {stage_number} of {unique_vehicle_name_list[i]}.")     
                else:
                    pass
                    print(f"Warning: Couldn't find any propellant for {stage_number} of {unique_vehicle_name_list[i]}.")
                    
            elif response.status_code == 429:
                message = f"On rocket {i+1} of {len(unique_vehicle_name_list)}."
                response_error_handler(response,message)
                continue
            else:
                response_error_handler(response,"")
            break
        
        temp_dict[f"{stage_number} Propellant Mass"] = temp_prop_mass
        
        # Get the dry stage mass.
        temp_dict[f"{stage_number} Stage Mass"] = stage["attributes"]["dryMass"]
        #print(temp_dict[f"{stage_number} Stage Mass"])
        
        return temp_dict
    
    def get_rocket_info(self):
        
        self.unique_rocket_list = []
        
        #Open file and get a list of unique vehicles.
        vehicle_name_list, vehicle_id_list = [], []
        
        ds = xr.open_dataset(f'./databases/launch_activity_data_{self.start_year}_{self.final_year}.nc', decode_times=False)
        vehicle_name_list = np.concatenate([vehicle_name_list,ds['Rocket_Name'].values])
        vehicle_id_list = np.concatenate([vehicle_id_list,ds['DISCOSweb_Rocket_ID'].values]).astype(int)
        ds.close()
        unique_vehicle_name_list = sorted(list(set(vehicle_name_list))) 
        
        #Loop over all rockets, and pull the information for each.
        for i in range(0,len(unique_vehicle_name_list)):
            
            temp_vehicle_id = vehicle_id_list[np.where(vehicle_name_list == unique_vehicle_name_list[i])[0][0]]
            temp_dict = {
                "id" : temp_vehicle_id,
                "name": unique_vehicle_name_list[i],
            }
            if temp_dict["name"] == "Astra Rocket 3":
                temp_dict["proxy"] = "Electron"
            elif temp_dict["name"] == "Ceres-1":
                temp_dict["proxy"] = "Shavit"
            elif temp_dict["name"] == "Long March (CZ) 11":
                temp_dict["proxy"] = "Minotaur 1"
            elif temp_dict["name"] == "Jielong-3":
                temp_dict["proxy"] = "Epsilon-2 CLPS"
            elif temp_dict["name"] == "Zhongke 1A":
                temp_dict["proxy"] = "Vega C"
            elif temp_dict["name"] == "Long March (CZ) 6A":
                temp_dict["proxy"] = "Long March (CZ) 7A"
            elif temp_dict["name"] == "Kuaizhou-11":
                temp_dict["proxy"] = "Long March (CZ) 6"
            elif temp_dict["name"] == "Zhuque-2":
                temp_dict["proxy"] = "Antares 230"
            else:
                temp_dict["proxy"] = ""
                    
            # Get the propellant mass info.
            while True:
                response = server_request({},f'/launch-vehicles/{temp_vehicle_id}/stages')
                if response.ok:
                    doc = response.json()
                    vehicle_info = doc['data']
                    # Usually 2-5 stages, so won't need to account for the maximum of 30 items in a list.
                
                    #Set up the arrays
                    temp_dict[f"Booster Propellant Mass"] = 0
                    temp_dict[f"Booster Fuel Type"] = ""
                    temp_dict[f"Booster Propellant Name"] = ""
                    temp_dict[f"Booster Stage Mass"] = 0
                    temp_dict[f"Booster Number"] = 0
                    temp_dict[f"Fairing Mass"] = 0
                    for j in range(0,4):
                        temp_dict[f"Stage{j+1} Propellant Mass"] = 0
                        temp_dict[f"Stage{j+1} Fuel Type"] = ""
                        temp_dict[f"Stage{j+1} Propellant Name"] = ""
                        temp_dict[f"Stage{j+1} Stage Mass"] = 0
                        
                    for count,stage in enumerate(vehicle_info):
                        temp_dict = self.handle_stage_info(count,stage,i,temp_dict,unique_vehicle_name_list) 
                                               
                elif response.status_code == 429:
                    message = f"On rocket {i+1} of {len(unique_vehicle_name_list)}. "
                    response_error_handler(response,message)
                    continue
                else:
                    response_error_handler(response,"")
                break      
            
            # Handle any incomplete/incorrect propellant/stage mass information.
            temp_dict = update_mass_info(i,temp_dict,unique_vehicle_name_list)
            # Update the rocket list.   
            self.unique_rocket_list.append(temp_dict)    
        
    def rocket_info_to_netcdf(self):
        """This saves the propellant information as a NetCDF file for later processing for GEOS-Chem.
        This is only needed if using GEOS-Chem. Recommended to ignore if you don't need this.
        """        
        #Set up the dimensions of the netcdf file.
        #launches = np.arange(len(self.Cospar_Id), dtype=np.int64)
        dims = ('rockets')
        
        #Set up the data
        rocket_name_list = []
        booster_propmass_list, booster_fuel_type_list, booster_stagemass_list, booster_number_list = [], [], [], []
        stage1_propmass_list,  stage1_fuel_type_list, stage1_stagemass_list  = [], [], []
        stage2_propmass_list,  stage2_fuel_type_list, stage2_stagemass_list  = [], [], []
        stage3_propmass_list,  stage3_fuel_type_list, stage3_stagemass_list  = [], [], []
        stage4_propmass_list,  stage4_fuel_type_list, stage4_stagemass_list  = [], [], []
        fairing_stagemass_list, proxy_list = [], []

        for rocket in self.unique_rocket_list:
            rocket_name_list.append(rocket["name"])
            booster_propmass_list.append(rocket["Booster Propellant Mass"])
            booster_fuel_type_list.append(rocket["Booster Fuel Type"])
            booster_stagemass_list.append(rocket["Booster Stage Mass"])
            booster_number_list.append(rocket["Booster Number"])
            stage1_propmass_list.append(rocket["Stage1 Propellant Mass"])
            stage1_fuel_type_list.append(rocket["Stage1 Fuel Type"])
            stage1_stagemass_list.append(rocket["Stage1 Stage Mass"])
            stage2_propmass_list.append(rocket["Stage2 Propellant Mass"])
            stage2_fuel_type_list.append(rocket["Stage2 Fuel Type"])
            stage2_stagemass_list.append(rocket["Stage2 Stage Mass"])
            stage3_propmass_list.append(rocket["Stage3 Propellant Mass"])
            stage3_fuel_type_list.append(rocket["Stage3 Fuel Type"])
            stage3_stagemass_list.append(rocket["Stage3 Stage Mass"])
            stage4_propmass_list.append(rocket["Stage4 Propellant Mass"])
            stage4_fuel_type_list.append(rocket["Stage4 Fuel Type"])
            stage4_stagemass_list.append(rocket["Stage4 Stage Mass"])
            fairing_stagemass_list.append(rocket["Fairing Mass"])
            proxy_list.append(rocket["proxy"])
            
        #Create the DataArrays.
        data_da_name                = xr.DataArray(rocket_name_list,       dims=dims, attrs=dict(long_name="Rocket Name"))
        data_da_booster_propmass    = xr.DataArray(booster_propmass_list,  dims=dims, attrs=dict(long_name="Booster Propellant Mass", units="kg"))
        data_da_booster_fuel_type   = xr.DataArray(booster_fuel_type_list, dims=dims, attrs=dict(long_name="Booster Fuel Type"))
        data_da_booster_stagemass   = xr.DataArray(booster_stagemass_list, dims=dims, attrs=dict(long_name="Booster Stage Mass", units="kg"))
        data_da_booster_number      = xr.DataArray(booster_number_list,    dims=dims, attrs=dict(long_name="Booster Number"))
        data_da_stage1_propmass     = xr.DataArray(stage1_propmass_list,   dims=dims, attrs=dict(long_name="Stage 1 Propellant Mass", units="kg"))
        data_da_stage1_fuel_type    = xr.DataArray(stage1_fuel_type_list,  dims=dims, attrs=dict(long_name="Stage 1 Fuel Type"))
        data_da_stage1_stagemass    = xr.DataArray(stage1_stagemass_list,  dims=dims, attrs=dict(long_name="Stage 1 Stage Mass", units="kg"))
        data_da_stage2_propmass     = xr.DataArray(stage2_propmass_list,   dims=dims, attrs=dict(long_name="Stage 2 Propellant Mass", units="kg"))
        data_da_stage2_fuel_type    = xr.DataArray(stage2_fuel_type_list,  dims=dims, attrs=dict(long_name="Stage 2 Fuel Type"))
        data_da_stage2_stagemass    = xr.DataArray(stage2_stagemass_list,  dims=dims, attrs=dict(long_name="Stage 2 Stage Mass", units="kg"))
        data_da_stage3_propmass     = xr.DataArray(stage3_propmass_list,   dims=dims, attrs=dict(long_name="Stage 3 Propellant Mass", units="kg"))
        data_da_stage3_fuel_type    = xr.DataArray(stage3_fuel_type_list,  dims=dims, attrs=dict(long_name="Stage 3 Fuel Type"))
        data_da_stage3_stagemass    = xr.DataArray(stage3_stagemass_list,  dims=dims, attrs=dict(long_name="Stage 3 Stage Mass", units="kg"))
        data_da_stage4_propmass     = xr.DataArray(stage4_propmass_list,   dims=dims, attrs=dict(long_name="Stage 4 Propellant Mass", units="kg"))
        data_da_stage4_fuel_type    = xr.DataArray(stage4_fuel_type_list,  dims=dims, attrs=dict(long_name="Stage 4 Fuel Type"))
        data_da_stage4_stagemass    = xr.DataArray(stage4_stagemass_list,  dims=dims, attrs=dict(long_name="Stage 4 Stage Mass", units="kg"))
        data_da_fairing_stagemass   = xr.DataArray(fairing_stagemass_list, dims=dims, attrs=dict(long_name="Fairing Mass", units="kg"))
        data_da_proxy               = xr.DataArray(proxy_list,             dims=dims, attrs=dict(long_name="Proxy Rocket"))
        
        # Create an xarray Dataset from the DataArrays.
        ds = xr.Dataset()
        ds['Rocket_Name']           = data_da_name
        ds['Booster_No']            = data_da_booster_number
        ds['Booster_PropMass']      = data_da_booster_propmass
        ds['Booster_Fuel_Type']     = data_da_booster_fuel_type
        ds['Booster_StageMass']     = data_da_booster_stagemass
        ds['Stage1_PropMass']       = data_da_stage1_propmass
        ds['Stage1_Fuel_Type']      = data_da_stage1_fuel_type
        ds['Stage1_StageMass']      = data_da_stage1_stagemass
        ds['Stage2_PropMass']       = data_da_stage2_propmass
        ds['Stage2_Fuel_Type']      = data_da_stage2_fuel_type
        ds['Stage2_StageMass']      = data_da_stage2_stagemass
        ds['Stage3_PropMass']       = data_da_stage3_propmass
        ds['Stage3_Fuel_Type']      = data_da_stage3_fuel_type
        ds['Stage3_StageMass']      = data_da_stage3_stagemass
        ds['Stage4_PropMass']       = data_da_stage4_propmass
        ds['Stage4_Fuel_Type']      = data_da_stage4_fuel_type
        ds['Stage4_StageMass']      = data_da_stage4_stagemass
        ds['Fairing_Mass']          = data_da_fairing_stagemass
        ds['Proxy_Rocket']          = data_da_proxy
             
        #Save to file and close the DataSet     
        #ds.to_netcdf(f'./databases/rocket_attributes_{self.start_year}-{self.final_year}_noupdate.nc')
        ds.to_netcdf(f'./databases/rocket_attributes_{self.start_year}-{self.final_year}.nc')
        ds.close()
        
    def save_discosweb_reentries(self):
        
        # Query DISCOSweb for all reentries.
        start_epoch = f'{str(self.start_year)}-01-01'
        end_epoch = f'{str(self.final_year+1)}-01-01'  
        full_reentry_list = []
        page_number = 1
        while True:
            #The filtering params currently get all reentries between start_epoch and end_epoch, sorted by epoch.
            params={
                    'filter': f"ge(epoch,epoch:'{start_epoch}')&lt(epoch,epoch:'{end_epoch}')",
                    'sort' : 'epoch',
                    'page[number]' : page_number
                }
            response = server_request(params, '/reentries')
            if response.ok:
                doc = response.json()
                temp_reentry_list = doc['data']
                full_reentry_list.extend(temp_reentry_list)
                #Each request gives a 'page' of results, and the maximum is 30 results per page. 
                #The script loops to request the next page.               
                if len(temp_reentry_list) == 30:
                    page_number += 1
                    continue
                else:
                    break
            elif response.status_code == 429:
                message = f"Found {len(full_reentry_list)} reentries."
                response_error_handler(response,message)
                continue
            else:
                response_error_handler(response,"")
            break  
        
        # Remove any duplicates.
        unique_DISCOSweb_reentry_list = []
        for item in full_reentry_list:
            current_id_list = []
            for unique_item in unique_DISCOSweb_reentry_list:
                current_id_list.append(unique_item["id"]) 
            if item["id"] not in current_id_list:
                unique_DISCOSweb_reentry_list.append(item)
        
        #We now have a list of reentries, however we need to query DISCOSweb for the object details.  
        name_list, mass_list, cosparId_list, objectClass_list, epoch_list = [], [], [], [], []      
        for count, reentry in enumerate(unique_DISCOSweb_reentry_list):
            while True:
                response = server_request({},f'/reentries/{reentry["id"]}/objects')
                if response.ok:
                    doc = response.json()
                    temp_object_list = doc['data']
                elif response.status_code == 429:
                    message = f"On reentry {count} of {len(unique_DISCOSweb_reentry_list)}."
                    response_error_handler(response,message)
                    continue
                else:
                    response_error_handler(response,"")
                break
            
            if len(temp_object_list) > 0:  
                name_list.append(temp_object_list[0]["attributes"]["name"])
                mass_list.append(temp_object_list[0]["attributes"]["mass"])
                cosparId_list.append(temp_object_list[0]["attributes"]["cosparId"])
                objectClass_list.append(temp_object_list[0]["attributes"]["objectClass"])
                epoch_list.append(reentry["attributes"]["epoch"])
            
        dims = ('discosweb_reentries')         
        #Create the DataArrays.
        data_da_name          = xr.DataArray(name_list,          dims=dims, attrs=dict(long_name="Object Name"))
        data_da_mass          = xr.DataArray(mass_list,          dims=dims, attrs=dict(long_name="Object Mass", units="kg"))
        data_da_cosparId      = xr.DataArray(cosparId_list,      dims=dims, attrs=dict(long_name="Object ID"))
        data_da_objectClass   = xr.DataArray(objectClass_list,   dims=dims, attrs=dict(long_name="Object Class"))
        data_da_epoch         = xr.DataArray(epoch_list,         dims=dims, attrs=dict(long_name="Reentry Epoch"))
    
        # Create an xarray Dataset from the DataArrays.
        ds = xr.Dataset()
        ds['DISCOSweb_Reentry_Name']     = data_da_name
        ds['DISCOSweb_Reentry_Mass']     = data_da_mass
        ds['DISCOSweb_Reentry_ID']       = data_da_cosparId
        ds['DISCOSweb_Reentry_Class']    = data_da_objectClass
        ds['DISCOSweb_Reentry_Epoch']    = data_da_epoch
             
        #Save to file and close the DataSet     
        ds.to_netcdf(f'./databases/reentry/DISCOSweb/discosweb_reentries_{self.start_year}-{self.final_year}.nc')
        ds.close()
                         
if __name__ == "__main__":
    """The main running of the program goes here. 
    All of the functions are inside the import_launches_from_discos class.
    Call each of the functions using import_launches_from_discos._function_
    e.g. import_launches_from_discos.launches_per_year(start_year,end_year)
    """         
    
    # Set up the arguments for each function.
    parser = argparse.ArgumentParser()
    parser.add_argument('-yl', "--yearly_launches", action='store_true', help="Get yearly launches.")
    parser.add_argument('-li', "--launch_info", action='store_true', help='Get launch info.')
    parser.add_argument('-sli', "--save_launch_info", action='store_true', help='Save launch info.')
    parser.add_argument('-ri', "--rocket_info", action='store_true', help='Get rocket info.')
    parser.add_argument('-sri', "--save_rocket_info", action='store_true', help='Save launch info.')
    parser.add_argument('-sdwr', "--save_discosweb_reentries", action='store_true', help='Save launch info.')
    parser.add_argument('-sy', "--start_year", default = "2020", choices=str(np.arange(1957,2024)), help='Start Year.')
    parser.add_argument('-fy', "--final_year", default = "2022", choices=str(np.arange(1957,2024)), help='Final Year.')
    args = parser.parse_args()
    
    # Sort out the year range.
    start_year = int(args.start_year)
    final_year = int(args.final_year)
    if start_year > final_year:
        final_year = start_year + 1
    print(f"Processing from year {start_year} to {final_year}.")  
    
    #Loop over all years and run functions depending on input arguments.
    LaunchData = import_launches_from_discos(start_year,final_year)    
    if args.yearly_launches == True:
        LaunchData.get_launch_list()
    if args.launch_info == True:
        LaunchData.get_launch_info()      
        if args.save_launch_info == True:
            LaunchData.launch_info_to_netcdf()               
    if args.save_discosweb_reentries == True:
        LaunchData.save_discosweb_reentries() 
    if args.rocket_info == True:
        LaunchData.get_rocket_info()      
        if args.save_rocket_info == True:
            LaunchData.rocket_info_to_netcdf()
