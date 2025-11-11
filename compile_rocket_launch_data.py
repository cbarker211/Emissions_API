from datetime import datetime
from time import sleep
import numpy as np
import xarray as xr
import argparse
import pandas as pd
import requests
from io import StringIO
import time
from tqdm import tqdm
from cProfile import Profile
from pstats import Stats

from python_modules.discosweb_api_func import server_request, response_error_handler
from update_rocket_launch_data import update_mass_info

""" Script to request data from the DISCOSweb database. 
NB: Each request is limited to 30 results, and there is a maximum limit of requests in a certain timeframe (20 requests per 60s).
To solve this, the script makes multiple requests over different filters, and automatically calculates the time required to wait."""

class import_launches:
    
    def __init__(self,start_year,final_year):
        """This function sets up the class.

        Args:
            year (_type_): The year.
        """        
        
        self.start_year = start_year
        self.final_year = final_year
        self.session = requests.Session()

    def server_loop(self,urlpath,params,message):

        while True:
            response = server_request(params,urlpath)
            if response.ok:
                return(response.json())
            elif response.status_code == 429:
                response_error_handler(response,message)
                continue
            else:
                response_error_handler(response,"")
                break
            
        return None
        
    def get_launch_list(self): 
        """
        Retrieve launch information for a given year range and store it in self.full_launch_list.

        This function fetches launch data between self.start_year and self.final_year (inclusive),
        handling pagination and rate limits automatically. It also removes irrelevant launches
        (e.g., Chang'e 5 ascender) and reports the number of successes and failures.

        Args:
            year (int): The year to get launches for.
        """        
        
        #Initialize variables.
        self.full_launch_list = []
        start_epoch = f'{str(self.start_year)}-01-01'
        end_epoch = f'{str(self.final_year+1)}-01-01'
        
        # Fetch all pages of launch activity data.
        page_number = 1
        while True:
            params={
                    'filter': f"ge(epoch,epoch:'{start_epoch}')&lt(epoch,epoch:'{end_epoch}')", # launches bteween start_epoch and end_epoch
                    'sort' : 'epoch', # sort by epoch
                    'page[number]' : page_number 
                }

            response = server_request(params, '/launches')
            
            # If the response is ok, extract the data and check if we need to request another page.
            if response.ok: 
                data = response.json().get("data", [])           
                self.full_launch_list.extend(data)
                
                # If the page has less than 30 items, its the last one.
                if len(data) < 30:
                    break
                page_number += 1
                    
            elif response.status_code == 429:
                message = ""
                response_error_handler(response,message)
                continue
            else:
                response_error_handler(response,"")
                    
        # Skip Chang'e 5 launch (this was from the moon and eventually 'crashed' back onto the moon to avoid space debris, so is irrelevant for this inventory).
        self.full_launch_list = [launch for launch in self.full_launch_list if launch["attributes"]["cosparLaunchNo"] != "2020-F11"]
        
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
        self.launches = []
        
        for count, launch in enumerate(self.full_launch_list):

            launch_info = {}
            launch_info["COSPAR_ID"] = launch["attributes"]["cosparLaunchNo"]
            message = f"On launch {count} of {len(self.full_launch_list)} in {launch_info['COSPAR_ID'][:4]}."
            
            if launch_info["COSPAR_ID"] == "2023-F07": # This is wrong in DISCOSweb, it puts F07 and F08 at the same date and time.
                launch_info["Date"] = "20230530"
                launch_info["Time(UTC)"] = 21.36
            else:
                temp_launch_epoch = launch["attributes"]["epoch"]
                dt = datetime.fromisoformat(temp_launch_epoch)
                launch_info["Date"] = dt.strftime("%Y%m%d")
                launch_info["Time(UTC)"] = dt.hour + dt.minute / 60 + dt.second / 3600

            # Mark as MCS any launches containing payloads with the name:
            #   - Starlink
            #   - Oneweb
            #   - Yinhe (and Lingxi)
            #   - Lynk
            #   - E-Space (and Protosat-1)
            #   - Tranche (covers all 2023-050, 1 2023-133, all 2024-028). 2023-133 also contains Wildfire 1-10 and BB 4 and BB 4. But this launch is already tagged.
            #   - Chinese National Constellation (Guowang/Xingwang/Guangwang/Hulianwang) - Ronghe, Digui, Hulianwang.
            #   - Qianfan
            params={'filter' : "(eq(objectClass,Payload))&(contains(name,Starlink)|contains(name,OneWeb)|contains(name,Oneweb)|contains(name,Yinhe)|contains(name,Lynk)|contains(name,E-Space)|contains(name,Lingxi)|contains(name,Protosat-1)|contains(name,Kuiper)|contains(name,Tranche)|contains(name,Ronghe)|contains(name,Digui)|contains(name,Hulianwang)|contains(name,Qianfan))",
                    'sort' : 'id'}
            
            doc = self.server_loop(f'/launches/{launch["id"]}/objects',params,message)
            
            if doc:
                # The Lynk 04 satellite is misassigned to 2020-011 instead of 2020-016 (confirmed used JSR). Fixed here.
                launch_info["Megaconstellation_Flag"] = bool(doc["data"]) and launch_info["COSPAR_ID"] != "2020-011" or launch_info["COSPAR_ID"] == "2020-016"
        
            # Get the launch site info.
            site_coords = {
                "Newquay, Spaceport Cornwall": (50.439240, -4.999055),
                "China Sea Launch": (35.4943, 123.7965),
                "Naro Space Center": (34.5, 127.5),
            }

            doc = self.server_loop(f'/launches/{launch["id"]}/site',{},message)
            if doc:
                site_name = doc['data']["attributes"]["name"]
                launch_info["Site"] = site_name

                lat_lon = site_coords.get(site_name)
                if lat_lon:
                    launch_info["Latitude"], launch_info["Longitude"] = lat_lon
                    print(f"{launch_info['COSPAR_ID']} set to {site_name}.")
                else:
                    launch_info["Latitude"] = np.float64(doc["data"]["attributes"]["latitude"])
                    launch_info["Longitude"] = np.float64(doc["data"]["attributes"]["longitude"])

            # Get the rocket info.
            doc = self.server_loop(f'/launches/{launch["id"]}/vehicle',{},message)
            if doc:
                rocket_name = doc['data']["attributes"]["name"]
                launch_date = datetime.strptime(launch_info["Date"], "%Y%m%d")
                launch_info["Rocket_Name"] = rocket_name

                # Sometimes the rocket has been upgraded, so we need to specify which version was used.
                # TODO: Need to check this doesn't mess things up pre-2020.
                if rocket_name.startswith("Atlas"):
                    if launch_date < datetime(2020,11,13):
                        launch_info["Rocket_Name"] = rocket_name
                    elif launch_info['COSPAR_ID'] in ["2021-042","2022-092"]:
                        launch_info["Rocket_Name"] = rocket_name + " v2021"
                    else:
                        launch_info["Rocket_Name"] = rocket_name + " v2020"
                elif rocket_name.startswith("Zhuque-2"):
                    if launch_date > datetime(2024,1,1):
                        launch_info["Rocket_Name"] = rocket_name + "E"

                launch_info["DISCOSweb_Rocket_ID"] = int(doc["data"]["id"])

            # Save to the dictionary list.
            self.launches.append(launch_info)

    def scrape_jsr(self,url):
        # Function to web scrape data and convert to a pandas DataFrame.
        response = self.session.get(url)
        if response.status_code == 200:
            # Convert the content to a file-like object for pandas
            tsv_data = StringIO(response.text)
            # Load the data into a pandas DataFrame
            df = pd.read_csv(tsv_data, delimiter="\t", dtype=object, low_memory=False)
        else:
            raise ImportError(f"Failed to fetch from JSR URL: {url}", response.status_code)
        return df
    
    def get_launch_info_jsr(self):
        
        ####################################################################
        # Web scrape the launch data from Jonathan McDowell's JSR website.
        ####################################################################

        df       = self.scrape_jsr("https://planet4589.org/space/gcat/tsv/launch/launch.tsv")
        df_sites = self.scrape_jsr("https://planet4589.org/space/gcat/tsv/tables/sites.tsv")

        jsr_data_dict = {}
        # satcat (main database), auxcat (should be in main but isn't), lcat (suborbital stages/objects). 
        # rcat (lower stages and fairings), ecat (capsules from ISS / crewed missions), ftocat (failed launches).
        files = ["satcat","auxcat","lcat","rcat","ecat","ftocat"]
        for file in files:
            url = "https://planet4589.org/space/gcat/tsv/cat/" + file + ".tsv"
            catalog = self.scrape_jsr(url)
            catalog = catalog[catalog["Type"].str[0].isin(["P","C","S","X","Z"])] # P = Payload, C = Component, S = Suborbital payload, X = Catalog entry that has been deleted, Z = Spurious catalog entry.
            catalog = catalog[catalog["Name"].str.lower().str.contains("starlink|oneweb|yinhe|lynk|e-space|protosat-1|kuiper|tranche|ronghe|digui|hulianwang|qianfan", na=False)]
            jsr_data_dict[file] = catalog

        # Filter the launches to remove unwanted entries.
        df.drop(0, inplace=True)
        df = df[df["#Launch_Tag"].str[:4].astype(int).between(self.start_year, self.final_year)]

        # Example dictionary mapping codes to orders
        code_to_order = {
            "ALP": 1, "BET": 2, "GAM": 3, "DEL": 4, "EPS": 5, "ZET": 6, "ETA": 7,
            "THE": 8, "IOT": 9, "KAP": 10, "LAM": 11, "MU": 12, "NU": 13, "XI": 14,
            "OMI": 15, "PI": 16, "RHO": 17, "SIG": 18, "TAU": 19, "UPS": 20,
            "PHI": 21, "CHI": 22, "PSI": 23, "OME": 24}
        
        def convert_launch_tag(tag):
            tag = tag.strip()
            if '-' in tag:
                return tag
            
            # Otherwise, assume format "YYYY (greek)"
            parts = tag.split()
            year = parts[0]
            mult = 0
            greek_letter = ""
            if len(parts) == 2:
                mult = 0
                greek_letter = parts[1]
            elif len(parts) == 3:
                greek_letter = parts[2]
                if parts[1] == "A":
                    mult = 1
                elif parts[1] == "B":
                    mult = 2
                else:
                    raise ValueError("Unexpected Launch Tag")

            number = code_to_order.get(greek_letter, 0) +mult*24
            return f"{year}-{str(number).zfill(3)}"          

        df['#Launch_Tag'] = df['#Launch_Tag'].apply(convert_launch_tag)

        df = df[~df["#Launch_Tag"].str.contains("A")]   # No endoatmospheric.
        df = df[~df["#Launch_Tag"].str.contains("C")]   # No 'not a launch'.
        df = df[~df["#Launch_Tag"].str.contains("E")]   # No pad explosions.
        df = df[~df["#Launch_Tag"].str.contains("S")]   # No suborbital.
        df = df[~df["#Launch_Tag"].str.contains("Y")]   # No suborbital (obselete old notation)
        df = df[~df["#Launch_Tag"].str.contains("M")]   # No mesosphere.
        df = df[~df["#Launch_Tag"].str.contains("W")]   # No mesosphere (obselete old notation)
        df = df[~df["LaunchCode"].str.startswith("X")]  # No launches from another planetary body or satellite.   
        df = df[~df["LaunchCode"].str.startswith("M")]  # No military launches. 
        df = df[~df["LaunchCode"].str.startswith("H")]  # No high-altitude sounding rockets.
        df = df[~df["LV_Type"].str.contains("Aerobee|R-UNK|Trailblazer", na=False)] # Skip sounding rockets.
        df = df[~df["#Launch_Tag"].str.contains("2013-U01", na=False)] # Skip sounding rockets.

        df = df.reset_index()
        
        # Filter for apogees above 50 km. Doesn't look like this is necessary.
        #df.loc[df["Apogee"].str.strip() == "-", "Apogee"] = 0
        #df = df[df["Apogee"].astype(int) >= 50] # No pad explosions
        
        print("Found",len(df),"launches from JSR for the years",self.start_year,"to",self.final_year)

        ####################################################################
        # Extract the required launch activity data and convert to a dict.
        ####################################################################

        mcs_tags = set()
        for file in files:
            catalog = jsr_data_dict[file]
            mcs_tags.update(catalog["Launch_Tag"].dropna().str.strip())

        def process_launches(df, df_sites):
            launches = []
            for n, row in tqdm(df.iterrows(), total=len(df), desc="Processing rows"):
                
                date = row["Launch_Date"].split()
                yearstr = str(date[0])
                monstr  = str(datetime.strptime(date[1].replace("?", ""), '%b').month).zfill(2) if len(date) > 1 else "01"
                daystr  = date[2].replace("?", "").zfill(2) if len(date) > 2 else "01"
                datestr = yearstr+monstr+daystr
                time_utc = np.float64(int(date[3].replace("?","")[0:2]) + int(date[3].replace("?","")[2:4]) / 60) if len(date) >= 4 else -1

                df_site = df_sites[df_sites["#Site"] == row["Launch_Site"].replace("?", "")]
                
                # Look in the payload catalogs to find which launches contain megaconstellations.
                mcs_flag = row["#Launch_Tag"].strip() in mcs_tags

                if row["Variant"] == "?":
                    variant = "-"
                else:
                    variant = row["Variant"]

                if "Chang Zheng" in row["LV_Type"]:
                    name = row["LV_Type"].replace("Chang Zheng","Long March (CZ)")
                else:
                    name = row["LV_Type"]

                launches.append({
                    "COSPAR_ID":              row["#Launch_Tag"].strip(),
                    "Time(UTC)":              round(time_utc, 2),
                    "Date":                   datestr,
                    "Site":                   df_site["Name"].values[0],
                    "Latitude":               float(df_site["Latitude"].values[0].strip()),
                    "Longitude":              float(df_site["Longitude"].values[0].strip()),
                    "Rocket_Name":            name,
                    "Rocket_Variant":         variant,
                    "DISCOSweb_Rocket_ID":    0,
                    "Megaconstellation_Flag": mcs_flag
                })

            return launches

        self.launches = process_launches(df, df_sites) 
   
    def launch_info_to_netcdf(self,source):
        """This saves the launch information as a NetCDF file for later processing for GEOS-Chem.
        """        
        
        if not self.launches:
            print("No launch data to save.")
            return

        dims = ('launches',)

        def get_field(field):
            return [launch.get(field, np.nan) for launch in self.launches]

        fields = ["COSPAR_ID","Time(UTC)","Date","Longitude","Latitude","Site","Rocket_Name","Rocket_Variant","DISCOSweb_Rocket_ID","Megaconstellation_Flag"]

        # Create DataArrays dynamically
        data_arrays = {}
        for field in fields:
            var_name = field.replace("_", " ").title()
            attrs = dict(long_name=var_name, short_name=var_name)
            data_arrays[field] = xr.DataArray(get_field(field), dims=dims, attrs=attrs)
        
        ds = xr.Dataset(data_arrays)

        # Add units only to specific variables
        ds['Longitude'].attrs['units'] = "Degrees"
        ds['Latitude'].attrs['units']  = "Degrees"
        ds['Date'].attrs['units']      = "YYYYMMDD"
        ds['Time(UTC)'].attrs['units'] = "UTC"
             
        #Save to file and close the dataset     
        ds.to_netcdf(f'./databases/launch_activity_data_{self.start_year}-{self.final_year}_{source}.nc')
        ds.close()
        
    def handle_stage_info(self, count, stage, i, temp_dict, unique_vehicle_name_list):
        
        ##########################
        # Manually rearrage stages
        ##########################
        
        stage_number = ""
        # Sometimes the stages are not in the right order on DISCOSweb, so we need to manually rearrange.
        if unique_vehicle_name_list[i] in ["Angara A5", "Angara A5 Persei", "Angara A5 Orion", "Ariane 5ECA", "Ariane 62", 
                                           "Atlas V N22 v2020", "Atlas V 401 v2020", "Atlas V 411", "Atlas V 421 v2021",
                                           "Atlas V 511 v2020", "Atlas V 531 v2020", 
                                           "Atlas V 541", "Atlas V 541 v2020", "Atlas V 551", "Atlas V 551 v2020",
                                           "Delta 4H", "Epsilon-2 CLPS", "H-IIA 202", "H-IIA 204", "H-IIB", "H-III 22" "GSLV Mk II ", "GSLV Mk III",
                                           "PSLV", "PSLV-DL", "PSLV-XL", "PSLV-CA", "Space Launch System - Block 1 Crew", "Falcon Heavy",
                                           "Long March (CZ) 11", "Long March (CZ) 2F", "Long March (CZ) 3B","Long March (CZ) 3B/YZ-1",
                                           "Long March (CZ) 3C", "Long March (CZ) 5", "Long March (CZ) 5B", "Long March (CZ) 5/YZ-2", 
                                           "Long March (CZ) 6", "Long March (CZ) 6A", "Long March (CZ) 7", "Long March (CZ) 7A", "Long March (CZ) 8",
                                           "Long March (CZ) 8A",
                                           "Minotaur 1", "Pegasus XL", "Vega C",
                                           "Soyuz-2-1A", "Soyuz-2-1A Fregat-M", "Soyuz-2-1B Fregat", "Soyuz-2-1B Fregat-M", "Soyuz-ST-A Fregat-M",
                                           "Soyuz-2-1A Fregat","Soyuz-2-1B","Soyuz-2-1V Volga","Soyuz-ST-B Fregat-MT",
                                           "NK Kerolox LV","Qaem-100", "Shavit 2"]:
            
            stage_name = stage["attributes"]["name"]

            if "booster" in stage_name.lower() or stage_name in ["EPS P241","Soyuz-FG Blok-B,V,G,D","K2-1","GSLV-SOM","Falcon 9 Merlin-1D+","S-200"]:
                stage_number = "Booster"
                
            # For stage ids 773 and 775, this is stages with the same name for the Delta 4H rocket.  
            # NK Kerolox LV has second and third stage literally just called "Stage 2" and "Stage 3", so to avoid complicating other rockets we hard code this. Same for third stage of Qaem-100.
            # For stage ids 114183 and 114190, this is stages with the same name for the Shavit 2 rocket.  
            elif (stage_name in ["URM-1","EPC H173", "Atlas V CCB", "H-II LE-7A", "H-II widebody LE-7A", "PS1", "Soyuz-FG Blok-A", 	"Soyuz-2.1V Blok-A",
                                 "Long March (CZ) 11 Stage 1", "L-186 (YF21B)", "L-172 (YF21C)", "L-186 (YF21C)", "CZ-5-500", "L-165 (H5-1)", "L-76 (YF100)",
                                 "GSLV-1", "Minuteman 1", "ORION 50SXL","SLS-B1 First Stage","Falcon Heavy Centre Core","L-110","JIE-3 first stage",
                                 "NK Kerolox LV","Rafe","LE-9"]) or stage["id"] in ["775","114183","117814"]:
                stage_number = "Stage1"
                
            elif (stage_name in ["URM-2","ESC-A","Centaur-5 SEC","Centaur-5 DEC", "H-II LE-5B", "PS2", "Blok-I", "Long March (CZ) 11 Stage 2", 
                                 "L-84 (YF26)", "L-45 (YF-24E)", "CZ-5-HO", "L-15 (YF115)", "L-15 (YF-115) (7)", "K3-2 short",
                                 "M-35","GSLV-2","SR19","SLS  iCPS","Zefiro 40C","Falcon Heavy second stage","JIE-3 second stage","C-25",
                                 "Salman","L-15 (YF-115) (6A)","LE-5B-3","Stage 2 - Upper Liquid Propulsion Module","WB"]) or stage["id"] in ["773","116151","114190"]:
                stage_number = "Stage2"
                
            elif (stage_name in ["Briz-M", "PS3", "Fregat-M", "Fregat", "Fregat-MT", "Long March (CZ) 11 Stage 3", "PBV",
                                "Blok-DM-3","KM-V2c","GS3 (GSLV Mk II)","Zefiro 9","Volga (Soyuz-2-1V Volga)","RSA-3-3","YZ-2"]) or stage["id"] in ["112791","116154","116721"]:
                stage_number = "Stage3"
                
            elif (stage_name in ["PSLV fourth stage (PS4)", "Long March (CZ) 11 Stage 4","CLPS","AVUM+","YZ-1"]) or stage["id"] == "112794":
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
            elif stage_name == "P120C":
                if unique_vehicle_name_list[i] == "Vega C":
                    stage_number = "Stage1"
                else:
                    stage_number = "Booster"

            
            # Long March often uses the same stage for different rockets, but not always at the same position.    
            elif stage_name == "H-18 (Long March (CZ) YF)":
                if unique_vehicle_name_list[i] in ["Long March (CZ) 3B", "Long March (CZ) 3B/YZ-1", "Long March (CZ) 3C", "Long March (CZ) 7A"]:
                    stage_number = "Stage3"
                elif unique_vehicle_name_list[i] == "Long March (CZ) 8":
                    stage_number = "Stage2"

            elif stage_name == "K3-1":
                if unique_vehicle_name_list[i] in ["Long March (CZ) 5", "Long March (CZ) 5B", "Long March (CZ) 5/YZ-2"]:
                    stage_number = "Booster"
                elif unique_vehicle_name_list[i] in ["Long March (CZ) 8", "Long March (CZ) 7", "Long March (CZ) 7A", "Long March (CZ) 8A"]:
                    stage_number = "Stage1"
                    
            else:
                print(f"Warning: Error manually assigning stages and boosters. Stage {stage_name} not in check lists.")
        else:
            stage_number = f"Stage{count+1}"
            
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
                    
                    # Check for solid/liquid hybrid fuel for a single stage.
                    if (propellant_info["fuel"] != None) and (propellant_info["oxidiser"] != None) and (propellant_info["solidPropellant"] != None):
                        print("HYBRID ROCKET DETECTED")
            
                    # Now check if its liquid or solid fuel.
                    if (propellant_info["fuel"] != None) and (propellant_info["oxidiser"] != None):
                        temp_dict[f"{stage_number} Propellant Name"] = propellant_info["fuel"] + "/" + propellant_info["oxidiser"]
                    elif propellant_info["solidPropellant"] != None:
                        temp_dict[f"{stage_number} Propellant Name"] = propellant_info["solidPropellant"]
                        
                    # Sort out rockets where the propellant type is missing but is still not equal to None.
                    elif unique_vehicle_name_list[i] == "GSLV Mk II":
                        if stage_number in ['Booster']:
                            temp_dict[f"{stage_number} Fuel Type"] = "Solid"
                        elif stage_number in ['Stage1','Stage2']:
                            temp_dict[f"{stage_number} Fuel Type"] = "Hypergolic"
                        else:
                            temp_dict[f"{stage_number} Fuel Type"] = "Hydrogen"
                    elif unique_vehicle_name_list[i] == "Long March (CZ) 11":
                        temp_dict[f"{stage_number} Fuel Type"] = "Solid"
                    else:
                        print(f"Warning 1: Missing propellant <{temp_dict[f'{stage_number} Propellant Name']}> for {stage_number} of {unique_vehicle_name_list[i]}.")
                        
                        
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
                    elif "Hydrazine" in temp_dict[f"{stage_number} Propellant Name"] or temp_dict[f"{stage_number} Propellant Name"] in ["UH25 (75%UDMH+25%N2H4)/N2O4","Liquid bi-propellant"] :
                        temp_dict[f"{stage_number} Fuel Type"] = "Hypergolic"  
                    elif temp_dict[f"{stage_number} Propellant Name"] in ['LH2 (Liquid Hydrogen)/LOX']:
                        temp_dict[f"{stage_number} Fuel Type"] = "Hydrogen"  
                    elif temp_dict[f"{stage_number} Propellant Name"] in ['LNG (Liquid Natural Gas)/LOX','LCH4 (Liquid Methan)/LOX','Methane/LOX']:
                        temp_dict[f"{stage_number} Fuel Type"] = "Methane"
                        
                # Sort out rockets where the propellant type is missing.  
                elif unique_vehicle_name_list[i] == "Ceres-1":
                    if stage_number in ['Stage1','Stage2','Stage3']:
                        temp_dict[f"{stage_number} Fuel Type"] = "Solid"
                    else:
                        temp_dict[f"{stage_number} Fuel Type"] = "Hypergolic"
                elif unique_vehicle_name_list[i] == "Electron (Curie)":
                    if stage_number in ['Stage2']:
                        temp_dict[f"{stage_number} Propellant Name"] = 'Kerosene/LOX'
                        temp_dict[f"{stage_number} Fuel Type"]       = 'Kerosene'
                    elif stage_number in ['Stage3']:
                        temp_dict[f"{stage_number} Propellant Name"] = "Liquid bi-propellant"
                        temp_dict[f"{stage_number} Fuel Type"]       = 'Kerosene'
                elif unique_vehicle_name_list[i] == "Goche Yeollyo Uju Balsache (GYUB) - TV2":
                    if stage_number in ['Stage1','Stage2']:
                        temp_dict[f"{stage_number} Fuel Type"] = "Solid"
                    elif stage_number in ['Stage3']:
                        temp_dict[f"{stage_number} Propellant Name"] = "N2O4/MMH" # From JSR https://planet4589.org/space/gcat/data/tables/engines.html GYUB4
                        temp_dict[f"{stage_number} Fuel Type"] = "Hypergolic"
                elif unique_vehicle_name_list[i] == "KAIROS" and stage_number == "Stage4":
                    temp_dict[f"{stage_number} Propellant Name"] = "Monopropellant hydrazine" # From JSR https://planet4589.org/space/gcat/data/tables/engines.html Kairos PBS
                    temp_dict[f"{stage_number} Fuel Type"] = "Hypergolic"
                elif unique_vehicle_name_list[i] in ["Kuaizhou-11","Qaem-100 ","Jielong-3","Zhongke 1A"]:
                    temp_dict[f"{stage_number} Fuel Type"] = "Solid"
                else:
                    print(f"Warning 2: Missing propellant <<{temp_dict[f'{stage_number} Propellant Name']}>> for {stage_number} of {unique_vehicle_name_list[i]}.")
                    
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

    def get_rocket_info(self,source):
        
        self.unique_rocket_list,vehicle_ids = [], []
        with xr.open_dataset(f'./databases/launch_activity_data_{self.start_year}-{self.final_year}_{source}.nc', decode_times=False) as ds:
            vehicle_names = ds['Rocket_Name'].values
            vehicle_variants = ds['Rocket_Variant'].values
            if source == "dw":
                vehicle_ids = ds['DISCOSweb_Rocket_ID'].values.astype(int)

        # Clean name list and ensure uniqueness
        if self.start_year >= 2007 and source == "dw": # DISCOSweb only starts calling it Shavit 2 from 2023, whereas JSR, SLR, and wiki say it was from 2007.
            vehicle_names = np.where(vehicle_names == "Shavit", "Shavit 2", vehicle_names)

        vehicles_combined = np.array(list(zip(vehicle_names, vehicle_variants)))
        unique_vehicle_names = np.unique(vehicles_combined, axis=0)

        # Map rocket names to proxy names
        proxy_map = {"dw": {"Astra Rocket 3":     "Electron",
                            "Ceres-1":            "Shavit 2",
                            "Long March (CZ) 11": "Minotaur 1",
                            "Jielong-3":          "Epsilon-2 CLPS",
                            "Zhongke 1A":         "Vega C",
                            "Long March (CZ) 6A": "Long March (CZ) 7A",
                            "Kuaizhou-11":        "Long March (CZ) 6",
                            "Zhuque-2":           "Antares 230",
                            },
                     "jsr": {}} # TODO: Add JSR proxy mappings.

        df_vehicles, df_stages, df_engines = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
        if source == "jsr":
            df_vehicles = self.scrape_jsr("https://planet4589.org/space/gcat/tsv/tables/lvs.tsv") 
            df_stages   = self.scrape_jsr("https://planet4589.org/space/gcat/tsv/tables/stages.tsv")
            df_engines  = self.scrape_jsr("https://planet4589.org/space/gcat/tsv/tables/engines.tsv")
        
        #Loop over all rockets, and pull the information for each.
        for i, (name, variant) in enumerate(unique_vehicle_names):

            #Set up the arrays to hold the rocket info.
            temp_dict = {
                "proxy":                   proxy_map[source].get(name, ""),
                "Booster Number":          0,
                "Fairing Mass":            0,
            }

            for j in range(0,5):
                temp_dict.update({
                    f"Stage{j} Fuel Type":       "",
                    f"Stage{j} Propellant Mass": 0,
                    f"Stage{j} Stage Mass":      0,
                })

            if source == "dw":
                temp_dict["name"] = name
                temp_dict["variant"] = ""
                vehicle_index = np.where(vehicle_names == name)[0][0]
                vehicle_id = vehicle_ids[vehicle_index]    
                # Get the propellant mass info.
                doc = self.server_loop(f'/launch-vehicles/{vehicle_id}/stages',{},f"On rocket {i+1} of {len(unique_vehicle_names)}.")
                if doc:
                    for count,stage in enumerate(doc['data']):
                        temp_dict = self.handle_stage_info(count,stage,i,temp_dict,unique_vehicle_names)

            elif source == "jsr":

                if "Long March (CZ)" in name:
                    name = name.replace("Long March (CZ)","Chang Zheng")

                # First find the vehicle and its stages in the vehicles databases.
                launch_vehicle, launch_variant = name, variant
                temp_dict["variant"] = launch_variant

                df_vehicle = df_vehicles[(df_vehicles["#LV_Name"] == name) & (df_vehicles["LV_Variant"] == launch_variant)].reset_index()

                # Then loop over all stages to extract mass information. 
                for i, row in df_vehicle.iterrows():

                    # Get the stage number and skip fairings and payloads as the mass isn't in JSR.
                    stage_number = row['Stage_No'].strip()
                    if stage_number in ["F","P"]:
                        continue
                    
                    # Get the number of boosters.
                    if stage_number == "0":
                        temp_dict["Booster Number"] = int(row["Multiplicity"])

                    # Skip adapters, ullage motors, kick motors (Start, Proton-M/DM-03 and Molniya 8K78).
                    if row["Stage_Name"] in ["Perekhodnik","SOZ","BOZ","DS"]:
                        continue
                    if "Molniya 8K78" in name and stage_number == "5":
                        stage_number = "4"
                    if (name.startswith("Proton-K/") or name.startswith("Proton-M/D") or name.startswith("UR-500K/")) and stage_number == "5":
                        stage_number = "4"

                    # Fix Saturn V.
                    if name.startswith("Saturn V"):
                        if stage_number == "2":
                            continue # Skip the interstage.
                        if stage_number in ["3","4"]:
                            stage_number = str(int(stage_number)-1)
                    
                    # Find the stage in the stages database.
                    df_stage = df_stages[df_stages["#Stage_Name"] == row["Stage_Name"]]
                    if len(df_stage) != 1:
                        print(f"Warning: Found {len(df_stage)} stages for stage {row['Stage_Name']} of {name}.")
                        continue

                    # Skip air-launched first stages - these are aircraft.
                    if df_stage["Stage_Family"].values[0].strip() == "Air":
                        continue

                    if int(stage_number) > 4: 
                        print("Too many stages for",name,variant,stage_number,"-",row["Stage_Name"])
                        pass
                    if int(stage_number) < 0:
                        if name == "H-II" and launch_variant == "(2S)":
                            continue # This rocket has two boosters but for some reason its duplicated in JSR.
                        print("Extra booster found for",name,variant,stage_number,"-",row["Stage_Name"])
                        pass

                    # TODO: Ariane 44LP and Atlas IIAS had two solid and two liquid boosters, so we need to handle this.
                    # Same for H-IIA 2022 and 2024, but these are all solid boosters.

                    # TODO: Minotaur IV, Minotaur V, NOTS EV1, Start all have five stages, so we need to deal with these.

                    # Get the dry mass and wet mass and handle.
                    dry_mass    = df_stage["Dry_Mass"].values[0]
                    launch_mass = df_stage["Launch_Mass"].values[0]
                    dry_mass    = None if dry_mass    == '-' or pd.isna(dry_mass)    else float(dry_mass)
                    launch_mass = None if launch_mass == '-' or pd.isna(launch_mass) else float(launch_mass)*1000
                    
                    if dry_mass is None:
                        print(f"Missing dry mass for Rocket: {name,variant}, Stage: {stage_number} - {row['Stage_Name']}")
                    if launch_mass is None:
                        print(f"Missing launch mass for Rocket: {name,variant}, Stage: {stage_number} - {row['Stage_Name']}")

                    temp_dict[f"Stage{stage_number} Stage Mass"] = dry_mass
                    if dry_mass is not None and launch_mass is not None:
                        temp_dict[f"Stage{stage_number} Propellant Mass"] = launch_mass - dry_mass

                    # Get the propellant name and type.
                    df_engine = df_engines[df_engines["#Name"] == df_stage["Engine"].values[0]]
                    if ("ZK-1A" not in df_stage["Engine"].values[0]) and ("YL-1" not in df_stage["Engine"].values[0]):
                        if len(df_engine) != 1:
                            print(f"Warning: Found {len(df_engine)} engines for stage {row['Stage_Name']} of {name,variant}.")
                            continue
                    
                    if df_engine["Group"].values[0] == "Solid":
                        fuel_type = "Solid"
                    elif df_engine["Group"].values[0] == "LOX/Methane":
                        fuel_type = "Methane"
                    elif df_engine["Group"].values[0] == "LOX/LH2":
                        fuel_type = "Hydrogen"
                    elif ("Kero" in df_engine["Group"].values[0]) or (df_engine["Group"].values[0] == "NA/Turps"): 
                        # Assuming that nitric acid + turpentine is similar to kerosene.
                        fuel_type = "Kerosene"
                    elif df_engine["Group"].values[0] in ["NTO/UDMH","MonoHyd","NA/UDMH","LOX/UDMH","NTO/Hyd","Green"]:
                        fuel_type = "Hypergolic"
                    else:
                        print("Missing fuel type", df_engine["Group"].values[0],name,variant,stage_number) 
                        fuel_type = None

                    temp_dict[f"Stage{stage_number} Fuel Type"] = fuel_type
            else:
                raise ValueError(f"Unknown source {source} for rocket data.")
            
            if "Chang Zheng" in name:
                temp_dict["name"] = name.replace("Chang Zheng","Long March (CZ)")
            else:
                temp_dict["name"] = name

            temp_dict_before = temp_dict.copy()
            # Handle any incomplete/incorrect propellant/stage mass information.
            temp_dict = update_mass_info(temp_dict, temp_dict["name"], variant)

            changed = False
            changed_key = ""
            for key in temp_dict:
                if key == "Fairing Mass":
                    continue
                if temp_dict_before.get(key) != temp_dict[key]:
                    changed = True
                    continue
            if changed:
                print(f"{name} {variant} changed.")

            # Update the rocket list.   
            self.unique_rocket_list.append(temp_dict) 
        
    def rocket_info_to_netcdf(self,source):
        """
        This saves the propellant information as a NetCDF file for later processing for GEOS-Chem.
        Recommended to ignore if you don't need this.
        """        
        #Set up the dimensions of the netcdf file.
        dims = ('rockets')

        def safe_mass(val):
            return float(val or 0.0)
        
        fields = {
            "Rocket_Name":       ("name", None),
            "Rocket_Variant":    ("variant", None),
            "Booster_No":        ("Booster Number", None),
            "Fairing_Mass":      ("Fairing Mass", "kg"),
            "Proxy_Rocket":      ("proxy", None),
        }

        # Add stages dynamically
        for stage in range(0, 5):  # Stage1 to Stage4
            fields.update({
                f"Stage{stage}_PropMass":  (f"Stage{stage} Propellant Mass", "kg"),
                f"Stage{stage}_Fuel_Type": (f"Stage{stage} Fuel Type", None),
                f"Stage{stage}_StageMass": (f"Stage{stage} Stage Mass", "kg"),
            })

        # Build data arrays in one pass
        data_vars = {}
        for field_name, (rocket_key, units) in fields.items():
            values = []
            for rocket in self.unique_rocket_list:
                val = rocket[rocket_key]
                if "Stage Mass" in rocket_key:
                    val = safe_mass(val)
                values.append(val)

            attrs = {"long_name": rocket_key}
            if units:
                attrs["units"] = units

            data_vars[field_name] = xr.DataArray(values, dims=dims, attrs=attrs)        
        
        # Create an xarray Dataset from the DataArrays.
        ds = xr.Dataset(data_vars)
        
        #Save to file and close the DataSet     
        #ds.to_netcdf(f'./databases/rocket_attributes_{self.start_year}-{self.final_year}_noupdate.nc')
        ds.to_netcdf(f'./databases/rocket_attributes_{self.start_year}-{self.final_year}_{source}.nc')
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
        name_list, mass_list, cosparId_list, objectClass_list, epoch_list, temp_object_list = [], [], [], [], [], []   
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
    All of the functions are inside the import_launches class.
    Call each of the functions using import_launches._function_
    e.g. import_launches.launches_per_year(start_year,end_year)
    """         
    
    # Set up the arguments for each function.
    parser = argparse.ArgumentParser()
    parser.add_argument('-so',   "--source",                   default = "jsr",                                     help="Source (DISCOSweb - dw or JSR - jsr).")
    parser.add_argument('-yl',   "--yearly_launches",          action='store_true',                                 help="Get yearly launches.")
    parser.add_argument('-li',   "--launch_info",              action='store_true',                                 help='Get launch info.')
    parser.add_argument('-sli',  "--save_launch_info",         action='store_true',                                 help='Save launch info.')
    parser.add_argument('-ri',   "--rocket_info",              action='store_true',                                 help='Get rocket info.')
    parser.add_argument('-sri',  "--save_rocket_info",         action='store_true',                                 help='Save launch info.')
    parser.add_argument('-sdwr', "--save_discosweb_reentries", action='store_true',                                 help='Save dw reentries.')
    parser.add_argument('-sy',   "--start_year",               default = "2023", choices=str(np.arange(1957,2026)), help='Start Year (1957-2025).')
    parser.add_argument('-fy',   "--final_year",               default = "2024", choices=str(np.arange(1942,2026)), help='Final Year (1957-2025).')
    args = parser.parse_args()
    
    # Sort out the year range.
    start_year = int(args.start_year)
    final_year = int(args.final_year)
    if start_year > final_year:
        final_year = start_year + 1
    print(f"Processing from year {start_year} to {final_year}.")  
    
    #Loop over all years and run functions depending on input arguments.
    LaunchData = import_launches(start_year,final_year)   
    if args.yearly_launches == True:
        LaunchData.get_launch_list()

    if args.launch_info == True:
        if args.source == "dw":
            LaunchData.get_launch_info()
        else:
            LaunchData.get_launch_info_jsr()
        if args.save_launch_info == True:
            LaunchData.launch_info_to_netcdf(args.source)  

    if args.save_discosweb_reentries == True:
        LaunchData.save_discosweb_reentries() 

    if args.rocket_info == True:
        LaunchData.get_rocket_info(args.source)      
        if args.save_rocket_info == True:
            LaunchData.rocket_info_to_netcdf(args.source)
