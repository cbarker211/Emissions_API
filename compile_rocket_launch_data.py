from datetime import datetime
import numpy as np
import xarray as xr
import argparse
import pandas as pd
import requests
from tqdm import tqdm
import sys

from python_modules.web_scrape_func import server_request, response_error_handler, scrape_jsr
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
        
        ####################################################################
        # Web scrape the launch data from Jonathan McDowell's JSR website.
        ####################################################################

        df       = scrape_jsr("https://planet4589.org/space/gcat/tsv/launch/launch.tsv",self.session)
        df_sites = scrape_jsr("https://planet4589.org/space/gcat/tsv/tables/sites.tsv",self.session)

        jsr_data_dict = {}
        # satcat (main database), auxcat (should be in main but isn't), lcat (suborbital stages/objects). 
        # rcat (lower stages and fairings), ecat (capsules from ISS / crewed missions), ftocat (failed launches).
        files = ["satcat","auxcat","lcat","rcat","ecat","ftocat"]
        for file in files:
            url = "https://planet4589.org/space/gcat/tsv/cat/" + file + ".tsv"
            catalog = scrape_jsr(url,self.session)
            catalog = catalog[catalog["Type"].str[0].isin(["P","C","S","X","Z"])] # P = Payload, C = Component, S = Suborbital payload, X = Catalog entry that has been deleted, Z = Spurious catalog entry.
            catalog = catalog[catalog["Name"].str.lower().str.contains("starlink|tintin|oneweb|yinhe|lynk|e-space|protosat-1|kuiper|tranche|ronghe|digui|hulianwang|qianfan|lingxi|whdw", na=False)]
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
        df = df[~df["LV_Type"].str.contains("NOTS EV1", na=False)] # These were failed launches from an expendable launch system / anti-satellite weapon.
        df = df[~df["#Launch_Tag"].str.contains("2013-U01", na=False)] # Skip sounding rockets.
        df = df[~df["#Launch_Tag"].str.contains("2021-U01", na=False)] # Skip military launches.

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
                date = [s.replace("?", "") for s in date]
                yearstr = str(date[0])
                monstr  = str(datetime.strptime(date[1], '%b').month).zfill(2) if len(date) > 1 else "01"
                daystr  = date[2].zfill(2) if len(date) > 2 else "01"
                if len(date) >= 4:
                    hrstr   = date[3][0:2].zfill(2)
                    minstr  = date[3][2:4].zfill(2)
                    secstr  = date[3][5:7].zfill(2) if ":" in date[3] else "00"
                else:
                    hrstr   = "00"
                    minstr  = "00"
                    secstr  = "00"
                datestr = f"{yearstr}-{monstr}-{daystr}T{hrstr}:{minstr}:{secstr}Z"
                date_obj = pd.to_datetime(datestr, errors="coerce", utc=True)

                df_site = df_sites[df_sites["#Site"] == row["Launch_Site"].replace("?", "")]

                def sort_sitenames(site):
                    site_map = {
                        "Centre Spatial Guyanais, Kourou, Guyane Francaise": "Guiana Space Center (Kourou)",
                        "Satish Dhawan Space Ctr, Sriharikota, Andhra Pradesh, India": "Sriharikota Space Center",
                        "Rocket Lab Launch Complex 1, Onenui Station, Mahia Peninsula": "Rocket Lab Launch Complex 1, Mahia Peninsula",
                        "Wenchang Space Center, Hainan": "Wenchang Satellite Launch Center",
                        "Xichang Space Center (Songlin), Sichuan, China": "Xichang Satellite Launch Center",
                        "Tanegashima Space Center, Tanegashima, Nippon": "Tanegashima Space Center",
                        "Naro Space Center (Naro Uju Senteo),Oenaro I,GoHeung, Jeollanam-do,Korea": "Naro Space Center",
                        "Israeli Air Force Test Range, Palmachim Beach, Israel": "Yavne Launch Facility (Palmachim)",
                        "Huang Hai CZ-11 launch zone": "China Sea Launch",
                        "NASA John F. Kennedy Space Center, Florida": "Kennedy Space Center (ETR)",
                        "Semnan missile launch site, Iran": "Semnan",
                        "Taiyuan weixing fashe zhongxin": "Taiyuan SLC (Wuzhai)",
                        "Sohae Launch Site, Tongch'ang-dong, Pyongang-bukdo (N Pyongan Prov), N Korea": "Sohae Satellite Launching Station",
                        "Jiuquan Space Center, Nei Monggol Zizhiqu, China": "Jiuquan SLC (Shuang Cheng Tzu)",
                        "Imam Reza Space Ctr, Damghan (Shahroud), Iran (Dasht-E-Kabir?)": "Shahroud Missile Test Site",
                    }

                    # Handle groups of equivalent site names
                    group_map = {
                        "Kodiak Launch Complex": [
                            "Kodiak Launch Complex, Kodiak Island, Alaska",
                            "Pacific Spaceport Complex Alaska, Kodiak Island, Alaska",
                        ],
                        "Spaceport Florida Authority": [
                            "Spaceport Florida, Cape Canaveral",
                            "Space Florida, Cape Canaveral",
                            "Cape Canaveral Air Station, Florida"
                        ],
                        "Kagoshima Space Center (Uchinoura)": [
                            "Kagoshima Space Center, Kagoshima, Kyushu, Nippon",
                            "Uchinoura Space Center, Kagoshima (formerly Kagoshima)"
                        ],
                        "Vandenberg Space Force Base": [
                            "South Base, Vandenberg Space Force Base, California",
                            "Vandenberg Space Force Base, California",
                            "Vandenberg AFB, California",
                            "South Vandenberg AFB, California",
                            "Naval Missile Facility, Point Arguello, California"
                        ],
                        "Mid-Atlantic Regional Spaceport": [
                            "Mid-Atlantic Regional Spaceport, Wallops Island, Virginia",
                            "Wallops Flight Facility, Wallops Island, Virginia",
                            "Wallops Island Main Base, NASA Wallops Flight Facility, Chincoteague, Virginia"
                        ],
                        "Baikonur Cosmodrome (Tyuratam)": [
                            "NIIP-5, Baykonur, Kazakstan",
                            "GIK-5, Baykonur, Kazakstan"
                        ],
                        "Kapustin Yar MSC": [
                            "GTsP-4, Kapustin Yar, Volgograd, Rossiya",
                            "GTsMP-4 MO RF, Znamensk (Kapustin Yar), Rossiya"
                        ],
                        "Vostochny Cosmodrome": [
                            "Vostochniy, Svobodniy, Amurskaya Oblast', Rossiya",
                            "GIK-2, Svobodniy, Amurskaya Oblast', Rossiya"
                        ],
                        "Plesetsk Cosmodrome": [
                            "GNIIP, Plesetsk, Rossiya",
                            "1-y Gosudarstvenniy Ispitatelniy Kosmodrom MO RF",
                            "GNIIP, VKS section, Plesetsk, Rossiya",
                            "53-y Nauchno-Issledovatelskiy Ispitatelniy Poligon,"
                        ]
                    }

                    # Direct lookup
                    if site in site_map:
                        return site_map[site]

                    # Group lookup
                    for normalized, variants in group_map.items():
                        if site in variants:
                            return normalized

                    # Default: return unchanged
                    return site
                
                # Look in the payload catalogs to find which launches contain megaconstellations.
                mcs_flag = row["#Launch_Tag"].strip() in mcs_tags

                if row["Variant"] == "?":
                    variant = "-"
                else:
                    variant = row["Variant"]

                name = row["LV_Type"]

                if "Proton-M" in name and date_obj >= pd.Timestamp("2007-07-07", tz="UTC"):
                    variant = "Enhanced"

                launches.append({
                    "COSPAR_ID":              row["#Launch_Tag"].strip(),
                    "Date":                   date_obj,
                    "Site":                   sort_sitenames(df_site["Name"].values[0]),
                    "Latitude":               float(df_site["Latitude"].values[0].strip()),
                    "Longitude":              float(df_site["Longitude"].values[0].strip()),
                    "Rocket_Name":            name,
                    "Rocket_Variant":         variant.replace("?",""),
                    "Megaconstellation_Flag": mcs_flag
                })

            return launches

        self.launches = process_launches(df, df_sites) 
   
    def launch_info_to_netcdf(self):
        """This saves the launch information as a NetCDF file for later processing for GEOS-Chem.
        """        
        
        if not self.launches:
            print("No launch data to save.")
            return

        dims = ('launches',)

        def get_field(field):
            values = [launch[field] for launch in self.launches]

            if field == "Date":
                values = np.array([v.tz_convert(None) for v in values], dtype="datetime64[ns]")
            
            elif field in ["COSPAR_ID", "Site", "Rocket_Name", "Rocket_Variant"]:
                values = np.array(values, dtype=str)

            return values

        fields = ["COSPAR_ID","Date","Longitude","Latitude","Site","Rocket_Name","Rocket_Variant","Megaconstellation_Flag"]

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
             
        #Save to file and close the dataset     
        ds.to_netcdf(f'./databases/launch_activity_data_{self.start_year}-{self.final_year}.nc')
        ds.close()

    def get_rocket_info(self):
        
        self.unique_rocket_list,vehicle_ids = [], []
        with xr.open_dataset(f'./databases/launch_activity_data_{self.start_year}-{self.final_year}.nc', decode_times=False) as ds:
            vehicle_names = ds['Rocket_Name'].values
            vehicle_variants = ds['Rocket_Variant'].values

        vehicles_combined = np.array(list(zip(vehicle_names, vehicle_variants)))
        unique_vehicle_names = np.unique(vehicles_combined, axis=0)

        # Map rocket names to proxy names
        proxy_map = {"Astra Rocket 3":     "Electron",
                     "Ceres-1":            "Shavit 2",
                     "Long March (CZ) 11": "Minotaur 1",
                     "Jielong-3":          "Epsilon-2 CLPS",
                     "Zhongke 1A":         "Vega C",
                     "Long March (CZ) 6A": "Long March (CZ) 7A",
                     "Kuaizhou-11":        "Long March (CZ) 6",
                     "Zhuque-1":           "Super Strypi",
                     "Zhuque-2":           "Antares 230"}

        # Load in the stage event altitude data from the file.
        stage_alt_dict = {}
        stage_alt_rockets = np.genfromtxt("./input_files/launch_event_altitudes.csv",dtype=str,skip_header=1,usecols=[0,1],delimiter=",")
        stage_alt_data = np.genfromtxt("./input_files/launch_event_altitudes.csv",dtype=np.float64,skip_header=1,usecols=[2,3,4,5],delimiter=",")

        stages = ["BECO", "MECO", "SEI1", "SECO"]
        for (name, variant), row in zip(stage_alt_rockets, stage_alt_data):
            for stage, value in zip(stages, row):
                stage_alt_dict[f"{name} {variant} {stage}"] = None if value == "" else np.float64(value)

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

        df_vehicles, df_stages, df_engines = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
        df_vehicles = scrape_jsr("https://planet4589.org/space/gcat/tsv/tables/lvs.tsv",self.session) 
        df_stages   = scrape_jsr("https://planet4589.org/space/gcat/tsv/tables/stages.tsv",self.session)
        df_engines  = scrape_jsr("https://planet4589.org/space/gcat/tsv/tables/engines.tsv",self.session)

        #Loop over all rockets, and pull the information for each.
        for i, (name, variant) in enumerate(unique_vehicle_names):

            #Set up the arrays to hold the rocket info.
            temp_dict = {
                "proxy":                   proxy_map.get(name, ""),
                "Booster Number":          0,
                "Fairing Mass":            0,
            }

            for j in range(-1,6):
                temp_dict.update({
                    f"Stage{j} Fuel Type":       "",
                    f"Stage{j} Propellant Mass": 0,
                    f"Stage{j} Stage Mass":      0,
                })

            # First find the vehicle and its stages in the vehicles databases.
            temp_dict["name"] = name
            temp_dict["variant"] = variant

            # This is because Proton-M Enhanced is only listed under the name Proton-M in JSR.
            if "Proton-M" in name and variant == "Enhanced":
                df_vehicle = df_vehicles[(df_vehicles["#LV_Name"] == name)].reset_index()
            else:
                df_vehicle = df_vehicles[(df_vehicles["#LV_Name"] == name) & (df_vehicles["LV_Variant"] == variant)].reset_index()

            # Then loop over all stages to extract mass information. 
            for i, row in df_vehicle.iterrows():

                # Get the stage number and skip fairings and payloads as the mass isn't in JSR.
                stage_number = row['Stage_No'].strip()
                if stage_number in ["F","P"]:
                    continue
                if name.startswith(("Soyuz", "Conestoga 1620", "Polyot", "Molniya 8K78", "Sputnik","Voskhod","Vostok", "Space Shuttle", "SLS Block 1")) and name != "Soyuz-2-1V":
                    stage_number = int(stage_number) - 1
                else:
                    stage_number = int(stage_number)

                # Skip adapters, ullage motors, kick motors (Start, Proton-M/DM-03 and Molniya 8K78).
                if row["Stage_Name"] in ["Perekhodnik","SOZ","BOZ","DS"]:
                    continue
                if "Molniya 8K78" in name and stage_number == 4:
                    stage_number = 3
                if (name.startswith("Proton-K/") or name.startswith("Proton-M/D") or name.startswith("UR-500K/")) and stage_number == 5:
                    stage_number = 4
                
                # Fix Saturn V.
                if name.startswith("Saturn V"):
                    if stage_number == 2:
                        continue # Skip the interstage.
                    if stage_number in [3,4]:
                        stage_number = stage_number-1
                
                # Find the stage in the stages database.
                df_stage = df_stages[df_stages["#Stage_Name"] == row["Stage_Name"]]
                if len(df_stage) != 1:
                    print(f"Warning: Found {len(df_stage)} stages for stage {row['Stage_Name']} of {name}.")
                    continue

                # Skip air-launched first stages - these are aircraft.
                if df_stage["Stage_Family"].values[0].strip() == "Air":
                    print("Air-launched - ",name,variant,stage_number,"-",row["Stage_Name"])
                    continue

                if stage_number > 5: 
                    print("Too many stages for",name,variant,stage_number,"-",row["Stage_Name"])
                    pass
                
                if name == "H-II" and variant == "(2S)" and stage_number < 0:
                    continue # This rocket has two boosters but for some reason its duplicated in JSR.

                # Get the dry mass and wet mass and handle.
                dry_mass    = df_stage["Dry_Mass"].values[0]
                launch_mass = df_stage["Launch_Mass"].values[0]
                dry_mass    = None if dry_mass    == '-' or pd.isna(dry_mass)    else float(dry_mass)
                launch_mass = None if launch_mass == '-' or pd.isna(launch_mass) else float(launch_mass)*1000
                if stage_number > 0 and int(row["Multiplicity"]) != 1 and name not in ["Conestoga 1620"]:
                    print(f'Potential booster found for {name} {variant} - stage {stage_number} x {int(row["Multiplicity"])}')
                if stage_number == 0:
                    temp_dict["Booster Number"] = int(row["Multiplicity"])
                if dry_mass != None:
                    dry_mass    = dry_mass * int(row["Multiplicity"])
                if launch_mass != None:
                    launch_mass = launch_mass * int(row["Multiplicity"])                           
                
                # Going to use a constant mass ratio when we have the launch mass but not the dry.
                # For the MG-18 stage (Thor MG-18, Scout X-2M, Scout X-3M), this is within the range of other scout 4th stages (10-33).
                # https://www.planet4589.org/space/book/lv/engines/kick/WIDELYUSEDMOTORS.html 'high-mass ratio'
                # Jonathan McDowell said this is a reasonable estimate.
                if launch_mass is not None and dry_mass is None:
                    #print(name,variant)
                    dry_mass = launch_mass / 10
                elif name == "Atlas D" and stage_number == 0:
                    launch_mass = dry_mass * 10 
                
                # These are fixed later in update_mass_info, so suppress warnings here.
                if name not in ["Angara-1.2","Diamant A","Diamant B","Diamant BP4","Electron","Kuaizhou","Kuaizhou-1A","Lambda 4S","Juno II","Jupiter C",
                                "Chang Zheng 11","Chang Zheng 2C/YZ-1S","Chang Zheng 2B/YZ-3","Chang Zheng 3B/YZ-1","Chang Zheng 3C/YZ-1","Chang Zheng 2D/YZ-3",
                                "Chang Zheng 5/YZ-2","Chang Zheng 5B/YZ-2","Chang Zheng 6","Chang Zheng 7/YZ-1A", "Chang Zheng 12", "Chang Zheng 12A",
                                "Minotaur-C 3210","Nuri","Safir","Shuang Quxian 1","Simorgh","Strela","Taurus 3110","Taurus 3210", "Scout F-1"]:
                    if not dry_mass:
                        print(f"Missing dry mass for Rocket: {name,variant}, Stage: {stage_number} - {row['Stage_Name']}")
                    if not launch_mass:
                        print(f"Missing launch mass for Rocket: {name,variant}, Stage: {stage_number} - {row['Stage_Name']}")

                temp_dict[f"Stage{stage_number} Stage Mass"] = dry_mass
                if dry_mass is not None and launch_mass is not None:
                    temp_dict[f"Stage{stage_number} Propellant Mass"] = launch_mass - dry_mass

                # Get the propellant info from the engine database.
                df_engine = df_engines[df_engines["#Name"] == df_stage["Engine"].values[0]]
                if len(df_engine) == 1:
                    
                    # Assign the propellant type.
                    if (df_engine["Group"].values[0] == "Solid") or (df_engine["Fuel"].values[0].strip() == "Polyethylene") or (df_engine["Fuel"].values[0].strip() == "Paraffin"):
                        # Assuming that hybrid H2O2/Polyethylene or Lox/Paraffin engines are close enough to solid fuel.
                        fuel_type = "Solid"
                    elif df_engine["Group"].values[0] == "LOX/Methane":
                        fuel_type = "Methane"
                    elif df_engine["Group"].values[0] == "LOX/LH2":
                        fuel_type = "Hydrogen"
                    elif ("Kero" in df_engine["Group"].values[0]) or (df_engine["Group"].values[0] in ["NA/Turps","LOX/Propane"]): 
                        # Assuming that nitric acid + turpentine is similar to kerosene.
                        fuel_type = "Kerosene"
                    elif df_engine["Group"].values[0] in ["NTO/UDMH","MonoHyd","NA/UDMH","LOX/UDMH","NTO/Hyd","Green"]:
                        fuel_type = "Hypergolic"
                    else:
                        print("Missing fuel type", df_engine["Group"].values[0],name,variant,stage_number) 
                        fuel_type = None
                elif name in ["Yinli-1","Jielong-3","Lijian-1"]:
                    fuel_type = "Solid"
                else:
                    print(f"Warning: Found {len(df_engine)} engines for stage {row['Stage_Name']} of {name,variant}.")
                    fuel_type = None

                temp_dict[f"Stage{stage_number} Fuel Type"] = fuel_type

            # Handle any incomplete/incorrect propellant/stage mass information.
            temp_dict = update_mass_info(temp_dict, temp_dict["name"], variant)

            # Find the event altitude information and look up configuration
            stage_alts = {key: stage_alt_dict.get(f"{name} {variant} {key}", np.nan) for key in ['BECO', 'MECO', 'SEI1', 'SECO']}

            stages = tuple(bool(temp_dict[f"Stage{i} Fuel Type"]) for i in range(6))
            rocket_config_type = config_map.get(stages) #type: ignore

            if rocket_config_type is None:
                raise IndexError(f"Incorrect rocket configuration for {name,variant} (stages={stages})")
            
            temp_dict["BECO"] = stage_alts['BECO']
            temp_dict["MECO"] = stage_alts['MECO']
            temp_dict["SEI1"] = stage_alts['SEI1']
            temp_dict["SECO"] = stage_alts['SECO']
            temp_dict["Rocket_Config"] = rocket_config_type

            # Update the rocket list.   
            self.unique_rocket_list.append(temp_dict) 
        
    def rocket_info_to_netcdf(self):
        """
        This saves the propellant information as a NetCDF file for later processing for GEOS-Chem.
        Recommended to ignore if you don't need this.
        """        
        #Set up the dimensions of the netcdf file.
        dims = ('rockets')

        def safe_mass(val):
            return float(val or 0.0)
        
        fields = {
            "Rocket_Name":    ("name", None),
            "Rocket_Variant": ("variant", None),
            "Booster_No":     ("Booster Number", None),
            "Fairing_Mass":   ("Fairing Mass", "kg"),
            "Proxy_Rocket":   ("proxy", None),
            "BECO":           ("BECO", None),
            "MECO":           ("MECO", None),
            "SEI1":           ("SEI1", None),
            "SECO":           ("SECO", None),
            "Rocket_Config":  ("Rocket_Config", None),
        }

        # Add stages dynamically
        for stage in range(-1, 6):  # Stage-1 to Stage5
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
        ds.to_netcdf(f'./databases/rocket_attributes_{self.start_year}-{self.final_year}.nc')
        ds.close()
                         
if __name__ == "__main__":
    """The main running of the program goes here. 
    All of the functions are inside the import_launches class.
    Call each of the functions using import_launches._function_
    e.g. import_launches.launches_per_year(start_year,end_year)
    """         
    
    # Set up the arguments for each function.
    parser = argparse.ArgumentParser()
    parser.add_argument('-yl',   "--yearly_launches",          action='store_true',                                 help="Get yearly launches.")
    parser.add_argument('-li',   "--launch_info",              action='store_true',                                 help='Get launch info.')
    parser.add_argument('-sli',  "--save_launch_info",         action='store_true',                                 help='Save launch info.')
    parser.add_argument('-ri',   "--rocket_info",              action='store_true',                                 help='Get rocket info.')
    parser.add_argument('-sri',  "--save_rocket_info",         action='store_true',                                 help='Save launch info.')
    parser.add_argument('-sy',   "--start_year",               default = "1957", choices=str(np.arange(1957,2026)), help='Start Year (1957-2025).')
    parser.add_argument('-fy',   "--final_year",               default = "2025", choices=str(np.arange(1942,2026)), help='Final Year (1957-2025).')
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
        LaunchData.get_launch_info()
        if args.save_launch_info == True:
            LaunchData.launch_info_to_netcdf()  

    if args.rocket_info == True:
        LaunchData.get_rocket_info()      
        if args.save_rocket_info == True:
            LaunchData.rocket_info_to_netcdf()
