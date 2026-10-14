from datetime import datetime
import numpy as np
import xarray as xr
import argparse
import pandas as pd
import requests
from tqdm import tqdm
import sys
from io import StringIO
from update_rocket_launch_data import update_mass_info

class import_launches:
    
    def __init__(self):
        """This function sets up the class."""
        self.session = requests.Session()

    def scrape_jsr(self,url,session):
        # Function to web scrape data and convert to a pandas DataFrame.
        response = session.get(url)
        if response.status_code == 200:
            
            filename = url.split("/")[-1]  # Use the last part of the URL as filename
            with open(f"gcat/{filename}", "w", encoding="utf-8") as f:
                f.write(response.text)
            
            # Convert the content to a file-like object for pandas
            tsv_data = StringIO(response.text)
            # Load the data into a pandas DataFrame
            df = pd.read_csv(tsv_data, delimiter="\t", dtype=object, low_memory=False)
        else:
            raise ImportError(f"Failed to fetch from JSR URL: {url}", response.status_code)
        return df

    def add_surface_altitude(self, demfile):
        """Sample ETOPO 2022 surface elevation (above mean sea level) at each launch site, in km."""
        lats = xr.DataArray([launch["Latitude"]  for launch in self.launches], dims="launches")
        lons = xr.DataArray([launch["Longitude"] for launch in self.launches], dims="launches")

        with xr.open_dataset(demfile) as dem:
            elev = dem["z"].sel(lat=lats, lon=lons, method="nearest").to_numpy()

        elev = np.clip(elev, 0, None) * 1e-3  # Make sure any negative values are clipped.
        for launch, z in zip(self.launches, elev):
            launch["Surface_Altitude"] = round(float(z),3)

    def get_launch_info(self):
        
        ####################################################################
        # Web scrape the launch data from Jonathan McDowell's JSR website.
        ####################################################################

        df       = self.scrape_jsr("https://planet4589.org/space/gcat/tsv/launch/launch.tsv",self.session)
        df_sites = self.scrape_jsr("https://planet4589.org/space/gcat/tsv/tables/sites.tsv",self.session)
        df_orgs  = self.scrape_jsr("https://planet4589.org/space/gcat/tsv/tables/orgs.tsv",self.session)
        df_lp    = self.scrape_jsr("https://planet4589.org/space/gcat/tsv/tables/lp.tsv", self.session)

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

        def convert_date(date):
            
            date = [s.replace("?", "") for s in date]
            yearstr = str(date[0])
            monstr  = str(datetime.strptime(date[1], '%b').month).zfill(2)
            try:
                daystr  = date[2].zfill(2)
            except:
                raise ValueError(f"Error converting date: {date}. Check the input format.")
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

            return date_obj
        
        df.drop(0, inplace=True)                                            # Remove the initial NaN row.
        df['#Launch_Tag'] = df['#Launch_Tag'].apply(convert_launch_tag)     # Convert the launch tag to a standard format.
        df = df[df["LaunchCode"].str.startswith(("O", "D", "N"), na=False)] # Include only orbital, near-orbital and deep space missions.
        df = df[~df["LaunchCode"].str[1].isin(["E", "A"])]                  # No pad explosions or pad aborts at the pad.
        df = df[~df["FailCode"].str.startswith("PE")]                       # No explosions near the pad.
        df = df[~df["#Launch_Tag"].str.contains("2021-U01", na=False)]      # Skip military launches.

        df["LaunchDate"] = df["Launch_Date"].apply(lambda x: convert_date(str(x).split()))
        df = df[df["LaunchDate"] < pd.Timestamp("2026-01-01", tz="UTC")]    # Only include launches before 2026.
        df = df.reset_index()
        print("Found",len(df),"launches from GCAT.")

        ####################################################################
        # Extract the required launch activity data and convert to a dict.
        ####################################################################

        def process_launches(df):
            launches = []
            for n, row in tqdm(df.iterrows(), total=len(df), desc="Processing rows"):
                
                df_site = df_sites[df_sites["#Site"] == row["Launch_Site"].replace("?", "")]
                site_row = df_site.iloc[0]
                    
                def sort_sitenames(site):

                    # Handle groups of equivalent site names
                    group_map = {
                        "Kodiak Launch Complex": [
                            "Kodiak Launch Complex, Kodiak Island, Alaska",
                            "Pacific Spaceport Complex Alaska, Kodiak Island, Alaska",
                        ],
                        "Cape Canaveral": [
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
                        "Wallops Island": [
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

                    # Group lookup
                    for normalized, variants in group_map.items():
                        if site in variants:
                            return normalized

                    # Default: return unchanged
                    return site
                
                # Use the group megaconstellation flag to classify megaconstellation launches.
                mcs_flag = False
                if " " in row["Group"]:
                    idx = row["Group"].index(" ")
                    char = row["Group"][idx + 1]
                    if char.lower() == "m":
                        mcs_flag = True

                if row["Variant"] == "?":
                    variant = "-"
                else:
                    variant = row["Variant"]

                if "/" in row["Agency"]:
                    agency_code = row["Agency"].split("/")[0].replace("?", "").strip()
                else:
                    agency_code = row["Agency"].replace("?", "").strip()

                df_org = df_orgs[df_orgs["#Code"] == agency_code]
                
                if len(df_org) != 1:
                    print(f"Warning: Found {len(df_org)} organizations for code {row['Agency']}.")
                    org_name = row["Agency"].strip()
                else:
                    org_name = df_org.iloc[0]["ShortName"].strip() 

                ascent_site = str(row["Ascent_Site"]).strip()
                if row["Launch_Site"].replace("?", "") == "TYSC":
                    lat = 38.85
                    lon = 111.61
                elif ascent_site != "-":
                    ascent_pad  = str(row["Ascent_Pad"]).strip()
                    lp_row = df_lp[(df_lp["#Site"].str.strip() == ascent_site) & (df_lp["Code"].str.strip() == ascent_pad)]
                    if len(lp_row) == 1:
                        lat = float(lp_row.iloc[0]["Latitude"].strip())
                        lon = float(lp_row.iloc[0]["Longitude"].strip())
                    else:
                        raise ValueError(f"Warning: Found {len(lp_row)} launch pads for ascent site {ascent_site} and ascent pad {ascent_pad}.")
                else:
                    lat = float(site_row['Latitude'].strip())
                    lon = float(site_row['Longitude'].strip())

                launches.append({
                    "COSPAR_ID":              row["#Launch_Tag"].strip(),
                    "Date":                   convert_date(row["Launch_Date"].split()),
                    "Site":                   sort_sitenames(df_site["Name"].values[0]),
                    "Latitude":               lat,
                    "Longitude":              lon,
                    "Rocket_Name":            row["LV_Type"],
                    "Rocket_Variant":         variant.replace("?",""),
                    "Megaconstellation_Flag": mcs_flag,
                    "Agency":                 org_name        
                })

            return launches

        self.launches = process_launches(df) 
        self.add_surface_altitude("./input_files/ETOPO_2022_v1_60s_N90W180_surface.nc")
   
    def launch_info_to_netcdf(self):
        """This saves the launch information as a NetCDF file for later processing for GEOS-Chem.
        """        
        
        if not self.launches:
            print("No launch data to save.")
            return

        dims = 'launches'

        def get_field(field):
            values = [launch[field] for launch in self.launches]

            if field == "Date":
                values = np.array([v.tz_convert(None) for v in values], dtype="datetime64[ns]")
            
            elif field in ["COSPAR_ID","Site","Rocket_Name","Rocket_Variant","Agency"]:
                values = np.array(values, dtype=str)

            return values

        fields = ["COSPAR_ID","Date","Longitude","Latitude","Surface_Altitude","Site","Rocket_Name","Rocket_Variant","Megaconstellation_Flag","Agency"]

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
        ds['Surface_Altitude'].attrs['units'] = "km"
             
        #Save to file and close the dataset     
        ds.to_netcdf("./databases/launch_activity_data.nc")
        df_out = ds.to_dataframe().sort_values("COSPAR_ID")
        df_out.to_csv("./databases/launch_activity_data.csv", index=False)
        ds.close()

    def get_rocket_info(self):
        
        self.unique_rocket_list = []
        with xr.open_dataset(f'./databases/launch_activity_data.nc', decode_times=False) as ds:
            vehicle_names = ds['Rocket_Name'].values
            vehicle_variants = ds['Rocket_Variant'].values

        vehicles_combined = np.array(list(zip(vehicle_names, vehicle_variants)))
        unique_vehicle_names = np.unique(vehicles_combined, axis=0)

        # Load in the stage event altitude data from the file.
        stage_alt_dict = {}
        df_alts = pd.read_csv("./input_files/launch_event_altitudes.csv", dtype={"Rocket": str, "Variant": str})
        for _, r in df_alts.iterrows():
            for stage in ["BECO", "MECO", "SEI1", "SECO"]:
                val = pd.to_numeric(r[stage], errors="coerce")
                stage_alt_dict[f"{r['Rocket']} {r['Variant']} {stage}"] = np.float64(val)

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
        _ = self.scrape_jsr("https://planet4589.org/space/gcat/tsv/tables/lv.tsv",self.session) 
        df_vehicles = self.scrape_jsr("https://planet4589.org/space/gcat/tsv/tables/lvs.tsv",self.session) 
        df_stages   = self.scrape_jsr("https://planet4589.org/space/gcat/tsv/tables/stages.tsv",self.session)
        df_engines  = self.scrape_jsr("https://planet4589.org/space/gcat/tsv/tables/engines.tsv",self.session)
        count,count2 = 0,0
        #Loop over all rockets, and pull the information for each.
        for i, (name, variant) in enumerate(unique_vehicle_names):

            #Set up the arrays to hold the rocket info.
            temp_dict = {
                "proxy":          "",
                "Booster Number": 0,
                "Fairing Mass":   0,
                "name":           name,
                "variant":        variant,
            }

            for j in range(-1,6):
                temp_dict.update({
                    f"Stage{j} Fuel Type":       "",
                    f"Stage{j} Propellant Mass": 0,
                    f"Stage{j} Stage Mass":      0,
                })

            # First find the vehicle and its stages in the vehicles databases.
            df_vehicle = df_vehicles[(df_vehicles["#LV_Name"] == name) & (df_vehicles["LV_Variant"] == variant)].reset_index()

            # Then loop over all stages to extract mass information. 
            for i, row in df_vehicle.iterrows():

                # Get the stage number.
                stage_number = row['Stage_No'].strip()
                if stage_number == "A":
                    print(f"{name} is air-launched.")
                if stage_number in ["F","P","C","A"]: # Skip fairings, payloads, components (interstages), air launches. 
                    continue
                else:
                    stage_number = int(stage_number)

                # Skip Post-Boost Propulsion System (Start rocket, uses N2).
                if row["Stage_Name"] in ["DS"]:
                    continue
                
                # Find the stage in the stages database.
                df_stage = df_stages[df_stages["#Stage_Name"] == row["Stage_Name"]]
                if len(df_stage) != 1:
                    raise ValueError(f"Warning: Found {len(df_stage)} stages for stage {row['Stage_Name']} of {name}.")

                if stage_number > 5: 
                    raise ValueError("Too many stages for",name,variant,stage_number,"-",row["Stage_Name"])
                elif stage_number < 0:
                    print("Negative stage number for",name,variant,stage_number,"-",row["Stage_Name"])

                # Get the dry mass and wet mass and handle.
                dry_mass    = df_stage["Dry_Mass"].values[0].replace("?","").strip()
                launch_mass = df_stage["Launch_Mass"].values[0].replace("?","").strip()
                dry_mass    = None if dry_mass    == '-' or pd.isna(dry_mass)    else float(dry_mass)
                launch_mass = None if launch_mass == '-' or pd.isna(launch_mass) else float(launch_mass)*1000
                if stage_number > 0 and int(row["Multiplicity"]) != 1 and name not in ["Conestoga 1620", "NOTS EV1"]:
                    # NOTS EV1 is air-launched, and just has multiple engines on first and second stages.
                    # Conestoga 1620 is modular so has multiples for 1st and 2nd stages.
                    print(f'Potential booster found for {name} {variant} - stage {stage_number} x {int(row["Multiplicity"])}')
                if stage_number == 0:
                    temp_dict["Booster Number"] = int(row["Multiplicity"])
                if dry_mass != None:
                    dry_mass    = dry_mass * int(row["Multiplicity"])
                if launch_mass != None:
                    launch_mass = launch_mass * int(row["Multiplicity"])   
                
                if launch_mass is not None and dry_mass is None:
                    count+=1
                else:
                    count2+=1 
                
                # Going to use a constant mass ratio when we have the launch mass but not the dry.
                # Jonathan McDowell said this is a reasonable estimate.
                if launch_mass is not None and dry_mass is None:
                    #print(name,variant)
                    dry_mass = launch_mass / 10
                elif row["Stage_Name"] in ["Atlas Booster","Atlas B Booster","Atlas D Booster","Atlas E/F Booster"]:
                    # Atlas boosters were fed propellant from the main stage, so they essentially have zero propellant mass (dry mass = launch mass in GCAT).
                    # We shall include them as a stage with zero propellant mass, as we need the dry mass for the re-entry emissions.
                    pass
                elif row["Stage_Name"] in ["NOTS 8in","NOTS 3in Sph."]:
                    # These have a zero launch mass.
                    launch_mass = dry_mass * 10 
                elif launch_mass == dry_mass and not (launch_mass is None or launch_mass == 0):
                    print(f'Warning: Launch Mass is equal to Dry Mass for {row["Stage_Name"]} - {launch_mass}')
                
                # These are fixed later in update_mass_info, so suppress warnings here.
                if (name not in ["Angara-1.2","Chang Zheng 2C/YZ-1S","Chang Zheng 12","Chang Zheng 12A","Strela","Nuri","Juno II","Jupiter C"] 
                    and row["Stage_Name"] not in [
                        "YZ-1",   # Using DW.
                        "YZ-1A",  # Using DW.
                        "YZ-1S",  # Using DW.
                        "YZ-2",   # Using DW.
                        "YZ-3",   # DW dry and wet are equal, so basing off YZ-1S
                        "CZ-K2-3" # Using NB.
                        ]): 
                    
                    if not dry_mass:
                        print(f"Missing dry mass for Rocket:    {name,variant}, Stage: {stage_number} - {row['Stage_Name']}")
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
                    # These catch when the engine is not in the database but the propellant type is known from other sources.
                    # Its the ZK-1A and YL-1 engines.
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
        print(f"Used launch/dry mass ratio for {count}/{count+count2} stages.")
        
    def rocket_info_to_netcdf(self):
        """
        This saves the propellant information as a NetCDF file for later processing for GEOS-Chem.
        Recommended to ignore if you don't need this.
        """        
        #Set up the dimensions of the netcdf file.
        dims = 'rockets'

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
        ds.to_netcdf(f'./databases/rocket_attributes.nc')
        ds.to_dataframe().sort_values(["Rocket_Name", "Rocket_Variant"]).to_csv("./databases/rocket_attributes.csv", index=False)
        ds.close()
                         
if __name__ == "__main__":
    """The main running of the program goes here. 
    All of the functions are inside the import_launches class.
    Call each of the functions using import_launches._function_
    e.g. import_launches.launches_per_year(start_year,end_year)
    """         

    current_year = (datetime.now().year)
    
    # Set up the arguments for each function.
    parser = argparse.ArgumentParser()
    parser.add_argument('-li',   "--launch_info",              action='store_true',                                 help='Get launch info.')
    parser.add_argument('-sli',  "--save_launch_info",         action='store_true',                                 help='Save launch info.')
    parser.add_argument('-ri',   "--rocket_info",              action='store_true',                                 help='Get rocket info.')
    parser.add_argument('-sri',  "--save_rocket_info",         action='store_true',                                 help='Save launch info.')
    args = parser.parse_args()
    
    print(f"Processing from 1957 to {current_year}.")  
    
    #Loop over all years and run functions depending on input arguments.
    LaunchData = import_launches()   

    if args.launch_info == True:
        LaunchData.get_launch_info()
        if args.save_launch_info == True:
            LaunchData.launch_info_to_netcdf()  

    if args.rocket_info == True:
        LaunchData.get_rocket_info()      
        if args.save_rocket_info == True:
            LaunchData.rocket_info_to_netcdf()
