import numpy as np
import argparse
import xarray as xr
import pandas as pd
from datetime import datetime, timedelta
import sys
import random
import geopandas as gpd
from shapely.geometry import Point, Polygon, shape
import fiona
import requests
from io import StringIO
from cProfile import Profile
from pstats import Stats
import time
from tqdm import tqdm

def convert_time(date):
    
    """Converts the "DDate" listing from GCAT to day, month, and time strings.

    Returns:
        datestr(str):  The date in YYYYMMDD format.
        time_utc(str): The decimal time in hours.
    """   
    
    yearstr=str(date[0]).replace("?","")
    if len(date) == 1:
        monstr = "01"   
        daystr = "01"
    elif date[1] == "Q3?": # 1995-005A
        monstr = "06"
        daystr = "30"  
    else:   
        date[1] = date[1].replace("?","")
        monstr = str(datetime.strptime(date[1], '%b').month).zfill(2)  
        if len(date) > 2: 
            daystr = str(date[2].replace("?","")).zfill(2)
        elif date == ['2022', 'Apr']:
            # DISCOSweb puts reentry on Mar 31st. Setting as March 31st.
            daystr = "31"
            monstr = "03"
        elif date == ['2022', 'Dec']:
            # DISCOSweb puts reentry on Nov 30th. Setting as Nov 30th.
            daystr = "30"
            monstr = "11"
        else:
            daystr = "01"
    
    datestr = yearstr+monstr+daystr
        
    if len(date) >= 4:    
        time_utc = np.float64(int(date[3].replace("?","")[0:2]) + int(date[3].replace("?","")[2:4]) / 60)
    else:
        time_utc = -1  
    return datestr, time_utc
    
class build_reentry_list:
    
    def __init__(self, start_year, final_year):
        """This function sets up the class and imports data.

        Args:
            start_year, final_year (int): The initial and final years of re-entry to cover.
        """        
        
        self.start_year = start_year
        self.final_year = final_year
        self.session = requests.Session()
        
        # Import the shapefiles.
        if use_gpd == True:

            shapefiles = {
                "Gulf of Mexico": "./databases/reentry/iho/gulf_of_mexico/iho.shp",
                "Indian Ocean":   "./databases/reentry/iho/indian_ocean/iho.shp",
                "Pacific":        "./databases/reentry/iho/pacific/iho.shp",
                "Atlantic":       "./databases/reentry/iho/atlantic/iho.shp",
                "Countries":      "./databases/reentry/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp",
                "States":         "./databases/reentry/ne_50m_admin_1_states_provinces/ne_50m_admin_1_states_provinces.shp"
            }

            self.loaded_shapefiles = {}
            for name, path in tqdm(shapefiles.items(), total=len(shapefiles.keys()), desc="Loading Shapefiles."):
                gdf = gpd.read_file(path, engine="pyogrio")
                self.loaded_shapefiles[name] = gdf
            
    def print_stats(self):
        """Print out statistics about the database at the end of the program.
        """        
        # TODO: This can be cleaned up and output as a pandas dataframe.
        # Set all variable counters to 0, then count how many items are missing location, time, and mass information.
        non_geo_count, non_time_count, non_mass_count, total_mass, geolocated_mass, starlink_count = 0, 0, 0, 0, 0, 0
        for reentry in self.unique_reentry_list:
            
            if reentry["lat"] == 0 and reentry["lon"] == 0:
                non_geo_count += 1
                if "starlink" in reentry["name"].lower():
                    starlink_count += 1      
            else:
                geolocated_mass += (reentry["abl_mass"] + reentry["other_mass"])
            if reentry["time"] == -1:
                non_time_count += 1
            if (reentry["abl_mass"] + reentry["other_mass"]) == 0:
                non_mass_count += 1
                   
            total_mass += (reentry["abl_mass"] + reentry["other_mass"])
        
        print(f"Total Reentries: {len(self.unique_reentry_list):.0f}.")
                
        geo_percent_all = (len(self.unique_reentry_list)-non_geo_count) / (len(self.unique_reentry_list)) * 100
        time_percent_all = (len(self.unique_reentry_list)-non_time_count) / (len(self.unique_reentry_list)) * 100
        mass_percent_all = (len(self.unique_reentry_list)-non_mass_count) / (len(self.unique_reentry_list)) * 100
        print(f"Total mass (Gg):   {total_mass*1e-6:.3f}")
        geolocated_mass_percent = geolocated_mass / total_mass * 100
        print(f"Geolocation (All): {geo_percent_all:.0f}%, {non_geo_count}.")
        print(f"Timed (All):       {time_percent_all:.0f}%, {non_time_count}.")    
        print(f"Mass (All):        {mass_percent_all:.0f}%, {non_mass_count}.")
        print(f"Geolocated Mass:   {geolocated_mass_percent:.0f}%.")
                
    def convert_launch_tag(self,tag):

        # Example dictionary mapping codes to orders
        code_to_order = {"ALP": 1, "BET": 2, "GAM": 3, "DEL": 4, "EPS": 5, "ZET": 6, "ETA": 7,
                            "THE": 8, "IOT": 9, "KAP": 10, "LAM": 11, "MU": 12, "NU": 13, "XI": 14,
                            "OMI": 15, "PI": 16, "RHO": 17, "SIG": 18, "TAU": 19, "UPS": 20,
                            "PHI": 21, "CHI": 22, "PSI": 23, "OME": 24}
        
        tag = str(tag).strip()
        if '-' in tag:
            return tag
        # Otherwise, assume format "YYYY (greek)"
        parts = tag.split()
        year = parts[0]
        mult = 0
        greek_letter = ""

        if len(parts[1]) > 1:
            mult = 0
            greek_letter = parts[1]
        else:
            greek_letter = parts[2]
            if parts[1] == "A":
                mult = 1
            elif parts[1] == "B":
                mult = 2
            else:
                raise ValueError("Unexpected Launch Tag")

        number = code_to_order.get(greek_letter, 0) + mult*24
        
        return f"{year}-{str(number).zfill(3)}" 
        
    def convert_lat_lon(self,latlonstr,inc,category,apogee,jsr_id):
        """Set the geolocation for all items that are not Falcon 9 fairings/1st stage.

        Args:
            latlonstr (str)     : The 
            inc (np.float64)    : The orbital inclination of the object.
            stage (int)         : The category of the object.
            apogee (int)        : The apogee of the object.
            dsl (xarray Dataset): A database of all rocket launches in the specified time range.
            jsr_id (str)        : The COSPAR ID of the object.

        Returns:
            lat(np.float64) : The latitude of the object.
            lon(np.float64) : The longitude of the object.
            location(int) : A toggle for whether the location was
                                1 - Lat/Lon Provided
                                2 - Launch Site / Named Location
                                3 - Political Region
                                4 - Ocean/Continent
                                5 - Falcon Reusable
                                6 - Inclination Bounded Random
                                7 - Non-bounded Random
        """        
        
        # Adjust the inclination so its useful for our purpose of bounding the lat.
        location = -1
        if inc > 90:
            inc = 180-inc
        elif inc == 0:
            inc=90
        #print(jsr_id,latlonstr)

        def sample_within_region(region, lat_range=(-inc, inc)):

            while True:
                coords = region.sample_points(1).get_coordinates()
                lon, lat = coords["x"].values[0], coords["y"].values[0]
                if lat_range[0] <= lat <= lat_range[1]:
                    break   
            return lat, lon             
        
        if use_gpd == True:    
            shapefile_map = {
                "Gujarat":      (self.loaded_shapefiles["States"], "name_en", "Gujarat", 3),
                "Kazakhstan":   (self.loaded_shapefiles["Countries"], "NAME", "Kazakhstan", 3),
                "Mexico":       (self.loaded_shapefiles["Countries"], "NAME", "Mexico", 3),
                "S Africa":     (self.loaded_shapefiles["Countries"], "NAME", "South Africa", 3),
                "N Zealand":    (self.loaded_shapefiles["Countries"], "NAME", "New Zealand", 3),
                "China":        (self.loaded_shapefiles["Countries"], "NAME", "China", 3),

                "Gulf":         (self.loaded_shapefiles["Gulf of Mexico"], "name", "Gulf of Mexico", 4),
                "Indian O":     (self.loaded_shapefiles["Indian Ocean"], "name", "Indian Ocean", 4),
                "Indian Ocean": (self.loaded_shapefiles["Indian Ocean"], "name", "Indian Ocean", 4),
                "Indian O.":    (self.loaded_shapefiles["Indian Ocean"], "name", "Indian Ocean", 4),
                "IOR":          (self.loaded_shapefiles["Indian Ocean"], "name", "Indian Ocean", 4), 
                "S POR":        (self.loaded_shapefiles["Pacific"], "name", "South Pacific Ocean", 4),
                "S Pacific":    (self.loaded_shapefiles["Pacific"], "name", "South Pacific Ocean", 4),
                "N Atlantic":   (self.loaded_shapefiles["Atlantic"], "name", "North Atlantic Ocean", 4),
                "S AOR":        (self.loaded_shapefiles["Atlantic"], "name", "South Atlantic Ocean", 4),
                "S Atl":        (self.loaded_shapefiles["Atlantic"], "name", "South Atlantic Ocean", 4),
            }
        else:
            shapefile_map = {}

        # If location is missing for lower stages for failed and successful launches, set as launch coordinates.    
        if latlonstr in ["-",""]:
            
            if (apogee <= 100) and (category in ["S0","S1"]): 
                for i, cpid in enumerate(self.dsl["COSPAR_ID"].values):
                    if jsr_id[:8] == cpid[:8]:
                        lat = self.dsl["Latitude"].values[i]
                        lon = self.dsl["Longitude"].values[i]
                        location = 2
                        break  
            else: 
                lat = round(random.uniform(-inc, inc),2)
                lon = round(random.uniform(-180, 180),2)
                location = 6  

        # Lat/Lon Provided
        elif latlonstr in ["GM 86.2W 29.7N", "32W 1E"]:
            location = 1
            if latlonstr == "GM 86.2W 29.7N": # Gulf of Mexico named coordinates.
                lat = 29.7
                lon = -86.2
            elif latlonstr == "32W 1E": # TODO: Check with Jonathan about this.
                lon = -32
                lat = 1
        
        # Launch Site / Named Location
        elif latlonstr in ["Alashan","Jiuquan","Dongfeng","LOPNOR RW05","Lop Nor","VSFB RW30/12","Koonibba","Jacklyn","Ocean","STB OLP1","WSSH","Splash","KSC SLF"]:
            location = 2
            if latlonstr == "Alashan": 
                # Landing expected east of Jiuquan launch site (this was actually launched from Wenchang).
                # Assuming that Alashan refers to Helen Shan Mountain (new name for Alashan Mountain).
                # https://space.skyrocket.de/doc_sdat/rcs-fc-sc.htm
                # https://planet4589.org/space/jsr/news.778
                # https://www.peakbagger.com/peak.aspx?pid=10693
                # https://www.nasaspaceflight.com/2020/05/china-next-generation-crew-capsule/
                lat = 38.833
                lon = 105.95
            elif latlonstr in ["Jiuquan","Dongfeng"]:
                # This refers to the Chinese launch center, sources say it landed 'in the desert near Jiuquan', and 'in China’s Inner Mongolia autonomous region'.
                # Using lat/lon for Jiuquan launch site.
                # https://discosweb.esoc.esa.int/launch-sites/21
                # https://nssdc.gsfc.nasa.gov/nmc/spacecraft/display.action?id=2020-027A
                # https://space.skyrocket.de/doc_sdat/xzf-sc.htm 
                # https://spaceflightnow.com/2020/05/08/chinas-next-generation-crew-spacecraft-lands-after-unpiloted-test-flight/ 
                lat = 41.3
                lon = 100.3               
            elif latlonstr in ["LOPNOR RW05","Lop Nor"]:
                # Reported test flight of a Chinese reusable experimental spacecraft.
                # Runway 05/22 at Lop Nor air base in Xinjiang (https://planet4589.org/space/gcat/data/tables/lp.html).
                # https://www.seradata.com/china-launches-own-mini-spaceplane-reusable-spacecraft-using-long-march-2f 
                # https://twitter.com/planet4589/status/1302486141090885632                
                lat = 40.78
                lon = 89.27
            elif latlonstr == "VSFB RW30/12":
                lat = 34.7
                lon = -120.6 
            elif latlonstr == "Koonibba":
                lat = -31.885558
                lon = 133.448686
            elif latlonstr == "Jacklyn": # https://x.com/spaceOFFSHORE/status/1878177979974775139/photo/1
                lat = 27.887984558404042
                lon = -74.15337755494835
            elif latlonstr == "Ocean": # TODO: CHeck this doesn't include others.
                # Electron Stage 1 (rcat) is occasionally recovered in the ocean.
                # From https://spaceflightnow.com/2020/11/05/rocket-lab-to-attempt-booster-recovery-on-next-mission/:
                #   Recovery vessels stationed near the booster’s splashdown zone around 250 miles (400 kilometers) 
                #   south of the launch site will move in to secure the first stage and hoist it onto a ship for return to New Zealand.
                # Geoloating 400km south of launch site.
                lat = -42.86 
                lon = 177.87
            elif latlonstr == "STB OLP1": # Starbase in Texas.
                lat = 25.996198
                lon = -97.154394
            elif latlonstr == "WSSH": # This is the White Sands Missile Range in New Mexico.
                
                lat = 33.238462 
                lon = -106.346383 
            elif latlonstr == "Splash":
                # This is the Electron failed helicopter stage 1 recovery in May 22. 
                # https://www.youtube.com/watch?v=BY0CXlOeWHI "Several hundred kilometers from the launch site."
                # Using inclination of 94 degrees (wiki), and distance of 300km.
                lat = -36.57
                lon = 177.63
            elif latlonstr == "KSC SLF": # https://www.world-airport-codes.com/united-states/nasa-shuttle-landing-facility-69738.html
                # https://www.spaceforce.mil/News/Article/3217077/x-37b-orbital-test-vehicle-concludes-sixth-successful-mission/
                lat = 28.61
                lon = -80.69
        
        # Ocean/Continent
        elif latlonstr in ["Antarctic", "Arctic", "S Ocean","S Oc.","S OCean"]:
            location = 4
            lon = round(random.uniform(-180, 180), 2)

            if latlonstr in ["S Ocean","S Oc.","S OCean"]:
                # For these objects in auxcat, the inclination is 51.65-74, outside of the general definition of the southern ocean (<60S).
                # Therefore we will just set the lat to the inclination.
                lat = -inc
            else:
                # Generate latitude until within target polar bounds
                while True:
                    lat = round(random.uniform(-inc, inc), 2)
                    if (latlonstr == "Antarctic" and lat <= -60) or (latlonstr == "Arctic" and lat >= 66):
                        break
            coordinate = Point(lon, lat)

        elif use_gpd == False and latlonstr in ["Pacific","PO","E Pacific","S Pacific","S POR","POR","SW POR","SE Pacific","W POR","NE RU/NW POR","SW POR","N Pacific",
                                                "Indian O","SE IOR","Indian Ocean","Indian O.","SW IOR","IOR",
                                                "Atlantic","N Atlantic","S AOR","SW AOR","S Atl","AOR","AO",
                                                "Kazakhstan","Gujarat","Gulf","Mexico","S Africa","E USA","N Zealand","SE Arkalyk","China","Africa","KRZ?",
                                                "S of Aus","S of Austr.","S of S Afr."
                                               ]:
            lat = 0
            lon = 0
            location = 0

        # Regions defined by shapefiles.
        elif latlonstr in shapefile_map:
            shp, col, name, location = shapefile_map[latlonstr]
            region = shp[shp[col].isin([name])]
            lat, lon = sample_within_region(region)
         
        elif latlonstr in ["Pacific", "PO", "POR"]:
            location = 4
            while True:
                random_location = self.loaded_shapefiles["Pacific"].sample_points(1).get_coordinates()
                lon = random_location["x"].values[0]
                lat = random_location["y"].values[0]
                if -inc <= lat <= inc:
                    break
        elif latlonstr in ["AO","Atlantic"]: 
            location = 4 
            while True:
                random_location = self.loaded_shapefiles["Atlantic"].sample_points(1).get_coordinates()
                lon = random_location["x"].values[0]
                lat = random_location["y"].values[0]
                if -inc <= lat <= inc:
                    break
        elif latlonstr in ["E Pacific"]:
            location = 4
            while True:
                random_location = self.loaded_shapefiles["Pacific"].sample_points(1).get_coordinates()
                lon = random_location["x"].values[0]
                lat = random_location["y"].values[0]
                if -inc <= lat <= inc and -180 <= lon <= -60: # Bounded to East Pacific only.
                    break 
        elif latlonstr in ["S of Tasmania"]:
            location = 4
            lon_lat_list = [[153.23, -30.02],[146.83, -43.64],[166.00, -50.92], [167.53, -47.29],
                            [168.14, -46.85], [168.85, -46.66], [174.28, -41.72], [175.28, -41.62], [173.01, -34.39]]
            tasman_sea_geom = Polygon(lon_lat_list)
            tasman_sea = gpd.GeoDataFrame(index=[0], crs='epsg:4326', geometry=[tasman_sea_geom])
            while True:
                lat = round(random.uniform(-inc, inc),2)
                lon = round(random.uniform(-180, 180),2)
                coordinate = Point(lon,lat)
                if tasman_sea.geometry.contains(coordinate).any():
                    break 
        elif latlonstr in ["SE IOR"]:
            location = 4
            indian_ocean = self.indian_ocean_shapefile[self.indian_ocean_shapefile.name.isin(['Indian Ocean'])] 
            while True:
                lat = round(random.uniform(-inc, -30),2)
                lon = round(random.uniform(77, 180),2)
                coordinate = Point(lon,lat)
                if indian_ocean.geometry.contains(coordinate).any():
                    break

        elif len(latlonstr.split()) == 4: # Coordinates in degrees-minutes
            geolocation = latlonstr.split()
            if geolocation[3][-1] in ["E","W","S","N"]: # Check its actually coordinates and not some long name.
                location = 1

                if geolocation[1][-1] in ["E","W"]: # If the longitude comes first.
                    lon_data = geolocation[:2]
                    lat_data = geolocation[2:]
                elif geolocation[3][-1] in ["E","W"]: # If the longitude comes last.
                    lon_data = geolocation[2:]
                    lat_data = geolocation[:2]  
                
                if "W" in lon_data[-1]:
                    lon = (np.float64(lon_data[0]) + np.float64(lon_data[1].replace("W",""))/60)*-1
                elif "E" in lon_data[-1]:
                    lon = (np.float64(lon_data[0]) + np.float64(lon_data[1].replace("E",""))/60)
                    
                if "N" in lat_data[-1]:
                    lat = (np.float64(lat_data[0]) + np.float64(lat_data[1].replace("N",""))/60)
                elif "S" in lat_data[-1]:
                    lat = (np.float64(lat_data[0]) + np.float64(lat_data[1].replace("S",""))/60)*-1
            
        elif latlonstr.split()[0][-1] in ["E","W","S","N"]: # Coordinates in degrees
            location = 1
            geolocation = latlonstr.split()
            if len(geolocation) == 1:
                if geolocation[0] == '155W50S':
                    lon_data = "155W"
                    lat_data = "50S"
            else:
                if geolocation[0][-1] in ["E","W"]: # If the longitude comes first.
                    lon_data = geolocation[0]
                    lat_data = geolocation[1]
                elif geolocation[1][-1] in ["E","W"]: # If the longitude comes last.
                    lon_data = geolocation[1]
                    lat_data = geolocation[0]  
                else:
                    raise ValueError(f"Can't find lon/lat coords for {jsr_id} - {latlonstr}.")

            if lon_data[-1] == "W":
                lon = np.float64(lon_data.replace("W",""))*-1
            elif lon_data[-1] == "E":
                lon = np.float64(lon_data.replace("E",""))
            else:
                print("Missing Lon Direction - Assuming W")
                lat = np.float64(lon_data.replace("W",""))*-1

            if lat_data[-1] == "N":
                lat = np.float64(lat_data.replace("N",""))
            elif lat_data[-1] == "S":
                lat = np.float64(lat_data.replace("S",""))*-1  
            else:
                print("Missing Lat Direction - Assuming N")
                lat = np.float64(lat_data) 

        else:
            lat = 0
            lon = 0
            #print(f"Need to sort out geolocation for {latlonstr}")
        
        try:
            lat = np.float64(lat) 
            lon = np.float64(lon)
        except:
            raise ValueError(f"Problem converting lat/lon for >{latlonstr}< - {jsr_id}.")

        if location == -1:
            pass
            #print(f"Missing location tag for {latlonstr}")
            
        return lat, lon, location
    
    def falcon_stage_lat_lon(self,datestr,jsr_id,jsr_dest):
        
        """Set the geolocation for all Falcon 9 1st stages.

        Returns:
            lat(np.float64) : The latitude of the object.
            lon(np.float64) : The longitude of the object..
        """  
        
        # TODO: Add more Falcon Heavy
        falcon_heavy_list = ["2022-144","2023-008","2023-108"]
        # Launches "2023-157","2023-210","2024-119" are also falcon heavy, but not in Raul's SpaceX map.
        
        # Find the ground landings in Raul's SpaceX Map
        if "LZ" in jsr_dest:
            matching = self.ground_landings[self.ground_landings["Date"].isin([f"{datestr[:4]}-{datestr[4:6]}-{datestr[6:]}"])].reset_index(drop=True)
            if matching.shape[0] == 1: # If its in the list of ground landings, then we can geolocate it.
                lat = np.array(matching["geometry"].y)[0]
                lon = np.array(matching["geometry"].x)[0]    
            elif jsr_id in falcon_heavy_list: # This is for Falcon Heavy where two boosters land on the ground.
                lat = np.array(matching["geometry"].y)[0]
                lon = np.array(matching["geometry"].x)[0]
            # If its not in the list, then as its a ground landing, we can set it to the launch site.
            elif matching.shape[0] == 0:
                found_launch = False
                for i in range(len(self.dsl["COSPAR_ID"])):
                    if self.dsl["COSPAR_ID"].values[i] == jsr_id:
                        found_launch = True
                        lat = self.dsl["Latitude"].values[i]
                        lon = self.dsl["Longitude"].values[i]
                if found_launch == False:
                    sys.exit(f"Problem geolocating Falcon Stage 1 recovery - {jsr_id}.")
            elif matching.shape[0] > 1:
                print(matching)
                sys.exit(f"Multiple ground entries for Falcon Stage 1 landing- {jsr_id}.")
            else:
                sys.exit(f"Problem geolocating Falcon Stage 1 landing- {jsr_id}.")
        else:
            # Find the ocean landings in Raul's SpaceX Map.
            matching = self.ocean_landings[self.ocean_landings["Date"].isin([f"{datestr[:4]}-{datestr[4:6]}-{datestr[6:]}"])].reset_index(drop=True)
            if (matching.shape[0] == 1) or (f"{datestr[:4]}-{datestr[4:6]}-{datestr[6:]}" == "2022-04-27"):
                # This is 27th April - typo in Raul's Space Map where second should be 2023. 
                lat = np.array(matching["geometry"].y)[0]
                lon = np.array(matching["geometry"].x)[0]
            # These next two are where there are two launches on the same day.
            elif jsr_id == "2022-124":
                lat = np.array(matching["geometry"].y)[0]
                lon = np.array(matching["geometry"].x)[0]
            elif jsr_id == "2022-125":
                lat = np.array(matching["geometry"].y)[1]
                lon = np.array(matching["geometry"].x)[1]
            # Again these are two on the same day.
            elif jsr_id == "2023-037":
                lat = np.array(matching["geometry"].y)[0]
                lon = np.array(matching["geometry"].x)[0]
            elif jsr_id == "2023-038":
                lat = np.array(matching["geometry"].y)[1]
                lon = np.array(matching["geometry"].x)[1]
            # If its not in the list, then set it to the ocean near the launch site. 
            # These coordinates are the most common grid box for reentries in 2020-2022.
            elif matching.shape[0] == 0:
                found_launch = False
                for i in range(len(self.dsl["COSPAR_ID"])):
                    if self.dsl["COSPAR_ID"].values[i] == jsr_id:
                        found_launch = True
                        # Kennedy Space Center (ETR)
                        if (np.round(self.dsl["Longitude"].values[i],0) == -81) and (np.round(self.dsl["Latitude"].values[i],0) in [28,29]):
                            lat = 34
                            lon = -75  
                        # Vandenberg
                        elif (np.round(self.dsl["Longitude"].values[i],0) == -121) and (np.round(self.dsl["Latitude"].values[i],0) == 35):
                            lat = 30
                            lon = -120 
                        else:
                            print(self.dsl["Longitude"].values[i],self.dsl["Latitude"].values[i])
                            sys.exit("Falcon launch site not found.")
            elif matching.shape[0] > 1:
                print(matching)
                sys.exit(f"Multiple ocean entries for Falcon Stage 1 landing- {jsr_id}.")
            else:
                sys.exit(f"Problem geolocating Falcon Stage 1 landing- {jsr_id}.") 
            
        return lat, lon 
    
    def falcon_fairing_lat_lon(self,datestr,jsr_id):   
        
        """Set the geolocation for all Falcon 9 fairings.

        Returns:
            lat(np.float64) : The latitude of the object.
            lon(np.float64) : The longitude of the object..
        """
        set_ocean = False
        if len(self.fairings) > 0:
            matching = self.fairings[self.fairings["Date"].isin([f"{datestr[:4]}-{datestr[4:6]}-{datestr[6:]}"])].reset_index(drop=True)
            if matching.shape[0] == 1:
                lat = np.array(matching["geometry"].y)[0]
                lon = np.array(matching["geometry"].x)[0]
            # These are two on the same day
            elif jsr_id == "2023-037":
                lat = np.array(matching["geometry"].y)[0]
                lon = np.array(matching["geometry"].x)[0]
            elif jsr_id == "2023-038":
                lat = np.array(matching["geometry"].y)[1]
                lon = np.array(matching["geometry"].x)[1]
            # If its not in the list, then set it to the ocean near the launch site. 
            # These coordinates are the most common grid box for reentries in 2020-2022.
            elif matching.shape[0] == 0:
                set_ocean = True
            elif matching.shape[0] > 1:
                print(matching)
                sys.exit(f"Multiple entries for Falcon fairing recovery - {jsr_id}.")
            else:
                sys.exit(f"Problem geolocating Falcon fairing recovery - {jsr_id}.")
        else:
            set_ocean = True

        if set_ocean == True:
            found_launch = False
            for i in range(len(self.dsl["COSPAR_ID"])):
                if self.dsl["COSPAR_ID"].values[i] == jsr_id:
                    found_launch = True
                    # Kennedy Space Center (ETR)
                    if (np.round(self.dsl["Longitude"].values[i],0) == -81) and (np.round(self.dsl["Latitude"].values[i],0) in [28,29]):
                        lat = 34
                        lon = -75  
                    # Vandenberg
                    elif (np.round(self.dsl["Longitude"].values[i],0) == -121) and (np.round(self.dsl["Latitude"].values[i],0) == 35):
                        lat = 30
                        lon = -120
                    else:
                        print(self.dsl["Longitude"].values[i],self.dsl["Latitude"].values[i])
                        sys.exit("Falcon launch site not found.")
            if found_launch == False:
                sys.exit(f"Problem geolocating Falcon fairing recovery - {jsr_id}.")
              
        return lat, lon
    
    def failed_launch_mass(self, jsr_id, jsr_name, reentry_category, drymass):
        
        """For all failed launches, set the masses manually as we need to differentiate wet vs dry mass reentry.

        Returns:
            mass(float)       : The mass of aluminium.
            other_mass(float): Extra non-aluminium mass.
        """        
            
        abl_mass = 0
        other_mass = 0

        rocket_ind = np.where(self.dsr["Rocket_Name"].values == self.dsl["Rocket_Name"].where(
                              self.dsl["COSPAR_ID"] == jsr_id, drop=True).values[0])[0][0]

        if reentry_category == "S0":
            abl_mass = self.dsr["Stage0_StageMass"].values[rocket_ind] / int(self.dsr["Booster_No"].values[rocket_ind])
        elif reentry_category in ["S1","S2","S3","S4"]:
            abl_mass = self.dsr[f"Stage{reentry_category[1]}_StageMass"].values[rocket_ind]
        elif "fairing" in jsr_name.lower():
            abl_mass = self.dsr["Fairing_Mass"].values[rocket_ind]  / 2    
        elif reentry_category in ["C","P"]:
            
            # 2020
            if jsr_name == "Zafar":
                abl_mass = 113
            elif jsr_name == "Xinjishu Yanzheng 6":
                abl_mass = 6550
            elif jsr_name == "Nusantara Dua":
                abl_mass = 5550
            elif jsr_name == "CE-SAT-I":
                abl_mass = 67
            elif "Flock 4e" in jsr_name:
                abl_mass = 5.8
            elif jsr_name == "Faraday-1":
                abl_mass = 10
            elif jsr_name == "Jilin-1 Gaofen 02E":
                abl_mass = 201    
            elif jsr_name == "Xiangrikui 2":
                abl_mass = 97
            elif jsr_name == "Neimonggol 1":
                abl_mass = 230
            elif jsr_name == "SEOSat-Ingenio":
                abl_mass = 750
            elif jsr_name == "Taranis":
                abl_mass = 175 
            elif jsr_name == "Astra Test Payload":
                abl_mass = 5

            # 2021   
            elif jsr_name == "Global-10":
                abl_mass = 60
            elif jsr_name == "Global-11":
                abl_mass = 60
            elif jsr_name == "Tolou-2?":
                abl_mass = 90
            elif jsr_name == "Jilin-1 Mofang 01":
                abl_mass = 18
            elif jsr_name == "GISAT-1":
                abl_mass = 2286
            elif jsr_name == "Wiseongmosache":
                abl_mass = 1500
            elif jsr_name == "GeeSat-1A":
                abl_mass = 130
            elif jsr_name == "GeeSat-1B":
                abl_mass = 130
            elif jsr_name == "Test payloads":
                abl_mass = 90
                
            # 2022
            elif jsr_name == "R5-S1":
                abl_mass = 1
            elif jsr_name == "INCA":
                abl_mass = 3.8
            elif jsr_name == "QubeSat":
                abl_mass = 4
            elif jsr_name == "BAMA-1":
                abl_mass = 4
            elif jsr_name == "Jilin-1 Mofang 01A/R":
                abl_mass = 18
            elif jsr_name == "TROPICS SV02":
                abl_mass = 5.34
            elif jsr_name == "TROPICS SV04":
                abl_mass = 5.34
            elif jsr_name == "EOS-02":
                abl_mass = 145
            elif jsr_name == "AzaadiSAT":
                abl_mass = 8
            elif jsr_name == "ESMS":
                abl_mass = 200
            elif jsr_name == "RAISE-3":
                abl_mass = 110
            elif jsr_name == "Amateru-I":
                abl_mass = 170
            elif jsr_name == "Amateru-II":
                abl_mass = 170
            elif jsr_name == "E-SSOD 1":
                abl_mass = 10
            elif jsr_name == "E-SSOD 2":
                abl_mass = 10
            elif jsr_name == "E-SSOD 3":
                abl_mass = 10
            elif jsr_name == "MITSUBA":
                abl_mass = 1.7
            elif jsr_name == "WASEDA-SAT-ZERO":
                abl_mass = 1.2
            elif jsr_name == "MAGNARO A":
                abl_mass = 3
            elif jsr_name == "MAGNARO B":
                abl_mass = 1.5
            elif jsr_name == "KOSEN-2":
                abl_mass = 2.7
            elif jsr_name == "FSI-SAT":
                abl_mass = 1.4
            elif jsr_name == "Zhixing 1B":
                abl_mass = 50       
            elif "Unknown payload" in jsr_name:
                abl_mass = 10
            elif jsr_name == "Pleiades Neo 5":
                abl_mass = 920
            elif jsr_name == "Pleiades Neo 6":
                abl_mass = 920
            else:
                abl_mass = np.float64(drymass)

        # Add second stage propellant mass to non aluminium mass.    
        if jsr_id in ["2020-F02","2020-F05",
                      "2021-F02",
                      "2022-F01","2022-F02","2022-F03",
                      "2023-F04","2023-F05","2023-F06","2023-F07","2023-F09",
                      "2024-F02","2024-F04"] and reentry_category == "S2":
            other_mass = self.dsr[f"Stage2_PropMass"].values[rocket_ind]
            if jsr_id == "2022-F03":
                other_mass = other_mass * 0.24
                
        # Add third stage propellant mass to non aluminium mass.    
        if jsr_id in ["2020-F02","2020-F03","2020-F05","2020-F06",
                      "2021-F02","2021-F06","2021-F09","2021-F10",
                      "2022-F02","2022-F05","2022-F07",
                      "2023-F07","2023-F08","2023-F09","2023-F10",
                      "2024-F02","2024-F04","2024-F05"] and reentry_category == "S3":   
            other_mass = self.dsr[f"Stage3_PropMass"].values[rocket_ind]
            if jsr_id == "2021-F09":
                other_mass = other_mass * 0.08
                    
        # Add fourth stage propellant mass to non aluminium mass.  
        if jsr_id in ["2020-F08","2020-F09",
                      "2021-F10",
                      "2022-F02","2022-F04","2022-F05","2022-F07",
                      "2023-F10",
                      "2024-F04"] and reentry_category == "S4":   
            other_mass = self.dsr[f"Stage4_PropMass"].values[rocket_ind]

        return abl_mass, other_mass
        
    def sort_inclination(self,jsr_inc,jsr_id):
        
        """When the inclination on file is zero (including failed), this function sets it appropriately.
        This is mainly the case for lower suborbital stages, which can then be set using the inclination of an upper stage or payload.

        Returns:
            jsr_inc(np.float64): The orbital inclination of the object.
        """        
        
        # Try to find an entry from satcat/auxcat that is already loaded.
        # NOTE: Occasionally, the inclinations are slightly different for each object.
        # The difference is usually <1 degree, but can be higher. As we don't know which inclination is 'correct',
        # we use the larger one to bound the geolocation.   

        for reentry in (self.unique_reentry_list):
            if reentry["id"][:8] == jsr_id[:8] and reentry["inc"] != 0:
                if reentry["inc"] > jsr_inc:
                    jsr_inc = reentry["inc"]
                    break            
            
        # Sometimes, there is no entry in the reentry list because it didn't meet the criteria before.
        # In this case we need to reload the satcat to look for any matching entries.
        if jsr_inc == 0:

            if "F" in jsr_id:
                jsr_data = self.jsr_data_dict["ftocat"]
            else:
                keys = ["satcat", "rcat", "auxcat", "ecat"]
                jsr_data = pd.concat([self.jsr_data_dict[k] for k in keys], ignore_index=True)

            jsr_data = jsr_data[jsr_data["Converted_Tag"].str[:8] == jsr_id[:8]]
            jsr_data = jsr_data[jsr_data["Inc"] != "-"]
            jsr_data = jsr_data[jsr_data["Inc"].astype(np.float64) != 0]
            jsr_data = jsr_data[jsr_data["Status"].isin(["O","R","DSO","DSA","AF","AS","F","S","GRP","AO","AR"])].reset_index(drop=True)
            
            if len(jsr_data) > 0:
                for i in range(len(jsr_data)):
                    if np.float64(jsr_data["Inc"][i]) > jsr_inc:
                        jsr_inc = np.float64(jsr_data["Inc"][i])  
                        break

        if jsr_inc == 0 and jsr_id[5] not in ["F","U"]:
            # TODO: Sort failed launches.
            print(f"inc still empty for {jsr_id}")
        return jsr_inc 
         
    def extract_jsr_info(self, jsr_data):
        """For each GCAT database, extract all the relevant information and output to the dictionary of all reentry objects.
        """        
        
        ################################################
        # Filter the database for the relevant entries.
        ################################################

        subset = jsr_data[jsr_data["Status"].isin(["R","R?","D","L","L?","S","F","AF","AS"])]           # Only include specific statuses (see https://planet4589.org/space/gcat/web/intro/phases.html)
        subset = subset[~subset["Launch_Tag"].str[5].isin(["A","C","E","S","Y","M","W"])]               # Exclude items from launches we ignore.
        subset = subset[~subset["Type"].str[0].isin(["Z","D"])]                                         # Ignore Z, this is a spurious entry according to JSR, and ignore D, this is debris objects. 
        subset = subset[subset["Primary"] == "Earth"]                                                   # Only objects whose primary body is Earth.
        subset = subset[subset["Piece"].str[5:6] != "S"]                                                # No suborbital launches (mainly military rockets).
        subset = subset[~subset["Piece"].isin(["2021-U01","2022-U01","2022-U03","2023-U01",             # Skipping military tests (mostly North Korea).
                                               "2023-U02","2023-U03","2023-U04","2024-U05"])]           # TODO: Will have to add more military launches here.
        subset = subset[~subset["Piece"].isin(["1961-U02","1962-U02","1964-U01","1964-U05",             # Skipping sounding rockets (Trailblazer).
                                               "1965-U03","1966-U05","1966-U07","1967-U01"])]                                                                
        subset["Apogee"] = subset["Apogee"].replace("-", 0).mask(subset["Apogee"].str.contains("Inf", na=False), 0).astype(int) # Sort apogees.
        subset = subset[~subset["#JCAT"].isin(["L80508","L80509","L80510","S57807"])]                   # Skipping 2021-F07 and H2A fairing debris object.
        subset = subset[subset["DDate"].str[0:4].ne("-")]                                               # Remove objects that haven't re-entered yet.
        subset = subset[subset["DDate"].str[0:4].astype(int).between(self.start_year, self.final_year)] # This time range only.
        subset = subset[(subset["Apogee"] >= 50) | (subset["Apogee"] == -54771)]                        # Apogee above 50 km, but include Hayabusa 2 Return Capsule.
        subset = subset[(subset["Status"] != "AS") | (subset["Piece"].str[5:6].isin(["F","U"]))]        # Only include AS re-entries if they are part of a failed or uncategorized launch.

        # Note: A negative apogee means that the object was in a hyperbolic orbit (velocity above escape velocity), so passes by Earth.
        # However, Hayabusa2 (apogee=-54771) did reenter on 5th Dec 2020, so including this.  

        # Work out how many launches there should be in the final year.
        # This is for in case we want to run an incomplete year and want to stop early.
        # TODO: Need to figure out a way this can work with the old notations before COSPAR.
        #self.dsl_id_list = [int(dsl_id[-3:]) for dsl_id in self.dsl["COSPAR_ID"].values if "F" not in dsl_id and "U" not in dsl_id]
            
        # Now loop over the list and format into a dictionary.
        for i, row in subset.iterrows():

            # First sort out COSPAR-IDs for pre-1962 launches.
            jsr_id = row["Piece"]
            jsr_id = self.convert_launch_tag(jsr_id)
            jsr_name   = str(row["Name"])
            jsr_inc    = np.float64(row["Inc"])
            jsr_dest   = row["Dest"]
            jsr_apogee = int(row["Apogee"])
            
            if str(jsr_name) == "Hayabusa 2 Return Capsule":
                jsr_apogee = jsr_apogee * -1
    
            # This is where we skip objects where the date is wrong in GCAT. 
            # TODO: Email Jonathan about this.
            # These are added manually from Aerospace Corp or DISCOSweb.
            if row["#JCAT"] in ["S46138","S44635","S40899"]:
                continue
            
            # TODO: Sort this when we want to do up to date launches.
            #if jsr_id[:4] == "2025" and jsr_id[5] != "F" and int(jsr_id[5:8]) > max(self.dsl_id_list):
            #    continue
            
            burnup = "Complete"
            # Set the burnup as partial for all objects landing or splashing down at the surface.
            if row["Status"] == "L":
                burnup = "Partial"
            # NOTE: Need to manually check for Electron booster recoveries.
            elif (jsr_name == "Electron Stage 1") and jsr_id in ["2019-084","2020-007","2020-085","2021-F02","2021-106",
                                                                 "2022-047","2022-147","2023-041","2023-100","2023-126","2024-022"]:
                    burnup = "Partial"
                 
            # Duplicates:
            #   2019-036H. This is TEPCE 1 and TEPCE 2, listed with same "Piece" in JSR. TEPCE 2 listed as 2019-036IA in DISCOSweb, so renaming to this.
            #   2020-027A. This is Xinyidai Zairen Feichuan (XZF) and XZF Service Module (not in DISCOSweb). Setting auxcat id to XZF 2020-027D.
            #   Remaining are where multiple objects in auxcat from same launch are just given the COSPAR ID of the launch.
            # TODO: Email Jonathan about these.
            if row["Name"] == "TEPCE 2":
                jsr_id = "2019-036IA" 
            elif row["Name"] == "XZF Service Module":
                jsr_id = "2020-027D"
                
            # Sort out the reentry time/date.
            datestr, time_utc = convert_time(row["DDate"].split()) 
            
            # Handle categories.
            if row["Type"][0] == "R":
                # Sort out the rocket stages. This is mainly Soyuz (Russian's just do the numbering differently), and errors in the GCAT.
                reentry_category = f"S{(row['Type'][1:2])}"

                # Setting Falcon Heavy Boosters as S0.
                if reentry_category == "S1" and "Falcon 9 Stage 1" in jsr_name:
                    if row["Parent"] != "-":
                        reentry_category = "S0"

                if row["#JCAT"] in ["R81714","R81715"]: # Falcon Heavy
                    reentry_category = "S0"
                elif row["Name"] in ["RSRMV-1L","RSRMV-1R"]: # SLS
                    reentry_category = "S0"   

                # Soyuz 
                elif "Blok-BVGD" in row["Bus"]:
                    reentry_category = "S0"
                elif "Blok-A" in row["Bus"]:
                    reentry_category = "S1"
                elif "Blok-I" in row["Name"]:
                    reentry_category = "S2"
                elif "Fregat" in row["Name"]:
                    reentry_category = "S3"

            elif row["Type"][0] in ["C","P"]:
                reentry_category = row["Type"][0]
            else:
                reentry_category = -1
                
            # Sort the inc.
            if jsr_inc == 0:
                jsr_inc = self.sort_inclination(jsr_inc,jsr_id)

            # Sort out the reentry lat/lon.
            if ("Falcon 9 Stage 1" in jsr_name):
                if row["Status"] == "L":
                    lat, lon = self.falcon_stage_lat_lon(datestr,jsr_id,jsr_dest) 
                    location = 5
                elif (row["Status"] in ["S","D"] or # Sometimes the 1st stage is expended to suit the mission.
                      row["Piece"] == "2025-021"):  # This one is listed as R for some reason.
                    if row["Status"] == "D":
                        print(f"Warning: Assuming missing landing of {row['#JCAT']}")
                    lat, lon, location = self.convert_lat_lon(jsr_dest.replace("?",""), jsr_inc, reentry_category, jsr_apogee, jsr_id)
                else:
                    sys.exit(f"Unexpected Falcon Stage 1 Status {row['Status']} {row['#JCAT']}")
            elif ("Falcon 9 Fairing" in jsr_name) or ("Falcon Heavy Fairing" in jsr_name):    
                lat, lon = self.falcon_fairing_lat_lon(datestr,jsr_id)
                burnup = "Partial"
                location = 5
            else:   
                lat, lon, location = self.convert_lat_lon(jsr_dest.replace("?",""), jsr_inc, reentry_category, jsr_apogee, jsr_id)
            
            # Sort the mass.    
            # # TODO: Sort failed launches back to 1957. 
            #if jsr_id[5:6] in ["F","U"]:
            #    abl_mass, other_mass = self.failed_launch_mass(jsr_id, jsr_name, reentry_category,row["DryMass"])
            #else:
            #    abl_mass = np.float64(row["DryMass"])
            #    other_mass = 0

            abl_mass = np.float64(row["DryMass"])
            other_mass = 0
            
            # Check for items with a missing geolocation (should have been dealt with already so this is just a sanity check).    
            if use_gpd == True and lat == 0 and lon == 0:
                print(f"No location for {jsr_id}, dest = {row['Dest']}")
            
            # Set up the dictionary.                      
            temp_reentry_dict = {
                "id"               : str(jsr_id),
                "jcat"             : str(row["#JCAT"]),
                "name"             : str(jsr_name),
                "category"         : reentry_category,
                "burnup"           : burnup,
                "time"             : time_utc,
                "datestr"          : datestr,
                "lat"              : lat,
                "lon"              : lon,
                "abl_mass"         : abl_mass,
                "other_mass"       : other_mass,
                "attached_abl_mass": 0,
                "location"         : location,
                "inc"              : jsr_inc,
                "apogee"           : jsr_apogee,
            }
            
            self.unique_reentry_list.append(temp_reentry_dict)   

    def add_cargo(self,jsr_data):
        
        """Short function to add the cargo mass to its parent object. jsr_data is ecat.
        """        
        
        jsr_data = jsr_data[jsr_data["Status"].isin(["AL IN","AR IN"])].reset_index(drop=True)
        
        for reentry_count in range(len(jsr_data)):           

            if jsr_data["#JCAT"][reentry_count] == "A11023": # This might be a mistake in GCAT, listed as AL IN but still in orbit?
                continue

            if (self.start_year <= int(jsr_data["DDate"][reentry_count][0:4]) <= self.final_year):
                if "*" in jsr_data["Parent"][reentry_count]:
                    print(f"Warning - asterix found in parent for {jsr_data['#JCAT'][reentry_count]}")
                jsr_parent = jsr_data["Parent"][reentry_count].replace("*","")
                cargo_found = False 
                for reentry in self.unique_reentry_list:
                    if reentry["jcat"] == jsr_parent:
                        reentry["other_mass"] += np.float64(jsr_data["DryMass"][reentry_count]) 
                        cargo_found = True    
                        if jsr_data["Piece"][reentry_count][5:6] in ["F","U","S"]:
                            sys.exit(f'Cargo in failed/suborbital object {jsr_data["Piece"][reentry_count]}.')
                    # These are where the parent object is also attached, so its not in the database.
                    # As these only happen very occasionally, just set it to the grandparent manually.
                    if reentry["jcat"] == "S51660" and jsr_data["#JCAT"][reentry_count] == "A09988":
                        reentry["other_mass"] += np.float64(jsr_data["DryMass"][reentry_count])
                        cargo_found = True
                    if reentry["jcat"] == "A06581" and jsr_data["#JCAT"][reentry_count] == "A09680":
                        reentry["other_mass"] += np.float64(jsr_data["DryMass"][reentry_count])
                        cargo_found = True 
                    if reentry["jcat"] == "S63628" and jsr_data["#JCAT"][reentry_count] == "A11468": # Cargo in Dragon CRS-32  
                        reentry["other_mass"] += np.float64(jsr_data["DryMass"][reentry_count])
                        cargo_found = True  
                if cargo_found == False:
                    print(f'Cargo missing for {jsr_data["#JCAT"][reentry_count]} - Parent {jsr_data["Parent"][reentry_count]}')                    
     
    def add_attached(self):
        
        for file in ["satcat","auxcat","lcat","rcat","lprcat","deepcat","ecat"]:
            attached_data = self.jsr_data_dict[file]
            attached_data_stripped = attached_data[attached_data["Status"].isin(["AR","AS","AL"])].reset_index(drop=True) # Attached objects only.
            attached_data_stripped = attached_data_stripped[~attached_data_stripped["Piece"].str[5:6].isin(["S","F","U"])].reset_index(drop=True) # Skip suborbital and failed launches (treated separately).
            attached_data_stripped = attached_data_stripped[
                (attached_data_stripped["DDate"].str[:4].astype(int) >= self.start_year) & 
                (attached_data_stripped["DDate"].str[:4].astype(int) <= self.final_year)
            ].reset_index(drop=True)
            attached_data_stripped = attached_data_stripped[attached_data_stripped["Primary"] == "Earth"].reset_index(drop=True)  
            
            # This should be listed as re-entering on the Moon, so ignore.
            attached_data_stripped = attached_data_stripped[attached_data_stripped["PLName"] != "Manfred Mem. Moon Mission"].reset_index(drop=True)  
            
            missing_list = []
            added = 0
            
            
            for reentry_count in range(len(attached_data_stripped)):
                
                # Skip objects from very recent launches.
                if (attached_data_stripped["Piece"][reentry_count][:4] == "2025" and 
                    int(attached_data_stripped["Piece"][reentry_count][5:8]) > max(self.dsl_id_list)):
                    continue

                # First see if the parent is already in the reentry list. Add the mass to the parent object if it is, or add the details to a list if missing.
                found_parent = False
                for reentry in self.unique_reentry_list:
                    if reentry["jcat"][:6] == attached_data_stripped["Parent"][reentry_count][:6]:
                        found_parent = True
                        added += 1
                        reentry["attached_abl_mass"] += np.float64(attached_data_stripped["DryMass"][reentry_count]) 

                # For A10426/7/8, the parent chain is A09939/42 > 38 > 37. So need to hard set if we don't want to search for great-grandparents.
                if attached_data_stripped['Parent'][reentry_count][:6] in ["A09942","A09939"]:
                    for reentry in self.unique_reentry_list:
                        if reentry["jcat"][:6] == "A09937":
                            reentry["attached_abl_mass"] += np.float64(attached_data_stripped["DryMass"][reentry_count]) 
                            added += 1
                            found_parent = True
                
                if found_parent == False and file not in ["lprcat.tsv"]:
                    missing_list.append([attached_data_stripped['Piece'][reentry_count],
                                            attached_data_stripped['Parent'][reentry_count][:6],
                                            attached_data_stripped["DryMass"][reentry_count],
                                            file,
                                            attached_data_stripped["#JCAT"][reentry_count]])
                            
            # Loopover all items where the parent couldn't be found (ie this means the parent is probably also attached.)
            for i in range(len(missing_list)):
                
                found_parent = False
                found_grandparent = False
                #print(f"Looking for parents of object {missing_list[i][0]} from {file}.")
                        
                # Look in all the files to find the parent object.
                for file_2 in ["satcat","auxcat","lcat","rcat","lprcat","deepcat","ecat"]:
                    file_data = self.jsr_data_dict[file_2]
                    file_data = file_data[file_data['#JCAT'] == missing_list[i][1]].reset_index(drop=True)
                    file_data = file_data[file_data['Status'].isin(["AR","AL"])].reset_index(drop=True)
                    if len(file_data) == 1:
                        found_parent = True
                        grandparent = file_data["Parent"][0] # Now we know what the grandparent is, its probably already in the reentry list.
                        break
                    elif len(file_data) > 1:
                        sys.exit("Multiple parents found.")
                if found_parent == True:
                    for reentry in self.unique_reentry_list:
                        if grandparent == reentry["jcat"][:6]:
                            #print(f'Found grandparent for {missing_list[i][0]} in {file}. {reentry["id"]}.')
                            reentry["attached_abl_mass"] += np.float64(missing_list[i][2])
                            added += 1
                            found_grandparent = True
                            break
                else:
                    print(f"Couldn't find parent {missing_list[i]}.")
                if found_grandparent == False:
                    print(f"Couldn't find grandparent {missing_list[i]}")
        
    def extract_aerospace_info(self, filepath):
        
        """Use the Aerospace Corp database to fill in missing time information and any other missing objects.
        """

        # First load the database and strip to only reentries in the year specified.        
        aerospace_corp_data = pd.read_csv(filepath)
        aerospace_corp_data_year_list = []
        for reentry_count in range(len(aerospace_corp_data)):
            if aerospace_corp_data["Aerospace Reentry Prediction (UTC)"][reentry_count][0:4] != "-":
                if self.start_year <= int(aerospace_corp_data["Aerospace Reentry Prediction (UTC)"][reentry_count][0:4]) <= self.final_year:
                    aerospace_corp_data_year_list.append(aerospace_corp_data.iloc[[reentry_count]])
        aerospace_corp_data_year = pd.concat(aerospace_corp_data_year_list, ignore_index=True)
        
        # Check what items are already in the list.
        jsr_id_list = []
        for reentry in (self.unique_reentry_list):
            jsr_id_list.append(reentry["id"])
        
        # Add time information from Aerospace Corp for any JSR entries missing time information.
        ac_time_update_mass, ac_time_update_count, ac_uncertainty_list = 0,0,[]
        for reentry in (self.unique_reentry_list):
            if reentry["time"] == -1:
                for reentry_count in range(len(aerospace_corp_data_year)): 
                    if aerospace_corp_data_year["International Designator"][reentry_count] == reentry["id"]:
                        ac_uncertainty_list.append(aerospace_corp_data_year["Aerospace Stated Uncertainty 20% Rule (+/- hrs)"][reentry_count])
                        reentry_time = aerospace_corp_data_year["Aerospace Reentry Prediction (UTC)"][reentry_count][11:16]
                        reentry["time"] = np.float64(int(reentry_time[0:2]) + int(reentry_time[3:5]) / 60)
                        ac_time_update_mass += (reentry["abl_mass"] +reentry["other_mass"])
                        ac_time_update_count += 1
                        
        print(f"Time set by Aerospace Corp: {ac_time_update_mass},{ac_time_update_count}")
        
        # Now add any new entries.
        for reentry_count in range(len(aerospace_corp_data_year)):  
            if "DEB" in aerospace_corp_data_year['Object Name'][reentry_count]:
                continue
            aerospace_corp_id = aerospace_corp_data_year['International Designator'][reentry_count].strip()

            if aerospace_corp_id[:4] == "2025" and int(aerospace_corp_id[5:8]) > max(self.dsl_id_list):
                continue

            if aerospace_corp_id in ["1010-012BF",  # Typo, should be 2020-012BF.
                                     "2021-041K",   # Already in GCAT with different COSPAR ID.
                                     "2019-024AG", # Already in GCAT with different COSPAR ID.
                                     "2020-001Z",  # Already in GCAT with different COSPAR ID.
                                     "2021-074Z",  # Already in GCAT with different COSPAR ID.
                                     "2021-073G",  # Already in GCAT with different COSPAR ID.
                                     "2021-073BG", # Already in GCAT with different COSPAR ID.
                                     #"2019-029Q",  # Starlink payload. DISCOSweb and JSR list this as reentering in July 2022. AC lists this as reentering in 2020. Ignoring.
                                     "2020-012B",  # Starlink payload. DISCOSweb has no re-entry data. AC lists this as reentering Mar 2020. JSR lists this as still in orbit. Ignoring.
                                     "2023-097A",  # Already in GCAT with different COSPAR ID.
                                     "1982-092A",  # Debris of Kosmos 1408.
                                     "2024-076Y",  # Still orbiting according to satellite tracker websites.
                                     "2019-018AC", # Orbital platform that re-entered attached to the 4th stage.
                                     "2023-046AZ", # Already in GCAT with different COSPAR ID.
                                     "2023-046AV", # Already in GCAT with different COSPAR ID.
                                     "2023-046BF", # Already in GCAT with different COSPAR ID.
                                     "2019-069B",  # GCAT has this reentering in 2022.
                                     "2020-074AM", # GCAT has this reentering in 2022.
                                     "1984-108B",  # GCAT has this reentering in 2022.
                                     "2025-082B",  # Already in GCAT with different COSPAR ID.
                                     "2024-253C",  # POEM 4, has definitley reentered but JSR lists the stage its attached to as still orbiting.
                                     "2025-037L",  # Still orbiting according to satellite tracker websites. 
                                     "2024-181B",  # Not yet reentered according to DW and GCAT.
                                     "2022-025AC", # Reentered in late Dec 2024.
            ]:
                continue

            if aerospace_corp_id not in jsr_id_list:
                print(f"Adding reentry from Aerospace Corp - {aerospace_corp_id}, {aerospace_corp_data_year['Object Name'][reentry_count]}") 
                # 1992-021C, Ariane 44L H10+ 3rd Rocket Stage. DISCOSweb and AC list this as reentering in Oct 2020. JSR lists this as exploding in Apr 1993.
                # 2020-057X, GCAT is wrong, should be 2025, not 2024.
                
                if aerospace_corp_data_year["International Designator"][reentry_count] == "1992-021C":
                    reentry_category = "S3"
                    abl_mass = 2080
                    inc = 3.98
                elif aerospace_corp_data_year["International Designator"][reentry_count] == "2020-057X":
                    reentry_category = "P"
                    abl_mass = 248
                    inc = 53.00
                else:
                    print("Added object not expected.")
                # Extract the UTC time.     
                reentry_time = aerospace_corp_data_year["Aerospace Reentry Prediction (UTC)"][reentry_count][11:16]
                time_utc = np.float64(int(reentry_time[0:2]) + int(reentry_time[3:5]) / 60)
                temp_reentry_dict = {
                    "id"               : aerospace_corp_data_year["International Designator"][reentry_count],
                    "jcat"             : "N/A",
                    "name"             : aerospace_corp_data_year["Object Name"][reentry_count],
                    "category"         : reentry_category,
                    "burnup"           : "Complete",
                    "time"             : time_utc,
                    "datestr"          : aerospace_corp_data_year["Aerospace Reentry Prediction (UTC)"][reentry_count][:10].replace("-",""),
                    "lat"              : round(random.uniform(-inc,inc),2),
                    "lon"              : round(random.uniform(-180, 180),2),
                    "abl_mass"         : abl_mass,
                    "other_mass"       : 0,
                    "attached_abl_mass": 0,
                    "location"         : 6,
                    "apogee"           : 100,
                }  
                
                self.unique_reentry_list.append(temp_reentry_dict)
    
    def extract_discosweb_info(self):
        
        """Look for and add any missing objects from DISCOSweb (contained within a netcdf file built in pull_from_discosweb.py).
        """        
        
        # Check what items are already in the list.
        jsr_ac_id_list = []
        for reentry in (self.unique_reentry_list):
            if reentry["id"] != "2020-029B":
                jsr_ac_id_list.append(reentry["id"]) 
        
        for reentry_count in range(len(self.ds_dw["DISCOSweb_Reentry_ID"].values)):

            discosweb_name = self.ds_dw["DISCOSweb_Reentry_Name"].values[reentry_count].item()
            discosweb_id   = self.ds_dw["DISCOSweb_Reentry_ID"].values[reentry_count].item()
            
            # Excluding Debris and some duplicated Starlinks with different IDs.
            if ("debris" in discosweb_name.lower() or "2024-129I" in discosweb_id or "2022-010I" in discosweb_id):
                continue
            if (discosweb_id not in jsr_ac_id_list):
                
                # Check if the stage has already been added.  
                found_stage = False               
                if self.ds_dw["DISCOSweb_Reentry_Class"].values[reentry_count] == "Rocket Body":
                    
                    if (discosweb_name in ["H-II LE-5B (H-IIB)", "L-53 (YF24B) (Long March (CZ) 2D)",
                                           "CZ-5-HO (Long March (CZ) 5)","Delta IV DCSS 5 (Delta 4H)",
                                           "LE-5B-3 (H-III 22)"]
                        or "Falcon 9 Merlin-V (1D" in discosweb_name or "Centaur-5" in discosweb_name
                        or "second" in discosweb_name.lower()):
                            reentry_category = "S2"
                    elif (("Fregat" in discosweb_name)
                        or (discosweb_name in ["PBV (Long March (CZ) 6)","YZ-1S (Long March (CZ) 2C/YZ-1S)",
                                               "H-18 (Long March (CZ) YF) (Long March (CZ) 3C/E)","Angara AM (Angara 1.2)",
                                               "Lunar Photon (Electron (Curie))","Stage-3 (Tianlong 2)","PS3 (PSLV-CA)",
                                               "KZ-11 Stage 3 (Kuaizhou-11)","YZ-3 (Long March (CZ) 2D/YZ-3)",
                                               "Blok-DM-3 (Angara A5 Orion)","YZ-2 (Long March (CZ) 5/YZ-2)"])
                        or "third" in discosweb_name.lower()):
                            reentry_category = "S3" 
                    elif (discosweb_name in ["AVUM (Vega)","Ceres-1 upperstage","YZ-1 (Long March (CZ) 3B/YZ-1)",
                                            "KZ-1 Stage 4 (Kuaizhou-1)","Long March (CZ) 11 Stage 4","AVUM+ (Vega C)",
                                            ""]
                        or "fourth" in discosweb_name.lower()):
                            reentry_category = "S4"
                    else:
                        print("Missing stage designation for : ",discosweb_id, discosweb_name)
                        reentry_category = "S0"
                    
                    for reentry in self.unique_reentry_list:
                        if reentry["id"][:8] == discosweb_id[:8] and str(reentry["category"]) == str(reentry_category):
                            found_stage = True
                            
                skipped_ids = ["2014-065A",  # Chang'e 5-T1, reentered in 2014.
                               "2021-011A",  # Progress MS-16 is listed in JSR as AR, so the mass is just added to the parent object and no object exists.
                               # Progress MS-16 has been added though, so we can ignore it.
                               "2020-001N",  # Starlink, still in orbit??
                               "2020-001T",  # Starlink, still in orbit??
                               "2020-001G",  # Starlink, still in orbit??
                               "2017-069C",  # Listed in JSR as still in orbit.
                               "2020-001AQ", # Starlink, still in orbit??
                               "2020-001BG", # Starlink, still in orbit??
                               "2020-001BF", # Starlink, still in orbit??
                               "2020-001BM", # Starlink, still in orbit??
                               "2020-001A",  # Starlink, still in orbit??
                               "2020-001AU", # Starlink, still in orbit??
                               "2015-007B",  # Debris.
                               "2021-110A", "2014-065B", "2023-098A", "2023-118A", "2022-168A", # Reentered or landed on another planet or moon.
                               "2024-030A", "2023-137D", "2025-010A", "2025-038A", # Reentered or landed on another planet or moon.
                               "2022-094B",  # Falcon Stage 2 # JSR has this has as deep space.
                               "1998-067NF", # Already in GCAT, tracable through GCAT.
                               "1998-067PU", # Already in GCAT, tracable through GCAT.
                               "1982-092A",  # Debris of Cosmos
                               "1967-048H", "2023-137B", "2023-137G", # GCAT lists as debris.
                               "2022-093J",  # Attached landed.
                               "2022-175M",  # Incorrect on DW, should be 2025. Confirmed using satellite tracker.
                               "2023-046AV", "2023-046AZ","2023-046BF", # Different COSPAR IDs in GCAT.
                               "2023-097A",  # Attached when reentering.
                               "2023-035A",  # Still in orbit, DW has probably just set re-entry of 2nd and 3rd stages together.
                               "2023-154C",  # Sent into deep space.
                               "2024-069IB", # Listed as in orbit in GCAT. Couldn't verify using satellite tracker.
                               "2024-110IA", # Listed as in orbit in GCAT. Couldn't verify using satellite tracker.
                               "2023-137H",  # Fairing Debris.
                               "2023-054AY", # Reentered in 2024.
                               "2023-193C",  # Still orbiting according to satellite tracker websites.
                               "2019-038K",  # Reentered in 2025, confirmed using satellite tracker.
                               "2020-029B",  # GCAT has DDate as "2023 Jan?", but DISCOSweb has 31st Dec 2022.
                               "2020-057X",  # GCAT / DW have this as 2024 but should be 2025.
                ]
                
                # Included IDs: 
                # 2019-069B Wrong in GCAT, should be 2023 not 2022. Confirmed using satellite tracker.
                # 2023-109A https://www.n2yo.com/satellite/?s=57481 Reentered in Jan 2025 
                # 2015-049A https://www.n2yo.com/satellite/?s=40899 Reentered in Apr 2025, not 2024 as in GCAT.

                if found_stage == False and self.ds_dw['DISCOSweb_Reentry_ID'].values[reentry_count] not in skipped_ids:
                    print(f"Adding object with COSPAR ID {self.ds_dw['DISCOSweb_Reentry_ID'].values[reentry_count]}")          
                
                    if self.ds_dw["DISCOSweb_Reentry_Class"].values[reentry_count] in ["Rocket Mission Related Object","Payload Mission Related Object"]:
                        reentry_category = "C"
                    elif self.ds_dw["DISCOSweb_Reentry_Class"].values[reentry_count] == "Payload":
                        reentry_category = "P" 
                    else:
                        print("Unexpected re-entry stage.") 
                    datestr = int(self.ds_dw["DISCOSweb_Reentry_Epoch"].values[reentry_count].item()[:10].replace("-",""))
                    time_utc = -1
                    
                    if discosweb_id == "2019-069B":
                        inc = 87.89
                    elif discosweb_id == "2023-109A":
                        inc = 5.00
                    elif discosweb_id == "2015-049A":
                        inc = 97.42
                    else: 
                        inc = 0
                        print(f"Added object not expected ({self.ds_dw['DISCOSweb_Reentry_ID'].values[reentry_count]}).")
                    lat = round(random.uniform(-inc, inc),2)
                    lon = round(random.uniform(-180, 180),2)
                    mass = self.ds_dw["DISCOSweb_Reentry_Mass"].values[reentry_count]
                    if np.isnan(mass):
                        mass = 0
    
                    temp_reentry_dict = {
                        "id"               : discosweb_id,
                        "jcat"             : "N/A",
                        "name"             : discosweb_name,
                        "category"         : reentry_category,
                        "burnup"           : "Complete",
                        "time"             : time_utc,
                        "datestr"          : datestr,
                        "lat"              : lat,
                        "lon"              : lon,
                        "abl_mass"         : mass,
                        "other_mass"       : 0, 
                        "attached_abl_mass": 0,
                        "location"         : 6, 
                        "apogee"           : 100,
                    }

                    self.unique_reentry_list.append(temp_reentry_dict)  
            
            else:
                # If the item already exists, check if the jsr entry is missing the mass information.
                # Only add if the jsr is missing the mass, and discosweb has the mass.
                # This is currently zero items, but leave this code in in case it is reused for future years.
                dict_ind = next((index for (index, d) in enumerate(self.unique_reentry_list) if d["id"] == self.ds_dw["DISCOSweb_Reentry_ID"].values[reentry_count]), False)
                if dict_ind != False:
                    if (self.unique_reentry_list[dict_ind]["category"] in ["C", "P"]) and ((self.unique_reentry_list[dict_ind]["abl_mass"] +self.unique_reentry_list[dict_ind]["other_mass"]) == 0):
                        if not np.isnan(self.ds_dw["DISCOSweb_Reentry_Mass"].values[reentry_count]):
                            self.unique_reentry_list[dict_ind]["payload_mass"] = self.ds_dw["DISCOSweb_Reentry_Mass"].values[reentry_count]
                            print("Added info for",self.ds_dw["DISCOSweb_Reentry_ID"].values[reentry_count])                       
         
    def add_missing_stages(self):
        
        """Add any missing stages for launches in the year.
        """        
                    
        stage_alt_dict = {}
        stage_alt_rockets = np.genfromtxt("./input_files/launch_event_altitudes.csv",dtype=str,skip_header=1,usecols=[0],delimiter=",")
        stage_alt_data = np.genfromtxt("./input_files/launch_event_altitudes.csv",dtype=np.float64,skip_header=1,usecols=[1,2,3,4],delimiter=",")
        for i in range(len(stage_alt_data)):
           
            if stage_alt_data[i][0] == '':
                stage_alt_dict[f"{stage_alt_rockets[i]} BECO"] = None
            else:
                stage_alt_dict[f"{stage_alt_rockets[i]} BECO"] = np.float64(stage_alt_data[i][0])
            if stage_alt_data[i][1] == '':
                stage_alt_dict[f"{stage_alt_rockets[i]} MECO"] = None
            else:
                stage_alt_dict[f"{stage_alt_rockets[i]} MECO"] = np.float64(stage_alt_data[i][1])
            if stage_alt_data[i][2] == '':
                stage_alt_dict[f"{stage_alt_rockets[i]} SEI1"] = None
            else:
                stage_alt_dict[f"{stage_alt_rockets[i]} SEI1"] =  np.float64(stage_alt_data[i][2])
        
        missing_boosters_count, missing_first_count = 0,0
        missing_boosters_mass, missing_first_mass = 0,0

        # Loop over each launch, locate the rocket and then filter for rocket configuration.
        for i in range(len(self.dsl["COSPAR_ID"])):
            for count, rocket_name in enumerate(self.dsr["Rocket_Name"].values):
                if rocket_name == self.dsl["Rocket_Name"].values[i] and self.dsl["COSPAR_ID"].values[i][5] != "F":

                    ################
                    # Boosters 
                    ################
                    
                    if int(self.dsr['Booster_No'].values[count]) > 0:
                        # Count how many boosters are currently added, and compare it to the number of boosters there should be.
                        booster_count = 0
                        for reentry in self.unique_reentry_list:
                            if reentry["id"][:8] == self.dsl["COSPAR_ID"].values[i] and str(reentry["category"])[0] == "B":
                                booster_count += 1
                        if booster_count != int(self.dsr['Booster_No'].values[count]):
                            #print(f"There are {booster_count} boosters for {self.dsl['COSPAR_ID'].values[i]} when there should be {self.dsr['Booster_No'].values[count]} boosters.")
                            # Adding boosters where BECO is above 50km. If average then its for B+1/2S and B+3S only.
                            beco = stage_alt_dict[f"{rocket_name} BECO"]
                            if (beco > 50) or (np.isnan(beco) and (self.dsr["Stage4_PropMass"].values[count] == 0)):
                                
                                # Set beco so we can add to the netcdf output.
                                if np.isnan(beco):
                                    if self.dsr["Stage3_PropMass"].values[count] == 0:
                                        beco = 66
                                    else:
                                        beco = 55
                                
                                # Now add the boosters.        
                                print(f"Adding {int(self.dsr['Booster_No'].values[count])-booster_count} boosters for {self.dsl['COSPAR_ID'].values[i]}")
                                for j in range(int(self.dsr["Booster_No"].values[count])-booster_count):
                                    abl_mass = float(self.dsr["Stage0_StageMass"].values[count]) / int(self.dsr["Booster_No"].values[count])
                                    if np.isnan(abl_mass):
                                        abl_mass = 0
                                    missing_boosters_count += 1
                                    missing_boosters_mass += abl_mass   
                                    temp_reentry_dict = {
                                        "id"               : f'{self.dsl["COSPAR_ID"].values[i]}',
                                        "jcat"             : "N/A",
                                        "name"             : f"{self.dsl['Rocket_Name'].values[i]} Booster",
                                        "burnup"           : "Complete",
                                        "category"         : "B" + str(j+1),
                                        "time"             : self.dsl["Time(UTC)"].values[i],
                                        "datestr"          : self.dsl["Date"].values[i],
                                        "lat"              : self.dsl["Latitude"].values[i],
                                        "lon"              : self.dsl["Longitude"].values[i],
                                        "abl_mass"         : abl_mass,
                                        "other_mass"       : 0,
                                        "attached_abl_mass": 0,
                                        "smc"              : False,
                                        "location"         : 2,
                                        "apogee"           : beco,
                                    }
                    
                                    self.unique_reentry_list.append(temp_reentry_dict) 
        
                    ################
                    # 1st Stage 
                    ################

                    meco = stage_alt_dict[f"{rocket_name} MECO"]
                    
                    add_first_stage = False
                    if (50 < meco <= 100):                 
                        add_first_stage = True
                    elif np.isnan(meco):
                        # 2S,3S,4S
                        if self.dsr["Stage0_PropMass"].values[count] == 0: 
                            add_first_stage = True
                            if ((self.dsr["Stage3_PropMass"].values[count] == 0) and (self.dsr["Stage4_PropMass"].values[count] == 0)):
                                meco = 90 #2S
                            elif ((self.dsr["Stage3_PropMass"].values[count] != 0) and (self.dsr["Stage4_PropMass"].values[count] == 0)):
                                meco = 56 #3S
                            elif ((self.dsr["Stage3_PropMass"].values[count] != 0) and (self.dsr["Stage4_PropMass"].values[count] != 0)):
                                meco = 52 #4S
                        # B+4S
                        if ((self.dsr["Stage0_PropMass"].values[count] != 0) and (self.dsr["Stage4_PropMass"].values[count] != 0)):    
                            add_first_stage = True
                            meco = 64
                    
                    if add_first_stage == True:    
                        # Check if there is a first stage, and add if not.
                        found_stage = False
                        for reentry in self.unique_reentry_list:
                            if reentry["id"][:8] == self.dsl["COSPAR_ID"].values[i] and reentry["category"] == "S1":
                                found_stage = True
                        if found_stage == False:
                            print(f"Adding 1st stage for {self.dsl['COSPAR_ID'].values[i]}")
                            abl_mass = self.dsr["Stage1_StageMass"].values[count]
                            if np.isnan(abl_mass):
                                abl_mass = 0
                            missing_first_count += 1
                            missing_first_mass += abl_mass
                            temp_reentry_dict = {
                                "id"               : f'{self.dsl["COSPAR_ID"].values[i]}S1',
                                "jcat"             : "N/A",
                                "name"             : f"{self.dsl['Rocket_Name'].values[i]} Stage 1",
                                "burnup"           : "Complete",
                                "category"         : "S1",
                                "time"             : self.dsl["Time(UTC)"].values[i],
                                "datestr"          : self.dsl["Date"].values[i],
                                "lat"              : self.dsl["Latitude"].values[i],
                                "lon"              : self.dsl["Longitude"].values[i],
                                "abl_mass"         : abl_mass,
                                "other_mass"       : 0,  
                                "attached_abl_mass": 0,                           
                                "smc"              : False,
                                "location"         : 2,
                                "apogee"           : meco,
                            }
        
                            self.unique_reentry_list.append(temp_reentry_dict)

                    # Check the number of stages.
                    stage_count = np.zeros((4))
                    for reentry in self.unique_reentry_list:
                        for j in range(1,5):
                            if reentry["id"][:8] == self.dsl["COSPAR_ID"].values[i] and reentry["category"] == f"S{j}":
                                stage_count[j-1] += 1
                    #print(stage_count)
                    for stage in stage_count:
                        if stage > 1:
                            print(f"Multiple stages for {self.dsl['COSPAR_ID'].values[i]}, {stage_count}")
        
        print(f"Missing Boosters:    {missing_boosters_mass},{missing_boosters_count}")
        print(f"Missing First Stage: {missing_first_mass},{missing_first_count}")
    
    def import_raul_spacex_map(self):
        """Import and set up the geolocation lists for falcon landings.
        """        
        
        # Import the data.
        fiona.drvsupport.supported_drivers['KML'] = 'rw'
        raul_data = gpd.read_file('./databases/reentry/General_SpaceX_Map_Raul.kml', driver='KML', layer =2)
        
        # Build the ocean landing list.
        ocean_landings = raul_data[(raul_data["Name"].str.contains("ASDS")) & (raul_data["Description"].str.contains("Landing -"))].reset_index(drop=True)
        ocean_landings = ocean_landings[~ocean_landings["Name"].str.contains("planned")].reset_index(drop=True)
        ocean_landings = ocean_landings[~ocean_landings["Name"].str.contains("planed")].reset_index(drop=True)
        date_col = []
        for landing in range(len(ocean_landings["Name"])):
            year_ind = ocean_landings["Description"][landing].index("Landing -")
            datestr = ocean_landings["Description"][landing][year_ind:].replace(" ","")[8:17].replace("-","")
            date_col.append(datestr)
        ocean_landings["Date"] = date_col
        ocean_landings['Date'] = pd.to_datetime(ocean_landings['Date'],format='%d%b%y')
        self.ocean_landings = ocean_landings[(ocean_landings["Date"] < f'{self.final_year+1}-01-01') & (ocean_landings["Date"] >= f'{self.start_year}-01-01')].reset_index(drop=True)
            
        # Build the ground landing list.
        ground_landings = raul_data[(raul_data["Name"].str.contains("landing")) & (raul_data["Name"].str.contains("LZ"))].reset_index(drop=True)
        ground_landings = ground_landings[~ground_landings["Name"].str.contains("planned")].reset_index(drop=True)
        date_col = []
        for landing in range(len(ground_landings["Name"])):
            year_ind = ground_landings["Description"][landing].index("Landing -")
            datestr = ground_landings["Description"][landing][year_ind:].replace(" ","")[8:17].replace("-","")
            date_col.append(datestr)
        ground_landings["Date"] = date_col
        ground_landings['Date'] = pd.to_datetime(ground_landings['Date'],format='%d%b%y')
        self.ground_landings = ground_landings[(ground_landings["Date"] < f'{self.final_year+1}-01-01') & (ground_landings["Date"] >= f'{self.start_year}-01-01')].reset_index(drop=True)
        
        # Build the fairing recovery list.
        ocean_missions, ocean_dates = self.ocean_landings["Name"].tolist(), self.ocean_landings["Date"].tolist()
        ground_missions, ground_dates = self.ground_landings["Name"].tolist(), self.ground_landings["Date"].tolist()

        if len(ocean_missions) + len(ground_missions) > 0:
            dataframe_size = 0
            for count,mission in enumerate(ocean_missions):
                mission_string = mission.replace("ASDS position","")
                extract_df = raul_data[(raul_data["Name"].str.contains(mission_string, regex=False)) & (raul_data["Name"].str.contains("fairing", regex=False))].reset_index(drop=True)
                extract_df["Date"] = ocean_dates[count]
                if dataframe_size == 0:
                    fairings = extract_df
                    dataframe_size += 1
                else: 
                    fairings = pd.concat([fairings,extract_df], ignore_index=True)
            for count,mission in enumerate(ground_missions):
                mission_string = mission.replace("ASDS position","")
                extract_df = raul_data[(raul_data["Name"].str.contains(mission_string, regex=False)) & (raul_data["Name"].str.contains("fairing", regex=False))].reset_index(drop=True)
                extract_df["Date"] = ground_dates[count]
                if dataframe_size == 0:
                    fairings = extract_df
                    dataframe_size += 1
                else: 
                    fairings = pd.concat([fairings,extract_df], ignore_index=True)
                self.fairings = fairings
        else:
            self.fairings = []
                         
    def get_reentry_info(self):

        """This is the main function of this class, and loops over all the key databases (GCAT, AC, DW).
        It also does small final adjustments:
            - checks for duplicates
            - sets smc info
            - sets rocket stage/fairing mass from database
            - fixes fairing geolocation
            - fixes timings where launch and reentry are on same day
            - sets any missing masses
            - fixes any geolocations on grid border (180E) 
        """        
        
        # Import the required files.
        #self.ds_dw = xr.open_dataset(f"./databases/reentry/DISCOSweb/discosweb_reentries_{self.start_year}-{self.final_year}.nc", decode_times=False)
        #class_mask = self.ds_dw["DISCOSweb_Reentry_Class"].isin(["Payload","Rocket Body","Rocket Mission Related Object","Payload Mission Related Object"])
        #self.ds_dw = self.ds_dw.where(class_mask, drop=True) # Filter for only the relevant classes.

        self.dsr = xr.open_dataset(f"./databases/rocket_attributes_{self.start_year}-{self.final_year}_jsr.nc", decode_times=False)
        self.dsl = xr.open_dataset(f"./databases/launch_activity_data_{self.start_year}-{self.final_year}_jsr.nc", decode_times=False) 
        self.import_raul_spacex_map() #(https://t.co/RAsQ9NDmEr)
             
        self.unique_reentry_list = []
        
        print("Loading JSR re-entries.") # (see https://planet4589.org/space/gcat/web/cat/cats.html)
        files = ["satcat","auxcat","lcat","rcat","lprcat","deepcat","ecat","ftocat"]
        # satcat (main database), auxcat (should be in main but isn't), lcat (suborbital stages/objects). 
        # rcat (lower stages and fairings), lprcat (several objects returning from space).
        # deepcat (several objects returning from space), ecat (capsules from ISS / crewed missions). 
        # ftocat (failed launches).
        self.jsr_data_dict = {}
        for file in files:
            url = "https://planet4589.org/space/gcat/tsv/cat/" + file + ".tsv"
            response = self.session.get(url)
            if response.status_code == 200:
                # Convert the content to a file-like object for pandas
                tsv_data = StringIO(response.text)
                # Load the data into a pandas DataFrame
                df = pd.read_csv(tsv_data, delimiter="\t", dtype=object)
                df = df[df["Piece"].notna()]
                df = df[df["Piece"] != "UNK"]
                df = df[~df["Type"].str[0].isin(["Z","D"])] 
                df["Converted_Tag"] = df["Piece"].apply(self.convert_launch_tag)
                self.jsr_data_dict[file] = df
            else:
                raise ImportError(f"Failed to fetch {file} from JSR", response.status_code)
            
        # Add each file sequentially. There is nothing relevant for this inventory in hcocat, tmpcat and csocat.
        for file in ["satcat","auxcat","lcat","rcat","lprcat","deepcat","ecat"]:
            print(file)
            self.extract_jsr_info(self.jsr_data_dict[file])
        print("ecat cargo")
        self.add_cargo(self.jsr_data_dict["ecat"])          
        print("ftocat")
        self.extract_jsr_info(self.jsr_data_dict["ftocat"])  
        print("attached")
        self.add_attached()

        # Check for duplicate reentries with the same jcat id.         
        self.jsr_jcat_list = []
        for reentry in self.unique_reentry_list:
            self.jsr_jcat_list.append(reentry["jcat"])
        duplicate_jcat_list = []
        for jsr_jcat in self.jsr_jcat_list:
            if self.jsr_jcat_list.count(jsr_jcat) > 1:
                duplicate_jcat_list.append(jsr_jcat)
        duplicate_jcat_list = sorted(list(set(duplicate_jcat_list)))
        if len(duplicate_jcat_list) > 0:
            print(f"Duplicates: {duplicate_jcat_list}")

        # Check for boosters and rename the category to B1-B4.
        booster_id_list = []
        for reentry in self.unique_reentry_list:
            if reentry["category"] == "S0":
                booster_id_list.append(reentry["id"])
        booster_id_list = sorted(list(set(booster_id_list)))
        for booster_id in booster_id_list:
            count = 1
            for reentry in self.unique_reentry_list:
                if reentry["id"] == booster_id and reentry["category"] == "S0":
                    reentry["category"] = "B" + str(count)
                    count += 1
        
        #Add missing stages and check the Aerospace Corp and DISCOSweb databases.
        print("Looking for missing rocket stages.")
        self.add_missing_stages()
        print("Searching for Aerospace Corp re-entries.")
        self.extract_aerospace_info("./databases/reentry/AerospaceCorp/AerospaceCorp_Reentries_07-07-25.csv")
        print("Searching for DISCOSweb re-entries.")
        self.extract_discosweb_info()
        
        #######################
        ## Final adjustments.
        #######################
        
        # Create a list of all smc-related launches (from launch database). 
        smc_dict = {}
        for i in range(len(self.dsl["COSPAR_ID"])): 
            smc_dict[self.dsl["COSPAR_ID"].values[i]] = self.dsl["Megaconstellation_Flag"].values[i]
            
        if self.start_year == 2023 or self.start_year == 2025:
            self.dsl2 = xr.open_dataset(f"./databases/launch_activity_data_2020-2022.nc", decode_times=False)
            for i in range(len(self.dsl2["COSPAR_ID"])): 
                smc_dict[self.dsl2["COSPAR_ID"].values[i]] = self.dsl2["Megaconstellation_Flag"].values[i]
            self.dsl2.close()
        if self.start_year == 2025:
            self.dsl3 = xr.open_dataset(f"./databases/launch_activity_data_2023-2024.nc", decode_times=False)
            for i in range(len(self.dsl3["COSPAR_ID"])): 
                smc_dict[self.dsl3["COSPAR_ID"].values[i]] = self.dsl3["Megaconstellation_Flag"].values[i]
            self.dsl3.close()

        for reentry in self.unique_reentry_list:
            try:
                reentry["smc"] = smc_dict[reentry["id"][:8]]
                if reentry["category"] == "P":
                    reentry["smc"] = False
                    if any(x.lower() in reentry["name"].lower() for x in ["Starlink", "OneWeb", "Yinhe", "Lingxi", "E-Space", "Lynk", "Kuiper"]):
                        reentry["smc"] = True
            except:
                if reentry["id"][:8] in ["2018-020","2019-029","2019-074"]:
                    reentry["smc"] = True
                else:
                    reentry["smc"] = False
                    
        # Alumina emissions are calculated using ablation and alumina content data from literature.
        # These values vary based on object class(core stage / upper stage / payload) and nature of launch (reusuable or not). 
        # https://www.sciencedirect.com/science/article/pii/B0122274105008887 "Fairings are typically made of aluminum or composite materials."
        # Therefore we treat fairings as a core stage, except for Falcon 9 which are recovered intact.
                    
        for reentry in self.unique_reentry_list:
            if ("fairing" in reentry["name"].lower()) or (reentry["category"][0] in ["B","S"]):
                reentry["alu_per"] = 0.7 # https://dspace.mit.edu/handle/1721.1/151443 (70% Al)
            elif reentry["category"] in ["C","P"]:
                reentry["alu_per"] = 0.4 # https://doi.org/10.1016/j.asr.2020.10.036 (40% Al)
            else:
                print(f"Couldn't assign aluminium mass information for {reentry['id']} with category {reentry['category']}")
            
            if reentry["burnup"] == "Complete":
                if ("fairing" in reentry["name"].lower()) or (reentry["category"][0] == "B") or (reentry["category"] == "S1"):
                    reentry["abl_deg"] = 0.3      # https://doi.org/10.1016/j.asr.2020.10.036 (70% survivability)
                elif reentry["category"] in ["C","P"]:
                    if reentry["smc"] == True:
                        reentry["abl_deg"] = 1.0  # https://doi.org/10.1016/j.asr.2020.10.036 (0% survivability)
                    elif reentry["smc"] == False:
                        reentry["abl_deg"] = 0.8  # https://doi.org/10.1016/j.asr.2020.10.036 (20% survivability)
                    else:
                        print(f"Couldn't assign smc information for {reentry['id']}")
                elif reentry["category"] in ["S2","S3","S4"]:
                    reentry["abl_deg"] = 0.65     # https://doi.org/10.1016/j.asr.2020.10.036 (35% survivability)
                else:
                    reentry["abl_deg"] = 0
                    print(f"Couldn't assign ablation information for complete burnup of {reentry['id']} with category {reentry['category']}")
            elif reentry["burnup"] == "Partial":
                reentry["abl_deg"] = 0
            else:
                print(f"Couldn't understand burnup information for {reentry['id']}")
        
        fairing_count_list = []
        for i in range(len(self.dsl["COSPAR_ID"])): 
             
            # Loop over all the fairings to check the number (should be two as it splits in half).
            # Also look for any difference in the inclination (should be zero).
            # Then set the geolocation to the location of the first fairing.  
            fairing_count = 0  
            fairing_inc = np.zeros((2))
            for count, reentry in enumerate(self.unique_reentry_list):  
                if (self.dsl["COSPAR_ID"].values[i][:8] == reentry["id"][:8] 
                and self.dsl["COSPAR_ID"].values[i][5] != "F"
                and "fairing" in reentry["name"].lower()):
                    fairing_count += 1
                    fairing_inc[fairing_count-1] = reentry["inc"]
                    if fairing_count == 1:
                        lat = reentry["lat"]
                        lon = reentry["lon"]
                    elif fairing_count == 2 and "Falcon 9" not in reentry["name"]:
                        reentry["lat"] = lat
                        reentry["lon"] = lon
            if fairing_inc[1]-fairing_inc[0] > 0:
                print(f"Warning: fairing inclinations differ for {self.dsl['COSPAR_ID'].values[i]} by {fairing_inc[1]-fairing_inc[0]}")      
            if fairing_count not in [0,2]:
                # NOTE: Apogee adjusted for 2022-150 fairing. One was 200km, one was 0km. Adjusted to make both 200km.
                sys.exit(f"Incorrect number of fairings ({fairing_count}) found for ID: {self.dsl['COSPAR_ID'].values[i]}")
            if fairing_count == 0 and self.dsl["COSPAR_ID"].values[i][5] != "F":
                pass # NOTE: If we want to manually add fairings, need to enable this to see which are missing.
                #print(f"No fairings found for {self.dsl['COSPAR_ID'].values[i]}")
            fairing_count_list.append(fairing_count)
            
        #################################################
        ## Update mass info using rocket_info databases. 
        ################################################# 
                    
        for reentry in self.unique_reentry_list:                
            for i in range(len(self.dsl["COSPAR_ID"])):
                if self.dsl["COSPAR_ID"].values[i][:8] == reentry["id"][:8]:
                    
                    # Update mass info for all rocket stages.
                    if reentry["category"] in ["B1","B2","B3","B4","B5","B6","S1","S2","S3","S4"] and self.dsl["COSPAR_ID"].values[i][5] != "F":
                        for count, rocket_name in enumerate(self.dsr["Rocket_Name"].values):
                            if rocket_name == self.dsl["Rocket_Name"].values[i]:
                                if reentry["category"] in ["B1","B2","B3","B4","B5","B6"]:  
                                    reentry["abl_mass"] = self.dsr["Stage0_StageMass"].values[count] / int(self.dsr["Booster_No"].values[count])
                                elif reentry["category"] in ["S1","S2","S3","S4"]:
                                    reentry["abl_mass"] = self.dsr[f"Stage{reentry['category'][1]}_StageMass"].values[count]
                    
                    # Update mass info for all fairings.
                    elif "fairing" in reentry["name"].lower() and self.dsl["COSPAR_ID"].values[i][5] != "F":
                        for count, rocket_name in enumerate(self.dsr["Rocket_Name"].values):
                            if rocket_name == self.dsl["Rocket_Name"].values[i]:
                                reentry["abl_mass"] = self.dsr["Fairing_Mass"].values[count] / fairing_count_list[i]
            
        time_update_mass_1, time_update_count_1 = 0,0     
        time_update_mass_2, time_update_count_2 = 0,0  
        missing_time_count = 0 
        for reentry in self.unique_reentry_list:                 
            # For West Ford dipoles, these are part of the West Ford Needles project. # https://space.skyrocket.de/doc_sdat/westford.htm
            # The needles each weigh 40 ng, and one 'clump' reentered in 2020. The needles are Copper, and so will not contribute to Al emissions.  
            if "sep motor cover" in reentry["name"]:
                # The same item is also listed as weighing 1 kg elsewhere.
                reentry["abl_mass"] = 1
            if reentry["id"] == "2020-086B": # Dest listed as 180E, this messes up other script so adjust to 179.9, will be in same grid square.
                reentry["lon"] = 179.9
            if reentry["name"] == "DLA-U": # Object 2013-009J (PSLV upper Dual Launch Adapter (DLA-U)) has this mass on DW.
                reentry["abl_mass"] = 100

            # Set time to midnight whenever the reentry is occurs on a different day than the launch and no time info is available.
            # Also set to midnight if the reentry is from a launch in a previous year.
            # When the reentry is on the same day, set it to the launch time.
            
            if reentry["time"] == -1:  
                missing_time_count +=1                
                for count, launch_id in enumerate(self.dsl["COSPAR_ID"].values):
                    if reentry["id"][:8] == launch_id:
                        if reentry["datestr"] == self.dsl["Date"].values[count]:
                            reentry["time"] = self.dsl["Time(UTC)"].values[count]
                            if np.isnan(reentry["abl_mass"]) or np.isnan(reentry["other_mass"]):
                                print(reentry["id"],reentry["abl_mass"],reentry["other_mass"])
                            time_update_mass_1 += (reentry["abl_mass"] +reentry["other_mass"])
                            time_update_count_1 += 1
                        else:    
                            reentry["time"] = 0
                            time_update_mass_2 += (reentry["abl_mass"] +reentry["other_mass"])
                            time_update_count_2 += 1
                            
            if reentry["time"] == -1:
                reentry["time"] = 0
                time_update_mass_2 += (reentry["abl_mass"] +reentry["other_mass"])
                time_update_count_2 += 1
    
        print(f"Time set to launch:   {int(time_update_mass_1)},{int(time_update_count_1)}")
        print(f"Time set to midnight: {int(time_update_mass_2)},{int(time_update_count_2)}")
        print(f"Time missing:         {missing_time_count}")

        self.dsl.close()
        self.dsr.close()
        self.ds_dw.close()
        
        self.print_stats()
        
    def reentry_info_to_netcdf(self):     
        """This saves the reentry information as a NetCDF file for later processing for GEOS-Chem.
        """        
        #Set up the dimensions of the netcdf file.
        dims = ('reentries')
        
        #Set up the data
        id_list, name_list, category_list, time_list = [], [], [], []
        datestr_list, lat_list, lon_list  = [], [], []
        abl_mass_list, abl_deg_list, abl_per_list, other_mass_list = [], [], [], []
        smc_list, location_list, apogee_list, burnup_list  = [], [], [], []
        
        for reentry in self.unique_reentry_list:
            id_list.append(reentry["id"])
            name_list.append(reentry["name"])
            category_list.append(reentry["category"])
            time_list.append(reentry["time"])
            datestr_list.append(reentry["datestr"])
            lat_list.append(reentry["lat"])
            lon_list.append(reentry["lon"])
            abl_mass_list.append(reentry["abl_mass"]+reentry["attached_abl_mass"])
            abl_deg_list.append(reentry["abl_deg"])
            abl_per_list.append(reentry["alu_per"])
            other_mass_list.append(reentry["other_mass"])
            smc_list.append(reentry["smc"])
            location_list.append(reentry["location"]) 
            apogee_list.append(reentry["apogee"])   
            burnup_list.append(reentry["burnup"])   
            
        #Create the DataArrays.
        data_da_id           = xr.DataArray(id_list,          dims=dims, attrs=dict(long_name="COSPAR_ID"))
        data_da_name         = xr.DataArray(name_list,        dims=dims, attrs=dict(long_name="Object Name"))
        data_da_category     = xr.DataArray(category_list,    dims=dims, attrs=dict(long_name="Category"))
        data_da_time         = xr.DataArray(time_list,        dims=dims, attrs=dict(long_name="Time (UTC)"))
        data_da_datestr      = xr.DataArray(datestr_list,     dims=dims, attrs=dict(long_name="Date"))
        data_da_lat          = xr.DataArray(lat_list,         dims=dims, attrs=dict(long_name="Latitude", units="Degrees"))
        data_da_lon          = xr.DataArray(lon_list,         dims=dims, attrs=dict(long_name="Longitude", units="Degrees"))
        data_da_abl_mass     = xr.DataArray(abl_mass_list,    dims=dims, attrs=dict(long_name="Ablatable Mass", units="kg"))
        data_da_abl_deg      = xr.DataArray(abl_deg_list,     dims=dims, attrs=dict(long_name="Ablation Degree"))
        data_da_per_alu      = xr.DataArray(abl_per_list,     dims=dims, attrs=dict(long_name="Percent Aluminium"))
        data_da_other_mass   = xr.DataArray(other_mass_list,  dims=dims, attrs=dict(long_name="Other Mass", units="kg"))
        data_da_smc          = xr.DataArray(smc_list,         dims=dims, attrs=dict(long_name="Megaconstellation_Flag"))
        data_da_location     = xr.DataArray(location_list,    dims=dims, attrs=dict(long_name="Location Constraint"))
        data_da_apogee       = xr.DataArray(apogee_list,      dims=dims, attrs=dict(long_name="Apogee", units="km"))
        data_da_burnup       = xr.DataArray(burnup_list,      dims=dims, attrs=dict(long_name="Burnup"))
    
        # Create an xarray Dataset from the DataArrays.
        ds = xr.Dataset()
        ds['COSPAR_ID']               = data_da_id
        ds['Object_Name']             = data_da_name
        ds['Category']                = data_da_category
        ds['Time (UTC)']              = data_da_time
        ds['Date']                    = data_da_datestr
        ds['Latitude']                = data_da_lat
        ds['Longitude']               = data_da_lon
        ds['Ablatable_Mass']          = data_da_abl_mass
        ds['Ablation_Degree']         = data_da_abl_deg
        ds['Percent_Aluminium']       = data_da_per_alu
        ds['Other_Mass']              = data_da_other_mass
        ds['Megaconstellation_Flag']  = data_da_smc
        ds['Location_Constraint']     = data_da_location
        ds['Apogee']                  = data_da_apogee
        ds['Burnup']                  = data_da_burnup
             
        #Save to file and close the DataSet  
        ds.to_netcdf(f'./databases/reentry_activity_data_{self.start_year}-{self.final_year}.nc')
        ds.close()
        
if __name__ == "__main__":  
    
    # Set up the arguments for each function.
    parser = argparse.ArgumentParser()
    parser.add_argument('-sv', "--save_reentry_info", action='store_true', help='Save reentry info.')
    parser.add_argument('-usegpd', "--use_geopandas", action='store_true', help='Load in geopandas dataframes.')
    parser.add_argument('-sy', "--start_year", default = "2023", choices=str(np.arange(1957,2026)), help='Start Year.')
    parser.add_argument('-fy', "--final_year", default = "2024", choices=str(np.arange(1957,2026)), help='Final Year.')
    args = parser.parse_args()
    use_gpd = args.use_geopandas
    
    # Sort out the year range.
    start_year = int(args.start_year)
    final_year = int(args.final_year)
    if start_year > final_year:
        final_year = start_year + 1
    print(f"Processing from year {start_year} to {final_year}.")
    
    # Compile the reentry data
    Data = build_reentry_list(start_year, final_year)
    Data.get_reentry_info()
    if args.save_reentry_info == True:
        Data.reentry_info_to_netcdf()
        
    # TODO: Check what Jonathan lists as the YYYY-UXX launches. These seem to be missing in DISCOSweb but include Starship launches above 50 km and a Chang'e-6.
    # Missing all 2024 starship failures (guess these are listed as suborbital).
    # Missing latest May 27th 2025 starship failure.
    # Sort out POEM 4
    # Add Tranche, Ronghe, Digui, Hulianwang, Qianfan eventually when they reenter.

    def test():
        return build_reentry_list(start_year, final_year)
    profiler = Profile()
    profiler.runcall(test)
    stats = Stats(profiler)
    stats.strip_dirs()
    stats.sort_stats('cumulative')
    stats.print_stats()