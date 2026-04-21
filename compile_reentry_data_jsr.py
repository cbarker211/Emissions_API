import numpy as np
import argparse
import xarray as xr
import pandas as pd
from datetime import datetime
import sys
import random
import geopandas as gpd
from shapely.geometry import Point, Polygon, shape, box
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
        yearstr = yearstr.replace("s","")
        monstr = "01"   
        daystr = "01"
    elif date[1] == "Q3?": # 1995-005A
        monstr = "07"
        daystr = "01"  
    elif date[1] == "Q4?": # 2015-025H
        monstr = "10"
        daystr = "01"
    else:   
        date[1] = date[1].replace("?","")
        try:
            monstr = str(datetime.strptime(date[1], '%b').month).zfill(2)  
        except:
            sys.exit(date)
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
        self.latlonerror = []
        
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

            self.shapefile_map = {
                "Gujarat":      (self.loaded_shapefiles["States"][self.loaded_shapefiles["States"]["name_en"] == "Gujarat"], 3),
                "Kazakhstan":   (self.loaded_shapefiles["Countries"][self.loaded_shapefiles["Countries"]["NAME"] == "Kazakhstan"], 3),
                "KRZ":          (self.loaded_shapefiles["Countries"][self.loaded_shapefiles["Countries"]["NAME"] == "Kazakhstan"], 3),
                "Kaz.":          (self.loaded_shapefiles["Countries"][self.loaded_shapefiles["Countries"]["NAME"] == "Kazakhstan"], 3),
                "Mexico":       (self.loaded_shapefiles["Countries"][self.loaded_shapefiles["Countries"]["NAME"] == "Mexico"], 3),
                "S Africa":     (self.loaded_shapefiles["Countries"][self.loaded_shapefiles["Countries"]["NAME"] == "South Africa"], 3),
                "N Zealand":    (self.loaded_shapefiles["Countries"][self.loaded_shapefiles["Countries"]["NAME"] == "New Zealand"], 3),
                "China":        (self.loaded_shapefiles["Countries"][self.loaded_shapefiles["Countries"]["NAME"] == "China"], 3),
                "Gulf":         (self.loaded_shapefiles["Gulf of Mexico"][self.loaded_shapefiles["Gulf of Mexico"]["name"] == "Gulf of Mexico"], 4),
                "Indian O":     (self.loaded_shapefiles["Indian Ocean"][self.loaded_shapefiles["Indian Ocean"]["name"] == "Indian Ocean"], 4),
                "Indian Ocean": (self.loaded_shapefiles["Indian Ocean"][self.loaded_shapefiles["Indian Ocean"]["name"] == "Indian Ocean"], 4),
                "Indian O.":    (self.loaded_shapefiles["Indian Ocean"][self.loaded_shapefiles["Indian Ocean"]["name"] == "Indian Ocean"], 4),
                "IOR":          (self.loaded_shapefiles["Indian Ocean"][self.loaded_shapefiles["Indian Ocean"]["name"] == "Indian Ocean"], 4), 
                "S POR":        (self.loaded_shapefiles["Pacific"][self.loaded_shapefiles["Pacific"]["name"] == "South Pacific Ocean"], 4),
                "S Pacific":    (self.loaded_shapefiles["Pacific"][self.loaded_shapefiles["Pacific"]["name"] == "South Pacific Ocean"], 4),
                "S Pac.":       (self.loaded_shapefiles["Pacific"][self.loaded_shapefiles["Pacific"]["name"] == "South Pacific Ocean"], 4),
                "N Atlantic":   (self.loaded_shapefiles["Atlantic"][self.loaded_shapefiles["Atlantic"]["name"] == "North Atlantic Ocean"], 4),
                "S AOR":        (self.loaded_shapefiles["Atlantic"][self.loaded_shapefiles["Atlantic"]["name"] == "South Atlantic Ocean"], 4),
                "S Atl":        (self.loaded_shapefiles["Atlantic"][self.loaded_shapefiles["Atlantic"]["name"] == "South Atlantic Ocean"], 4),
            }
        else:
            self.shapefile_map = dict.fromkeys(["Gujarat","Kazakhstan","KRZ","Kaz.","Mexico","S Africa","N Zealand","China","Gulf",
                                                "Indian O""Indian Ocean","Indian O.","IOR","S POR","S Pacific","S Pac.",
                                                "N Atlantic","S AOR","S Atl"])

        self.provided_coords = {
            "GM 86.2W 29.7N": (29.7, -86.2, 1),     # Gulf of Mexico named coordinates.
            "32W 1E": (1, -32, 1),                  # Mistake in GCAT.
            "Alashan": (38.833, 105.95, 2),         # Landing expected east of Jiuquan launch site (this was actually launched from Wenchang).
                                                        # Assuming that Alashan refers to Helen Shan Mountain (new name for Alashan Mountain).
                                                        # https://space.skyrocket.de/doc_sdat/rcs-fc-sc.htm
                                                        # https://planet4589.org/space/jsr/news.778
                                                        # https://www.peakbagger.com/peak.aspx?pid=10693
                                                        # https://www.nasaspaceflight.com/2020/05/china-next-generation-crew-capsule/
            "Jiuquan": (41.3, 100.3, 2),            # This refers to the Chinese launch center, sources say it landed 'in the desert near Jiuquan', and 'in China’s Inner Mongolia autonomous region'.
                                                        # Using lat/lon for Jiuquan launch site.
                                                        # https://discosweb.esoc.esa.int/launch-sites/21
                                                        # https://nssdc.gsfc.nasa.gov/nmc/spacecraft/display.action?id=2020-027A
                                                        # https://space.skyrocket.de/doc_sdat/xzf-sc.htm 
                                                        # https://spaceflightnow.com/2020/05/08/chinas-next-generation-crew-spacecraft-lands-after-unpiloted-test-flight/ 
            "Dongfeng": (41.3, 100.3, 2),           # Same as Jiuquan.
            "LOPNOR RW05": (40.78, 89.27, 2),       # Reported test flight of a Chinese reusable experimental spacecraft.
                                                        # Runway 05/22 at Lop Nor air base in Xinjiang (https://planet4589.org/space/gcat/data/tables/lp.html).
                                                        # https://www.seradata.com/china-launches-own-mini-spaceplane-reusable-spacecraft-using-long-march-2f 
                                                        # https://twitter.com/planet4589/status/1302486141090885632
            "Lop Nor": (40.78, 89.27, 2),           # Same as LORNOR. 
            "VSFB RW30/12": (34.7, -120.6, 2),      # Vandenberg Space Force Base, California, USA
            "Koonibba": (-31.886, 133.449, 2),      # Australian Test Range
            "Jacklyn": (27.888, -74.153, 2),        # https://x.com/spaceOFFSHORE/status/1878177979974775139/photo/1
            "Ocean": (-42.86, 177.87, 2),           # TODO: Check this doesn't include others
                                                        # Electron Stage 1 (rcat) is occasionally recovered in the ocean.
                                                        # From https://spaceflightnow.com/2020/11/05/rocket-lab-to-attempt-booster-recovery-on-next-mission/:
                                                        #   Recovery vessels stationed near the booster’s splashdown zone around 250 miles (400 kilometers) 
                                                        #   south of the launch site will move in to secure the first stage and hoist it onto a ship for return to New Zealand.
                                                        # Geoloating 400km south of launch site.
            "STB OLP1": (25.996, -97.154, 2),       # Starbase in Texas.
            "WSSH": (33.238, -106.346, 2),          # White Sands Missile Range in New Mexico.
            "Splash": (-36.57, 177.63, 2),          # This is the Electron failed helicopter stage 1 recovery in May 22.
                                                        # https://www.youtube.com/watch?v=BY0CXlOeWHI "Several hundred kilometers from the launch site."
                                                        # Using inclination of 94 degrees (wiki), and distance of 300km.
            "KSC RW33": (28.61, -80.69, 2),         # Kennedy Space Center Shuttle Landing Facility  
            "KSC RW15": (28.61, -80.69, 2),         # Kennedy Space Center Shuttle Landing Facility                                 
            "KSC SLF": (28.61, -80.69, 2),          # https://www.world-airport-codes.com/united-states/nasa-shuttle-landing-facility-69738.html
                                                        # https://www.spaceforce.mil/News/Article/3217077/x-37b-orbital-test-vehicle-concludes-sixth-successful-mission/
            "S Pole": (-90.0, 0.0, 2),              # South Pole (obviously...)
            "Africa": (16.4, 18.6, 2),              # Pioneer 3 re-entered over Africa
                                                        # https://web.archive.org/web/20200410070534/https://nssdc.gsfc.nasa.gov/nmc/spacecraft/display.action?id=1958-008A
            "Akmolinsk SW": (50.218, 69.910, 2),    # Kosmos 10 "Landed 150 km SW of Akmolinsk." https://www.orbitalfocus.uk/Diaries/Zenit/Zenit-2.php
            "EAFB RW22":    (34.924, -117.892, 2),  # Edwards Air Force Base
            "EAFB RW05R":   (34.924, -117.892, 2),  # Edwards Air Force Base
            "EAFB RW04L":   (34.924, -117.892, 2),  # Edwards Air Force Base
            "EAFB RW17":    (34.924, -117.892, 2),  # Edwards Air Force Base
            "EAFB RW23":    (34.924, -117.892, 2),  # Edwards Air Force Base
            "EAFB RW04":    (34.924, -117.892, 2),  # Edwards Air Force Base
            "EAFB RW33":    (34.924, -117.892, 2),  # Edwards Air Force Base
            "Arkalyk":      (50.249, 66.902, 2),    # Arkalyk
            "Wisconsin":    (44.099, -87.658, 2),   # https://www.smithsonianmag.com/air-space-magazine/when-sputnik-crashed-wisconsin-180952388/
            "Uralsk W 130": (51.224, 51.373, 2),    # Uralsk
            "Kustanai":     (53.2, 63.62, 2),       # Kostanay
        }
            
    def print_stats(self):
        """Print out statistics about the database at the end of the program.
        """        

        total_reentries = len(self.df_reentry)
        non_geo_mask = self.df_reentry["lat"].isna() & self.df_reentry["lon"].isna()
        non_geo_count = non_geo_mask.sum()
        geolocated_mass = ((self.df_reentry.loc[~non_geo_mask, "abl_mass"] + self.df_reentry.loc[~non_geo_mask, "other_mass"]).sum())

        non_time_mask = self.df_reentry["time"] == -1
        non_time_count = non_time_mask.sum()

        total_mass_array = self.df_reentry["abl_mass"] + self.df_reentry["other_mass"]
        total_mass = total_mass_array.sum()
        non_mass_mask = total_mass_array == 0
        non_mass_count = non_mass_mask.sum()

        geo_percent_all =  non_geo_count  / total_reentries * 100
        time_percent_all = non_time_count / total_reentries * 100
        mass_percent_all = non_mass_count / total_reentries * 100
        geolocated_mass_percent = geolocated_mass / total_mass * 100

        data = {'Property':  ['Total Reentries', 'Total mass [Gg]','Non-geolocated Objects','Non-geolocated Objects [%]',
                              'Non-timed Objects','Non-timed Objects [%]','Objects with no Mass','Objects with no Mass [%]','Geolocated Mass'],
                'Value':     [int(total_reentries), total_mass*1e-6, int(non_geo_count), geo_percent_all,
                              int(non_time_count), time_percent_all, int(non_mass_count), mass_percent_all, geolocated_mass_percent]}
        
        df = pd.DataFrame(data)
        df_rounded = df.round(1)
        print(df_rounded)  
        df_rounded.to_csv(f"./out_files/reentry_stats.csv",sep=',',index=False)                    
                
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
        
        lat, lon, location = None, None, -1

        def sample_within_region(region, lat_range=(-inc, inc)):

            while True:
                coords = region.sample_points(1).get_coordinates()
                lon, lat = coords["x"].values[0], coords["y"].values[0]
                if lat_range[0] <= lat <= lat_range[1]:
                    break   
            return lat, lon

        # If location is missing for lower stages for failed and successful launches, set as launch coordinates.    
        if latlonstr in ["-",""]:
            
            # If this object didn't go above 100 km and is a lower stage, then set the re-entry coordinates to the launch site.
            # This assumes the objects haven't travelled far from the launch site before re-entering (usually <400 km).
            if (apogee <= 100) and (category in ["S0","S1"]): 
                for i, cpid in enumerate(self.dsl["COSPAR_ID"].values):
                    if jsr_id[:8] == cpid[:8]:
                        lat = self.dsl["Latitude"].values[i]
                        lon = self.dsl["Longitude"].values[i]
                        location = 2
                        break  
            else: 
                lat = round(np.random.uniform(-inc, inc),2)
                lon = round(np.random.uniform(-180, 180),2)
                location = 6  
    
        # Lat/Lon Provided, Launch Site, or Named Location 
        elif latlonstr in self.provided_coords:
            lat, lon, location = self.provided_coords[latlonstr]
 
        # Ocean/Continent
        elif latlonstr in ["Antarctic", "Arctic", "S Ocean","S Oc.","S OCean"]:
            location = 4
            lon = round(np.random.uniform(-180, 180), 2)

            if latlonstr in ["S Ocean","S Oc.","S OCean"]:
                # For these objects in auxcat, the inclination is 51.65-74, outside of the general definition of the southern ocean (<60S).
                # Therefore we will just set the lat to the inclination.
                lat = -inc
            elif latlonstr == "Antarctic":
                lat = round(np.random.uniform(-90, min(-60, -inc)), 2)
            else:  # Arctic
                lat = round(np.random.uniform(max(66, inc), 90), 2)
            coordinate = Point(lon, lat)

        # Regions defined by shapefiles.
        elif use_gpd == False and (latlonstr in list(self.shapefile_map.keys()) or latlonstr in [
            "Pacific","PO","E Pacific","POR","SW POR","SE Pacific","W POR","NE RU/NW POR","SW POR","N Pacific",
            "SE IOR","SW IOR",
            "Atlantic","SW AOR","AOR","AO",
            "E USA","SE Arkalyk","S of Aus","S of Austr.","S of S Afr.", "SE Hawaii"]):
        
            lat, lon, location = 0, 0, -1

        elif latlonstr in self.shapefile_map:
            region, location = self.shapefile_map[latlonstr]
            lat, lon = sample_within_region(region, lat_range=(-inc, inc))
         
        elif latlonstr in ["Pacific", "PO", "POR","E Pacific"]:

            PACIFIC_BOUNDS = {
                "Pacific":      (-180,  180, None, None),
                "PO":           (-180,  180, None, None),
                "POR":          (-180,  180, None, None),
                "E Pacific":    (-180,  -60, None, None),
                "SW POR":       (-180,    0, None,    0),  # southern: lat_max = 0
                "SE Pacific":   (-180,  -60, None,    0),  # southern: lat_max = 0
                "W POR":        (-180,    0, None, None),
                "N Pacific":    (-180,  180,    0, None),  # northern: lat_min = 0
                "NE RU/NW POR": (-180,    0,    0, None),  # northern: lat_min = 0
            }
                    
            location = 4
            lon_min, lon_max, elat_min, elat_max = PACIFIC_BOUNDS[latlonstr]
            lat_min = max(-inc, elat_min) if elat_min is not None else -inc
            lat_max = min( inc, elat_max) if elat_max is not None else  inc
            coords = self.loaded_shapefiles["Pacific"].clip(box(lon_min, lat_min, lon_max, lat_max)).sample_points(1).get_coordinates()
            lon, lat = coords["x"].values[0], coords["y"].values[0]

        elif latlonstr in ["AO","Atlantic"]: 
            location = 4 
            while True:
                random_location = self.loaded_shapefiles["Atlantic"].sample_points(1).get_coordinates()
                lon = random_location["x"].values[0]
                lat = random_location["y"].values[0]
                if -inc <= lat <= inc:
                    break
        elif latlonstr in ["S of Tasmania"]:
            location = 4
            lon_lat_list = [[153.23, -30.02],[146.83, -43.64],[166.00, -50.92], [167.53, -47.29],
                            [168.14, -46.85], [168.85, -46.66], [174.28, -41.72], [175.28, -41.62], [173.01, -34.39]]
            tasman_sea_geom = Polygon(lon_lat_list)
            tasman_sea = gpd.GeoDataFrame(index=[0], crs='epsg:4326', geometry=[tasman_sea_geom])
            while True:
                lat = round(np.random.uniform(-inc, inc),2)
                lon = round(np.random.uniform(-180, 180),2)
                coordinate = Point(lon,lat)
                if tasman_sea.geometry.contains(coordinate).any():
                    break 
        elif latlonstr in ["SE IOR","SW IOR"]:
            location = 4
            region, location = self.shapefile_map[latlonstr]
            while True:
                coords = region.sample_points(1).get_coordinates()
                lon, lat = coords["x"].values[0], coords["y"].values[0]
                if (-inc <= lat <= 77) and (-30 <= lon <= 180):
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
                else:
                    raise ValueError(f"Can't find lon/lat coords for {jsr_id} - {latlonstr}.")
                
                if "W" in lon_data[-1]:
                    lon = np.float64((lon_data[0]) + str(np.float64(lon_data[1].replace("W",""))/60))*-1
                elif "E" in lon_data[-1]:
                    lon = np.float64((lon_data[0]) + str(np.float64(lon_data[1].replace("E",""))/60))
                    
                if "N" in lat_data[-1]:
                    lat = np.float64((lat_data[0]) + str(np.float64(lat_data[1].replace("N",""))/60))
                elif "S" in lat_data[-1]:
                    lat = np.float64((lat_data[0]) + str(np.float64(lat_data[1].replace("S",""))/60))*-1
            
        elif latlonstr.split()[0][-1] in ["E","W","S","N"]: # Coordinates in degrees
            location = 1
            geolocation = latlonstr.split()
            if len(geolocation) == 1:
                if geolocation[0] == '155W50S':
                    lon_data = "155W"
                    lat_data = "50S"
                else:
                    raise ValueError(f"Can't find lon/lat coords for {jsr_id} - {latlonstr}.")
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
                print(f"Missing lon Direction for {jsr_id} - {latlonstr} - Assuming E")
                lon = lon_data

            if lat_data[-1] == "N":
                lat = np.float64(lat_data.replace("N",""))
            elif lat_data[-1] == "S":
                lat = np.float64(lat_data.replace("S",""))*-1  
            else:
                print(f"Missing lat direction for {jsr_id} - {latlonstr} - Assuming N")
                lat = lat_data 

        if location is None or lat is None or lon is None:
            self.latlonerror.append(latlonstr)
        else:
            lat = np.float64(lat)
            lon = np.float64(lon)

        return lat, lon, location
    
    def falcon_stage_lat_lon(self,datestr,jsr_id,dest,plname):
        
        """Set the geolocation for all Falcon 9 1st stages.

        Returns:
            lat(np.float64) : The latitude of the object.
            lon(np.float64) : The longitude of the object..
        """ 

        lat, lon = None, None

        # Find the ground landings in Raul's SpaceX Map
        if "LZ" in dest:
            matching = self.ground_landings[self.ground_landings["Date"] == datestr].reset_index(drop=True)

            # If its in the list of ground landings, then we can geolocate it.
            if matching.shape[0] == 1: 
                lat = matching["geometry"].y.iloc[0]
                lon = matching["geometry"].x.iloc[0]   

            # If its not in the list, then as its a ground landing, we can set it to the launch site.
            elif matching.shape[0] == 0:
                idx = np.where(self.dsl["COSPAR_ID"] == jsr_id)[0]
                if len(idx) == 0:
                    raise ValueError(f"Problem geolocating Falcon Stage 1 recovery - {jsr_id}.")
                lat = np.round(self.dsl["Latitude"].values[idx[0]]) 
                lon = np.round(self.dsl["Longitude"].values[idx[0]])
                
            elif matching.shape[0] > 1:
                # If its a Falcon Heavy launch with two landings, then we set both to the same location.
                if plname.startswith("FH-"):
                    lat = matching["geometry"].y.iloc[0]
                    lon = matching["geometry"].x.iloc[0]
                # Otherwise, raise an error.
                else:
                    raise ValueError(f"Multiple ground entries for Falcon Stage 1 landing- {jsr_id}.")
            else:
                raise ValueError(f"Problem geolocating Falcon Stage 1 landing- {jsr_id}.")
        else:
            # Find the ocean landings in Raul's SpaceX Map.
            matching = self.ocean_landings[self.ocean_landings["Date"] == (datestr)].reset_index(drop=True)
            if (matching.shape[0] == 1) or (datestr == "20220427"):
                # This is 27th April - typo in Raul's Space Map where second should be 2023. 
                lat = matching["geometry"].y.iloc[0]
                lon = matching["geometry"].x.iloc[0]
            # These next two are where there are two launches on the same day.
            elif jsr_id == "2022-124":
                lat = matching["geometry"].y.iloc[0]
                lon = matching["geometry"].x.iloc[0]
            elif jsr_id == "2022-125":
                lat = matching["geometry"].y.iloc[1]
                lon = matching["geometry"].x.iloc[1]
            # Again these are two on the same day.
            elif jsr_id == "2023-037":
                lat = matching["geometry"].y.iloc[0]
                lon = matching["geometry"].x.iloc[0]
            elif jsr_id == "2023-038":
                lat = matching["geometry"].y.iloc[1]
                lon = matching["geometry"].x.iloc[1]

            # If its not in the list, then set it to the ocean near the launch site. 
            # These coordinates are the most common grid box for reentries in 2020-2022.
            elif matching.shape[0] == 0:

                idx = np.where(self.dsl["COSPAR_ID"] == jsr_id)[0]
                if len(idx) == 0:
                    raise ValueError(f"Problem geolocating Falcon fairing recovery - {jsr_id}.")
                lon_val = np.round(self.dsl["Longitude"].values[idx[0]])
                lat_val = np.round(self.dsl["Latitude"].values[idx[0]])

                # Kennedy Space Center (ETR)
                if (lon_val == -81) and (lat_val in [28,29]):
                    lat = 34
                    lon = -75  
                # Vandenberg
                elif (lon_val == -121) and (lat_val == 35):
                    lat = 30
                    lon = -120 
                else:
                    raise ValueError(f'Falcon launch site not found - {lon_val} {lat_val}.')
            elif matching.shape[0] > 1:
                print(matching,datestr)
                raise ValueError(f"Multiple ocean entries for Falcon Stage 1 landing- {jsr_id}.")
            else:
                raise ValueError(f"Problem geolocating Falcon Stage 1 landing- {jsr_id}.") 
            
        return lat, lon 
    
    def falcon_fairing_lat_lon(self,datestr,jsr_id):   
        
        """Set the geolocation for all Falcon 9 fairings.

        Returns:
            lat(np.float64) : The latitude of the object.
            lon(np.float64) : The longitude of the object..
        """

        lat, lon = None, None

        set_ocean = False
        if len(self.fairings) > 0:
            matching = self.fairings[self.fairings["Date"] == datestr].reset_index(drop=True)
            if matching.shape[0] == 1:
                lat = matching["geometry"].y.iloc[0]
                lon = matching["geometry"].x.iloc[0]
            # These are two on the same day
            elif jsr_id == "2023-037":
                lat = matching["geometry"].y.iloc[0]
                lon = matching["geometry"].x.iloc[0]
            elif jsr_id == "2023-038":
                lat = matching["geometry"].y.iloc[0]
                lon = matching["geometry"].x.iloc[0]
            # If its not in the list, then set it to the ocean near the launch site. 
            # These coordinates are the most common grid box for reentries in 2020-2022.
            elif matching.shape[0] == 0:
                set_ocean = True
            elif matching.shape[0] > 1:
                print(matching)
                raise ValueError(f"Multiple entries for Falcon fairing recovery - {jsr_id}.")
            else:
                raise ValueError(f"Problem geolocating Falcon fairing recovery - {jsr_id}.")
        else:
            set_ocean = True

        if set_ocean:
            idx = np.where(self.dsl["COSPAR_ID"] == jsr_id)[0]
            if len(idx) == 0:
                raise ValueError(f"Problem geolocating Falcon fairing recovery - {jsr_id}.")
            
            lon_val = np.round(self.dsl["Longitude"].values[idx[0]])
            lat_val = np.round(self.dsl["Latitude"].values[idx[0]])

            # Kennedy Space Center (ETR)
            if (lon_val == -81) and (lat_val in [28,29]):
                lat = 34
                lon = -75  
            # Vandenberg
            elif (lon_val == -121) and (lat_val == 35):
                lat = 30
                lon = -120
            else:
                raise ValueError(f'Falcon launch site not found - {lon_val} {lat_val}.')
              
        return lat, lon
    
    def failed_launch_mass(self, jsr_id, jsr_name, reentry_category, drymass):
        
        """For all failed launches, set the masses manually as we need to differentiate wet vs dry mass reentry.

        Returns:
            mass(float)       : The mass of aluminium.
            other_mass(float): Extra non-aluminium mass.
        """        
            
        abl_mass = 0
        other_mass = 0

        rocket_name = self.dsl.loc[self.dsl["COSPAR_ID"] == jsr_id, "Rocket_Name"].iloc[0]
        rocket_ind = np.where(self.dsr["Rocket_Name"].values == rocket_name)[0][0]

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
        
    def sort_inclination(self,inc_eff,jsr_id):
        
        """When the inclination on file is zero (including failed), this function sets it appropriately.
        This is mainly the case for lower suborbital stages, which can then be set using the inclination of an upper stage or payload.

        Returns:
            inc_eff(np.float64): The orbital inclination of the object.
        """        
        
        # Use the maximum inclination for all objects with the same launch prefix.
        inc_eff = max(inc_eff, self.max_inc_by_prefix.get(jsr_id[:8], 0))       
            
        # Sometimes, there is no entry in the reentry list because it didn't meet the criteria before.
        # In this case we need to look in the other inventories to look for any matching entries.
        if inc_eff == 0:
            if "F" in jsr_id:
                inc_eff = self.inc_by_prefix_failed.get(jsr_id[:8], 0)
            else:
                inc_eff = self.inc_by_prefix_all.get(jsr_id[:8], 0)

        if inc_eff == 0 and jsr_id[5] not in ["F","U"]:
            # TODO: Sort failed launches.
            raise ValueError(f"Inclination still empty for {jsr_id}.")
        return inc_eff 
         
    def extract_jsr_info(self, file):
        """For each GCAT database, extract all the relevant information and output to the dictionary of all reentry objects.
        """        
        
        jsr_data = self.jsr_data_dict[file]
        self.max_inc_by_prefix = {}
        
        ################################################
        # Filter the database for the relevant entries.
        ################################################

        subset = jsr_data[
            (jsr_data["Status"].isin(["R","R?","D","L","L?","S","F","AF","AS"])) &            # Only include specific statuses (see https://planet4589.org/space/gcat/web/intro/phases.html)       
            (~jsr_data["#JCAT"].isin(["L80508","L80509","L80510"])) &                         # Skipping 2021-F07, GCAT lists apogee as 100km, but it actually failed around 30 km and reached a height of 50 km.
            (jsr_data["DYear"].between(self.start_year, self.final_year)) &                   # This time range only.
            ((jsr_data["Apogee"] >= 50) | (jsr_data["Apogee"] < 0)) &                         # Apogee above 50 km, or below 0. This makes sure to include objects in lprcat returning from space.
            ((jsr_data["Status"] != "AS") | (jsr_data["Piece"].str[5:6].isin(["F","U"]))) &   # Only include AS re-entries if they are part of a failed or uncategorized launch.
            (~jsr_data["#JCAT"].isin(["S46138","S44635","S40899"]))                           # Skip objects where the date is wrong in GCAT.
        ]

        # Ignore failed launches pre-2020.
        mask = ~(subset["Piece"].str[5].isin(["F", "U"]) & (subset["DYear"] < 2020))          
        subset = subset[mask]

        # Extract the arrays before looping.
        cols = ["DDate","Piece","Name","Type","Inc","Dest","Apogee","Status","DryMass","#JCAT","Bus","Parent","PLName"]
        arrays = {col: subset[col].to_numpy() for col in cols}

        # Now loop over the list and format into a dictionary.
        for i in tqdm(range(len(arrays["Piece"])), total=len(subset), desc=f"Processing {file}"):

            # Sort out the reentry time/date.
            datestr, time_utc = convert_time(arrays["DDate"][i].split()) 

            # Create some variables for easier access.
            jsr_id      = self.convert_launch_tag(arrays["Piece"][i])
            jsr_name    = str(arrays["Name"][i])
            jsr_type    = str(arrays["Type"][i]).strip()
            jsr_inc     = arrays["Inc"][i]
            jsr_dest    = arrays["Dest"][i]
            jsr_apogee  = arrays["Apogee"][i]
            jsr_status  = arrays["Status"][i]
            jsr_drymass = arrays["DryMass"][i]
            jsr_jcat    = arrays["#JCAT"][i]
            jsr_parent  = arrays["Parent"][i]
            jsr_bus     = arrays["Bus"][i]
            jsr_plname  = arrays["PLName"][i]
            
            burnup = "Complete"
            # Set the burnup as partial for all objects landing or splashing down at the surface.
            if jsr_status == "L":
                burnup = "Partial"
            # NOTE: Need to manually check for Electron booster recoveries if they ever restart recoveries.
            elif (jsr_name == "Electron Stage 1") and jsr_id in ["2019-084","2020-007","2020-085","2021-F02","2021-106",
                                                                 "2022-047","2022-147","2023-041","2023-100","2023-126","2024-022"]:
                    burnup = "Partial"
                 
            # Duplicates:
            #   2019-036H. This is TEPCE 1 and TEPCE 2, listed with same "Piece" in JSR. TEPCE 2 listed as 2019-036IA in DISCOSweb, so renaming to this.
            #   2020-027A. This is Xinyidai Zairen Feichuan (XZF) and XZF Service Module (not in DISCOSweb). Setting auxcat id to XZF 2020-027D.
            #   Remaining are where multiple objects in auxcat from same launch are just given the COSPAR ID of the launch.
            # TODO: Email Jonathan about these.
            if jsr_name == "TEPCE 2":
                jsr_id = "2019-036IA" 
            elif jsr_name == "XZF Service Module":
                jsr_id = "2020-027D"
            
            # Handle categories.
            if jsr_type[0] == "R":
                rocket_name = self.cospar_to_rocket[jsr_id[:8]]
                
                # Sort out the rocket stages. This is mainly Soyuz (Russian's just do the numbering differently), and errors in the GCAT.
                reentry_category = f"S{(jsr_type[1:2])}"

                # Setting Falcon Heavy Boosters as S0.
                if rocket_name == "Falcon Heavy" and jsr_name.startswith("Falcon 9 Stage 1") and jsr_parent != "-":
                    reentry_category = "S0"

                # Shift the stage numbers down by one.
                # Saturn V has the interstage listed as a separate stage in the GCAT launch vehicle list, but seems to not include it in the re-entry list.
                # So we don't need to shift for Saturn V.
                if rocket_name in ["Space Shuttle","SLS Block 1","Conestoga 1620"]:
                    stage = int(jsr_type[1:2]) - 1
                    reentry_category = f"S{stage}"

                # Soyuz
                if "Blok-BVGD" in jsr_bus:
                    reentry_category = "S0"
                elif "Blok-A" in jsr_bus:
                    reentry_category = "S1"
                elif "Blok-I" in jsr_name:
                    reentry_category = "S2"
                elif any(name in jsr_name for name in ["Blok-L", "Fregat"]):
                    reentry_category = "S3"

                if len(reentry_category) == 1:
                    print(jsr_type)
            else:
                reentry_category = jsr_type[0]

            # Adjust the inclination so its useful for our purpose of bounding the lat.
            if jsr_inc > 90:
                inc_eff = 180 - jsr_inc
            else:
                inc_eff = jsr_inc
                
            # Sort the inc.
            if jsr_inc == 0:
                inc_eff = self.sort_inclination(inc_eff,jsr_id)

            # Keep track of the maximum inclination for each object prefix.
            if inc_eff != 0:
                self.max_inc_by_prefix[jsr_id[:8]] = max(inc_eff, self.max_inc_by_prefix.get(jsr_id[:8], 0))

            # Sort out the reentry lat/lon for Falcon 9 first stages and fairings.
            if ("Falcon 9 Stage 1" in jsr_name):

                # If the first stage landed on a drone ship or landing zone, use the known location from Raul's SpaceX map.
                if "LZ" in jsr_dest or jsr_status in ["L","LF"] or jsr_dest == "OCISLY":
                    lat, lon = self.falcon_stage_lat_lon(datestr,jsr_id,jsr_dest, jsr_plname) 
                    location = 5

                # If the first stage was expended, then use the normal lat/lon conversion.
                elif jsr_status in ["S","R"]:
                    lat, lon, location = self.convert_lat_lon(jsr_dest.replace("?",""), inc_eff, reentry_category, jsr_apogee, jsr_id)
                else:
                    raise ValueError(f"Unexpected Falcon Stage 1 Status {jsr_status} {jsr_jcat}")
            
            elif ("Falcon 9 Fairing" in jsr_name) or ("Falcon Heavy Fairing" in jsr_name):    
                lat, lon = self.falcon_fairing_lat_lon(datestr,jsr_id)
                burnup = "Partial"
                location = 5

            else:   
                lat, lon, location = self.convert_lat_lon(jsr_dest.replace("?",""), inc_eff, reentry_category, jsr_apogee, jsr_id)
            
            # Sort the mass.    
            # # TODO: Sort failed launches back to 1957. 
            if jsr_id[5:6] in ["F","U"]:
                abl_mass, other_mass = self.failed_launch_mass(jsr_id, jsr_name, reentry_category,jsr_drymass)
            else:
                abl_mass = np.float64(jsr_drymass)
                other_mass = 0
            
            # Check for items with a missing geolocation (should have been dealt with already so this is just a sanity check).    
            if use_gpd == True and lat is None and lon is None:
                print(f"No location for {jsr_id}, dest = {jsr_dest}")
            
            # Set up the dictionary.                      
            self.reentry_rows.append({
                "id"               : jsr_id,
                "jcat"             : jsr_jcat,
                "name"             : jsr_name,
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
                "parent"           : jsr_parent
            })
            
    def add_attached(self):
        
        """
        Function to add cargo or attached mass to its parent object.
        Skipping objects attached or inside an object that has exploded, as we can't accurately estimate the emissions.
        """

        # Combine the catalogues, trimming to only include objects which have landed or re-entered inside or attached to another object within the date range.
        keys = ["satcat","auxcat","lcat","rcat","lprcat","deepcat","ecat"]
        attached_cat  = pd.concat([self.jsr_data_dict[k] for k in keys], ignore_index=True)
        attached_cat  = attached_cat[
            (attached_cat["Primary"] == "Earth") &                               # Primary body is Earth only.
            (attached_cat["Status"].isin(["AR","AL","AL IN","AR IN","E"])) &     # Attached objects only.
            (~attached_cat["Piece"].str[5:6].isin(["S","F","U"])) &              # Skip suborbital and failed launches (treated separately).
            (attached_cat["DYear"].between(self.start_year, self.final_year)) &  # Only the selected timeframe.
            (attached_cat["PLName"] != "Manfred Mem. Moon Mission")              # This should be listed as re-entering on the Moon, so ignore.
        ].copy()

        attached_cat["parent_clean"] = attached_cat["Parent"].str.replace("*","", regex=False)            

        # A function to find the parent of any jcat.
        def find_parent(jcat):
            parent = None
            match = attached_cat[attached_cat['#JCAT'] == jcat]
            if len(match) > 0:
                if match["Status"].values[0] == "E":
                    parent = f"Explosion {match['#JCAT'].values[0]} {match['Name'].values[0].strip()}"
                else:
                    parent = match["Parent"].values[0].split()[0]
            
            return parent

        # Loop over the objects and add the mass to the ancestor.  
        failed_mass, failed_count = 0, 0          
        for idx, row in attached_cat.iterrows():

            if row["Status"] == "E":
                continue

            parent = row["parent_clean"].split()[0]
            if parent == "S16482": # This parent object is Satcom K1, but should be Columbia (STS 61-C) (S16481). 
                parent = "S16481"
            found = False

            # If parent is already in the full reentry list.
            if parent in self.df_reentry.index:
                if "IN" in row["Status"]:
                    self.df_reentry.loc[parent, "other_mass"] += float(row["DryMass"]) # type: ignore
                else:
                    self.df_reentry.loc[parent, "attached_abl_mass"] += float(row["DryMass"]) # type: ignore
                found = True

            # If not, then need to find the grandparent.
            if not found and parent is not None:
                ancestor = parent
                max_searches = 5

                #  It might be that the grandparent is also attached or inside, so find the great-grandparent.
                for search in range(max_searches):
                    ancestor = find_parent(ancestor)
                    if ancestor is None:
                        break
                    elif ancestor.startswith("Explosion"):
                        found = True
                        break
                    if ancestor in self.df_reentry.index:
                        if "IN" in row["Status"]:
                            self.df_reentry.loc[ancestor, "other_mass"] += float(row["DryMass"]) # type: ignore
                        else:
                            self.df_reentry.loc[ancestor, "attached_abl_mass"] += float(row["DryMass"]) # type: ignore
                        found = True
                        break
                        
            if not found:
                failed_count += 1
                failed_mass += float(row["DryMass"])
        
        print(f'Failed to add mass for {failed_count} objects totalling {failed_mass} kg.')
         
    def add_missing_stages(self):
        
        """Add any missing stages for launches in the year.
        """   
        apogee_limit = 100  # km

        # Load in the launch list and only choose successful launches.     
        success_dsl  = self.dsl[~self.dsl["COSPAR_ID"].str[5].isin(["F", "U"])]
        success_dsl = success_dsl.rename(columns={"Time(UTC)": "Time_UTC"})

        self.df_reentry = self.df_reentry.copy()
        self.df_reentry["_cospar_prefix"] = self.df_reentry["cospar_id_new"].str[:8]
        reentry_grouped = self.df_reentry.groupby(["_cospar_prefix", "category"])
        rocket_index = self.dsr.set_index(["Rocket_Name", "Rocket_Variant"])

        # BECO, MECO, SEI1, SECO
        # B+1/2S, B+3S, B+4S, 2S, 3S, 4S 
        self.event_alts = [[66,55,29,0,0,0],
                           [220,120,64,90,56,52],
                           [229,120,64,103,61,59],
                           [356,232,216,312,176,149]]

        missing_stages_list = []
        missing_boosters_count, missing_first_count = 0,0
        missing_boosters_mass, missing_first_mass = 0,0

        # Loop over each launch, locate the rocket and then filter for rocket configuration.
        for launch_row in success_dsl.itertuples():

            # Locate the rocket.
            key = (launch_row.Rocket_Name, launch_row.Rocket_Variant)
            try:
                rocket_row = rocket_index.loc[key]          # keep as DataFrame # type: ignore
            except KeyError:
                raise ValueError(f"Rocket not found: {key}")

            rocket_config_type = rocket_row["Rocket_Config"]

            temp_dict = {
                "id"               : f'{launch_row.COSPAR_ID}',
                "jcat"             : "N/A",
                "burnup"           : "Complete",
                "category"         : "S0",
                "time"             : launch_row.Time_UTC,
                "datestr"          : launch_row.Date,
                "lat"              : launch_row.Latitude,
                "lon"              : launch_row.Longitude,
                "other_mass"       : 0,
                "attached_abl_mass": 0,
                "location"         : 2,
            }

            ################s
            # Boosters 
            ################
            boosters = int(rocket_row.Booster_No)
            if boosters > 0:
                
                # Find the event altitudes (use default if NaN) 
                beco = rocket_row.BECO
                if np.isnan(beco):
                    beco = self.event_alts[0][rocket_config_type]

                # Count how many boosters are currently added, and compare it to the number of boosters there should be.
                try:
                    booster_count = len(reentry_grouped.get_group((launch_row.COSPAR_ID, "S0")))
                except KeyError:
                    booster_count = 0

                # Adding boosters where BECO is above 50km.
                if (booster_count != boosters) and (beco > apogee_limit):
                    print(f"{boosters-booster_count} missing boosters for {launch_row.COSPAR_ID} {launch_row.Rocket_Name} {launch_row.Rocket_Variant}.")
                    missing_boosters_count += (boosters-booster_count)

                    for i in range(boosters-booster_count):
                        abl_mass = float(rocket_row.Stage0_StageMass) / boosters
                        if np.isnan(abl_mass):
                            raise ValueError(f"Missing booster mass for {launch_row.Rocket_Name} {launch_row.Rocket_Variant}.")
                        missing_boosters_mass += abl_mass   
                        missing_stages_list.append({
                            **temp_dict,
                            "name"             : f"{launch_row.Rocket_Name} {launch_row.Rocket_Variant} Booster",
                            "category"         : "S0",
                            "apogee"           : beco,
                            "abl_mass"         : abl_mass,
                        })
            
            ################
            # 1st Stage 
            ################

            # Find the event altitudes (use default if NaN) 
            meco = rocket_row.MECO
            if np.isnan(meco):
                meco = self.event_alts[1][rocket_config_type]
            
            # See if there is a matching first stage already in the dictionary.
            try:
                first_count = len(reentry_grouped.get_group((launch_row.COSPAR_ID, "S1")))
            except KeyError:
                first_count = 0

            if meco > apogee_limit and (first_count < 1):
                # Check if there is a first stage, and add if not.
                print(f"Missing 1st stage for {launch_row.COSPAR_ID} {launch_row.Rocket_Name} {launch_row.Rocket_Variant} {rocket_row.Stage1_StageMass} kg")
                abl_mass = rocket_row.Stage1_StageMass
                if np.isnan(abl_mass):
                    raise ValueError(f"Missing first stage mass for {launch_row.Rocket_Name} {launch_row.Rocket_Variant}.")
                missing_first_count += 1
                missing_first_mass += abl_mass
                missing_stages_list.append({
                    **temp_dict,
                    "name"             : f"{launch_row.Rocket_Name} {launch_row.Rocket_Variant} Stage 1",
                    "category"         : "S1",
                    "apogee"           : meco,
                    "abl_mass"         : abl_mass,
                })

        # Check the number of stages.
        stage_counts = reentry_grouped.size()
        stage_mask = (stage_counts.index.get_level_values('category').str.startswith('S') & 
                        (stage_counts.index.get_level_values('category') != 'S0'))
        stage_counts = stage_counts[stage_mask]
        duplicates   = stage_counts[stage_counts > 1]
        if len(duplicates) > 0:
            print("Duplicate stages found:")
            for (cospar, category), count in duplicates.items():
                print(f"  {cospar} has {count} × {category}")
        
        print(f"Missing Boosters:    {missing_boosters_mass},{missing_boosters_count}")
        print(f"Missing First Stage: {missing_first_mass},{missing_first_count}")
    
    def import_raul_spacex_map(self):
        """Import and set up the geolocation lists for falcon landings.
        """        
        
        # Import the data.
        fiona.drvsupport.supported_drivers['KML'] = 'rw' # type: ignore
        raul_data = gpd.read_file('./databases/reentry/General_SpaceX_Map_Raul.kml', driver='KML', layer=2)

        # Subset to only landings, excluding planned landings.
        raul_data_subset = raul_data[
            (~raul_data["Name"].str.contains("planned", regex=False)) &
            (~raul_data["Name"].str.contains("planed", regex=False))
        ].copy()
        
        # Build the ocean and ground landing lists.
        ocean_landings = raul_data_subset[
            (raul_data_subset["Name"].str.contains("ASDS", regex=False)) & 
            (raul_data_subset["Description"].str.contains("Landing -", regex=False))
        ].copy()

        ground_landings = raul_data_subset[
            (raul_data_subset["Name"].str.contains("landing", regex=False)) & 
            (raul_data_subset["Name"].str.contains("LZ", regex=False))
        ].copy()

        # Make a date column and filter to only the years of interest.
        def extract_landing_date(desc):
            date_ind = desc.index("Landing -")
            return desc[date_ind:].replace(" ", "")[8:17].replace("-","")
        
        ocean_landings["Date"] = ocean_landings["Description"].map(extract_landing_date)
        ocean_landings['Date'] = pd.to_datetime(ocean_landings['Date'],format='%d%b%y')
        ground_landings["Date"] = ground_landings["Description"].map(extract_landing_date)
        ground_landings['Date'] = pd.to_datetime(ground_landings['Date'],format='%d%b%y')

        self.ocean_landings = ocean_landings[
            ocean_landings["Date"].between(f"{self.start_year}-01-01",f"{self.final_year}-12-31",inclusive="both")
            ].reset_index(drop=True)

        self.ground_landings = ground_landings[
            ground_landings["Date"].between(f"{self.start_year}-01-01",f"{self.final_year}-12-31",inclusive="both")
            ].reset_index(drop=True)
        
        # Build the fairing recovery list.
        ocean_missions, ocean_dates = self.ocean_landings["Name"].tolist(), self.ocean_landings["Date"].tolist()
        ground_missions, ground_dates = self.ground_landings["Name"].tolist(), self.ground_landings["Date"].tolist()
        missions = ocean_missions + ground_missions
        dates    = ocean_dates    + ground_dates
        fairings = []

        if len(ocean_missions) + len(ground_missions) > 0:

            fairings_all = raul_data[raul_data["Name"].str.contains("fairing", regex=False, na=False)].copy()

            for mission, date in zip(missions, dates):
                mission_string = mission.replace("ASDS position","")
                mask = fairings_all["Name"].str.contains(mission_string, regex=False)
                if mask.any():
                    tmp = fairings_all.loc[mask].copy()
                    tmp["Date"] = date
                    fairings.append(tmp)

            self.fairings = pd.concat(fairings, ignore_index=True)

        else:
            self.fairings = pd.DataFrame(columns=raul_data.columns)

        for df in [self.ocean_landings, self.ground_landings, self.fairings]:
            df["Date"] = df["Date"].dt.strftime("%Y%m%d")
                         
    def get_reentry_info(self):

        """This is the main function of this class, and loops over all the key databases (GCAT, AC).
        It also does small final adjustments:
            - checks for duplicates
            - sets smc info
            - sets rocket stage/fairing mass from database
            - fixes fairing geolocation
            - fixes timings where launch and reentry are on same day
            - sets any missing masses
            - fixes any geolocations on grid border (180E) 
        """

        self.dsr = xr.open_dataset(f"./databases/rocket_attributes_{self.start_year}-{self.final_year}.nc", decode_times=False).to_dataframe().reset_index()
        self.dsl = xr.open_dataset(f"./databases/launch_activity_data_{self.start_year}-{self.final_year}.nc", decode_times=False).to_dataframe().reset_index()
        self.import_raul_spacex_map() #(https://t.co/RAsQ9NDmEr)

        self.cospar_to_rocket = (self.dsl.set_index("COSPAR_ID")["Rocket_Name"].to_dict())
        
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
                
                # TODO: Will have to add more military launches here.
                military_launches = ["2021-U01","2022-U01","2022-U03","2023-U01","2023-U02","2023-U03","2023-U04","2024-U05"]
                sounding_rockets  = ["1961-U02","1962-U02","1964-U01","1964-U05","1965-U03","1966-U05","1966-U07","1967-U01","2013-U01"]

                # Cut down the size of the dataframe.
                df = df[
                      (df["Primary"] == "Earth")                                     # Only objects whose primary body is Earth.
                    & (df["Piece"].notna())                                          # No NaN pieces.
                    & (df["Piece"] != "UNK")                                         # No unknown pieces.
                    & (~df["Type"].str[0].isin(["Z","D", "S"]))                      # No spurious, suborbital or debris objects.
                    & (df["Piece"].str[5:6] != "S")                                  # No suborbital launches (mainly military rockets).
                    & (~df["Piece"].isin(military_launches))                         # Skipping military tests (mostly North Korea).
                    & (~df["Piece"].isin(sounding_rockets))                          # Skipping sounding rockets (Trailblazer).
                ]

                # Convert the tags with the old notation and exclude items from launches we ignore.
                df["Converted_Tag"] = df["Piece"].apply(self.convert_launch_tag)
                df = df[~df["Converted_Tag"].str[5].isin(["A","C","E","S","Y","M","W"])]

                # Set the Apogee and Inc column to numeric, replacing errors with 0.
                df["Apogee"] = pd.to_numeric(df["Apogee"].astype(str).str.rstrip("?"), errors="coerce").fillna(0).astype(int)
                df["Inc"]    = pd.to_numeric(df["Inc"], errors="coerce").fillna(0).astype(np.float64)

                # Create a year column, with errors as NaN.
                df["DYear"]  = pd.to_numeric(df["DDate"].astype(str).str[0:4],errors="coerce")
                                                                                    
                self.jsr_data_dict[file] = df
            else:
                raise ImportError(f"Failed to fetch {file} from JSR", response.status_code)
        
        #######################################################################################
        # Make the combined lists for when we need to search all JSR data for inclination data.
        #######################################################################################

        # This cuts down a lot of time when looking for missing inclination data.
        def groupbyinc(keys):
            combined_jsr  = pd.concat([self.jsr_data_dict[k] for k in keys], ignore_index=True)
            combined_jsr  = combined_jsr[
                (combined_jsr["Inc"] != "-") & 
                (combined_jsr["Inc"].astype(np.float64) != 0) &
                (combined_jsr["Status"].isin(["O","R","DSO","DSA","AF","AS","F","S","GRP","AO","AR","R?","AL"]))
            ].copy()
            combined_jsr["Inc"] = combined_jsr["Inc"].astype(np.float64)
            combined_jsr["Prefix"] = combined_jsr["Converted_Tag"].str[:8]
            combined_jsr["Inc_eff"] = np.where(combined_jsr["Inc"] > 90,180.0 - combined_jsr["Inc"],combined_jsr["Inc"])
            combined_jsr = combined_jsr.groupby("Prefix")["Inc"].max().to_dict()
            
            return combined_jsr

        self.inc_by_prefix_all    = groupbyinc(["satcat", "rcat", "auxcat", "ecat"])
        self.inc_by_prefix_failed = groupbyinc(["ftocat"])

        ##############################
        # Add each file sequentially. 
        ##############################
        
        # There is nothing relevant for this inventory in hcocat, tmpcat and csocat.
        self.reentry_rows = []
        for file in ["satcat","auxcat","lcat","rcat","lprcat","deepcat","ecat","ftocat"]:
            self.extract_jsr_info(file)

        column_types = {
            **dict.fromkeys(["id", "jcat", "name", "category", "burnup", "datestr", "parent"], str),
            **dict.fromkeys(["time", "lat", "lon", "abl_mass", "other_mass", "attached_abl_mass", "inc", "apogee"], float),
            "location": int
        }
        
        self.df_reentry = pd.DataFrame(self.reentry_rows).astype(column_types).set_index("jcat")
        self.df_reentry["cospar_id_new"] = self.df_reentry["id"].apply(self.convert_launch_tag)
        print("cargo")
        self.add_attached()        

        # Check for duplicate reentries with the same jcat id.
        jcat_counts = self.df_reentry.index.value_counts()
        duplicate_jcat_list = jcat_counts[jcat_counts > 1].index.sort_values().tolist()
        if duplicate_jcat_list:
            print(f"Duplicates: {duplicate_jcat_list}")            
        
        # Add missing stages and check the Aerospace Corp and DISCOSweb databases.
        print("Looking for missing rocket stages.")
        self.add_missing_stages()

        print(f"Problem converting lat/lon for the following locations:")
        print(list(set(self.latlonerror)))

        ########################
        # Final adjustments.
        ########################
        
        # Create a list of all smc-related launches (from launch database). 
        smc_map = dict(zip(self.dsl["COSPAR_ID"], self.dsl["Megaconstellation_Flag"]))
        self.df_reentry["smc"] = self.df_reentry["id"].str[:8].map(smc_map)  
                    
        # Alumina emissions are calculated using ablation and alumina content data from literature.
        # These values vary based on object class(core stage / upper stage / payload) and nature of launch (reusuable or not). 
        # https://www.sciencedirect.com/science/article/pii/B0122274105008887 "Fairings are typically made of aluminum or composite materials."
        # Therefore we treat fairings as a core stage, except for Falcon 9 which are recovered intact.
        
        print("Setting ablation information.")
        for idx, row in self.df_reentry.iterrows():
            if ("fairing" in row["name"].lower()) or (row["category"][0]  == "S"):
                self.df_reentry.loc[idx, "alu_per"] = 0.7 # https://dspace.mit.edu/handle/1721.1/151443 (70% Al)
            elif row["category"] in ["C","P"]:
                self.df_reentry.loc[idx, "alu_per"] = 0.4 # https://doi.org/10.1016/j.asr.2020.10.036 (40% Al)
            else:
                print(f"Couldn't assign aluminium mass information for \n{row}")
            
            if row["burnup"] == "Complete":
                if ("fairing" in row["name"].lower()) or (row["category"] in ["S0", "S1"]):
                     self.df_reentry.loc[idx, "abl_deg"] = 0.3 # https://doi.org/10.1016/j.asr.2020.10.036 (70% survivability)

                elif row["category"] in ["C","P"]:
                    if row["smc"] == True:
                        self.df_reentry.loc[idx, "abl_deg"] = 1.0  # https://doi.org/10.1016/j.asr.2020.10.036 (0% survivability)
                    elif row["smc"] == False:
                        self.df_reentry.loc[idx, "abl_deg"] = 0.8  # https://doi.org/10.1016/j.asr.2020.10.036 (20% survivability)
                    else:
                        print(f"Couldn't assign smc information for \n{row}")

                elif row["category"] in ["S2","S3","S4","S5"]:
                    self.df_reentry.loc[idx, "abl_deg"] = 0.65     # https://doi.org/10.1016/j.asr.2020.10.036 (35% survivability)

                else:
                    self.df_reentry.loc[idx, "abl_deg"] = 0
                    print(f"Couldn't assign ablation information for complete burnup for \n{row}")

            elif row["burnup"] == "Partial":
                self.df_reentry.loc[idx, "abl_deg"] = 0
            else:
                print(f"Couldn't understand burnup information for \n{row}")

        # TODO: Add missing fairings.
        #fairing_count_list = []
        #for i in range(len(self.dsl["COSPAR_ID"])): 
        #     
        #    # Loop over all the fairings to check the number (should be two as it splits in half).
        #    # Also look for any difference in the inclination (should be zero).
        #    # Then set the geolocation to the location of the first fairing.  
        #    fairing_count = 0  
        #    fairing_inc = np.zeros((2))
        #    for count, reentry in enumerate(self.unique_reentry_list):  
        #        if (self.dsl["COSPAR_ID"].values[i][:8] == reentry["id"][:8] 
        #        and self.dsl["COSPAR_ID"].values[i][5] != "F"
        #        and "fairing" in reentry["name"].lower()):
        #            fairing_count += 1
        #            fairing_inc[fairing_count-1] = reentry["inc"]
        #            if fairing_count == 1:
        #                lat = reentry["lat"]
        #                lon = reentry["lon"]
        #            elif fairing_count == 2 and "Falcon 9" not in reentry["name"]:
        #                reentry["lat"] = lat
        #                reentry["lon"] = lon
        #    if fairing_inc[1]-fairing_inc[0] > 0:
        #        print(f"Warning: fairing inclinations differ for {self.dsl['COSPAR_ID'].values[i]} by {fairing_inc[1]-fairing_inc[0]}")      
        #    if fairing_count not in [0,2]:
        #        # NOTE: Apogee adjusted for 2022-150 fairing. One was 200km, one was 0km. Adjusted to make both 200km.
        #        sys.exit(f"Incorrect number of fairings ({fairing_count}) found for ID: {self.dsl['COSPAR_ID'].values[i]}")
        #    if fairing_count == 0 and self.dsl["COSPAR_ID"].values[i][5] != "F":
        #        pass # NOTE: If we want to manually add fairings, need to enable this to see which are missing.
        #        #print(f"No fairings found for {self.dsl['COSPAR_ID'].values[i]}")
        #    fairing_count_list.append(fairing_count)
            
        #################################################
        ## Update mass info using rocket_info databases. 
        #################################################

        print("Updating rocket mass.")

        # Precompute launch lookup: COSPAR prefix -> (Rocket_Name, Rocket_Variant)
        dsl_lookup = (
            self.dsl.assign(cospar_prefix=self.dsl["COSPAR_ID"].str[:8])
            .drop_duplicates("cospar_prefix")
            .set_index("cospar_prefix")[["Rocket_Name", "Rocket_Variant"]]
            .apply(tuple, axis=1)
            .to_dict()
        )

        # Precompute rocket lookup: (Rocket_Name, Rocket_Variant) -> row
        dsr_lookup = self.dsr.set_index(["Rocket_Name", "Rocket_Variant"]).to_dict("index")

        for idx, row in self.df_reentry.iterrows():
            
            # Is this a stage or fairing from a successful launch?
            if (row["category"].startswith("S") or "fairing" in row["name"].lower()) and row["id"][5] != "F":

                rocket_tuple = dsl_lookup.get(row["id"][:8])
                rocket_row = dsr_lookup.get(rocket_tuple)

                # Update mass info for all rocket stages.
                if row["category"] == "S0":
                    self.df_reentry.loc[idx, "abl_mass"] = rocket_row[f"Stage{row['category'][1]}_StageMass"] / int(rocket_row["Booster_No"])
                elif row["category"] in ["S1", "S2", "S3", "S4", "S5"]:
                    self.df_reentry.loc[idx, "abl_mass"] = rocket_row[f"Stage{row['category'][1]}_StageMass"]
                # TODO: Add fairing mass.
                elif "fairing" in row["name"].lower():
                   self.df_reentry.loc[idx, "abl_mass"] = rocket_row["Fairing_Mass"] / 2

        ############################
        # Clean up individual cases.
        ############################ 

        # The same item is also listed as weighing 1 kg elsewhere.
        mask_sep = self.df_reentry["name"].str.contains("sep motor cover", case=False, na=False)
        self.df_reentry.loc[mask_sep, "abl_mass"] = 1

        # Dest listed as 180E, this messes up other script so adjust to 179.9, will be in same grid square.
        self.df_reentry.loc[self.df_reentry["id"].eq("2020-086B"), "lon"] = 179.9
        
        # Object 2013-009J (PSLV upper Dual Launch Adapter (DLA-U)) has this mass on DW.
        self.df_reentry.loc[self.df_reentry["name"].eq("DLA-U"), "abl_mass"] = 100 

        # Set time to midnight whenever the reentry is occurs on a different day than the launch and no time info is available.
        # Also set to midnight if the reentry is from a launch in a previous year.
        # When the reentry is on the same day, set it to the launch time.

        # TODO: Sort out timings.
        #time_update_mass_1, time_update_count_1 = 0,0     
        #time_update_mass_2, time_update_count_2 = 0,0  
        #missing_time_count = 0 
        #
        #for idx, row in self.df_reentry.iterrows():
        #    
        #    if row["time"] == -1:  
        #        missing_time_count +=1                
        #        for count, launch_id in enumerate(self.dsl["COSPAR_ID"].values):
        #            if reentry["id"][:8] == launch_id:
        #                if reentry["datestr"] == self.dsl["Date"].values[count]:
        #                    reentry["time"] = self.dsl["Time(UTC)"].values[count]
        #                    if np.isnan(reentry["abl_mass"]) or np.isnan(reentry["other_mass"]):
        #                        print(reentry["id"],reentry["abl_mass"],reentry["other_mass"])
        #                    time_update_mass_1 += (reentry["abl_mass"] +reentry["other_mass"])
        #                    time_update_count_1 += 1
        #                else:    
        #                    reentry["time"] = 0
        #                    time_update_mass_2 += (reentry["abl_mass"] +reentry["other_mass"])
        #                    time_update_count_2 += 1
        #                    
        #    if reentry["time"] == -1:
        #        reentry["time"] = 0
        #        time_update_mass_2 += (reentry["abl_mass"] +reentry["other_mass"])
        #        time_update_count_2 += 1
        #
        #print(f"Time set to launch:   {int(time_update_mass_1)},{int(time_update_count_1)}")
        #print(f"Time set to midnight: {int(time_update_mass_2)},{int(time_update_count_2)}")
        #print(f"Time missing:         {missing_time_count}")
        
        self.print_stats()
        
    def reentry_info_to_netcdf(self):     
        """This saves the reentry information as a NetCDF file for later processing for GEOS-Chem.
        """        
        #Set up the dimensions of the netcdf file.
        dims = ('reentries')
    
        self.df_reentry["Ablatable_Mass"] = self.df_reentry["abl_mass"] + self.df_reentry["attached_abl_mass"]
            
        #Create the DataArrays.

        ds = xr.Dataset({
            'COSPAR_ID'              : xr.DataArray(self.df_reentry["id"].to_numpy(),  dims=dims, attrs=dict(long_name="COSPAR_ID")),
            'Object_Name'            : xr.DataArray(self.df_reentry["name"].to_numpy(),  dims=dims, attrs=dict(long_name="Object Name")),
            'Category'               : xr.DataArray(self.df_reentry["category"].to_numpy(),  dims=dims, attrs=dict(long_name="Category")),
            'Time (UTC)'             : xr.DataArray(self.df_reentry["time"].to_numpy(),  dims=dims, attrs=dict(long_name="Time (UTC)")),
            'Date'                   : xr.DataArray(self.df_reentry["datestr"].to_numpy(),  dims=dims, attrs=dict(long_name="Date")),
            'Latitude'               : xr.DataArray(self.df_reentry["lat"].to_numpy(),  dims=dims, attrs=dict(long_name="Latitude", units="Degrees")),
            'Longitude'              : xr.DataArray(self.df_reentry["lon"].to_numpy(),  dims=dims, attrs=dict(long_name="Longitude", units="Degrees")),
            'Ablatable_Mass'         : xr.DataArray(self.df_reentry["Ablatable_Mass"].to_numpy(),  dims=dims, attrs=dict(long_name="Ablatable Mass", units="kg")),
            'Ablation_Degree'        : xr.DataArray(self.df_reentry["abl_deg"].to_numpy(),  dims=dims, attrs=dict(long_name="Ablation Degree")),
            'Percent_Aluminium'      : xr.DataArray(self.df_reentry["alu_per"].to_numpy(),  dims=dims, attrs=dict(long_name="Percent Aluminium")),
            'Other_Mass'             : xr.DataArray(self.df_reentry["other_mass"].to_numpy(),  dims=dims, attrs=dict(long_name="Other Mass", units="kg")),
            'Megaconstellation_Flag' : xr.DataArray(self.df_reentry["smc"].to_numpy(),  dims=dims, attrs=dict(long_name="Megaconstellation_Flag")),
            'Location_Constraint'    : xr.DataArray(self.df_reentry["location"].to_numpy(),  dims=dims, attrs=dict(long_name="Location Constraint")),
            'Apogee'                 : xr.DataArray(self.df_reentry["apogee"].to_numpy(),  dims=dims, attrs=dict(long_name="Apogee", units="km")),
            'Burnup'                 : xr.DataArray(self.df_reentry["burnup"].to_numpy(),  dims=dims, attrs=dict(long_name="Burnup")),
        })
             
        #Save to file and close the DataSet  
        ds.to_netcdf(f'./databases/reentry_activity_data_{self.start_year}-{self.final_year}.nc')
        ds.close()
        
if __name__ == "__main__":  
    
    # Set up the arguments for each function.
    parser = argparse.ArgumentParser()
    parser.add_argument('-sv', "--save_reentry_info", action='store_true', help='Save reentry info.')
    parser.add_argument('-usegpd', "--use_geopandas", action='store_true', help='Load in geopandas dataframes.')
    parser.add_argument('-sy', "--start_year", default = "1957", choices=str(np.arange(1957,2026)), help='Start Year.')
    parser.add_argument('-fy', "--final_year", default = "2025", choices=str(np.arange(1957,2026)), help='Final Year.')
    args = parser.parse_args()
    use_gpd = args.use_geopandas
    
    # Sort out the year range.
    start_year = int(args.start_year)
    final_year = int(args.final_year)
    if start_year > final_year:
        final_year = start_year + 1
    print(f"Processing from year {start_year} to {final_year}.")
    
    # Compile the reentry data
    def run_full_pipeline():
        Data = build_reentry_list(start_year, final_year)
        Data.get_reentry_info()
        if args.save_reentry_info == True:
            Data.reentry_info_to_netcdf()
        
    # Missing latest May 27th 2025 starship failure.
    # Add Tranche, Ronghe, Digui, Hulianwang, Qianfan eventually when they reenter.
    
    profiler = Profile()
    profiler.enable()
    run_full_pipeline()
    profiler.disable()
    stats = Stats(profiler)
    stats.strip_dirs()
    stats.sort_stats('cumulative')
    #stats.print_stats(30)