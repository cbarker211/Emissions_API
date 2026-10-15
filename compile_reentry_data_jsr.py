import numpy as np
import argparse
import xarray as xr
import pandas as pd
from datetime import datetime
import sys
import geopandas as gpd
import requests
from io import StringIO
from cProfile import Profile
from pstats import Stats
from tqdm import tqdm
from collections import Counter
import warnings
warnings.filterwarnings("ignore", message="Geometry is in a geographic CRS")

from python_modules.distribute_emis_func import get_failed_launch_dict

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
        self.failed_launch_dict = get_failed_launch_dict()

        self.region_specs = {
            "Gujarat":      ("States","name_en","Gujarat", 3),
            "Orenburg":     ("States","name_en","Orenburg", 3),
            "Orenburg RZ":  ("States","name_en","Orenburg", 3),
            "B.Columbia":   ("States","name_en","British Columbia", 3),
            "Amur":         ("States","name_en","Amur", 3),
            "Amur, Russia": ("States","name_en","Amur", 3),
            "Utah":         ("States","name_en","Utah", 3),
            "Nova Scotia":  ("States","name_en","Nova Scotia", 3),
            "Texas":        ("States","name_en","Texas", 3),
            "Tasmania":     ("States","name_en","Tasmania", 3),
            "Oregon":       ("States","name_en","Oregon", 3),
            "Yukon":        ("States","name_en","Yukon", 3),
            "Sichuan":      ("States","name_en","Sichuan", 3),
            "I. Mongolia":  ("States","name_en","Inner Mongolia", 3),
            "Hawaii":       ("States","name_en","Hawaii", 3),
            "SE Hawaii":    ("States","name_en","Hawaii", 3),
            "Tibet":        ("States","name_en","Tibet", 3),
            "Alaska":       ("States","name_en","Alaska", 3),
            "Baja":         ("Mexico","st_name",["Baja California Norte","Baja California Sur"], 3),
            "Kazakhstan":   ("Countries","NAME","Kazakhstan", 3),
            "KRZ":          ("Countries","NAME","Kazakhstan", 3),
            "Kaz.":         ("Countries","NAME","Kazakhstan", 3),
            "Kaz":          ("Countries","NAME","Kazakhstan", 3),
            "Mexico":       ("Countries","NAME","Mexico", 3),
            "S Africa":     ("Countries","NAME","South Africa", 3),
            "N Zealand":    ("Countries","NAME","New Zealand", 3),
            "China":        ("Countries","NAME","China", 3),
            "Zhongguo":     ("Countries","NAME","China", 3),
            "Greenland":    ("Countries","NAME","Greenland", 3),
            "Saudi Ar.":    ("Countries","NAME","Saudi Arabia", 3),
            "Saudi":        ("Countries","NAME","Saudi Arabia", 3),
            "Saudi Arabia": ("Countries","NAME","Saudi Arabia", 3),
            "Zimbabwe":     ("Countries","NAME","Zimbabwe", 3),
            "Indonesia":    ("Countries","NAME","Indonesia", 3),
            "Russia":       ("Countries","NAME","Russia", 3),
            "RF":           ("Countries","NAME","Russia", 3),
            "USA":          ("Countries","NAME","United States of America", 3),
            "Morocco":      ("Countries","NAME","Morocco", 3),
            "Uruguay":      ("Countries","NAME","Uruguay", 3),
            "Bolivia":      ("Countries","NAME","Bolivia", 3),
            "Brazil":       ("Countries","NAME","Brazil", 3),
            "Mongolia":     ("Countries","NAME","Mongolia", 3),
            "Australia":    ("Countries","NAME","Australia", 3),
            "Spain":        ("Countries","NAME","Spain", 3),
            "UK":           ("Countries","NAME","United Kingdom", 3),
            "Argentina":    ("Countries","NAME","Argentina", 3),  
            "Eastern US":   ("Countries", "NAME", "United States of America", 3, (-90, 24, -66, 50)),
            "E USA":        ("Countries", "NAME", "United States of America", 3, (-90, 24, -66, 50)), 
            "New Guinea":   ("Countries", "NAME", ["Papua New Guinea", "Indonesia"], 4, (130.5, -11, 151, 0)), 
            "Spitzbergen":  ("Countries", "NAME", "Norway", 3, (10, 74, 35, 81)),
            "C. Africa":    ("Countries", "SUBREGION", "Middle Africa", 3),
            "Africa":       ("Countries", "CONTINENT", "Africa", 4),
            "Europe":       ("Countries", "CONTINENT", "Europe", 4, (-25, 34, 45, 72)),
            "Gulf":         ("Gulf of Mexico","name","Gulf of Mexico", 4),
            "Indian O":     ("Indian Ocean","name","Indian Ocean", 4),
            "Indian Ocean": ("Indian Ocean","name","Indian Ocean", 4),
            "Indian O.":    ("Indian Ocean","name","Indian Ocean", 4),
            "IOR":          ("Indian Ocean","name","Indian Ocean", 4), 
            "Bay of Bengal":("Indian Ocean","name","Bay of Bengal", 4),
            "Pacific":      ("Pacific",  None, None, 4),
            "PO":           ("Pacific",  None, None, 4),
            "POR":          ("Pacific",  None, None, 4),
            "POR SEA":      ("Pacific",  None, None, 4),
            "POR MAR":      ("Pacific",  None, None, 4),
            "POR (Sank)":   ("Pacific",  None, None, 4),
            "POR (Sea)":    ("Pacific",  None, None, 4),
            "N Pacific":    ("Pacific","name","North Pacific Ocean", 4),
            "S POR":        ("Pacific","name","South Pacific Ocean", 4),
            "S Pacific":    ("Pacific","name","South Pacific Ocean", 4),
            "S Pac.":       ("Pacific","name","South Pacific Ocean", 4),
            "Philipp.Sea":  ("Pacific","name","Philippine Sea", 4),
            "S of Tasmania":("Pacific","name","Tasman Sea", 4),
            "Atlantic":     ("Atlantic", None, None, 4),
            "AO":           ("Atlantic", None, None, 4),
            "AOR":          ("Atlantic", None, None, 4),
            "N Atlantic":   ("Atlantic","name","North Atlantic Ocean", 4),
            "S AOR":        ("Atlantic","name","South Atlantic Ocean", 4),
            "S Atl":        ("Atlantic","name","South Atlantic Ocean", 4),
            "G of Guinea":  ("Atlantic","name","Gulf of Guinea", 4),
            "Caribbean":    ("Atlantic","name","Caribbean Sea", 4),
            "Black Sea":    ("Black Sea","name","Black Sea", 4),
            "Aral Sea":     ("Aral Sea","name","Aral Sea", 4),
            "Barents":      ("Barentz Sea","name","Barentsz Sea", 4),
            "SW AOR":       ("Atlantic",     None, None, 4, (-180, -90,   0,   0)),
            "E Pacific":    ("Pacific",      None, None, 4, (-180, -90,   0,  90)),  # lon boundary is at 180W/E
            "SE Pacific":   ("Pacific",      None, None, 4, (-180, -90,   0,   0)),
            "W POR":        ("Pacific",      None, None, 4, (   0, -90, 180,  90)),
            "SW POR":       ("Pacific",      None, None, 4, (   0, -90, 180,   0)),
            "NE RU/NW POR": ("Pacific",      None, None, 4, (   0,   0, 180,  90)),
            "SE IOR":       ("Indian Ocean", None, None, 4, (  80, -90, 180,   0)),  # lon boundary is at 80E
            "SW IOR":       ("Indian Ocean", None, None, 4, (   0, -90,  80,   0)),  # lon boundary is at 80E
            "S of Aus":     ("Indian Ocean", None, None, 4, ( 113, -90, 147, -35)),
            "S of Austr.":  ("Indian Ocean", None, None, 4, ( 113, -90, 147, -35)),
            "S of S Afr.":  ("Indian Ocean", None, None, 4, (  20, -90,  35, -35)),
            "IOR/SO":       ("Indian Ocean", None, None, 4, (  20, -90, 147, -40)),
        }
        
        # Import the shapefiles.
        if use_gpd == True:

            shapefiles = {
                "Pacific":        "./databases/reentry/shapefiles/pacific/iho.shp",
                "Atlantic":       "./databases/reentry/shapefiles/atlantic/iho.shp",
                "Gulf of Mexico": "./databases/reentry/shapefiles/gulf_of_mexico/iho.shp",
                "Indian Ocean":   "./databases/reentry/shapefiles/indian_ocean/iho.shp",
                "Southern":       "./databases/reentry/shapefiles/southern/iho.shp",
                "Black Sea":      "./databases/reentry/shapefiles/black_sea/iho.shp",
                "Aral Sea":       "./databases/reentry/shapefiles/aral_sea/worldlakes.shp",
                "Barentz Sea":    "./databases/reentry/shapefiles/barentsz_sea/iho.shp",
                "Countries":      "./databases/reentry/shapefiles/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp",
                "States":         "./databases/reentry/shapefiles/ne_50m_admin_1_states_provinces/ne_50m_admin_1_states_provinces.shp",
                "Mexico":         "./databases/reentry/shapefiles/mexico_states/mexico_states.shp",
            }

            self.loaded_shapefiles = {}
            for name, path in tqdm(shapefiles.items(), total=len(shapefiles.keys()), desc="Loading Shapefiles."):
                gdf = gpd.read_file(path, engine="pyogrio")
                self.loaded_shapefiles[name] = gdf

            def build_geom(layer, col, val, _loc, bbox=None):
                gdf = self.loaded_shapefiles[layer]
                if col is not None:
                    gdf = gdf[gdf[col].isin(val if isinstance(val, list) else [val])]
                if bbox:
                    gdf = gdf.clip(bbox)
                return gdf

            self.region_geoms = {k: build_geom(*spec) for k, spec in self.region_specs.items()}

            empty = [k for k, g in self.region_geoms.items() if g.empty or g.geometry.is_empty.all()]
            if empty:
                raise ValueError(f"No geometry matched for: {empty}")
        
        self.provided_coords = {
            "GM 86.2W 29.7N": (29.7, -86.2, 1),     # Gulf of Mexico named coordinates.
            "Alashan":      (38.833, 105.95, 2),    # Landing expected east of Jiuquan launch site (this was actually launched from Wenchang).
                                                        # Assuming that Alashan refers to Helen Shan Mountain (new name for Alashan Mountain).
                                                        # https://space.skyrocket.de/doc_sdat/rcs-fc-sc.htm
                                                        # https://planet4589.org/space/jsr/news.778
                                                        # https://www.peakbagger.com/peak.aspx?pid=10693
                                                        # https://www.nasaspaceflight.com/2020/05/china-next-generation-crew-capsule/
            "Jiuquan":      (41.3, 100.3, 2),       # This refers to the Chinese launch center, sources say it landed 'in the desert near Jiuquan', and 'in China’s Inner Mongolia autonomous region'.
                                                        # Using lat/lon for Jiuquan launch site.
                                                        # https://discosweb.esoc.esa.int/launch-sites/21
                                                        # https://nssdc.gsfc.nasa.gov/nmc/spacecraft/display.action?id=2020-027A
                                                        # https://space.skyrocket.de/doc_sdat/xzf-sc.htm 
                                                        # https://spaceflightnow.com/2020/05/08/chinas-next-generation-crew-spacecraft-lands-after-unpiloted-test-flight/ 
            "Dongfeng":     (41.3, 100.3, 2),       # Same as Jiuquan.
            "LOPNOR RW05":  (40.78, 89.27, 2),      # Reported test flight of a Chinese reusable experimental spacecraft.
                                                        # Runway 05/22 at Lop Nor air base in Xinjiang (https://planet4589.org/space/gcat/data/tables/lp.html).
                                                        # https://www.seradata.com/china-launches-own-mini-spaceplane-reusable-spacecraft-using-long-march-2f 
                                                        # https://twitter.com/planet4589/status/1302486141090885632
            "Lop Nor":      (40.78, 89.27, 2),      # Same as LOPNOR. 
            "Koonibba":     (-31.886, 133.449, 2),  # Australian Test Range
            "Jacklyn":      (27.888, -74.153, 2),   # https://x.com/spaceOFFSHORE/status/1878177979974775139/photo/1
            "WSSH":         (32.95, -106.41, 2),    
            "WSMR RW17":    (32.40, -106.53, 2),    # White Sands Missile Range from sites.tsv
            "Splash":       (-36.57, 177.63, 2),    # This is the Electron failed helicopter stage 1 recovery in May 22.
                                                        # https://www.youtube.com/watch?v=BY0CXlOeWHI "Several hundred kilometers from the launch site."
                                                        # Using inclination of 94 degrees (wiki), and distance of 300km.
            "KSC RW33":     (28.61, -80.69, 2),     # Kennedy Space Center Shuttle Landing Facility  
            "KSC RW15":     (28.61, -80.69, 2),     # Kennedy Space Center Shuttle Landing Facility                                 
            "KSC SLF":      (28.61, -80.69, 2),     # https://www.world-airport-codes.com/united-states/nasa-shuttle-landing-facility-69738.html
                                                        # https://www.spaceforce.mil/News/Article/3217077/x-37b-orbital-test-vehicle-concludes-sixth-successful-mission/
            "S Pole":       (-90.0, 0.0, 2),        # South Pole (obviously...)
            "Akmolinsk SW": (50.218, 69.910, 2),    # Kosmos 10 "Landed 150 km SW of Akmolinsk." https://www.orbitalfocus.uk/Diaries/Zenit/Zenit-2.php
            "Arkalyk":      (50.249, 66.902, 2),    # Arkalyk
            "SE Arkalyk":   (49.61, 67.90, 2),      # 100 km SE of Arkalyk
            "Wisconsin":    (44.099, -87.658, 2),   # https://www.smithsonianmag.com/air-space-magazine/when-sputnik-crashed-wisconsin-180952388/
            "Uralsk W 130": (51.224, 51.373, 2),    # Uralsk
            "Kustanai":     (53.2, 63.62, 2),       # Kostanay
            "Nr. Kustanai": (53.2, 63.62, 2),       # Kostanay
            "110 km Orsk":  (52.28, 61.14, 2),      # Landing of German Mirka satellite.
                                                        # https://www.esa.int/ESA_Multimedia/Images/1998/01/Foton-11_interface_tests9, 200 km south-east of Kustenay.
                                                        # Distance between Kustenay to Orsk is ~415 km. 
                                                        # http://www.astronautix.com/m/mirka.html, landed in Kazakhstan. 
                                                        # Setting to 200 km SE of Kustenay.
            "Woomera":      (-30.93, 136.50, 2),   # Woomera Instrumented Range from sites.tsv   
            "Irkutsk":      (52.29, 104.29, 2),
            "UTTR":         (41, -113.25, 2),      # Utah Test and Training Range
            "Tashkent":     (41.31, 69.28, 2),     # Tashkent, Uzbekistan
            "GTsP-4":       (48.50, 45.80, 2),     # Kapustin Yar, coords from sites.tsv
            "Eniwetok":     (11.50 , 162.25, 2),   # Enewetak Atoll, coords from sites.tsv
            "Canary Is.":   (28.3, -15.6, 2),
            "KB":           (46.0, 63.0, 2),       # Not sure why this is KB, but this is the group for GIK-5 / NIIP-5 at Baikonur.
        }

        self.site_codes = {
            "VAFB": "V",
            "VSFB": "V",
            "EAFB": "EAFB",
            "STB": "STB",
        }

        url = "https://planet4589.org/space/gcat/tsv/tables/lp.tsv"
        response = self.session.get(url)
        if response.status_code == 200:

            # Save the file so we can search and track updates.
            filename = url.split("/")[-1]  # Use the last part of the URL as filename
            with open(f"gcat/{filename}", "w", encoding="utf-8") as f:
                f.write(response.text)

            # Convert the content to a file-like object for pandas
            tsv_data = StringIO(response.text)
            # Load the data into a pandas DataFrame
            self.lp_data = pd.read_csv(tsv_data, delimiter="\t", dtype=object, skiprows=[1])

    def convert_time(self, date):
    
        """Converts the "DDate" listing from GCAT to day, month, and time strings.

        Returns:
            date_obj(pandas datetime):  The date in YYYY-MM-DDTHH:MM:SSZ format.
        """   
        date = [s.replace("?", "") for s in date]
        yearstr=str(date[0])
        if len(date) == 1:
            yearstr = yearstr.replace("s","")
            monstr = "01"   
            daystr = "01"
        elif date[1] == "Q3": # 1995-005A
            monstr = "07"
            daystr = "01"  
        elif date[1] == "Q4": # 2015-025H
            monstr = "10"
            daystr = "01"
        else:   
            try:
                monstr = str(datetime.strptime(date[1], '%b').month).zfill(2)  
            except:
                sys.exit(date)
            if len(date) > 2: 
                daystr = str(date[2]).zfill(2)
            else:
                daystr = "01"
                
        if len(date) >= 4:    
            hrstr   = date[3][0:2].zfill(2)
            minstr  = date[3][2:4].zfill(2)
            secstr  = date[3][5:7].zfill(2) if ":" in date[3] else "00"
        else:
            hrstr   = "00"
            minstr  = "00"
            secstr  = "00"          
        datestr = f"{yearstr}-{monstr}-{daystr}T{hrstr}:{minstr}:{secstr}Z"
        date_obj = pd.to_datetime(datestr, errors="coerce", utc=True, format="%Y-%m-%dT%H:%M:%SZ")

        return date_obj       

    def print_stats(self):
        """Print out statistics about the database at the end of the program.
        """        

        total_reentries = len(self.df_reentry)
        non_geo_count = len(self.latlonerror)
        geo_percent_all =  non_geo_count  / total_reentries * 100

        total_mass_array = self.df_reentry["abl_mass"] + self.df_reentry["other_mass"]
        total_mass = total_mass_array.sum()
        non_mass_mask = total_mass_array == 0
        non_mass_count = non_mass_mask.sum()
        mass_percent_all = non_mass_count / total_reentries * 100

        data = {'Property':  ['Total Reentries', 'Total mass [Gg]','Non-geolocated Objects','Non-geolocated Objects [%]',
                              'Objects with no Mass','Objects with no Mass [%]'],
                'Value':     [int(total_reentries), total_mass*1e-6, int(non_geo_count), geo_percent_all,
                              int(non_mass_count), mass_percent_all]}
        
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
        
    def convert_lat_lon(self,latlonstr,inc,category,apogee,jsr_id, jsr_name):
        """Set the geolocation for all items that are not Falcon 9 fairings/1st stage.

        Args:
            latlonstr (str)     : The dest field from GCAT.
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
                                3 - Political Region (Country, State, Province)
                                4 - Natural Region (Ocean, Landmass, Continent)
                                5 - Falcon Reusable
                                6 - Inclination Bounded Random
        """        
        
        lat, lon, location = None, None, -1

        def sample_within_region(region, lat_range=(-inc, inc), lon_range=(-180, 180)):

            band = region.clip((lon_range[0], lat_range[0], lon_range[1], lat_range[1]))
            if band.empty or band.geometry.is_empty.all():
                print(f"Inc band misses region: {jsr_id} {jsr_name} {latlonstr} inc={inc}")
                band = region

            row = band.sample(1, weights=band.area)
            coords = row.sample_points(1).get_coordinates()
            lon, lat = coords["x"].values[0], coords["y"].values[0]
            return lat, lon
        
        if latlonstr in self.region_specs:
            location = self.region_specs[latlonstr][3] # Get the location tag.
            if use_gpd:
                lat, lon = sample_within_region(self.region_geoms[latlonstr], lat_range=(-inc, inc))

        # Missing location - 6
        elif latlonstr in ["-","","LEO", "CLO"]:
            
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

        # Natural Region - 4
        elif latlonstr in ["S Ocean", "S Oc.", "S OCean", "Antarctic", "Arctic", "Ocean"]:
            location = 4

            # Polar
            if latlonstr in ["Antarctic", "Arctic", "S Ocean","S Oc.","S OCean"]:
                lon = round(np.random.uniform(-180, 180), 2)
            
                if latlonstr in ["S Ocean","S Oc.","S OCean"]:
                    # For these objects in auxcat, the inclination is 51.65-74, outside of the general definition of the southern ocean (<60S).
                    # Therefore we will just set the lat to the inclination.
                    lat = -inc
                elif latlonstr == "Antarctic":
                    lat = round(np.random.uniform(-inc, -60), 2)
                elif latlonstr == "Arctic":
                    lat = round(np.random.uniform(66, inc), 2)

            # Ocean landings (Electron and Falcon 9)
            elif latlonstr == "Ocean":

                if jsr_name == "Electron Stage 1":
                    # Electron Stage 1 (rcat) is occasionally recovered in the ocean.
                    # From https://spaceflightnow.com/2020/11/05/rocket-lab-to-attempt-booster-recovery-on-next-mission/:
                    # Recovery vessels stationed near the booster’s splashdown zone around 250 miles (400 kilometers) 
                    # south of the launch site will move in to secure the first stage and hoist it onto a ship for return to New Zealand.
                    # Geolocating 400km south of launch site.
                    lat, lon = -42.86, 177.87

                elif jsr_name == "Falcon 9 Stage 1 V1.1":
                    # Setting to the same coords as used in the other falcon functions.
                    site_idx = np.where(self.dsl["COSPAR_ID"] == jsr_id)[0]
                    site_name = self.dsl["Site"].values[site_idx[0]]
                    # Kennedy Space Center (ETR)
                    if site_name == "Cape Canaveral":
                        lat = 34
                        lon = -75  
                    # Vandenberg
                    elif site_name == "Vandenberg Space Force Base":
                        lat = 30
                        lon = -120

        # Lat/Lon Provided, Launch Site, or Named Location - 1-2
        elif latlonstr in self.provided_coords:
            lat, lon, location = self.provided_coords[latlonstr]

        elif latlonstr.split()[0] in self.site_codes:
            location = 2
            match = self.lp_data[
                (self.lp_data["#Site"] == self.site_codes[latlonstr.split()[0]]) &
                (self.lp_data["Code"] == latlonstr.split()[1])
            ]
            lat = match["Latitude"].item()
            lon = match["Longitude"].item()

        elif any(char.isdigit() for char in latlonstr):

            def degree_coords(latlonstr):
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

                return lat, lon, location

            if len(latlonstr.split()) == 3 and len(latlonstr.split()[-1]) == 1:
                latlonstr = latlonstr.split()[0] + " " + latlonstr.split()[1]
                lat, lon, location = degree_coords(latlonstr)
            
            elif latlonstr == "IOR 90E 10S":
                lat, lon, location = degree_coords("90E 10S")

            elif latlonstr.split()[0][-1] in ["E","W","S","N"]:
                lat, lon, location = degree_coords(latlonstr)

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
                        lon = -(np.float64(lon_data[0]) + np.float64(lon_data[1].replace("W",""))/60)
                    elif "E" in lon_data[-1]:
                        lon =   np.float64(lon_data[0]) + np.float64(lon_data[1].replace("E",""))/60
                        
                    if "N" in lat_data[-1]:
                        lat =   np.float64(lat_data[0]) + np.float64(lat_data[1].replace("N",""))/60
                    elif "S" in lat_data[-1]:
                        lat = -(np.float64(lat_data[0]) + np.float64(lat_data[1].replace("S",""))/60)

        if lat == 0.0 and lon == 0.0:
            print("Re-entry location is 0.0")
        if lat is None or lon is None or location == -1:
            self.latlonerror.append(latlonstr)
            lat = 0.0
            lon = 0.0
        else:
            lat = np.float64(lat)
            lon = np.float64(lon)

        return lat, lon, location
    
    def falcon_stage_lat_lon(self,date_obj,jsr_id,dest,plname):
        
        """Set the geolocation for all Falcon 9 1st stages.

        Returns:
            lat(np.float64) : The latitude of the object.
            lon(np.float64) : The longitude of the object..
        """ 

        lat, lon = None, None

        # Find the ground landings in Raul's SpaceX Map
        if "LZ" in dest:
            matching = self.ground_landings[self.ground_landings["Date"].dt.date == date_obj.date()].reset_index(drop=True)

            # If its in the list of ground landings, then we can geolocate it.
            if matching.shape[0] == 1: 
                lat = matching["Latitude"].iloc[0]
                lon = matching["Longitude"].iloc[0]   

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
                    lat = matching["Latitude"].iloc[0]
                    lon = matching["Longitude"].iloc[0]
                # Otherwise, raise an error.
                else:
                    raise ValueError(f"Multiple ground entries for Falcon Stage 1 landing- {jsr_id}.")
            else:
                raise ValueError(f"Problem geolocating Falcon Stage 1 landing- {jsr_id}.")
        else:
            # Find the ocean landings in Raul's SpaceX Map.
            matching = self.ocean_landings[self.ocean_landings["Date"].dt.date == date_obj.date()].reset_index(drop=True)
            if (matching.shape[0] == 1) or (date_obj.tz_convert(None).normalize() == pd.Timestamp("2022-04-27")):
                # This is 27th April - typo in Raul's Space Map where second should be 2023. 
                lat = matching["Latitude"].iloc[0]
                lon = matching["Longitude"].iloc[0]
            # These next two are where there are two launches on the same day.
            elif jsr_id == "2022-124":
                lat = matching["Latitude"].iloc[0]
                lon = matching["Longitude"].iloc[0]
            elif jsr_id == "2022-125":
                lat = matching["Latitude"].iloc[1]
                lon = matching["Longitude"].iloc[1]
            # Again these are two on the same day.
            elif jsr_id == "2023-037":
                lat = matching["Latitude"].iloc[0]
                lon = matching["Longitude"].iloc[0]
            elif jsr_id == "2023-038":
                lat = matching["Latitude"].iloc[1]
                lon = matching["Longitude"].iloc[1]

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
                    lat = 32
                    lon = -75  
                # Vandenberg
                elif (lon_val == -121) and (lat_val == 35):
                    lat = 30
                    lon = -122.5 
                else:
                    raise ValueError(f'Falcon launch site not found - {lon_val} {lat_val}.')
            elif matching.shape[0] > 1:
                print(matching,date_obj)
                raise ValueError(f"Multiple ocean entries for Falcon Stage 1 landing- {jsr_id}.")
            else:
                raise ValueError(f"Problem geolocating Falcon Stage 1 landing- {jsr_id}.") 
            
        return lat, lon 
    
    def falcon_fairing_lat_lon(self,date_obj,jsr_id):   
        
        """Set the geolocation for all Falcon 9 fairings.

        Returns:
            lat(np.float64) : The latitude of the object.
            lon(np.float64) : The longitude of the object..
        """

        lat, lon = None, None

        set_ocean = False
        if len(self.fairings) > 0:
            matching = self.fairings[self.fairings["Date"].dt.date == date_obj.date()].reset_index(drop=True)
            if matching.shape[0] == 1:
                lat = matching["Latitude"].iloc[0]
                lon = matching["Longitude"].iloc[0]
            # These are two on the same day
            elif jsr_id == "2023-037":
                lat = matching["Latitude"].iloc[0]
                lon = matching["Longitude"].iloc[0]
            elif jsr_id == "2023-038":
                lat = matching["Latitude"].iloc[1]
                lon = matching["Longitude"].iloc[1]
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
    
    def failed_launch_mass(self, jsr_id, reentry_category, drymass, mass):
        
        """Split failed-launch reentry mass into ablatable (dry) and non-aluminium (unburnt propellant) mass.

        Returns:
            abl_mass (float)  : Ablatable mass (GCAT dry mass; stage/fairing masses are overwritten later from rocket_attributes).
            other_mass (float): Unburnt propellant mass.
        """        
            
        abl_mass   = np.float64(drymass)
        other_mass = 0.0

        # Payloads that fail to reach orbit re-enter with their propellant.
        if reentry_category == "P":
            other_mass = max(np.float64(mass) - abl_mass, 0.0)

        failure = self.failed_launch_dict.get(jsr_id[:8])
        if failure is not None and reentry_category in ["S1","S2","S3","S4","S5"]:
            
            stage_no = int(reentry_category[1])
            launch = self.dsl.loc[self.dsl["COSPAR_ID"] == jsr_id[:8]].iloc[0]
            match = (self.dsr["Rocket_Name"] == launch["Rocket_Name"]) & (self.dsr["Rocket_Variant"] == launch["Rocket_Variant"])
            rocket_ind = np.where(match)[0][0]
            prop_mass = self.dsr[f"Stage{stage_no}_PropMass"].values[rocket_ind]

            if stage_no > failure["failed_stage"]:
                other_mass = prop_mass
            elif stage_no == failure["failed_stage"]:
                other_mass = round((1 - failure["burned"]) * prop_mass, 1)

        # TODO: Email Jonathan about these payload masses.
        # ID       Name            DryMass  Mass    Lit    Source
        # 2020-F01 Zafar           90.0     90.0    113.0  https://space.skyrocket.de/doc_sdat/zafar.htm https://spaceflightnow.com/2020/02/10/iran-fails-in-satellite-launch-attempt/ https://www.rferl.org/a/iran-says-satellite-launch-imminent-in-next-few-hours-/30424791.html
        # 2020-F05 CE-SAT-I        65.0     65.0    67.0   https://spaceflightnow.com/2020/07/04/rocket-lab-satellite-launch-fails-before-reaching-orbit/ https://spacenews.com/rocket-lab-electron-launch-fails/
        # 2020-F05 Flock 4e-1      5.7      5.7     5.8    https://pubs.usgs.gov/of/2021/1030/f/ofr20211030f.pdf
        # 2020-F05 Flock 4e-2      5.7      5.7     5.8    https://pubs.usgs.gov/of/2021/1030/f/ofr20211030f.pdf
        # 2020-F05 Flock 4e-3      5.7      5.7     5.8    https://pubs.usgs.gov/of/2021/1030/f/ofr20211030f.pdf
        # 2020-F05 Flock 4e-4      5.7      5.7     5.8    https://pubs.usgs.gov/of/2021/1030/f/ofr20211030f.pdf
        # 2020-F05 Flock 4e-5      5.7      5.7     5.8    https://pubs.usgs.gov/of/2021/1030/f/ofr20211030f.pdf
        # 2022-F01 INCA            4.0      4.0     3.8    https://www.nasaspaceflight.com/2022/02/astra-cape-canaveral-elana-41/ https://space.skyrocket.de/doc_sdat/inca.htm
        # 2022-F03 TROPICS SV02    4.0      4.0     5.34   https://web.archive.org/web/20250719051151/https://tropics.ll.mit.edu/CMS/tropics/The-MicroMAS-2-Cubesat
        # 2022-F03 TROPICS SV04    4.0      4.0     5.34   https://web.archive.org/web/20250719051151/https://tropics.ll.mit.edu/CMS/tropics/The-MicroMAS-2-Cubesat
        # 2022-F06 Pleiades Neo 5  895.0    922.0   920.0  https://earth.esa.int/eogateway/missions/pleiades-neo/description
        # 2022-F06 Pleiades Neo 6  895.0    922.0   920.0  https://earth.esa.int/eogateway/missions/pleiades-neo/description 
        # 2022-U02 EOS-02          130.0    135.0   145.0  https://spaceflightnow.com/2022/08/07/indias-new-small-satellite-launcher-fails-to-put-satellites-into-correct-orbit/ https://web.archive.org/web/20220803133923/https://www.isro.gov.in/sites/default/files/sslvd1-eos2_-_mission_brochure.pdf https://space.skyrocket.de/doc_lau/sslv.htm

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
            inc_eff = self.inc_by_prefix.get(jsr_id[:8], 0)

        if inc_eff == 0:
            print(f"Warning: Inclination still empty for {jsr_id}.")
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
            (jsr_data["Status"].isin(["R","R?","D","L","L?","LF","S","F","AF","AS"])) &             # Only include specific statuses (see https://planet4589.org/space/gcat/web/intro/phases.html)       
            (jsr_data["DYear"].between(self.start_year, self.final_year)) &                         # This time range only.
            ((jsr_data["Apogee"] >= 50) | (jsr_data["Apogee"] < 0)) &                               # Apogee above 50 km, or below 0. This makes sure to include objects in lprcat returning from space.
            ((jsr_data["Status"] != "AS") | (jsr_data["Converted_Tag"].str[5:6].isin(["F","U"]))) & # Only include AS re-entries if they are part of a failed or uncategorized launch.
            (jsr_data["Dest"] != "Z") &                                                             # Skip duplicate objects.
            (jsr_data["Dest"].str.lower() != "lost")                                                # Jonathan said this means the same as the OX status tag.
        ]

        #s1 = len(subset)
        subset = subset[(subset["Dest"] != "Exploded")] # Skip objects that exploded in orbit, as we can't accurately estimate the emissions.
        #s2 = len(subset)
        #print(f"{s1-s2} exploded objects removed.")

        # Ignore failed launches pre-2020.
        mask = ~(subset["Converted_Tag"].str[5].isin(["F", "U"]) & (subset["Converted_Tag"].str[:4].astype(int) < 2020))          
        subset = subset[mask]

        # Extract the arrays before looping.
        cols = ["DDate","Converted_Tag","Name","Type","Inc","Dest","Apogee","Status","DryMass","#JCAT","Bus","Parent","PLName","Mass"]
        arrays = {col: subset[col].to_numpy() for col in cols}

        # Now loop over the list and format into a dictionary.
        for i in tqdm(range(len(arrays["Converted_Tag"])), total=len(subset), desc=f"Processing {file}"):

            # Sort out the reentry time/date.
            ddate_parts = arrays["DDate"][i].split()
            date_obj = self.convert_time(ddate_parts)
            time_known = len(ddate_parts) >= 4   # has HHMM
            day_known  = len(ddate_parts) >= 3   # has DD

            # Create some variables for easier access.
            jsr_id      = arrays["Converted_Tag"][i]
            jsr_name    = str(arrays["Name"][i])
            jsr_type    = str(arrays["Type"][i]).strip()
            jsr_inc     = arrays["Inc"][i]
            jsr_dest    = arrays["Dest"][i]
            jsr_apogee  = arrays["Apogee"][i]
            jsr_status  = arrays["Status"][i]
            jsr_drymass = arrays["DryMass"][i]
            jsr_mass    = arrays["Mass"][i]
            jsr_jcat    = arrays["#JCAT"][i]
            jsr_parent  = arrays["Parent"][i]
            jsr_bus     = arrays["Bus"][i]
            jsr_plname  = arrays["PLName"][i]

            # If the rocket failed at the first stage, it doesn't reach re-entry altitude.
            failure = self.failed_launch_dict.get(jsr_id[:8])
            if failure is not None and failure["failed_stage"] == 1:
                continue
            
            burnup = "Complete"
            # Set the burnup as partial for all objects landing or splashing down at the surface.
            if jsr_status in ["L","L?","LF"]:
                burnup = "Partial"
            elif (jsr_name == "Electron Stage 1") and jsr_id in ["2019-084","2020-007","2020-085","2021-F02","2021-106",
                                                                 "2022-047","2022-147","2023-041","2023-100","2023-126","2024-022"]:
                burnup = "Partial"
            elif (jsr_name == "Starship Super Heavy") and jsr_id in ["2024-U03","2024-U06","2025-U03","2025-U04"]:
                # Controlled soft splashdowns in the Gulf.
                burnup = "Partial"
            
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
                if rocket_name.startswith(("Space Shuttle", "SLS Block 1", "Conestoga 1620", "Polyot", "Molniya", "Sputnik", "Voskhod", "Vostok")):
                    stage = int(jsr_type[1:2]) - 1
                    reentry_category = f"S{stage}"

                if rocket_name.startswith("Space Shuttle") and reentry_category == "S0":
                    burnup = "Partial"

                # Soyuz
                if "Blok-BVGD" in jsr_bus:
                    reentry_category = "S0"
                elif "Blok-A" in jsr_bus:
                    reentry_category = "S1"
                elif ("Blok-I" in jsr_name) or ("Yuanzheng-2" in jsr_name) or (jsr_name.startswith("Soyuz-U") and jsr_type[1:2] == "3"):
                    reentry_category = "S2"
                elif any(name in jsr_name for name in ["Blok-L", "Fregat", "Volga", "Blok-ML", "Blok SO-L", "Blok-MVL", "Blok-2BL"]):
                    reentry_category = "S3"

                if len(reentry_category) == 1:
                    raise ValueError(f"Unknown reentry category for {jsr_id} - {jsr_type}")
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
                    lat, lon = self.falcon_stage_lat_lon(date_obj,jsr_id,jsr_dest, jsr_plname) 
                    location = 5

                # If the first stage was expended, then use the normal lat/lon conversion.
                elif jsr_status in ["S","R"]:
                    lat, lon, location = self.convert_lat_lon(jsr_dest.replace("?",""), inc_eff, reentry_category, jsr_apogee, jsr_id, jsr_name)
                else:
                    raise ValueError(f"Unexpected Falcon Stage 1 Status {jsr_status} {jsr_jcat}")
            
            elif ("Falcon 9 Fairing" in jsr_name) or ("Falcon Heavy Fairing" in jsr_name):    
                lat, lon = self.falcon_fairing_lat_lon(date_obj,jsr_id)
                burnup = "Partial"
                location = 5

            else:   
                lat, lon, location = self.convert_lat_lon(jsr_dest.replace("?",""), inc_eff, reentry_category, jsr_apogee, jsr_id, jsr_name)
            
            # Sort the mass.    
            if jsr_id[5:6] in ["F","U"]:
                abl_mass, other_mass = self.failed_launch_mass(jsr_id, reentry_category,jsr_drymass,jsr_mass)
            else:
                abl_mass = np.float64(jsr_drymass)
                other_mass = 0
            
            # Set up the dictionary.                      
            self.reentry_rows.append({
                "id"               : jsr_id,
                "jcat"             : jsr_jcat,
                "name"             : jsr_name,
                "category"         : reentry_category,
                "burnup"           : burnup,
                "date"             : date_obj,
                "lat"              : lat,
                "lon"              : lon,
                "abl_mass"         : abl_mass,
                "other_mass"       : other_mass,
                "attached_abl_mass": 0,
                "location"         : location,
                "inc"              : jsr_inc,
                "apogee"           : jsr_apogee,
                "parent"           : jsr_parent,
                "time_known"       : time_known,
                "day_known"        : day_known,
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
            (~attached_cat["Converted_Tag"].str[5:6].isin(["S","F","U"])) &      # Skip suborbital and failed launches (treated separately).
            (attached_cat["DYear"].between(self.start_year, self.final_year)) &  # Only the selected timeframe.
            (attached_cat["PLName"] != "Manfred Mem. Moon Mission") &            # This should be listed as re-entering on the Moon, so ignore.
            (attached_cat["Parent"] != "A04313") &                               # Parent collided, producing debris.
            (attached_cat["Parent"] != "S29402") &                               # Parent exploded.
            (attached_cat["Parent"] != "S27542") &                               # Parent went in deep space.
            (attached_cat["#JCAT"] != "A11023")                                  # AL IN is wrong.
        ].copy()

        attached_cat["parent_clean"] = attached_cat["Parent"].str.replace("*","", regex=False)            

        # A function to find the parent of any jcat.
        def find_parent(jcat):
            parent = None
            match = attached_cat[attached_cat['#JCAT'] == jcat]
            if len(match) == 1:
                if match["Status"].values[0] == "E":
                    parent = f"Explosion {match['#JCAT'].values[0]} {match['Name'].values[0].strip()}"
                else:
                    parent = match["Parent"].values[0].replace("*", "").split()[0]
            elif len(match) > 1:
                raise ValueError(f"Multiple matches for parent object {jcat}.")
            
            return parent

        # Loop over the objects and add the mass to the ancestor.  
        failed_mass, failed_count = 0, 0          
        for idx, row in attached_cat.iterrows():

            if row["Status"] == "E":
                continue

            parent = row["parent_clean"].split()[0]
            child_parent_override = {
                "A04252": "S16481",  # Landed inside Columbia (STS 61-C)
                "A04253": "S16481",  # Landed inside Columbia (STS 61-C)
                "A06730": "S33441",  # Landed inside Endeavour (STS-126)
                "A07376": "S32960",  # Landed inside Discovery (STS-124)
                "A07705": "S36507",  # Landed inside Discovery (STS-131)
            }
            parent = child_parent_override.get(row["#JCAT"], parent)
            found = False

            def add_mass(target, row):
                # If the payload is inside the spacecraft or in the payload bay, then add to other mass.
                in_payload_bay = "(STS" in self.df_reentry.loc[target, "name"]
                if "IN" in row["Status"] or in_payload_bay:
                    self.df_reentry.loc[target, "other_mass"] += float(row["DryMass"])

                # Otherwise add to ablatable mass.
                else:
                    self.df_reentry.loc[target, "attached_abl_mass"] += float(row["DryMass"])

            # If parent is already in the full reentry list.
            if parent in self.df_reentry.index:
                add_mass(parent, row)
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
                        add_mass(ancestor, row)
                        found = True
                        break
                        
            if not found:
                print(row["#JCAT"], row["Name"], row["Status"], row["Parent"], row["DryMass"], "Mass not added.")
                failed_count += 1
                failed_mass += float(row["DryMass"])
        
        print(f'Failed to add mass for {failed_count} objects totalling {failed_mass} kg.')
         
    def add_missing_stages(self):
        
        """Add any missing stages for launches in the year.
        """   
        apogee_limit = 100  # km

        # Load in the launch list and only choose successful launches.     
        success_dsl  = self.dsl[~self.dsl["COSPAR_ID"].str[5].isin(["F", "U"])]

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
                "date"             : launch_row.Date,
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
        raul_data = pd.read_parquet('./databases/reentry/raul_map.parquet')
        raul_data['Name'] = raul_data['Name'].str.rstrip('\n')

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

        self.dsr = xr.open_dataset(f"./databases/rocket_attributes.nc", decode_times=False).to_dataframe().reset_index()
        self.dsl = xr.open_dataset(f"./databases/launch_activity_data.nc").to_dataframe().reset_index()
        self.import_raul_spacex_map() #(https://t.co/RAsQ9NDmEr)

        self.cospar_to_rocket = (self.dsl.set_index("COSPAR_ID")["Rocket_Name"].to_dict())
        
        print("Loading JSR re-entries.") # (see https://planet4589.org/space/gcat/web/cat/cats.html)
        files = ["satcat","auxcat","lcat","rcat","lprcat","deepcat","ecat","ftocat"]
        # satcat (main database), auxcat (should be in main but isn't), lcat (suborbital stages/objects). 
        # rcat (lower stages and fairings), lprcat (several objects returning from space).
        # deepcat (several objects returning from space), ecat (capsules from ISS / crewed missions). 
        # ftocat (failed launches).
        # There is nothing relevant for this inventory in hcocat, tmpcat and csocat.
        self.jsr_data_dict = {}
        for file in files:
            url = "https://planet4589.org/space/gcat/tsv/cat/" + file + ".tsv"
            response = self.session.get(url)
            if response.status_code == 200:

                # Save the file so we can search and track updates.
                filename = url.split("/")[-1]  # Use the last part of the URL as filename
                with open(f"gcat/{filename}", "w", encoding="utf-8") as f:
                    f.write(response.text)

                # Convert the content to a file-like object for pandas
                tsv_data = StringIO(response.text)
                # Load the data into a pandas DataFrame
                df = pd.read_csv(tsv_data, delimiter="\t", dtype=object, skiprows=[1])

                # Strip whitespace from the '#jcat' column to ensure consistent formatting
                df["#JCAT"] = df["#JCAT"].str.strip()

                # Convert the tags with the old notation and exclude items from launches we ignore.
                df = df[df["Piece"] != "UNK"] # Remove UNK before converting.
                df["Converted_Tag"] = df["Piece"].apply(self.convert_launch_tag)
                df = df[df["Converted_Tag"].str[:8].isin(self.dsl["COSPAR_ID"])]

                # Cut down the size of the dataframe.
                named_debris = df["Name"].str.lower().str.contains("debris") | df["Name"].str.contains("deb")
                remove_debris = df["Name"].str.startswith("DebrisSat") | (df["Name"] == "RemoveDebris")

                df = df[(df["Primary"] == "Earth")                       # Only objects whose primary body is Earth.
                      & (~df["Type"].str[0].isin(["Z","D"]))             # No spurious or debris objects.
                      & (~named_debris | remove_debris)                  # Skipping debris objects.
                ]

                # Sort out the Apogee, Inc, and DYear columns. 
                df["Apogee"] = pd.to_numeric(df["Apogee"].astype(str).str.rstrip("?"), errors="coerce").fillna(0).astype(int)
                df["Inc"]    = pd.to_numeric(df["Inc"], errors="coerce").fillna(0).astype(np.float64)
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
            combined_jsr = combined_jsr.groupby("Prefix")["Inc_eff"].max().to_dict()
            
            return combined_jsr

        self.inc_by_prefix = groupbyinc(["satcat", "rcat", "auxcat", "ecat","ftocat"])

        ##############################
        # Add each file sequentially. 
        ##############################
        
        self.reentry_rows = []
        for file in files:
            self.extract_jsr_info(file)

        column_types = {
            **dict.fromkeys(["id", "jcat", "name", "category", "burnup", "parent"], str),
            **dict.fromkeys(["lat", "lon", "abl_mass", "other_mass", "attached_abl_mass", "inc", "apogee"], float),
            "location": int,
            "time_known": bool,
            "day_known": bool,
        }
        
        self.df_reentry = pd.DataFrame(self.reentry_rows).astype(column_types).set_index("jcat")
        self.df_reentry["date"] = pd.to_datetime(self.df_reentry["date"], utc=True).dt.tz_convert(None)
        self.df_reentry["cospar_id_new"] = self.df_reentry["id"].apply(self.convert_launch_tag)
        print("cargo")
        self.add_attached()        

        # Check for duplicate reentries with the same jcat id.
        jcat_counts = self.df_reentry.index.value_counts()
        duplicate_jcat_list = jcat_counts[jcat_counts > 1].index.sort_values().tolist()
        if duplicate_jcat_list:
            raise ValueError(f"Duplicates: {duplicate_jcat_list}")            
        
        # Add missing stages.
        #print("Looking for missing rocket stages.")
        #self.add_missing_stages()

        print(f"Need to convert lat/lon for the following {len(self.latlonerror)} locations:")
        counts = Counter(self.latlonerror)
        for loc, n in counts.most_common():
            print(f"{n:5d}  {loc}")

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

        NOT_LV_FAIRING = ["s-12 experiment fairing", # Gemini experiment fairing
                    "nose fairing",            # Gemini spacecraft fairing
                    "horizon sensor fairing",  # Gemini spacecraft fairing
                    "kosmos-316",              # Spacecraft object, not part of the 11K69 fairing.
                    "kosmos-1867"              # Spacecraft object, not part of the 11K69 fairing.
                    ]

        def is_lv_fairing(name):
            name = name.lower()
            return ("fairing" in name or "plf" in name) and not any(s in name for s in NOT_LV_FAIRING)

        # Precompute reusable masks
        is_fairing      = self.df_reentry["name"].map(is_lv_fairing)
        is_payload_comp = self.df_reentry["category"].isin(["C", "P", "S"])
        is_lower_stage  = self.df_reentry["category"].isin(["S0", "S1"])
        is_upper_stage  = self.df_reentry["category"].isin(["S2", "S3", "S4", "S5"])
        is_complete     = self.df_reentry["burnup"] == "Complete"
        is_partial      = self.df_reentry["burnup"] == "Partial"
        smc = self.df_reentry["smc"]

        # Initialize columns if needed
        self.df_reentry["alu_per"] = np.nan
        self.df_reentry["abl_deg"] = np.nan

        # Aluminium percentage
        self.df_reentry.loc[is_payload_comp, "alu_per"] = 0.4 # https://doi.org/10.1016/j.asr.2020.10.036 (40% Al)
        self.df_reentry.loc[is_fairing | is_lower_stage | is_upper_stage, "alu_per"]  = 0.7 # https://dspace.mit.edu/handle/1721.1/151443 (70% Al)
        
        # Ablation degree
        self.df_reentry.loc[is_complete & is_payload_comp & smc.eq(True), "abl_deg"]  = 1.0  # https://doi.org/10.1016/j.asr.2020.10.036 (0% survivability)
        self.df_reentry.loc[is_complete & is_payload_comp & smc.eq(False), "abl_deg"] = 0.8  # https://doi.org/10.1016/j.asr.2020.10.036 (20% survivability)

        # This has to go after the previous two so its not reset for fairings.
        self.df_reentry.loc[is_complete & (is_fairing | is_lower_stage), "abl_deg"]   = 0.3  # https://doi.org/10.1016/j.asr.2020.10.036 (70% survivability)
        self.df_reentry.loc[is_complete & is_upper_stage, "abl_deg"]                  = 0.65 # https://doi.org/10.1016/j.asr.2020.10.036 (35% survivability)
        self.df_reentry.loc[is_partial, "abl_deg"]                                    = 0.0

        # Optional: report rows that were not classified
        unknown_alu = self.df_reentry["alu_per"].isna()
        unknown_abl = self.df_reentry["abl_deg"].isna()

        if unknown_alu.any():
            print("Couldn't assign aluminium mass information for:")
            print(self.df_reentry.loc[unknown_alu])

        if unknown_abl.any():
            print("Couldn't assign ablation information for:")
            print(self.df_reentry.loc[unknown_abl])
            
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

        # Number of fairing pieces per launch, counted from the full GCAT catalogues.
        all_cat   = pd.concat(self.jsr_data_dict.values(), ignore_index=True).drop_duplicates("#JCAT")
        all_cat = all_cat[all_cat["Dest"] != "Z"]
        fair_cat  = all_cat[all_cat["Name"].map(is_lv_fairing)].copy()
        fair_cat["prefix"]  = fair_cat["Converted_Tag"].str[:8]
        fair_cat["DryMass"] = pd.to_numeric(fair_cat["DryMass"], errors="coerce").fillna(0)
        n_fairing = fair_cat.groupby("prefix").size().to_dict()
        n_fairing_no_gcat = fair_cat[fair_cat["DryMass"] == 0].groupby("prefix").size().to_dict()

        for idx, row in self.df_reentry.iterrows():
            
            # Is this a stage or fairing from a successful launch?
            if (row["category"].startswith("S") or is_lv_fairing(row["name"])):

                rocket_tuple = dsl_lookup.get(row["id"][:8])
                rocket_row = dsr_lookup.get(rocket_tuple)

                # Update mass info for all rocket stages.
                if row["category"] == "S0":
                    self.df_reentry.loc[idx, "abl_mass"] = rocket_row[f"Stage{row['category'][1]}_StageMass"] / int(rocket_row["Booster_No"])
                elif row["category"] in ["S1", "S2", "S3", "S4", "S5"]:
                    self.df_reentry.loc[idx, "abl_mass"] = rocket_row[f"Stage{row['category'][1]}_StageMass"]

                # For fairings, use values set in update_rocket_launch_data, then use GCAT, then use default.    
                elif is_lv_fairing(row["name"]):
                    prefix = row["id"][:8]
                    if rocket_row["Fairing_Mass"] != 1756:
                        self.df_reentry.loc[idx, "abl_mass"] = rocket_row["Fairing_Mass"] / n_fairing.get(prefix, 2)
                    elif row["abl_mass"] == 0:
                        self.df_reentry.loc[idx, "abl_mass"] = 1756 / n_fairing_no_gcat.get(prefix, 2)

        ############################
        # Clean up individual cases.
        ############################ 

        # The same item is also listed as weighing 1 kg elsewhere.
        mask_sep = self.df_reentry["name"].str.contains("sep motor cover", case=False, na=False)
        self.df_reentry.loc[mask_sep, "abl_mass"] = 1
        
        # Object 2013-009J (PSLV upper Dual Launch Adapter (DLA-U)) has this mass on DW.
        self.df_reentry.loc[self.df_reentry["name"].eq("DLA-U"), "abl_mass"] = 100 

        # Fixing skylab fairing section with no mass.
        self.df_reentry.loc[self.df_reentry["name"].eq("Skylab Fairing section"), "abl_mass"] = 2950.0

        # 1965-047A is a lead ballast dummy payload inside a fairing (GCAT PLName "Lead Ballast"), not aluminium.
        lead = self.df_reentry["id"].eq("1965-047A")
        self.df_reentry.loc[lead, "other_mass"] = self.df_reentry.loc[lead, "abl_mass"]
        self.df_reentry.loc[lead, "abl_mass"] = 0

        print("Setting missing reentry times.")

        launch_dt = self.dsl.set_index("COSPAR_ID")["Date"]
        launch_time = self.df_reentry["id"].str[:8].map(launch_dt)
        total_mass  = self.df_reentry["abl_mass"] + self.df_reentry["other_mass"]
        missing     = ~self.df_reentry["time_known"]
        same_day    = (~self.df_reentry["time_known"] &
                       self.df_reentry["day_known"] &
                       (launch_time.dt.normalize() == self.df_reentry["date"].dt.normalize()))
        
        # If the re-entry is on a different day to the launch, the time is set to midnight (done in convert_time).
        midnight    = missing & ~same_day

        # If the re-entry is on the same day as the launch, set the time to the launch time. 
        self.df_reentry.loc[same_day, "date"] = launch_time[same_day]

        print(f"Time missing:         {int(missing.sum())}")
        print(f"Time set to launch:   {int(total_mass[same_day].sum())},{int(same_day.sum())}")
        print(f"Time set to midnight: {int(total_mass[midnight].sum())},{int(midnight.sum())}")
        
        self.print_stats()
        
    def reentry_info_to_netcdf(self):     
        """This saves the reentry information as a NetCDF file for later processing for GEOS-Chem.
        """        
        #Set up the dimensions of the netcdf file.
        dims = 'reentries'
    
        self.df_reentry["Ablatable_Mass"] = self.df_reentry["abl_mass"] + self.df_reentry["attached_abl_mass"]
            
        #Create the DataArrays.

        ds = xr.Dataset({
            'Date'                   : xr.DataArray(self.df_reentry["date"].to_numpy(),  dims=dims, attrs=dict(long_name="Date")),
            'COSPAR_ID'              : xr.DataArray(self.df_reentry["id"].to_numpy(),  dims=dims, attrs=dict(long_name="COSPAR_ID")),
            'Category'               : xr.DataArray(self.df_reentry["category"].to_numpy(),  dims=dims, attrs=dict(long_name="Category")),
            'Object_Name'            : xr.DataArray(self.df_reentry["name"].to_numpy(),  dims=dims, attrs=dict(long_name="Object Name")),
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
        ds.to_netcdf(f'./databases/reentry_activity_data.nc')
        ds.to_dataframe().drop(columns=["Latitude", "Longitude"]).sort_values(["Date","COSPAR_ID","Category","Object_Name"]).to_csv("./databases/reentry_activity_data.csv", index=False)
        ds.close()
        
if __name__ == "__main__":  
    
    # Set up the arguments for each function.
    parser = argparse.ArgumentParser()
    parser.add_argument('-sv', "--save_reentry_info", action='store_true', help='Save reentry info.')
    parser.add_argument('-usegpd', "--use_geopandas", action='store_true', help='Load in geopandas dataframes.')
    args = parser.parse_args()
    use_gpd = args.use_geopandas

    current_year = (datetime.now().year)
    print(f"Processing from 1957 to {current_year}.")
    
    # Compile the reentry data
    def run_full_pipeline():
        Data = build_reentry_list(1957, current_year-1)
        Data.get_reentry_info()
        if args.save_reentry_info == True:
            Data.reentry_info_to_netcdf()
    
    profiler = Profile()
    profiler.enable()
    run_full_pipeline()
    profiler.disable()
    stats = Stats(profiler)
    stats.strip_dirs()
    stats.sort_stats('cumulative')
    #stats.print_stats(30)