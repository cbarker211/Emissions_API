import numpy as np
import argparse
import xarray as xr
import pandas as pd
from datetime import datetime
import sys
import random
import geopandas as gpd
from shapely.geometry import Point, Polygon
import fiona

class build_reentry_list:
    
    def __init__(self, start_year, final_year):
        """This function sets up the class and imports data.

        Args:
            start_year, final_year (int): The initial and final years of re-entry to cover.
        """        
        
        self.start_year = start_year
        self.final_year = final_year
        
        # Import the shapefiles and create unions for the main oceans.
        ocean_shapefile     = gpd.read_file("./databases/reentry/GOaS_v1_20211214/goas_v01.shp")
        pacific             = ocean_shapefile[ocean_shapefile.name.isin(['North Pacific Ocean','South Pacific Ocean'])]
        pacific_union       = gpd.GeoDataFrame(geometry=[pacific.unary_union]) 
        atlantic            = ocean_shapefile[ocean_shapefile.name.isin(['North Atlantic Ocean','South Atlantic Ocean'])]
        atlantic_union      = gpd.GeoDataFrame(geometry=[atlantic.unary_union])
        country_shapefile   = gpd.read_file("./databases/reentry/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp") 
        states_shapefile    = gpd.read_file("./databases/reentry/ne_10m_admin_1_states_provinces/ne_10m_admin_1_states_provinces.shp")
        self.ocean_shapefile   = ocean_shapefile
        self.pacific_union     = pacific_union
        self.atlantic_union    = atlantic_union
        self.country_shapefile = country_shapefile
        self.states_shapefile  = states_shapefile
            
    def print_stats(self):
        """Print out statistics about the database at the end of the program.
        """        
        
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
            print(f"0 Inc for {jsr_id}")
            inc=90
        
        if latlonstr == "-":
            # If location is missing for lower stages for failed and successful launches, set as launch coordinates.
            if (apogee <= 100) and (category in ["S0","S1"]): 
                for i in range(len(self.dsl["COSPAR_ID"])):
                    if jsr_id[:8] == self.dsl["COSPAR_ID"].values[i][:8]:
                        lat = self.dsl["Latitude"].values[i]
                        lon = self.dsl["Longitude"].values[i]
                        location = 2  
            else: 
                lat = round(random.uniform(-inc, inc),2)
                lon = round(random.uniform(-180, 180),2)
                location = 6  
        elif latlonstr == "Alashan": 
            # Landing expected east of Jiuquan launch site (this was actually launched from Wenchang).
            # Assuming that Alashan refers to Helen Shan Mountain (new name for Alashan Mountain).
            # https://space.skyrocket.de/doc_sdat/rcs-fc-sc.htm
            # https://planet4589.org/space/jsr/news.778
            # https://www.peakbagger.com/peak.aspx?pid=10693
            # https://www.nasaspaceflight.com/2020/05/china-next-generation-crew-capsule/
            lat = 38.833
            lon = 105.95
            location = 2
        elif latlonstr == "Jiuquan":
            # This refers to the Chinese launch center, sources say it landed 'in the desert near Jiuquan', and 'in China’s Inner Mongolia autonomous region'.
            # Using lat/lon for Jiuquan launch site.
            # https://discosweb.esoc.esa.int/launch-sites/21
            # https://nssdc.gsfc.nasa.gov/nmc/spacecraft/display.action?id=2020-027A
            # https://space.skyrocket.de/doc_sdat/xzf-sc.htm 
            # https://spaceflightnow.com/2020/05/08/chinas-next-generation-crew-spacecraft-lands-after-unpiloted-test-flight/ 
            lat = 41.3
            lon = 100.3               
            location = 2              
        elif latlonstr == "LOPNOR RW05":
            # Reported test flight of a Chinese reusable experimental spacecraft.
            # Runway 05/22 at Lop Nor air base in Xinjiang (https://planet4589.org/space/gcat/data/tables/lp.html).
            # https://www.seradata.com/china-launches-own-mini-spaceplane-reusable-spacecraft-using-long-march-2f 
            # https://twitter.com/planet4589/status/1302486141090885632                
            lat = 40.78
            lon = 89.27
            location = 2
        elif latlonstr == "Ocean":
            # Electron Stage 1 (rcat) is occasionally recovered in the ocean.
            # From https://spaceflightnow.com/2020/11/05/rocket-lab-to-attempt-booster-recovery-on-next-mission/:
            #   Recovery vessels stationed near the booster’s splashdown zone around 250 miles (400 kilometers) 
            #   south of the launch site will move in to secure the first stage and hoist it onto a ship for return to New Zealand.
            # Geoloating 400km south of launch site.
            lat = -42.86 
            lon = 177.87
            location = 2
        elif latlonstr == "Antarctic?":
            location = 4
            while True:
                lat = round(random.uniform(-inc, inc),2)
                lon = round(random.uniform(-180, 180),2)
                coordinate = Point(lon,lat)
                if lat <= -60:
                    break 
        elif latlonstr == "Arctic": 
            location = 4
            while True:
                lat = round(random.uniform(-inc, inc),2)
                lon = round(random.uniform(-180, 180),2)
                coordinate = Point(lon,lat)
                if lat >= 66:
                    break
        elif latlonstr in ["S Ocean", "S Ocean?"]:
            location = 4
            # For these three objects in auxcat, the inclination is 51.65-53.
            # This is outside of the general definition of the southern ocean (<60S).
            # Therefore we will just set the lat to the inclination.
            lat = -inc
            lon = round(random.uniform(-180, 180),2)
        #elif latlonstr in ["Pacific","Pacific?","PO","AO","E Pacific","E Pacific?","S Pacific",
        #                   "S POR","Indian O?","SE IOR?","SE IOR","Atlantic","POR","Kazakhstan",
        #                   "Mexico","S Africa","Gujarat"]:
        #    lat = 0
        #    lon = 0
        elif latlonstr == "WSSH":
            location = 2
            # This is the White Sands Missile Range in New Mexico.
            lat = 33.238462 
            lon = -106.346383 
        elif latlonstr == "GM 86.2W 29.7N":
            location = 1
            # Gulf of Mexico
            lat = 29.7
            lon = -86.2
        elif latlonstr == "Splash?":
            location = 2
            # This is the Electron failed helicopter stage 1 recovery in May 22. 
            # https://www.youtube.com/watch?v=BY0CXlOeWHI "Several hundred kilometers from the launch site."
            # Using inclination of 94 degrees (wiki), and distance of 300km.
            lat = -36.57
            lon = 177.63
        elif latlonstr == "Gujarat":
            location = 3
            gujurat     = self.states_shapefile[self.states_shapefile.name_en.isin(['Gujarat'])] 
            while True:
                lat = round(random.uniform(-inc, inc),2)
                lon = round(random.uniform(-180, 180),2)
                coordinate = Point(lon,lat)
                if gujurat.geometry.contains(coordinate).any():
                    break
        elif latlonstr == "KSC SLF": # https://www.world-airport-codes.com/united-states/nasa-shuttle-landing-facility-69738.html
            # https://www.spaceforce.mil/News/Article/3217077/x-37b-orbital-test-vehicle-concludes-sixth-successful-mission/
            location = 2
            lat = 28.61
            lon = -80.69
        elif latlonstr == "Kazakhstan":
            location = 3
            kazakhstan     = self.country_shapefile[self.country_shapefile.NAME.isin(['Kazakhstan'])]
            while True:
                lat = round(random.uniform(-inc, inc),2)
                lon = round(random.uniform(-180, 180),2)
                coordinate = Point(lon,lat)
                if kazakhstan.geometry.contains(coordinate).any():
                    break
        elif latlonstr == "Mexico":
            location = 3
            mexico     = self.country_shapefile[self.country_shapefile.NAME.isin(['Mexico'])]
            while True:
                lat = round(random.uniform(-inc, inc),2)
                lon = round(random.uniform(-180, 180),2)
                coordinate = Point(lon,lat)
                if mexico.geometry.contains(coordinate).any():
                    break
        elif latlonstr == "S Africa":
            location = 3
            s_africa     = self.country_shapefile[self.country_shapefile.NAME.isin(['South Africa'])]
            while True:
                lat = round(random.uniform(-inc, inc),2)
                lon = round(random.uniform(-180, 180),2)
                coordinate = Point(lon,lat)
                if s_africa.geometry.contains(coordinate).any():
                    break
        elif latlonstr == "100 03E 41 39N":
            location = 1
            lon = 100.05
            lat = 41.65
        elif latlonstr == "47 20N 69 34E":     
            location = 1     
            lon = 69.57
            lat = 47.33
        elif latlonstr == "41 39N 100 09E":     
            location = 1     
            lon = 100.15
            lat = 41.65            
        elif latlonstr == "47 20N 69 35E?":     
            location = 1     
            lon = 69.58
            lat = 47.33      
        elif latlonstr == "100 04E 41 37N":
            location = 1
            lon = 100.07
            lat = 41.62
        elif latlonstr == "69 37E 47 21N":
            location = 1
            lon = 69.62
            lat = 47.35
        elif latlonstr in ["Pacific", "Pacific?", "PO", "POR"]:
            location = 4
            while True:
                random_location = self.pacific_union.sample_points(1).get_coordinates()
                lon = random_location["x"].values[0]
                lat = random_location["y"].values[0]
                if -inc <= lat <= inc:
                    break
        elif latlonstr in ["AO","Atlantic"]: 
            location = 4 
            while True:
                random_location = self.atlantic_union.sample_points(1).get_coordinates()
                lon = random_location["x"].values[0]
                lat = random_location["y"].values[0]
                if -inc <= lat <= inc:
                    break
        elif latlonstr in ["E Pacific", "E Pacific?"]:
            location = 4
            while True:
                random_location = self.pacific_union.sample_points(1).get_coordinates()
                lon = random_location["x"].values[0]
                lat = random_location["y"].values[0]
                if -inc <= lat <= inc and -180 <= lon <= -60: # Bounded to East Pacific only.
                    break 
        elif latlonstr in ["S Pacific", "S POR"]:
            location = 4
            while True:
                random_location = self.pacific_union.sample_points(1).get_coordinates()
                lon = random_location["x"].values[0]
                lat = random_location["y"].values[0]
                if -inc <= lat <= inc and lat <= 0: # Bounded to South Pacific only.
                    break 
        elif latlonstr in ["Indian O?"]:
            location = 4
            indian_ocean = self.ocean_shapefile[self.ocean_shapefile.name.isin(['Indian Ocean'])] 
            while True:
                lat = round(random.uniform(-inc, inc),2)
                lon = round(random.uniform(-180, 180),2)
                coordinate = Point(lon,lat)
                if indian_ocean.geometry.contains(coordinate).any():
                    break 
        elif latlonstr in ["S of Tasmania?"]:
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
        elif latlonstr in ["SE IOR?","SE IOR"]:
            location = 4
            indian_ocean = self.ocean_shapefile[self.ocean_shapefile.name.isin(['Indian Ocean'])] 
            while True:
                lat = round(random.uniform(-inc, -30),2)
                lon = round(random.uniform(77, 180),2)
                coordinate = Point(lon,lat)
                if indian_ocean.geometry.contains(coordinate).any():
                    break 
        elif latlonstr.replace("?","").split()[0][-1] in ["E","W","S","N"]:
            location = 1
            geolocation = latlonstr.replace("?","").split()
            if geolocation[0][-1] in ["E","W"]: # If the longitude comes first.
                lon_data = geolocation[0]
                lat_data = geolocation[1]
            elif geolocation[1][-1] in ["E","W"]: # If the longitude comes last.
                lon_data = geolocation[1]
                lat_data = geolocation[0]  
            if lon_data[-1] == "W":
                lon = np.float64(lon_data.replace("W",""))*-1
            elif lon_data[-1] == "E":
                lon = np.float64(lon_data.replace("E",""))
            if lat_data[-1] == "N":
                lat = np.float64(lat_data.replace("N",""))
            elif lat_data[-1] == "S":
                lat = np.float64(lat_data.replace("S",""))*-1   
        else:
            lat = 0
            lon = 0
            print(f"Need to sort out geolocation for {latlonstr}")
        
        lat = np.float64(lat) 
        lon = np.float64(lon) 
            
        return lat, lon, location
    
    def falcon_stage_lat_lon(self,datestr,jsr_id):
        
        """Set the geolocation for all Falcon 9 1st stages.

        Returns:
            lat(np.float64) : The latitude of the object.
            lon(np.float64) : The longitude of the object..
        """        
        
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
        elif matching.shape[0] == 0:
            matching = self.ground_landings[self.ground_landings["Date"].isin([f"{datestr[:4]}-{datestr[4:6]}-{datestr[6:]}"])].reset_index(drop=True)
            if matching.shape[0] == 1:
                lat = np.array(matching["geometry"].y)[0]
                lon = np.array(matching["geometry"].x)[0]
            elif jsr_id == "2022-144":
                # This is 1st Nov - Falcon Heavy where two boosters land. 
                lat = np.array(matching["geometry"].y)[0]
                lon = np.array(matching["geometry"].x)[0]
            elif matching.shape[0] == 0:
                found_launch = False
                for i in range(len(self.dsl["COSPAR_ID"])):
                    if self.dsl["COSPAR_ID"].values[i] == jsr_id:
                        found_launch = True
                        # Kennedy Space Center (ETR)
                        if (self.dsl["Longitude"].values[i] == -81) and (self.dsl["Latitude"].values[i] == 28.5):
                            lat = 34
                            lon = -75  
                        # Vandenberg
                        elif (self.dsl["Longitude"].values[i] == -120.6) and (self.dsl["Latitude"].values[i] == 34.7):
                            lat = 30
                            lon = -120 
                        else:
                            sys.exit("Falcon launch site not found.")
                if found_launch == False:
                    sys.exit(f"Problem geolocating Falcon Stage 1 recovery - {jsr_id}.")
            elif matching.shape[0] > 1:
                sys.exit(f"Multiple ground entries for Falcon Stage 1 landing- {jsr_id}.")
            else:
                sys.exit(f"Problem geolocating Falcon Stage 1 landing- {jsr_id}.")
        elif matching.shape[0] > 1:
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
        
        matching = self.fairings[self.fairings["Date"].isin([f"{datestr[:4]}-{datestr[4:6]}-{datestr[6:]}"])].reset_index(drop=True)
        if matching.shape[0] == 1:
            lat = np.array(matching["geometry"].y)[0]
            lon = np.array(matching["geometry"].x)[0]
        elif matching.shape[0] == 0:
            found_launch = False
            for i in range(len(self.dsl["COSPAR_ID"])):
                if self.dsl["COSPAR_ID"].values[i] == jsr_id:
                    found_launch = True
                    # Kennedy Space Center (ETR)
                    if (self.dsl["Longitude"].values[i] == -81) and (self.dsl["Latitude"].values[i] == 28.5):
                        lat = 34
                        lon = -75   
                    # Vandenberg
                    elif (self.dsl["Longitude"].values[i] == -120.6) and (self.dsl["Latitude"].values[i] == 34.7):
                        lat = 30
                        lon = -120 
                    else:
                        sys.exit("Falcon launch site not found.")
            if found_launch == False:
                sys.exit(f"Problem geolocating Falcon fairing recovery - {jsr_id}.")
        elif matching.shape[0] > 1:
            sys.exit(f"Multiple entries for Falcon fairing recovery - {jsr_id}.")
        else:
            sys.exit(f"Problem geolocating Falcon fairing recovery - {jsr_id}.")
              
        return lat, lon
        
    def convert_time(self,date):
        
        """Converts the "DDate" listing from GCAT to day, month, and time strings.

        Returns:
            datestr(str):  The date in YYYYMMDD format.
            time_utc(str): The decimal time in hours.
        """   
        yearstr=str(date[0])     
        date[1] = date[1].replace("?","")
        mon = datetime.strptime(date[1], '%b').month  
        if mon in np.arange(10):
            monstr = "0" + str(mon)
        else:
            monstr = str(mon)
        if len(date) > 2: 
            day = date[2].replace("?","")
            if int(day) in np.arange(10):
                daystr = "0" + day
            else:
                daystr = day
        elif date == ['2022', 'Apr']:
            # DISCOSweb puts reentry on Mar 31st. Setting as March 31st.
            daystr = "31"
            monstr = "03"
        elif date == ['2022', 'Dec']:
            # DISCOSweb puts reentry on Nov 30th. Setting as Nov 30th.
            daystr = "30"
            monstr = "11"
        else:
            print(date)
            daystr = "01"
        datestr = yearstr+monstr+daystr
        if len(date) >= 4:    
            time_utc = np.float64(int(date[3].replace("?","")[0:2]) + int(date[3].replace("?","")[2:4]) / 60)
        else:
            time_utc = -1  
        return datestr, time_utc 
    
    def failed_launch_mass(self, jsr_id, jsr_name, reentry_category):
        
        """For all failed launches, set the masses manually as we need to differentiate wet vs dry mass reentry.

        Returns:
            mass(float)       : The mass of aluminium.
            other_mass(float): Extra non-aluminium mass.
        """        
            
        abl_mass = 0
        other_mass = 0
        
        for i in range(len(self.dsl["COSPAR_ID"])):
            if self.dsl["COSPAR_ID"].values[i] == jsr_id:
                for count, rocket_name in enumerate(self.dsr["Rocket_Name"].values):
                    if rocket_name == self.dsl["Rocket_Name"].values[i]:
                        rocket_ind = count
                        
        if reentry_category == "S0":
            abl_mass = self.dsr["Booster_StageMass"].values[rocket_ind] / int(self.dsr["Booster_No"].values[rocket_ind])
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
                     
        # Add second stage propellant mass to non aluminium mass.    
        if jsr_id in ["2020-F02","2020-F05",
                      "2021-F02",
                      "2022-F01","2022-F02","2022-F03"] and reentry_category == "S2":
            other_mass = self.dsr[f"Stage2_PropMass"].values[rocket_ind]
            if jsr_id == "2022-F03":
                other_mass = other_mass * 0.24
                
        # Add third stage propellant mass to non aluminium mass.    
        if jsr_id in ["2020-F02","2020-F03","2020-F05","2020-F06",
                      "2021-F02","2021-F06","2021-F09","2021-F10",
                      "2022-F02","2022-F05","2022-F07"] and reentry_category == "S3":   
            other_mass = self.dsr[f"Stage3_PropMass"].values[rocket_ind]
            if jsr_id == "2021-F09":
                other_mass = other_mass * 0.08
                    
        # Add fourth stage propellant mass to non aluminium mass.  
        if jsr_id in ["2020-F08","2020-F09",
                      "2021-F10",
                      "2022-F02","2022-F04","2022-F05","2022-F07"] and reentry_category == "S4":   
            other_mass = self.dsr[f"Stage4_PropMass"].values[rocket_ind]

        return abl_mass, other_mass
        
    def sort_inclination(self,jsr_inc,jsr_id):      
        
        """When the inclination on file is zero (including failed), this function sets it appropriately.
        This is mainly the case for lower suborbital stages, which can then be set using the inclination of an upper stage or payload.

        Returns:
            jsr_inc(np.float64): The orbital inclination of the object.
        """        
        
        # Try to find an entry from satcat/auxcat that is already loaded.
        # NOTE: Occasionally, the inclinations are slightly different for each object
        # The difference is usually <1 degree, but can be higher. As we don't know which inclination is 'correct',
        # we use the larger one to bound the geolocation.
        count = 0
        for reentry in (self.unique_reentry_list):
            if reentry["id"][:8] == jsr_id[:8] and reentry["inc"] != 0:
                if count == 0:
                    jsr_inc = reentry["inc"] 
                elif count > 0:
                    if reentry["inc"] > jsr_inc:
                        #if (reentry["inc"] - jsr_inc) > 2:
                            #print(f"Inclination difference of {reentry['inc'] - jsr_inc}, {reentry['inc']}, {jsr_inc}, {jsr_id}")
                        jsr_inc == reentry["inc"]
                count += 1
            
        # Sometimes, there is no entry in the reentry list because it didn't meet the criteria before.
        # In this case we need to reload the satcat to look for any matching entries.
        if jsr_inc == 0:
            for filename in ["./databases/reentry/GCAT/satcat_2ndNov.tsv","./databases/reentry/GCAT/ftocat_2ndNov.tsv"]:
                jsr_data = pd.read_csv(filename, delimiter='\t', dtype=object)
                satcat_list = []
                for reentry_count in range(1,len(jsr_data)):
            
                    if jsr_data["Piece"][reentry_count][:8] == "2020-F09":
                        new_id = "2020-F10"
                    elif jsr_data["Piece"][reentry_count][:8] == "2020-U01":
                        new_id = "2020-F09"
                    elif jsr_data["Piece"][reentry_count][:8] == "2022-U02":
                        new_id = "2022-F04"
                    elif jsr_data["Piece"][reentry_count][:8] == "2022-F04":
                        new_id = "2022-F05"
                    elif jsr_data["Piece"][reentry_count][:8] == "2022-F05":
                        new_id = "2022-F06"
                    elif jsr_data["Piece"][reentry_count][:8] == "2022-F06":
                        new_id = "2022-F07"
                    else:
                        new_id = jsr_data["Piece"][reentry_count][:8]
                    
                    if (
                        jsr_data["Primary"][reentry_count] == "Earth" 
                        and new_id == jsr_id[:8]
                        and np.float64(jsr_data["Inc"][reentry_count]) != 0
                    ):
                        satcat_list.append(jsr_data.iloc[[reentry_count]])
                if len(satcat_list) > 0:
                    satcat = pd.concat(satcat_list, ignore_index=True)
                    count = 0
                    for i in range(len(satcat)):
                        if satcat["Status"][i] in ["O","R","DSO","DSA","AF","AS","F","S","GRP","AO","AR"]:
                            if count == 0:
                                jsr_inc = np.float64(satcat["Inc"][i])
                            elif count > 0:
                                if np.float64(satcat["Inc"][i]) > jsr_inc:
                                    #if (np.float64(satcat["Inc"][i]) - jsr_inc) > 2:
                                        #print(f"Inclination difference of {np.float64(satcat['Inc'][i]) - jsr_inc}, {np.float64(satcat['Inc'][i])}, {jsr_inc}, {jsr_id}")
                                    jsr_inc == np.float64(satcat["Inc"][i])
                            count += 1
                        else:
                            print(f'Extra item loaded for inclination, {jsr_id},{satcat["Status"][i]}')
                
        if jsr_inc == 0:
            print(f"inc still empty for {jsr_id}")
            
        return jsr_inc 
         
    def extract_jsr_info(self, filepath):
        
        """For each GCAT database, extract all the relevant information and output to the dictionary of all reentry objects.
        """        
        
        # Load in the database.
        jsr_data = pd.read_csv(filepath, delimiter='\t', dtype=object)
        # Only include specific statuses (see https://planet4589.org/space/gcat/web/intro/phases.html)
        jsr_data_stripped = jsr_data[jsr_data["Status"].isin(["R","R?","D","L","L?","S","F","AF","AS"])].reset_index(drop=True)
        
        jsr_data_stripped_list = []
        for reentry_count in range(len(jsr_data_stripped)):
            
            # For some reason, there are 4 CZ-8 boosters listed in JSR rcat for 2020-102 (R80052-54), when there should be two (Jonathan has since updated his database to correct this).
            if jsr_data_stripped["DDate"][reentry_count][0:4] != "-" and jsr_data_stripped["#JCAT"][reentry_count] not in ["R80054","R80055", # 2020-102 boosters.
                                                                                                                           "L80508","L80509","L80510"  # Skipping 2021-F07
                                                                                                                           ]:
                # Note: A negative apogee means that the object was in a hyperbolic orbit (velocity above escape velocity), so passes by Earth.
                # However, Hayabusa2 (apogee=-54771) did reenter on 5th Dec, so including this.

                if (jsr_data_stripped["Apogee"][reentry_count] == "-") or ("Inf" in jsr_data_stripped["Apogee"][reentry_count]):
                    jsr_data_stripped["Apogee"][reentry_count] = 0
                
                if (
                    self.start_year <= int(jsr_data_stripped["DDate"][reentry_count][0:4]) <= self.final_year     # This year only
                    and (int(jsr_data_stripped["Apogee"][reentry_count]) >= 50 or int(jsr_data_stripped["Apogee"][reentry_count]) == -54771) # Apogee above 50 km.           
                    and jsr_data_stripped["Piece"][reentry_count][5:6] != "S"                                     # No suborbital launches (mainly military rockets.)
                    and jsr_data_stripped["Piece"][reentry_count][:8] not in ["2021-U01","2022-U01","2022-U03"]   # Skipping 2021-U01, 2022-U01 and 2022-U03, all military tests.
                    and jsr_data_stripped["Primary"][reentry_count] == "Earth"                                    # Only objects whose primary body is Earth.
                    and jsr_data_stripped["Type"][reentry_count][0] not in ["Z", "D"]                             # Ignore Z, this is a spurious entry according to JSR, and ignore D, this is debris objects. 
                ): 
                    if jsr_data_stripped["Status"][reentry_count] == "AS":
                        if jsr_data_stripped["Piece"][reentry_count][5:6] in ["F","U"]: 
                            jsr_data_stripped_list.append(jsr_data_stripped.iloc[[reentry_count]])
                    else:
                        jsr_data_stripped_list.append(jsr_data_stripped.iloc[[reentry_count]])
        # Rebuild the database.
        if len(jsr_data_stripped_list) > 0:          
            jsr_data_stripped_range = pd.concat(jsr_data_stripped_list, ignore_index=True)
        else:
            jsr_data_stripped_range = []
        # Now loop over the list and format into a dictionary.
        for reentry_count in range(len(jsr_data_stripped_range)):
            
            jsr_id     = jsr_data_stripped_range["Piece"][reentry_count] 
            jsr_name   = jsr_data_stripped_range["Name"][reentry_count]
            jsr_inc    = np.float64(jsr_data_stripped_range["Inc"][reentry_count])
            jsr_dest   = jsr_data_stripped_range["Dest"][reentry_count]
            jsr_apogee = int(jsr_data_stripped_range["Apogee"][reentry_count])
            
            if str(jsr_name) == "Hayabusa 2 Return Capsule":
                jsr_apogee = jsr_apogee * -1
         
            # The failure numbers are labelled differently in DW and JSR. To make sure the main script runs OK, adjust the COSPAR ID here.
            if jsr_id == "2020-F09":
                jsr_id = "2020-F10"
            elif jsr_id == "2020-U01":
                jsr_id = "2020-F09"
            elif jsr_id == "2022-U02":
                jsr_id = "2022-F04"
            elif jsr_id == "2022-F04":
                jsr_id = "2022-F05"
            elif jsr_id == "2022-F05":
                jsr_id = "2022-F06"
            elif jsr_id == "2022-F06":
                jsr_id = "2022-F07"
            
            # Set the burnup as partial for all re-entering capsules.
            if jsr_name in ["Soyuz MS-13","Soyuz MS-15","Soyuz MS-16","Soyuz MS-17",
                            "Soyuz MS-18","Soyuz MS-19","Soyuz MS-20","Soyuz MS-21",
                            "Shenzhou 12", "Shenzhou 13", "Shenzhou 14",
                            "Endeavour","Endurance","Freedom","Resilience",
                            "Axiom Ax-1","Inspiration4","XZF-SC","X-37B OTV-6",
                            "Dragon CRS-18","Dragon CRS-19", "Dragon CRS-20",
                            "Dragon CRS-21","Dragon CRS-22", "Dragon CRS-23",
                            "Dragon CRS-24","Dragon CRS-25", "Dragon CRS-20",
                            "Starliner OFT-2","Artemis I","Chonfu Shiyong Shiyan HQ",
                            "Bernard Kutter LOFTID","Hayabusa 2 Return Capsule","Chang'e-5 Fanhui Qi",
                            ]:
                burnup = "Partial"
            elif jsr_name == "Electron Stage 1":
                if jsr_id in ["2020-007","2020-085","2021-F02","2021-106","2022-047","2022-147"]:
                    burnup = "Partial"
                else:
                    burnup = "Complete"   
            else:
                burnup = "Complete"
                
            # Duplicates:
            #   2019-036H. This is TEPCE 1 and TEPCE 2, listed with same "Piece" in JSR. TEPCE 2 listed as 2019-036IA in DISCOSweb, so renaming to this.
            #   2020-027A. This is Xinyidai Zairen Feichuan (XZF) and XZF Service Module (not in DISCOSweb). Setting auxcat id to XZF 2020-027D.
            #   Remaining are where multiple objects in auxcat from same launch are just given the COSPAR ID of the launch.
            if jsr_data_stripped_range["Name"][reentry_count] == "TEPCE 2":
                jsr_id = "2019-036IA" 
            elif jsr_data_stripped_range["Name"][reentry_count] == "XZF Service Module":
                jsr_id = "2020-027D"  
            elif jsr_data_stripped_range["Piece"][reentry_count] == "2020-078B": # This is just a mistake in JSR.
                jsr_name = "Falcon 9-098 Stage 2"
                
            # Sort out the reentry time/date.
            datestr, time_utc = self.convert_time(jsr_data_stripped_range["DDate"][reentry_count].split()) 
            
            # Handle categories.
            if jsr_data_stripped_range["Type"][reentry_count][0] == "R":
                # Sort out the rocket stages. This is mainly Soyuz (Russian's just do the numbering differently), and errors in the GCAT.
                reentry_category = f"S{(jsr_data_stripped_range['Type'][reentry_count][1:2])}"
                if jsr_data_stripped_range["#JCAT"][reentry_count] in ["R81714","R81715"]:
                    reentry_category = "S0"
                if jsr_data_stripped_range["Name"][reentry_count] in ["RSRMV-1L","RSRMV-1R"]:
                    reentry_category = "S0"   
                if "Blok-BVGD" in jsr_data_stripped_range["Bus"][reentry_count]:
                    reentry_category = "S0"
                elif "Blok-A" in jsr_data_stripped_range["Bus"][reentry_count]:
                    reentry_category = "S1"
                elif "Blok-I" in jsr_data_stripped_range["Name"][reentry_count]:
                    reentry_category = "S2"
                elif "Fregat" in jsr_data_stripped_range["Name"][reentry_count]:
                    reentry_category = "S3"
                elif jsr_data_stripped_range["Piece"][reentry_count] in ["2020-028C","2021-086B","2021-097B","2021-112B","2022-066B","2022-102B"]:
                    reentry_category = "S4"  
            elif jsr_data_stripped_range["Type"][reentry_count][0] in ["C","P"]:
                reentry_category = jsr_data_stripped_range["Type"][reentry_count][0]
            else:
                reentry_category = -1
                
            # Sort the inc.
            if jsr_inc == 0:
                jsr_inc = self.sort_inclination(jsr_inc,jsr_id)
                
            # Sort out the reentry lat/lon.
            if ("Falcon 9 Stage 1" in jsr_name) and (jsr_data_stripped_range["#JCAT"][reentry_count] != "R81716"):
                lat, lon = self.falcon_stage_lat_lon(datestr,jsr_id)
                burnup = "Partial"
                location = 5
            elif ("Falcon 9 Fairing" in jsr_name) or ("Falcon Heavy Fairing" in jsr_name):    
                lat, lon = self.falcon_fairing_lat_lon(datestr,jsr_id)
                burnup = "Partial"
                location = 5
            else:   
                lat, lon, location = self.convert_lat_lon(jsr_dest, jsr_inc, reentry_category, jsr_apogee, jsr_id)
            
            # Sort the mass.      
            if jsr_id[5:6] in ["F","U"]:
                abl_mass, other_mass = self.failed_launch_mass(jsr_id, jsr_name, reentry_category)
            else:
                abl_mass = np.float64(jsr_data_stripped_range["DryMass"][reentry_count])
                other_mass = 0
            
            # Check for items with a missing geolocation (should have been dealt with already so this is just a sanity check).    
            if lat == 0 and lon == 0:
                print(f"No location for {jsr_id}, dest = {jsr_data_stripped_range['Dest'][reentry_count]}")
            
            # Set up the dictionary.                      
            temp_reentry_dict = {
                "id"               : str(jsr_id),
                "jcat"             : str(jsr_data_stripped_range["#JCAT"][reentry_count]),
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
    
    def add_cargo(self,filepath):
        
        """Short function to add the cargo mass to its parent object.
        """        
        
        jsr_data = pd.read_csv(filepath, delimiter='\t', dtype=object)        
        for reentry_count in range(len(jsr_data)):           
            if jsr_data["Status"][reentry_count] in ["AL IN","AR IN"] and (self.start_year <= int(jsr_data["DDate"][reentry_count][0:4]) <= self.final_year):
                cargo_found = False 
                for reentry in self.unique_reentry_list:
                    if reentry["jcat"] == jsr_data["Parent"][reentry_count]:
                        reentry["other_mass"] += np.float64(jsr_data["DryMass"][reentry_count]) 
                        cargo_found = True    
                        if jsr_data["Piece"][reentry_count][5:6] in ["F","U","S"]:
                            sys.exit(f'Cargo in failed/suborbital object {jsr_data["Piece"][reentry_count]}.')
                    if reentry["jcat"] == "S51660" and jsr_data["#JCAT"][reentry_count] == "A09988":
                        reentry["other_mass"] += np.float64(jsr_data["DryMass"][reentry_count])
                        cargo_found = True
                    if reentry["jcat"] == "A06581" and jsr_data["#JCAT"][reentry_count] == "A09680":
                        reentry["other_mass"] += np.float64(jsr_data["DryMass"][reentry_count])
                        cargo_found = True    
                if cargo_found == False:
                    print(f'Cargo missing for {jsr_data["Piece"][reentry_count]}')                    
     
    def add_attached(self):
        
        for filepath in ["./databases/reentry/GCAT/satcat_2ndNov.tsv","./databases/reentry/GCAT/auxcat_2ndNov.tsv",
                            "./databases/reentry/GCAT/lcat_2ndNov.tsv",  "./databases/reentry/GCAT/rcat_2ndNov.tsv",
                            "./databases/reentry/GCAT/lprcat_2ndNov.tsv","./databases/reentry/GCAT/deepcat_7thNov.tsv",
                            "./databases/reentry/GCAT/ecat_8thNov.tsv"]:
            attached_data = pd.read_csv(filepath, delimiter='\t', dtype=object)
            attached_data_stripped = attached_data[attached_data["Status"].isin(["AR","AS","AL"])].reset_index(drop=True)   
            missing_list = []
            added = 0
            # First see if the parent is already in the reentry list. Add the mass to the parent object if it is, or add the details to a list if missing.
            
            for reentry_count in range(len(attached_data_stripped)): 
                if ((self.start_year <= int(attached_data_stripped["DDate"][reentry_count][0:4]) <= self.final_year) and 
                    attached_data_stripped["Piece"][reentry_count][5:6] not in ["S","F","U"]):
                    
                    found_parent = False
                    for reentry in self.unique_reentry_list:
                        if reentry["jcat"][:6] == attached_data_stripped["Parent"][reentry_count][:6]:
                            found_parent = True
                            added += 1
                            if reentry["id"][:8] != attached_data_stripped["Piece"][reentry_count][:8]:
                                pass
                                #print(f'Found missing item in original list. {reentry["id"]},{attached_data_stripped["Piece"][reentry_count]},{filepath}.')
                            reentry["attached_abl_mass"] += np.float64(attached_data_stripped["DryMass"][reentry_count]) 
                            

                    if found_parent == False and filepath not in ["./databases/reentry/GCAT/lprcat_2ndNov.tsv"]:
                        missing_list.append([attached_data_stripped['Piece'][reentry_count],
                                             attached_data_stripped['Parent'][reentry_count][:6],
                                             attached_data_stripped["DryMass"][reentry_count],
                                             filepath])          
            # Loopover all items where the parent couldn't be found. The parent is always starts with A, but is usually in the ecat.
            for i in range(len(missing_list)):
                found_parent = False
                if missing_list[i][1][0] != "A":
                    sys.exit("Item has a missing parent from another database than auxcat.")
                ecat_data = pd.read_csv("./databases/reentry/GCAT/ecat_8thNov.tsv", delimiter='\t', dtype=object)  
                # Now loop over the ecat and look for the parent.
                for reentry_count in range(len(ecat_data)): 
                    if (missing_list[i][1] == ecat_data['#JCAT'][reentry_count][:6]
                    and ecat_data['Status'][reentry_count] in ["AR","AL"]):
                        for reentry_count_2 in range(len(ecat_data)): 
                            if (ecat_data['Parent'][reentry_count][:6] == ecat_data['#JCAT'][reentry_count_2][:6]
                            and ecat_data['Status'][reentry_count_2] in ["D","L"]):
                                # If we find the grandparent, then hopefully the object is already in the reentry list.
                                for reentry in self.unique_reentry_list:
                                    if ecat_data['#JCAT'][reentry_count_2][:6] == reentry["jcat"][:6]:
                                        if reentry["id"][:8] != missing_list[i][0][:8]:
                                            pass
                                            #print(f'Found missing item in ecat. {reentry["id"]},{missing_list[i][0]},{filepath}.')
                                        reentry["attached_abl_mass"] += np.float64(missing_list[i][2])
                                        added += 1
                                        found_parent = True
                        
                if found_parent == False:
                    aux_data = pd.read_csv("./databases/reentry/GCAT/auxcat_2ndNov.tsv", delimiter='\t', dtype=object)
                    for reentry_count in range(len(aux_data)): 
                        if (missing_list[i][1] == aux_data['#JCAT'][reentry_count][:6]
                        and aux_data['Status'][reentry_count] in ["AR","AL"]):    
                               
                            # If the grandparent is also in the auxcat, then look for the grandparent.
                            if aux_data['Parent'][reentry_count][0] == "A":             
                                for reentry_count_2 in range(len(aux_data)): 
                                    if (aux_data['Parent'][reentry_count][:6] == aux_data['#JCAT'][reentry_count_2][:6]
                                        and aux_data['Status'][reentry_count_2] in ["D","L"]):
                                        # If we find the grandparent, then hopefully the object is already in the reentry list.
                                        for reentry in self.unique_reentry_list:
                                            if aux_data['#JCAT'][reentry_count_2][:6] == reentry["jcat"][:6]:
                                                if reentry["id"][:8] != missing_list[i][0][:8]:
                                                    pass
                                                    print(f'Found missing item in existing reentries. {reentry["id"]},{missing_list[i][0]},{filepath}.')
                                                reentry["attached_abl_mass"] += np.float64(missing_list[i][2])
                                                added += 1
                                                found_parent = True
                            # If the grandparent starts with S, its either in the satcat or ecat, and has probably already been added to the list.
                            elif aux_data['Parent'][reentry_count][0] == "S":
                                for reentry in self.unique_reentry_list:
                                    if aux_data['Parent'][reentry_count][:6] == reentry["jcat"][:6]:
                                        if reentry["id"][:8] != missing_list[i][0][:8]:
                                            pass
                                            print(f'Found missing item in existing reentries. {reentry["id"]},{missing_list[i][0]},{filepath}.')
                                        reentry["attached_abl_mass"] += np.float64(missing_list[i][2])
                                        added += 1
                                        found_parent = True
                            else:
                                print("Grandparent does not start with S or A.")
                # For 2022-002, for A10426/7/8, the parent chain is A09939/42 > 38 > 37. So just hard setting to avoid a huge indented block.            
                if missing_list[i][1] in ["A09942","A09939"]:
                    for reentry in self.unique_reentry_list:
                        if reentry["jcat"][:6] == "A09937":
                            reentry["attached_abl_mass"] += np.float64(missing_list[i][2])
                            added += 1
                            found_parent = True
                                        
                if found_parent == False:
                    print(missing_list[i])
                    
            #print(added)
        
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
        print(ac_uncertainty_list)
        
        # Now add any new entries.
        for reentry_count in range(len(aerospace_corp_data_year)):  
            if aerospace_corp_data_year["International Designator"][reentry_count] not in jsr_id_list:
                #print(f"Found new reentry from Aerospace Corp - {aerospace_corp_data_year['International Designator'][reentry_count]}, {aerospace_corp_data_year['Object Name'][reentry_count]}") 
                # 1992-021C, Ariane 44L H10+ 3rd Rocket Stage. DISCOSweb and AC list this as reentering in Oct 2020. JSR lists this as exploding in Apr 1993.
                # 2019-029Q, Starlink payload. DISCOSweb and JSR list this as reentering in July 2022. AC lists this as reentering in 2020. Ignoring.
                # 2020-012B, Starlink payload. DISCOSweb has no re-entry data. AC lists this as reentering Mar 2020. JSR lists this as still in orbit. Ignoring.
                # 2019-024AG, Starlink 2259, this should be 2021-024AG. Modified AC csv file.
                # 2020-001Z, Starlink 1072, listed as 2020-001X in JSR. Modified AC csv file.
                # 2021-074Z, Starlink 1901, this should be 2020-074Z.   Modified AC csv file.
                # 2021-073G, STARLINK-1731, this should be 2020-073G.   Modified AC csv file.
                # 2021-073BG, STARLINK-1827, this should be 2020-073BG. Modified AC csv file.
                # 2020-050B, this is wildly wrong, should be 2019-073B. Modified AC csv file. SSN seems to be correct each time.
                # 1977-010A , also wrong, should be 2020-064B. Modified AC csv file.
                # 1010-012BF, typo, should be 2020-012BF. Modified AC csv file.
                # 1998-067PU, this is traceable in GCAT using the JCAT. Should have correct COSPAR in GCAT now.
                # 2021-024AS, extra spaces in AC csv file, fixed.    
                # 2021-041K, STARLINK-2173, should be 2021-041H. Modified AC csv file.
                # 1982-092AJF, debris on DISCOSweb, JSR and AC, ignoring.  
                # 2021-125AE, this is already included in the 2021 reentry list, as JSR says it reentered on 28/12/21.
                if aerospace_corp_data_year["International Designator"][reentry_count] in ["1992-021C"]:
                    print(f"Adding reentry from Aerospace Corp - {aerospace_corp_data_year['International Designator'][reentry_count]}") 
                    if aerospace_corp_data_year["International Designator"][reentry_count] == "1992-021C":
                        reentry_category = "S3"
                        abl_mass = 2080
                        inc = 3.98
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
            
            # Excluding Debris
            valid_types = ["Rocket Body","Rocket Mission Related Object","Payload Mission Related Object","Payload"]
            if (self.ds_dw["DISCOSweb_Reentry_ID"].values[reentry_count] not in jsr_ac_id_list) and self.ds_dw["DISCOSweb_Reentry_Class"].values[reentry_count] in valid_types:
                
                # Check if the stage has already been added.  
                found_stage = False               
                if self.ds_dw["DISCOSweb_Reentry_Class"].values[reentry_count] == "Rocket Body":
                    if (self.ds_dw["DISCOSweb_Reentry_Name"].values[reentry_count] in ["H-II LE-5B (H-IIB)", "L-53 (YF24B) (Long March (CZ) 2D)",
                                                                                       "CZ-5-HO (Long March (CZ) 5)","Delta IV DCSS 5 (Delta 4H)"]
                        or "Falcon 9 Merlin-V (1D" in self.ds_dw["DISCOSweb_Reentry_Name"].values[reentry_count]
                        or "Centaur-5" in self.ds_dw["DISCOSweb_Reentry_Name"].values[reentry_count]):
                        reentry_category = "S2"
                    elif (("Fregat" in self.ds_dw["DISCOSweb_Reentry_Name"].values[reentry_count])
                        or (self.ds_dw["DISCOSweb_Reentry_Name"].values[reentry_count] in ["PBV (Long March (CZ) 6)","YZ-1S (Long March (CZ) 2C/YZ-1S)",
                                                                                      "H-18 (Long March (CZ) YF) (Long March (CZ) 3C/E)","Angara AM (Angara 1.2)"])):
                        reentry_category = "S3" 
                    elif self.ds_dw["DISCOSweb_Reentry_Name"].values[reentry_count] in ["AVUM (Vega)","Ceres-1 upperstage","YZ-1 (Long March (CZ) 3B/YZ-1)",
                                                                                   "KZ-1 Stage 4 (Kuaizhou-1)","Long March (CZ) 11 Stage 4"]:
                        reentry_category = "S4"
                    else:
                        print("Missing stage designation for : ",self.ds_dw["DISCOSweb_Reentry_ID"].values[reentry_count], self.ds_dw["DISCOSweb_Reentry_Name"].values[reentry_count])
                        reentry_category = "S0"
                    
                    for reentry in self.unique_reentry_list:
                        if reentry["id"][:8] == self.ds_dw["DISCOSweb_Reentry_ID"].values[reentry_count][:8] and str(reentry["category"]) == str(reentry_category):
                            found_stage = True
                            
                skipped_ids = ["2014-065A", # Chang'e 5-T1, reentered in 2014.
                               "2021-011A", # Progress MS-16 is listed in JSR as AR, so the mass is just added to the parent object and no object exists.
                               # Progress MS-16 has been added though, so we can ignore it.
                               "2020-001N", # Starlink, still in orbit??
                               "2020-001T", # Starlink, still in orbit??
                               "2020-001G", # Starlink, still in orbit??
                               "2017-069C", # Listed in JSR as still in orbit.
                               "2020-001AQ",# Starlink, still in orbit??
                               "2020-001BG",# Starlink, still in orbit??
                               "2020-001BF",# Starlink, still in orbit??
                               "2020-001BM",# Starlink, still in orbit??
                               "2020-001A", # Starlink, still in orbit??
                               "2021-122E", # Debris.
                               "2020-001AU",# Starlink, still in orbit??
                               "1967-014M", # Debris.
                               "2022-010ID", # https://www.spacex.com/launches/sl4-7/ 49 starlink satellites launched.
                               "2022-010IE","2022-010IS", # https://en.wikipedia.org/wiki/List_of_Starlink_and_Starshield_launches 
                               "2022-010IF", # Wiki says 38 satellites reentered by 12th Feb.
                               "2022-010IT", # 38 already included in GCAT, so IDs must be different. Ignoring.
                               "2022-010IG", "2022-010IQ", "2022-010IR", "2022-010II", "2022-010IM", "2022-010IAD",
                               "2022-010IO", "2022-010IB", "2022-010IY", "2022-010IJ", "2022-010IU", "2022-010IV",
                               "2022-010IAF","2022-010IH", "2022-010IL", "2022-010IX", "2022-010IAE","2022-010IW",       
                               "2022-010IK", "2022-010IC", "2022-010IZ", "2022-010IP", "2022-010IAA",
                               "2015-007B", # Debris.
                               "2014-065B", # Reentered on the moon.
                               "1990-010C", # Debris.
                               "2022-085C", # Debris.
                               "2022-094B", # Falcon Stage 2 # JSR has this has as deep space.
                               "2014-047H", # Debris.
                               "2021-110A", # Deorbited on Didymos.
                               "1998-067NF",# Already in GCAT, tracable through JCAT.
                               "1998-067PU",# Already in GCAT, tracable through JCAT.
                               ]
                
                # Included IDs: 
                #               2020-022A https://www.n2yo.com/satellite/?s=45464     
                #               2020-029B DDate = "2023 Jan?", DISCOSweb has 31st Dec, so it has reentered and we can include.   
                if found_stage == False and self.ds_dw['DISCOSweb_Reentry_ID'].values[reentry_count] not in skipped_ids:
                    print(f"Adding object with COSPAR ID {self.ds_dw['DISCOSweb_Reentry_ID'].values[reentry_count]}")          
                
                    if self.ds_dw["DISCOSweb_Reentry_Class"].values[reentry_count] in ["Rocket Mission Related Object","Payload Mission Related Object"]:
                        reentry_category = "C"
                    elif self.ds_dw["DISCOSweb_Reentry_Class"].values[reentry_count] == "Payload":
                        reentry_category = "P" 
                    else:
                        print("Unexpected re-entry stage.") 
                    datestr = int(self.ds_dw["DISCOSweb_Reentry_Epoch"].values[reentry_count][:10].replace("-",""))
                    time_utc = -1
                    
                    if self.ds_dw["DISCOSweb_Reentry_ID"].values[reentry_count] == "2020-029B":
                        inc = 45.01
                    elif self.ds_dw["DISCOSweb_Reentry_ID"].values[reentry_count] == "2020-022A":
                        inc = 26.50
                    else: 
                        inc = 0
                        print(f"Added object not expected ({self.ds_dw['DISCOSweb_Reentry_ID'].values[reentry_count]}).")
                    lat = round(random.uniform(-inc, inc),2)
                    lon = round(random.uniform(-180, 180),2)
                    mass = self.ds_dw["DISCOSweb_Reentry_Mass"].values[reentry_count]
                    if np.isnan(mass):
                        mass = 0
    
                    temp_reentry_dict = {
                        "id"               : self.ds_dw["DISCOSweb_Reentry_ID"].values[reentry_count],
                        "jcat"             : "N/A",
                        "name"             : self.ds_dw["DISCOSweb_Reentry_Name"].values[reentry_count],
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
        # Loop over each launch in 2020, locate the rocket and then filter for rocket configuration.
        for i in range(len(self.dsl["COSPAR_ID"])):
            for count, rocket_name in enumerate(self.dsr["Rocket_Name"].values):
                if rocket_name == self.dsl["Rocket_Name"].values[i] and self.dsl["COSPAR_ID"].values[i][5] != "F":
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
                                    abl_mass = self.dsr["Booster_StageMass"].values[count] / int(self.dsr["Booster_No"].values[count])
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
        
                    # Check the number of stages.
                    stage_count = np.zeros((4))
                    for reentry in self.unique_reentry_list:
                        for j in range(1,5):
                            if reentry["id"][:8] == self.dsl["COSPAR_ID"].values[i] and reentry["category"] == f"S{j}":
                                stage_count[j-1] += 1
                    #print(stage_count)
                    for stage in stage_count:
                        if stage > 1:
                            print(f"Multiple stages for {self.dsl['COSPAR_ID'].values[i]}")
                    
                    meco = stage_alt_dict[f"{rocket_name} MECO"]
                    add_first_stage = False
                    if (50 < meco <= 100):                 
                        add_first_stage = True
                    elif np.isnan(meco):
                        # 2S,3S,4S
                        if self.dsr["Booster_PropMass"].values[count] == 0: 
                            add_first_stage = True
                            if ((self.dsr["Stage3_PropMass"].values[count] == 0) and (self.dsr["Stage4_PropMass"].values[count] == 0)):
                                meco = 90 #2S
                            elif ((self.dsr["Stage3_PropMass"].values[count] != 0) and (self.dsr["Stage4_PropMass"].values[count] == 0)):
                                meco = 56 #3S
                            elif ((self.dsr["Stage3_PropMass"].values[count] != 0) and (self.dsr["Stage4_PropMass"].values[count] != 0)):
                                meco = 52 #4S
                        # B+4S
                        if ((self.dsr["Booster_PropMass"].values[count] != 0) and (self.dsr["Stage4_PropMass"].values[count] != 0)):    
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
        self.ds_dw = xr.open_dataset(f"./databases/reentry/DISCOSweb/discosweb_reentries_{self.start_year}-{self.final_year}.nc", decode_times=False)
        self.dsr = xr.open_dataset(f"./databases/rocket_attributes_{self.start_year}-{self.final_year}.nc", decode_times=False)
        self.dsl = xr.open_dataset(f"./databases/launch_activity_data_{self.start_year}-{self.final_year}.nc", decode_times=False) 
        self.import_raul_spacex_map() #(https://t.co/RAsQ9NDmEr)
             
        self.unique_reentry_list = []
        
        print("Adding JSR re-entries.") # (see https://planet4589.org/space/gcat/web/cat/cats.html)
        print("satcat")
        self.extract_jsr_info("./databases/reentry/GCAT/satcat_2ndNov.tsv")  # Load in the JSR satcat database (main database).
        print("auxcat") 
        self.extract_jsr_info("./databases/reentry/GCAT/auxcat_2ndNov.tsv")  # Load in the JSR auxcat database (should be in main but isn't).
        print("lcat") 
        self.extract_jsr_info("./databases/reentry/GCAT/lcat_2ndNov.tsv")    # Load in the JSR lcat database (suborbital stages/objects).
        print("rcat") 
        self.extract_jsr_info("./databases/reentry/GCAT/rcat_2ndNov.tsv")    # Load in the JSR rcat database (lower stages and fairings).
        print("lprcat") 
        self.extract_jsr_info("./databases/reentry/GCAT/lprcat_2ndNov.tsv")  # Load in the JSR lprcat database (several objects returning from space).
        print("deepcat")
        self.extract_jsr_info("./databases/reentry/GCAT/deepcat_7thNov.tsv") # Load in the JSR deepcat database (several objects returning from space).
        print("ecat")
        self.extract_jsr_info("./databases/reentry/GCAT/ecat_8thNov.tsv")    # Load in the JSR ecat database (capsules from ISS / crewed missions).
        print("ecat cargo")
        self.add_cargo("./databases/reentry/GCAT/ecat_8thNov.tsv")           # Load in the JSR ecat database.
        print("ftocat")
        self.extract_jsr_info("./databases/reentry/GCAT/ftocat_2ndNov.tsv")  # Load in the JSR ftocat database.
        print("attached")
        self.add_attached()

        # There is nothing relevant for this inventory in hcocat, tmpcat and csocat.   
        
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
        
        # Add missing stages and check the Aerospace Corp and DISCOSweb databases.
        print("Looking for missing rocket stages.")
        self.add_missing_stages()
        print("Searching for Aerospace Corp re-entries.")
        self.extract_aerospace_info("./databases/reentry/AerospaceCorp/AerospaceCorp_Reentries_14thDec23.csv")
        print("Searching for DISCOSweb re-entries.")
        self.extract_discosweb_info()
        
        #######################
        ## Final adjustments.
        #######################
        
        # Create a list of all smc-related launches (from launch database). 
        smc_dict = {}
        for i in range(len(self.dsl["COSPAR_ID"])): 
            smc_dict[self.dsl["COSPAR_ID"].values[i]] = self.dsl["Megaconstellation_Flag"].values[i]
            
        for reentry in self.unique_reentry_list:
            try:
                reentry["smc"] = smc_dict[reentry["id"][:8]]
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
                    print(f"Couldn't assign ablation information for complete burnup of {reentry['id']}")
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
            
            
        ################################################
        # Update mass info using rocket_info databases. 
        ################################################ 
                    
        for reentry in self.unique_reentry_list:                
            for i in range(len(self.dsl["COSPAR_ID"])):
                if self.dsl["COSPAR_ID"].values[i][:8] == reentry["id"][:8]:
                    
                    # Update mass info for all rocket stages.
                    if reentry["category"] in ["B1","B2","B3","B4","B5","B6","S1","S2","S3","S4"] and self.dsl["COSPAR_ID"].values[i][5] != "F":
                        for count, rocket_name in enumerate(self.dsr["Rocket_Name"].values):
                            if rocket_name == self.dsl["Rocket_Name"].values[i]:
                                if reentry["category"] in ["B1","B2","B3","B4","B5","B6"]:  
                                    reentry["abl_mass"] = self.dsr["Booster_StageMass"].values[count] / int(self.dsr["Booster_No"].values[count])
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
            if reentry["name"] in ["CZ-2C sep motor cover","CZ-4C sep motor cover?","CZ-2D sep motor cover",
                                   "CZ-2D sep motor cover"]:
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
        
        print(f"Time set to launch:   {time_update_mass_1},{time_update_count_1}")
        print(f"Time set to midnight: {time_update_mass_2},{time_update_count_2}")
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
        smc_list, location_list, apogee_list  = [], [], []
        
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
             
        #Save to file and close the DataSet  
        ds.to_netcdf(f'./databases/reentry_activity_data_{self.start_year}-{self.final_year}.nc')
        ds.close()
        
if __name__ == "__main__":  
    
    # Set up the arguments for each function.
    parser = argparse.ArgumentParser()
    parser.add_argument('-sv', "--save_reentry_info", action='store_true', help='Save reentry info.')
    parser.add_argument('-sy', "--start_year", default = "2020", choices=str(np.arange(1957,2023)), help='Start Year.')
    parser.add_argument('-fy', "--final_year", default = "2022", choices=str(np.arange(1957,2023)), help='Final Year.')
    args = parser.parse_args()
    
    # Sort out the year range.
    start_year = int(args.start_year)
    final_year = int(args.final_year)
    if start_year > final_year:
        final_year = start_year + 1
    print(f"Processing from year {start_year} to {final_year}.")
    
    # Compile the reentry data
    Data = build_reentry_list(start_year, final_year,)
    Data.get_reentry_info()
    if args.save_reentry_info == True:
        Data.reentry_info_to_netcdf()