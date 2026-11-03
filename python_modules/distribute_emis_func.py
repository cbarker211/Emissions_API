import numpy as np
from netCDF4 import Dataset

def make_grid_LL(llres, in_extent=[-180, 180, -90, 90], out_extent=[]):
    """
    Creates a lat/lon grid description.

    Args:
        llres: str
            lat/lon resolution in 'latxlon' format (e.g. '4x5')

    Keyword Args (optional):
        in_extent: list[float, float, float, float]
            Describes minimum and maximum latitude and longitude of initial grid
            in the format [minlon, maxlon, minlat, maxlat]
            Default value: [-180, 180, -90, 90]
        out_extent: list[float, float, float, float]
            Describes minimum and maximum latitude and longitude of target grid
            in the format [minlon, maxlon, minlat, maxlat]. Needed when intending
            to use grid to trim extent of input data
            Default value: [] (assumes value of in_extent)

    Returns:
        llgrid: dict
            dict grid description of format {'lat'   : lat midpoints,
                                            'lon'   : lon midpoints,
                                            'lat_b' : lat edges,
                                            'lon_b' : lon edges}
    """

    # get initial bounds of grid
    [minlon, maxlon, minlat, maxlat] = in_extent
    [dlat, dlon] = list(map(float, llres.split('x')))
    lon_b = np.linspace(minlon - dlon / 2, maxlon - dlon /
                        2, int((maxlon - minlon) / dlon) + 1)
    lat_b = np.linspace(minlat - dlat / 2, maxlat + dlat / 2,
                        int((maxlat - minlat) / dlat) + 2)
    if minlat <= -90:
        lat_b = lat_b.clip(-90, None)
    if maxlat >= 90:
        lat_b = lat_b.clip(None, 90)
    lat = (lat_b[1:] + lat_b[:-1]) / 2
    lon = (lon_b[1:] + lon_b[:-1]) / 2
    

    # Trim grid bounds when your desired extent is not the same as your initial grid extent
    if out_extent == []:
        out_extent = in_extent
    if out_extent != in_extent:
        [minlon, maxlon, minlat, maxlat] = out_extent
        minlon_ind = np.nonzero(lon >= minlon)
        maxlon_ind = np.nonzero(lon <= maxlon)
        lon_inds = np.intersect1d(minlon_ind, maxlon_ind)
        lon = lon[lon_inds]
        # make sure to get edges of grid correctly
        lon_inds = np.append(lon_inds, np.max(lon_inds) + 1)
        lon_b = lon_b[lon_inds]

        minlat_ind = np.nonzero(lat >= minlat)
        maxlat_ind = np.nonzero(lat <= maxlat)
        lat_inds = np.intersect1d(minlat_ind, maxlat_ind)
        lat = lat[lat_inds]
        # make sure to get edges of grid correctly
        lat_inds = np.append(lat_inds, np.max(lat_inds) + 1)
        lat_b = lat_b[lat_inds]

    llgrid = {'lat': lat,
            'lon': lon,
            'lat_b': lat_b,
            'lon_b': lon_b}
    return llgrid

def read_gc_box_height(gc_file,nlev):
    ''' This function reads GEOS-Chem box height from StateMet files. '''
    
    #Open Dataset.
    fk = Dataset(gc_file, mod='r')
    box_height = fk.variables['Met_BXHEIGHT'][:]
    box_height = box_height[0,:,:,:]
    area = fk.variables['AREA'][:]
    fk.close()

    # Calculate the top, mid, and bottom altitudes of each box.
    top_alt = np.zeros_like(box_height)
    bot_alt = np.zeros_like(box_height)
    mid_alt = np.zeros_like(box_height)
    for l in range(nlev):
        if l==0:
            top_alt[l,:,:] = box_height[l,:,:]
            mid_alt[l,:,:] = 0.5*np.add(bot_alt[l,:,:], top_alt[l,:,:])
            continue
        bot_alt[l,:,:] = top_alt[l-1,:,:]
        top_alt[l,:,:] = np.add(top_alt[l-1,:,:], box_height[l,:,:])
        mid_alt[l,:,:] = 0.5*np.add(bot_alt[l,:,:], top_alt[l,:,:])
    
    return box_height, area, bot_alt, mid_alt, top_alt

def get_ross_profiles():
    """Get Ross profiles provided in 5-km increments in Figure 1 of
        https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013EF000160
    """
    # Centre altitude (km):
    # These are in 5-km increments, so 80-85, 75-80 and so on.
    ross_alt = np.arange(2.5, 107.5, 5)
    ross_alt_edge =  np.arange(0, 110, 5)
    
    # New values using webplotdigitizer. Percentages add to 79% so 21% unaccounted for, and must be emitted above 100 km.
    # So adding an extra bin at midpoint 102.5 containing 21% of emissions.    
    ross_prop_mass = np.asarray([21.0, 2.0, 2.1, 2.3, 2.4, 2.6, 2.7, 2.9, 3.5, 4.0, 3.2, 2.3, 2.5, 2.7, 3.0, 3.3, 3.9, 4.4, 5.0, 6.9, 17.3])
    ross_prop_mass = ross_prop_mass[::-1] / np.sum(ross_prop_mass) * 100
    
    # Create cumulative distribution of propellant mass. 
    # This is used in the interp_prop_mass function to work out how much propellant is emitted above the gc limits.
    ross_cumulative_mass = np.zeros((len(ross_prop_mass)+1))
    for i in range(1,len(ross_prop_mass)+1):
        ross_cumulative_mass[i] = ross_cumulative_mass[i-1] + ross_prop_mass[i-1]
            
    return ross_alt_edge, ross_prop_mass, ross_cumulative_mass

def interp_prop_mass(bot_alt, mid_alt, top_alt, ross_alt_edge, ross_cumulative_mass):
    
    """Interpolate proportional propellant mass burned from Ross and Sheaffer
        5-km incremental grid to desired vertical grid.
        Input  (i) and  (j) 
        
        Vertical grids are ordered [alt, lat, lon].

    Args:
        i (numpy.int64): longitude index of model gridbox where launch occurs.
        j (numpy.int64): latitude index of model gridbox where launch occurs.
    """        
    # Interpolate the cumulative distribution to work out how much propellant is emitted above the top altitude.
    propellant_in_model = np.interp(top_alt[-1]*1e-3, ross_alt_edge, ross_cumulative_mass)
    
    # Define array of percent distribution of propellant mass on model grid:
    # This sets up an array over number of model levels.
    model_relative_mass = np.zeros_like(mid_alt)
    # Loop over model layers and interpolate the top and bottom cumulative mass of each bin:
    for l, h in enumerate(mid_alt):
        top_mass = np.interp(top_alt[l]*1e-3,ross_alt_edge, ross_cumulative_mass)
        bot_mass = np.interp(bot_alt[l]*1e-3,ross_alt_edge, ross_cumulative_mass)
        model_relative_mass[l] = top_mass - bot_mass
    
    return propellant_in_model, model_relative_mass

def get_failed_launch_dict(id, fail_code):

    if id[5] not in ["F","U"]:
        return None
    
    failed_launch_dict = {
        # Launches which did full burns of all stages.
        '2017-F02': {"full_burn": True}, #                                                (AD)  https://rocketlabcorp.com/missions/launches/its-a-test/
        '2020-F01': {"full_burn": True}, # Didn't quite reach orbital velocity.           (S2U?)
        '2020-F09': {"full_burn": True}, # Early shutdown                                 (S2P)
        '2021-F03': {"full_burn": True}, # Underperformance                               (U)
        '2021-F05': {"full_burn": True}, # Fairing separation issue                       (FF)
        '2021-F11': {"full_burn": True}, # Didn't reach orbital velocity.                 (S2U?)
        '2022-F05': {"full_burn": True}, # Early shutdown                                 (S2P)
        '2023-F01': {"full_burn": True}, # Underperformance                               (S2P)
        '2023-F11': {"full_burn": True}, # Went through burns                             (S2U)
        '2024-F03': {"full_burn": True}, # Fairing separation issue                       (S4I?)
        '2024-U01': {"full_burn": True}, # Starship V1.0 reached SECO                     (-)
        '2024-U03': {"full_burn": True}, # Starship V1.0 reached SECO                     (-)
        '2024-U04': {"full_burn": True}, # Starship V1.0 reached SECO                     (-)
        '2024-U06': {"full_burn": True}, # Starship V1.0 reached SECO                     (-)
        '2025-F01': {"full_burn": True}, # Starship V2 Early shutdown                     (S2P)
        '2025-F03': {"full_burn": True}, # Starship V2 Early shutdown                     (S2P)
        '2025-U01': {"full_burn": True}, # Underperformance                               (S1S)
        '2025-U02': {"full_burn": True}, # Starship V2 reached SECO                       (-)
        '2025-U03': {"full_burn": True}, # Starship V2 reached SECO                       (-)
        '2025-U04': {"full_burn": True}, # Starship V2 reached SECO                       (-)
        '2025-U05': {"full_burn": True}, # Early shutdown of 4th stage                    (S4P)
        '2025-U06': {"full_burn": True}, # Fairing separation issue                       (FF)

        # Specific fractions.
        '1963-F02': {"failed_stage": 1, "cutoff_alt": 44000}, # Detonated at 44 km        (S0I) https://kirkmcd.princeton.edu/JEMcDonald/mcdonald_weatherwise_16_99_63.pdf
        '2015-F02': {"failed_stage": 1, "cutoff_alt": 44900}, #                           (S2P) https://www.youtube.com/watch?v=PuNymhcTtSQ
        '2020-F07': {"failed_stage": 1, "cutoff_alt": 900},   #                           (S1G)
        '2021-F01': {"failed_stage": 1, "cutoff_alt": 10700}, #                           (S1G)
        '2021-F07': {"failed_stage": 1, "cutoff_alt": 31000}, #                           (S1P)  
        '2025-F08': {"failed_stage": 1, "cutoff_alt": 3000},  #                           (S1)
        '2022-F03': {"failed_stage": 2, "burned": 240/315},   #                           (S2P)
        '2021-F09': {"failed_stage": 3, "burned": 475/521},   #                           (S3U)
        '2022-F06': {"failed_stage": 2, "burned": 2/3},       #                           (S2P)

        # Corrections.
        '1963-F06': {"failed_stage": 1, "cutoff_alt": 0}, # No info, so treating as U     (AD)
        '1964-F08': {"failed_stage": 2, "burned": 0},     # Second stage exploded         (AD)  https://astronautix.com/t/thewrongstuhiclefailures.html
        '1975-F01': {"failed_stage": 2, "burned": 0},     # Reached S2 ignition           (S0S) https://en.wikipedia.org/wiki/Atlas-Centaur#AC-33_failure
        '1987-F02': {"failed_stage": 1, "cutoff_alt": 0}, # First stage hit by lightning  (WX)  https://astronautix.com/t/thewrongstuhiclefailures.html
        '2003-F02': {"failed_stage": 2, "burned": 0},     # Destroyed sometime during S2  (S0S) 
        '2010-F02': {"failed_stage": 1, "cutoff_alt": 0}, # First stage failure           (AD)  https://astronautix.com/2/2010chronology.html
        '2021-F08': {"failed_stage": 2, "burned": 0}, # Exploded only 12s before MECO     (S1P)
        '2022-F01': {"failed_stage": 2, "burned": 0}, # Fairing failed before S2 ignition (FF)
        '2022-F04': {"failed_stage": 3, "burned": 0}, # Failed just before 2/3 separation (S2G)
        '2023-F05': {"failed_stage": 2, "burned": 0}, # Exploded near end of booster burn (S1S)
        '2023-F07': {"failed_stage": 3, "burned": 0}, # Accidental destruct               (AD)
        '2023-F09': {"failed_stage": 3, "burned": 0}, # Stage 3 attitude failure          (S1P)
        '2023-F10': {"failed_stage": 1, "cutoff_alt": 0}, # Unknown. so ignoring          (U)
        '2024-F01': {"failed_stage": 1, "cutoff_alt": 0}, # Accidental destruct           (AD)
        '2024-F02': {"failed_stage": 2, "burned": 0}, # Failed/exploded after 2 mins      (S1P)
        '2024-F04': {"failed_stage": 3, "burned": 0}, # S1 error but still did S2 burns   (S1S)
    }
    
    failure = failed_launch_dict.get(id)
    
    # If the failure isn't in the manual list, then calculate automatically.
    if failure is None:
        failure_gcat = fail_code.replace("?","")

        # Sometimes the fail code has two parts, so we use the stage part.
        if "/" in failure_gcat:
            for part in failure_gcat.split("/"):
                if part.startswith("S") and part[1].isdigit():
                    failure_gcat = part

        # This is where the stage number is given by the fail_code.
        if failure_gcat.startswith("S"):
            failed_stage_gcat = int(failure_gcat[1])
            if failed_stage_gcat in [0,1]:
                failure = {'failed_stage': 1, 'cutoff_alt': 0}
            else:
                failure = {'failed_stage': failed_stage_gcat, 'burned': 0}

        # Unknown failure codes, ignoring.
        elif failure_gcat in ["U",   # Failure, all details unknown
                              "UX"]: # Unconfirmed failure (may not have launched, or reaching orbit may not have been intended).
            failure = {"failed_stage": 1, "cutoff_alt": 0}

        # Failure during payload deployment, so treat as normal.
        elif failure_gcat in [
            "-",   # No FailCode
            "FF",  # Fairing related failure (no fairing sep, or MaxQ failure, fairing destroyed)
            "FS",  # Structure failure other than fairing and other than due to previous loss of control
            "PI",  # Payload orbit insertion system, ignition problem.
            "PG"]: # Payload orbit insertion system, attitude or guidance problem
            return None
        else:
            raise ValueError(f"Unhandled GCAT fail code {fail_code} for {id}.")
        
    # Full burn but didn't reach orbit, so treat as normal.
    if failure is not None and failure.get("full_burn"):
        return None

    return failure