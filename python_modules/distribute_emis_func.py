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