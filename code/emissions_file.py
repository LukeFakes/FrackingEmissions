# Functions and code used to create emissions point source files for fracking work

###########################################################################################
# 					FIND NEAREST SITE TO POINT
###########################################################################################
# emis_site = (54.2042,-0.892848)
# site_dict = {}
# site_list = defra['site'].unique()
# for site in site_list:
# 	if not np.isnan(site):
# 		site_dict[site] = defra.loc[defra.site==site][['latitude','longitude']].iloc[0].values


def find_nearest_site(site_dict, location, return_distances=False):
    """ 
	Given a site_dict and a location, return either the name of the nearest site to the 
	location or a dictionary of the distance of each site to the locations
	
	Arguments
	---------
	- site_dict (dict): Site locations to look through, in format {'sitename':(lat,lon)}
	- location (tuple, (lat,lon)): Location to refer to.
	- return_distance (bool): Default False.
	
	Returns
	-------
	either: 
		-  closest_site (str) , distance (float)
		-  distance_dict (dict). Distance between location and each site in site_dict.
		
	Example
	-------
	nearest_site, distance = find_nearest_site(site_dict,(54.2042,-0.892848))  
	"""
    import numpy as np

    distance_dict = {}
    closest_site = None
    closest_distance = np.inf
    for site in site_dict:
        lat_dist = site_dict[site][0] - location[0]
        lon_dist = site_dict[site][1] - location[1]
        distance = np.sqrt((lat_dist ** 2 + lon_dist ** 2))
        distance_dict[site] = distance
        if closest_site is None:
            closest_site = site
        if closest_distance > distance:
            closest_site = site
            closest_distance = distance
    if not return_distances:
        return closest_site, closest_distance
    else:
        print(f"Nearest site: {closest_site}, Distance: {closest_distance}")
        return distance_dict


###########################################################################################
# 			MAKING emissions point source files:
###########################################################################################
# check emissions for site location from HEMCO diagnostics
# hemco = xr.open_mfdataset(hemco_files)
# #trial site: Kirby Misperton (54.204185, -0.892848) or Lancashire Preston New Road (53.787282,2.951474)
# KirMis = hemco['EmisNO_Anthro'].sel(time='2019-01-01').isel(lev=0).sel(lat=54.204185,lon=-0.892848,method='nearest')
# for k in frack_emis.keys():
# 	NO_emissions = frack_emis[k].isel(lev=0,time=0).sel(lat=54.204185,lon=-0.892848,method='nearest')['EmisNO_Anthro'].values
# 	print(k, NO_emissions)
# Lanc_Frack = hemco['EmisNO_Anthro'].isel(lev=0).sel(lat=53.787282,lon=-2.951474,method='nearest').isel(lev=0).values


def make_emission_file(
    value, lat, lon, pct, savename, savedir, varname="FRACKING", **kwargs
):
    """
	Create a COARDS-compliant emissions netcdf for fracking work.
	
	Arguments
	--------- 
	value (float): value to emission at the site for a month
	lat,lon (floats): position of site
	pct: the percentage of that value to set emission at that point to
	savename, savedir (str):
	varname (str): Default 'FRACKING'
	**kwargs:
		- title (str): Title for netcdf file
	"""
    import datetime
    from netCDF4 import Dataset, date2num
    import numpy as np
    import xarray as xr

    # arbitrary domain size, just making sure its bigger than the UK nested runs ((45,65),(-15,5))
    lats = np.arange(45, 65.5, 0.5)
    lons = np.arange(-15, 5.625, 0.625)
    if len(value) == 1:
        times = [datetime.datetime(2019, 1, 1, 0, 0, 0, 0)]
    elif len(value) == 12:
        times = [datetime.datetime(2019, m, 1, 0, 0, 0, 0) for m in np.arange(1, 13)]
    else:
        print("length of values is not 1 or 12")
    times = date2num(
        times, units="hours since 2019-01-01 00:00:00 00:00", calendar="standard"
    )
    array = np.full((len(times), len(lats), len(lons)), fill_value=1e-32)
    lat_ix = np.abs(lats - lat).argmin()
    lon_ix = np.abs(lons - lon).argmin()

    fill_value = value * (pct / 100)
    print(f"Fill value at location is {fill_value}")
    for i, t in enumerate(times):
        if len(times > 1):
            array[i, lat_ix, lon_ix] = fill_value[i]
        else:
            array[i, lat_ix, lon_ix] = fill_value
    print(f"array max is {np.max(array)}")

    New_f = Dataset(savedir + f"{savename}.nc", mode="w", format="NETCDF4")
    New_f.createDimension("time", len(times))
    New_f.createDimension("lat", len(lats))
    New_f.createDimension("lon", len(lons))
    New_f.Conventions = "COARDS"
    if "title" not in kwargs:
        New_f.title = "Fracking NOx emissions 2020-01-01"
    else:
        New_f.title = kwargs["title"]
    New_f.history = f'created {datetime.datetime.today().strftime("%Y-%m-%d")} by LF'
    New_f.SpatialCoverage = "UK nested grid domain"
    New_f.Resolution = "0.1x0.1"

    New_f.createVariable("time", "f8", ("time"))
    New_f["time"].long_name = "Time"
    New_f["time"].units = f"hours since 2019-01-01 00:00:00 00:00"
    New_f["time"].calendar = "standard"
    New_f["time"][:] = times
    New_f["time"].axis = "T"

    New_f.createVariable("lat", "f4", ("lat"))
    New_f["lat"].long_name = "Latitude"
    New_f["lat"].units = "degrees_north"
    New_f["lat"][:] = lats  # _list
    New_f["lat"].axis = "Y"

    New_f.createVariable("lon", "f4", ("lon"))
    New_f["lon"].long_name = "Longitude"
    New_f["lon"].units = "degrees_east"
    New_f["lon"][:] = lons  # _list
    New_f["lon"].axis = "X"

    New_f.createVariable(varname, "f8", ("time", "lat", "lon"))
    New_f[varname][:] = array
    New_f[varname].units = "kg/m2/s"
    New_f[varname].missing_value = 1e15
    New_f.add_offset = 0
    New_f[varname].long_name = f"{varname}"

    New_f.set_auto_mask(False)
    New_f.close()
    print(f"done")


# savedir = '.'
# files = {}
# for P in [5,10,15,20]:
# 	savename = f'Fracking_KirbyMisperton_{P}pct_emis'
# 	make_emission_file(value,lat=54.204185,lon=-0.892848,pct=P,savename=savename,savedir='')
# 	files[f'NO_{P}pct'] = xr.open_dataset(f'{savename}.nc')['FRACKING']
# 	print(P, files[f'NO_{P}pct'].max(),files[f'NO_{P}pct'].max()/value )
#
# for P in [5,10,15,20]:
# 	savename = f'Fracking_Lancashire_{P}pct_emis'
# 	make_emission_file(v,lat=53.787282,lon=-2.951474,pct=P,savename=savename,savedir='')
# 	files[f'NO_{P}pct'] = xr.open_dataset(f'{savename}.nc')['FRACKING']
# 	print(P, files[f'NO_{P}pct'].max(),files[f'NO_{P}pct'].max()/v )
