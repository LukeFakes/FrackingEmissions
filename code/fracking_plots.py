#!/usr/bin/env python

# code for plotting up fracking results
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from matplotlib import colors
import sys
from defra_geos_plotting import plot_diff_between_outputs
sys.path.append('../DEFRA_network/code')

def diverge_map(high=(0.565, 0.392, 0.173), low=(0.094, 0.310, 0.635)): 
	'''Make a custom diverging cmap. Will go from high, to white, to low.
	Arguments
	---------
	- high, low (tuple) : RGB values of the colorbar limits
	Returns
	-------
	- LinearSegmentedColormap
	'''
	import matplotlib.colors as colors
	def make_colormap(seq):
		seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3] 
		cdict = {'red': [], 'green': [], 'blue': []} 
		for i, item in enumerate(seq): 
			if isinstance(item, float): 
				r1, g1, b1 = seq[i - 1] 
				r2, g2, b2 = seq[i + 1] 
				cdict['red'].append([item, r1, r2]) 
				cdict['green'].append([item, g1, g2]) 
				cdict['blue'].append([item, b1, b2]) 
		return colors.LinearSegmentedColormap('CustomMap', cdict) 
	c = colors.ColorConverter().to_rgb 
	return make_colormap([low, c('white'), 0.5, c('white'), high]) 

if __name__ == __main___:
	
	# vik_data = np.loadtxt('../vik.txt')
	# vik2 = diverge_map(vik_data[-1],vik_data[0])
	# Reds_cmap = colors.LinearSegmentedColormap.from_list('custom_reds',['white',(0.5,0.0,0.0)])
	
	# # work out how many kg per month for each month
	# months = np.arange(1,13)
	# emis = xr.open_dataset('Fracking_Lancashire_10pct_emis.nc') 
	# hemco = xr.open_dataset('geosfp_0.5x0.625_fracking_test_0pct/OutputDir/HEMCO_diagnostics.201902010000.nc')
	# area = hemco['AREA']
	# # for each month, check position of the max value (in kg/m2/s) is the same
	# seconds_in_month = 2.628e6
	# scale_by = []
	# for month in months: 
	# 	data = emis.sel(time=f'2019-{month:02}')['FRACKING']
	# 	coords_max = data.where(data==data.max(),drop=True).squeeze().drop('time').coords
	# 	#print(f'{month}:{data.values.max()}, lat {coords_max["lat"].values}, lon {coords_max["lon"].values}') 
	# 	area_size = area.sel(lat=coords_max['lat'],lon=coords_max['lon']).values
	# 	mass_emitted = data.values.max()*area_size*seconds_in_month
	# 	scale_by.append(mass_emitted)
	# scale_by_inverse = [1/i for i in scale_by]
	# scale_by_year = np.sum(scale_by)
	# scale_by_year_inverse = 1/scale_by_year
	
	#loading data
	files_0pct = glob.glob('geosfp_0.5x0.625_fracking_test_0pct/OutputDir/GEOSChem.SpeciesConcSubset.2019*')
	files_0pct.sort()
	files_0pct
	data_0pct = xr.concat((xr.open_mfdataset(files_0pct[:11]).sel(time=slice('2019-01-07','2019-11-08')),xr.open_mfdataset(files_0pct[11:])),dim='time')
	files_10pct = glob.glob('geosfp_0.5x0.625_fracking_test_10pct_weekqueue/OutputDir/GEOSChem.SpeciesConcSubset.2019*')
	files_10pct.sort()
	data_10pct = xr.open_mfdataset(files_10pct).sel(time=slice('2019-01-07','2019-12-31'))
	fracking2 = {'0pct':data_0pct,'10pct':data_10pct}
	 
	pm_files_0pct = glob.glob('geosfp_0.5x0.625_fracking_test_0pct/OutputDir/GEOSChem.AerosolMass.2019*')
	pm_files_0pct.sort()
	pm_ds_0pct_1 = xr.open_mfdataset(pm_files_0pct[:11]).sel(time=slice('2019-01-08','2019-11-08')) 
	pm_ds_0pct_2 = xr.open_mfdataset(pm_files_0pct[11:])    
	pm_ds_0pct = xr.concat((pm_ds_0pct_1,pm_ds_0pct_2),dim='time')
	pm_files_10pct = glob.glob('geosfp_0.5x0.625_fracking_test_10pct_weekqueue/OutputDir/GEOSChem.AerosolMass.2019*')
	pm_files_10pct.sort()
	pm_ds_10pct = xr.open_mfdataset(pm_files_10pct)  
	pm_ds_10pct = pm_ds_10pct.sel(time=slice('2019-01-08','2019-12-31'))
	fracking_pm2 = {'0pct':pm_ds_0pct,'10pct':pm_ds_10pct}
	for k in fracking_pm2.keys():
		fracking_pm2[k] = fracking_pm2[k].sel(time=slice('2019-01-08','2019-12-31'))


	# plotting
	######################################################
	#   	FRACKING (KIRBY MISPERTON/HIGH MUFFLES)
	######################################################
	fig,ax=plt.subplots()
	plot_DEFRA_v_GEOS_multi_case(defra.loc[defra.site=='High Muffles'],geos_concs,'2019-01-01','2019-01-31','no',colors={'0pct':'r','5pct':'c','10pct':'g', '15pct':'b','20pct':'m'},site_type=None,scale_factor=1e9,title='Modelled NO at High Muffles with additional fracking emissions',fill=False,ax=ax,label='High Muffles obs',legend=True)   
	fig.savefig('HighMuffles_NO_frackingcases.png') 
	
	fig,ax = plt.subplots()
	plot_DEFRA_v_GEOS_base_case(defra.loc[defra.site=='High Muffles'],geos_concs,'2019-01-01','2019-01-31','no',site_type=None,title='Modelled NO at High Muffles with additional fracking emissions',fill=False,ax=ax,label='High Muffles obs',legend=True)   
	plot_DEFRA_v_GEOS_base_case(defra.loc[defra.site=='High Muffles'],fracking_5pct,'2019-01-01','2019-01-31','no',site_type=None,fill=False,ax=ax,label='High Muffles obs',legend=True,color='blue')   
	ax.legend()
	fig.savefig('HighMuffles_NO_fracking5pct.png')
	fig,ax=plt.subplots()
	plot_DEFRA_v_GEOS_diff(defra.loc[defra['code']=='CHBO'],geos_concs,'2019-01-01','2019-12-31','o3',site_type=None,title='Chilbolton % difference in O3',ax=ax,fill=False)
	fig.savefig('plots/individual_sites/Chilbolton/CHBO_o3_pctdiff.png',dpi=320)
	fig,axes = plt.subplots(nrows=2,sharex=True)
	for i, ax in enumerate(axes.ravel()):
		if i == 0:
			plot_GEOS_at_site(lat=54,lon=-0.625,geos_data=all_fracking,start_date='2019-01-01',
					end_date='2019-01-31',species='no',title='Kirby Misperton',time_freq=None,time_fmt='%d',
					geos_label={'0pct':'+0 %','5pct':'+5 %','10pct':'+10 %','15pct':'+15 %','20pct':'+20 %'},
					ax=ax,color={'0pct':'red','5pct':'blue','10pct':'green','15pct':'cyan','20pct':'magenta'})
		else:
		#plot_DEFRA_v_GEOS_base_case(defra.loc[defra['site']=='High Muffles'],{'0pct':geos_concs['0pct']},'2019-01-01','2019-01-31','no',site_type=None,geos_label='0pct',ax=ax,fill=False,time_freq=None)
			plot_DEFRA_v_GEOS_multi_case(DEFRA_data=defra.loc[defra['site']=='High Muffles'],geos_data=all_fracking,start_date='2019-01-01',end_date='2019-01-31',species='no',colors={'0pct':'red','5pct':'blue','10pct':'green','15pct':'cyan','20pct':'magenta'},site_type=None,title='High Muffles',scale_factor=1e9,time_fmt='%d',ax=ax,label='Observations',fill=False,time_freq=None)
			ax.legend()
	plt.subplots_adjust(hspace=0.3)
	fig.savefig('KirbyMisperton_test.png') 
	plot_DEFRA_v_GEOS_base_case(DEFRA_data=defra,geos_data=BT_tower,start_date='2019-01-07',end_date='2019-12-30',species='o3',avg_diurnal=True,site_type='Rural Background',title='Average O3 diurnal at Rural Background sites',geos_label='BT tower diurnal',ax=ax,label='Observations',fill=True,time_fmt=None,time_freq='hour',plot_DEFRA=True,ax_settings=False)  

	#####################################################################################
	#								ONE DAY- CONTOURS									#
	#####################################################################################
	# NO2
	plot_diff_between_outputs(geos_dict=fracking2,basecase='0pct',latmin=45,latmax=65,lonmin=-15,
			lonmax=5,species=f'SpeciesConc_NO2',savename=f'NO2_coolwarm_Contour_diff_20190710.png',
			title='Change in NO2 with fracking NOx emissions, 2019-07-10',
			diff='mass_per_kg',scale_by=scale_by_inverse[6]*1000,cmap='coolwarm',
			cbar_label=r'$\Delta$ NO2 $\left(\frac{\mu gm^{-3}}{Tonne NOx emitted month^{-1}}\right)$',
			monthly=False,time='2019-07-10',contour=True,levels=27)

	# O3
	plot_diff_between_outputs(geos_dict=fracking2,basecase='0pct',latmin=45,latmax=65,lonmin=-15,
			lonmax=5,species=f'SpeciesConc_O3',savename=f'O3_cReds_Contour_diff_20190710.png',
			title='Change in O3 with fracking NOx emissions, 2019-07-10',
			diff='mass_per_kg',scale_by=scale_by_inverse*1000,cmap=Reds_cmap,vmax=0,
			cbar_label=r'$\Delta$ NO2 $\left(\frac{\mu gm^{-3}}{Tonne NOx emitted month^{-1}}\right)$',
			monthly=False,time='2019-07-10',contour=True,levels=27)
	# PM25
	plot_diff_between_outputs(geos_dict=fracking_pm2,basecase='0pct',latmin=45,latmax=65,lonmin=-15,
			lonmax=5,species=f'PM25',savename=f'PM25_coolwarm_Contour_diff_20190710.png',
			title='Change in PM2.5 with fracking NOx emissions, 2019-07-10',
			diff='mass_per_kg',scale_by=scale_by_inverse[6]*1000,cmap='coolwarm',
			cbar_label=r'$\Delta$ PM2.5 $\left(\frac{\mu g m^-3}{Tonne NOx emitted month^{-1}}\right)$',
			monthly=False,time='2019-07-10',contour=True,levels=27)

	####################################################################################
	#									 DAILY										   #
	####################################################################################
	# NO2
	plot_diff_between_outputs(geos_dict=fracking2,basecase='0pct',latmin=45,latmax=65,lonmin=-15,
			lonmax=5,species=f'SpeciesConc_NO2',savename=f'NO2_cReds_20190710.png',
			title='Change in NO2 with fracking NOx emissions, 2019-07-10',
			diff='mass_per_kg',scale_by=scale_by_inverse[6]*1000,cmap=Reds_cmap,vmin=0,
			cbar_label=r'$\Delta$ NO2 $\left(\frac{\mu gm^{-3}}{Tonne NOx emitted month^{-1}}\right)$',
			monthly=False,time='2019-07-10',contour=True,levels=30)
	# O3
	plot_diff_between_outputs(geos_dict=fracking2,basecase='0pct',latmin=45,latmax=65,lonmin=-15,
			lonmax=5,species=f'SpeciesConc_O3',savename=f'O3_vik2_20190710.png',
			title='Change in O3 with fracking NOx emissions, 2019-07-10',
			diff='mass_per_kg',scale_by=scale_by_inverse[6]*1000,cmap=vik2,
			cbar_label=r'$\Delta$ O3 $\left(\frac{\mu gm^{-3}}{Tonne NOx emitted month^{-1}}\right)$',
			monthly=False,time='2019-07-10',contour=True,levels=31)
	# PM2.5
	plot_diff_between_outputs(geos_dict=fracking_pm2,basecase='0pct',latmin=45,latmax=65,lonmin=-15,
			lonmax=5,species=f'PM25',savename=f'PM25_cReds_20190710.png',
			title='Change in PM2.5 with fracking NOx emissions, 2019-07-10',
			diff='mass_per_kg',scale_by=scale_by_inverse[6]*1000,cmap=Reds_cmap,vmin=0,
			cbar_label=r'$\Delta$ PM2.5 $\left(\frac{\mu g m^-3}{Tonne NOx emitted month^{-1}}\right)$',
			monthly=False,time='2019-07-10',contour=True,levels=30)

	#for each day in july, plot contours for PM and NO2
	for day in np.arange(1,31):
		plot_diff_between_outputs(geos_dict=fracking_pm2,basecase='0pct',latmin=45,latmax=65,lonmin=-15,
				lonmax=5,species=f'PM25',savename=f'JulyPM/PM25_201906{day:02}.png',
				title=f'Change in PM2.5 with fracking NOx emissions,201905{day:02}',
				diff='mass_per_kg',scale_by=scale_by_inverse[6]*1000,cmap=vik2,
				cbar_label=r'$\Delta$ PM2.5 $\left(\frac{\mu g m^-3}{Tonne NOx emitted month^{-1}}\right)$',
				monthly=False,time=f'2019-05-{day:02}',contour=True,levels=30)
		plt.close('all')	
	for day in np.arange(1,31):
		plot_diff_between_outputs(geos_dict=fracking2,basecase='0pct',latmin=45,latmax=65,lonmin=-15,
				lonmax=5,species=f'SpeciesConc_NO2',savename=f'JulyNO2/NO2_201907{day:02}.png',
				title=f'Change in NO2 with fracking NOx emissions,201905{day:02}',
				diff='mass_per_kg',scale_by=scale_by_inverse[6]*1000,cmap=vik2,
				cbar_label=r'$\Delta$ NO2 $\left(\frac{\mu g m^-3}{Tonne NOx emitted month^{-1}}\right)$',
				monthly=False,time=f'2019-07-{day:02}',contour=True,levels=31)
		plt.close('all')	

	####################################################################################
	#									 ANNUAL										   #
	####################################################################################
	#NO2
	plot_diff_between_outputs(geos_dict=fracking2,basecase='0pct',latmin=45,latmax=65,lonmin=-15,
			lonmax=5,species=f'SpeciesConc_NO2',savename=f'NO2_Annual_cReds.png',
			title='Annual change in NO2 with fracking NOx emissions',
			diff='mass_per_kg',scale_by=scale_by_inverse[6]*1000,cmap=Reds_cmap,vmin=0,
			cbar_label=r'$\Delta$ NO2 $\left(\frac{\mu gm^{-3}}{Tonne NOx emitted month^{-1}}\right)$',
			monthly=False,time='2019',contour=True,levels=30)
	# O3
	plot_diff_between_outputs(geos_dict=fracking2,basecase='0pct',latmin=45,latmax=65,lonmin=-15,
			lonmax=5,species=f'SpeciesConc_O3',savename=f'O3_Annual_vik2.png',
			title='Annual change in O3 with fracking NOx emissions',
			diff='mass_per_kg',scale_by=scale_by_inverse[6]*1000,cmap=vik2,
			cbar_label=r'$\Delta$ O3 $\left(\frac{\mu gm^{-3}}{Tonne NOx emitted month^{-1}}\right)$',
			monthly=False,time='2019',contour=True,levels=31)
	# PM2.5
	plot_diff_between_outputs(geos_dict=fracking_pm2,basecase='0pct',latmin=45,latmax=65,lonmin=-15,
			lonmax=5,species=f'PM25',savename=f'PM25_Annual_cReds.png',
			title='Annual change in PM2.5 with fracking NOx emissions',
			diff='mass_per_kg',scale_by=scale_by_inverse[6]*1000,cmap=Reds_cmap,vmin=0,
			cbar_label=r'$\Delta$ PM2.5 $\left(\frac{\mu g m^-3}{Tonne NOx emitted month^{-1}}\right)$',
		monthly=False,time='2019-07',contour=True,levels=30)

	####################################################################################
	#									 ANNUAL										   #
	####################################################################################
	plot_diff_between_outputs(geos_dict=fracking2,basecase='0pct',latmin=49,latmax=60,lonmin=-12,
			lonmax=3,species=f'SpeciesConc_NO2',savename=f'NO2_cReds_diff_2019.png',
			title='Annual average change in NO2 with fracking NOx emissions',
			diff='mass_per_kg',scale_by=scale_by_year_inverse*1000,cmap=Reds_cmap,vmin=0,vmax=.0003,
			cbar_label=r'$\Delta$ NO2 $\left(\frac{\mu gm^{-3}}{Tonne NOx emitted annually}\right)$',
			monthly=False,time='2019',contour=True,levels=30)

	plot_diff_between_outputs(geos_dict=fracking3,basecase='0pct',latmin=49,latmax=60,lonmin=-12,
			lonmax=3,species=f'SpeciesConc_O3',savename=f'O3_vik2_diff_2019.png',
			title='Annual average change in O3 with fracking NOx emissions',
			diff='mass_per_kg',scale_by=scale_by_year_inverse*1000,cmap=vik2,
			cbar_label=r'$\Delta$ NO2 $\left(\frac{\mu gm^{-3}}{Tonne NOx emitted annually}\right)$',
			monthly=False,time='2019',contour=True,levels=31)

	plot_diff_between_outputs(geos_dict=fracking_pm2,basecase='0pct',latmin=49,latmax=60,lonmin=-12,
			lonmax=3,species=f'PM25',savename=f'PM25_cReds_diff_2019.png',
			title='Annual average change in PM2.5 with fracking NOx emissions',
			diff='mass_per_kg',scale_by=scale_by_year_inverse*1000,cmap=Reds_cmap,vmin=0,
			cbar_label=r'$\Delta$ PM2.5 $\left(\frac{\mu g m^-3}{Tonne NOx emitted annually}\right)$',
			monthly=False,time='2019',contour=True,levels=30)

	plot_diff_between_outputs(geos_dict=fracking2,basecase='0pct',latmin=45,latmax=65,lonmin=-15,
			lonmax=5,species=f'SpeciesConc_NO2',savename=f'NO2_test_annual.png',
			title='Change in NO2 with fracking NOx emissions, 2019-07-10',
			diff='mass_per_kg',scale_by=scale_by_year_inverse*1000,cmap=Reds_cmap,vmin=0,
			cbar_label=r'$\Delta$ NO2 $\left(\frac{\mu gm^{-3}}{Tonne NOx emitted month^{-1}}\right)$',
			monthly=False,time='2019',contour=True,levels=30)

	####################################################################################
	#									 MONTHLY									   #
	####################################################################################
	# NO2
	plot_diff_between_outputs(geos_dict=fracking,basecase='0pct_weekqueue',latmin=45,latmax=65,lonmin=-15,
			lonmax=5,species=f'SpeciesConc_NO2',savename=f'v2_NO2_Monthly_cReds.png',
			title='Monthly change in NO2 with fracking NOx emissions',
			diff='mass_per_kg',scale_by=scale_by_inverse*1000,cmap=Reds_cmap,vmin=0,
			cbar_label=r'$\Delta$ NO2 $\left(\frac{\mu gm^{-3}}{Tonne NOx emitted month^{-1}}\right)$',
			monthly=True,contour=True,levels=30)
	# O3
	plot_diff_between_outputs(geos_dict=fracking,basecase='0pct_weekqueue',latmin=45,latmax=65,lonmin=-15,
			lonmax=5,species=f'SpeciesConc_O3',savename=f'v2_O3_Monthly_vik2.png',
			title='Monthly change in O3 with fracking NOx emissions',
			diff='mass_per_kg',scale_by=scale_by_inverse*1000,cmap=vik2,
			cbar_label=r'$\Delta$ O3 $\left(\frac{\mu gm^{-3}}{Tonne NOx emitted month^{-1}}\right)$',
			monthly=True,contour=True,levels=31)

	# PM2.5
	plot_diff_between_outputs(geos_dict=frackingPM,basecase='0pct_weekqueue',latmin=45,latmax=65,lonmin=-15,
			lonmax=5,species=f'PM25',savename=f'v2_PM25_Monthly.png',
			title='Monthly change in PM2.5 with fracking NOx emissions',
			diff='mass_per_kg',scale_by=scale_by_inverse*1000,cmap=vik2,
			cbar_label=r'$\Delta$ PM2.5 $\left(\frac{\mu g m^-3}{Tonne NOx emitted month^{-1}}\right)$',
			monthly=True,time='2019-07',contour=True,levels=30)
