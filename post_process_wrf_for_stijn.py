"""
Neil Berg, February 2018

Post process wrfout files from the Annenberg simulation for Stijn 

Creates monthly files of 3-hourly output
"""
import os
import sys
import datetime
import calendar
import netCDF4
import numpy as np
import fnmatch
import subprocess

from dateutil.rrule import rrule,MONTHLY	

# Create dictionairy listing Fengpeng's directory names and the real model names
mod_dict = {
			'MPI-ESM-LR': 'MPIESMLR',
			'CNRM-CM5': 'CNRMCM5',
			'IPSL-CM5A-LR': 'IPSLCM5ALR',
			'GFDL-CM3': 'GFDLCM3', 
			'inmcm4': 'INMCM4'
			}

def monthly_file_maker(strt_dt, end_dt, dom, mod=None):
	""" Create monthly files of 3-hourly WRF output """

	# Define the data paths based on historical or future run (i.e if mod==None or not)
	if not mod: # historical 
		print('Defining paths for baseline run')
		in_path1 = '/data/modelOutput/WRF/Annenberg/Baseline_1980_1990/' # 1980-01-01 - 1990-12-31
		in_path2 = '/data/modelOutput/WRF/Annenberg/Baseline/' # 1991-01-01 - 2014-12-05
		in_path3 = '/data/modelOutput/WRF/Annenberg/Baseline_2007_d1d2/' # 2014-01-01 - 2015-07-31
		out_path = '/work/nberg/data_transfers/stijn_hantson/output/historical/'

	else: # future
		print('Defining paths for future run using model {0}'.format(mod))
		in_path = '/data/modelOutput/WRF/Annenberg/Future_'+mod_dict[mod]+'/'
		out_path = '/work/nberg/data_transfers/stijn_hantson/output/future/'+mod+'/'

	# Read in the invariant file for this domain.
	inv_fl = '/data/modelOutput/WRF/Annenberg/invariant/invariant_'+str(dom)+'.nc'
	try:
		inv_nc = netCDF4.Dataset(inv_fl, 'r')
	except IOError as e:
		sys.exit(e)
	hgt		= inv_nc.variables['HGT'][0,:,:]
	xlat	= inv_nc.variables['XLAT'][0,:,:]
	xlon 	= inv_nc.variables['XLONG'][0,:,:]
	xlon_u	= inv_nc.variables['XLONG_U'][0,:,:]
	xlat_u 	= inv_nc.variables['XLAT_U'][0,:,:]
	xlon_v 	= inv_nc.variables['XLONG_V'][0,:,:]
	xlat_v 	= inv_nc.variables['XLAT_V'][0,:,:]
	sina	= inv_nc.variables['SINALPHA'][0,:,:]
	cosa	= inv_nc.variables['COSALPHA'][0,:,:]
	inv_nc.close()

	# Define dimension sizes for unstaggering and incremental variables.
	nlat = xlat.shape[0]
	nlat_stag = xlat_v.shape[0]
	nlon = xlon.shape[1]
	nlon_stag = xlon_u.shape[1]
	nlev = 43
	nlev_stag = 44
	if dom in ['d02', 'd03']:
		ntim = 7 #8, but zero based indexing
	else:
		ntim = 3 #4, but zero based indexing

	# Loop over months and output monthly files of 3-hourly variables.
	for dt in rrule(MONTHLY, dtstart=strt_dt, until=end_dt):
		# Determine number of days in this month.
		ndays = calendar.monthrange(int(dt.year),int(dt.month))[1]
		nhrs = ndays*8
		
		# Create empty arrays for variables in netCDF file.
		times_arr = np.zeros([nhrs], dtype='int')
		tracker = 0
		for d in range(1,ndays+1):
			for h in range(0,22,3):
				date_str = '{0}{1:02d}{2:02d}{3:02d}'.format(
							dt.year,dt.month,d,h) 
				times_arr[tracker] = int(date_str)
				tracker += 1
		
		# 3D variables needed by Stijn: Precip, T2, T skin, U10, V10 
		pr_arr	= np.zeros([nhrs,nlat,nlon], dtype='f4')
		tsk_arr	= np.zeros([nhrs,nlat,nlon], dtype='f4')
		t2_arr 	= np.zeros([nhrs,nlat,nlon], dtype='f4')
		u10_arr = np.zeros([nhrs,nlat,nlon], dtype='f4')
		v10_arr = np.zeros([nhrs,nlat,nlon], dtype='f4')

		# 4D variables [time,lev,lat,lon]
		# 4D variables needed to create RH for Stijn: QVAPOR, P, PB, and T
		qv_arr	= np.zeros([nhrs,nlev,nlat,nlon], dtype='f4')
		p_arr	= np.zeros([nhrs,nlev,nlat,nlon], dtype='f4')
		pb_arr	= np.zeros([nhrs,nlev,nlat,nlon], dtype='f4')
		ta_arr	= np.zeros([nhrs,nlev,nlat,nlon], dtype='f4')

		# Read in today's and yesterday's file for incremental calculations.
		for dy in np.arange(0,ndays,1):
			curr_dt = datetime.datetime(int(dt.year),int(dt.month),int(dy+1))
			yest_dt	= curr_dt + datetime.timedelta(days=-1)
			
			# Switch between the three historical directories based on the date
			if 1980 <= curr_dt.year <= 1990:
				src_dir = in_path1
			elif 1991 <= curr_dt.year < 2014:
				src_dir = in_path2
			else:
				src_dir = in_path3

			curr_fl = src_dir+'wrfout_'+dom+'_'+str(curr_dt).replace(" ","_")
			yest_fl = src_dir+'wrfout_'+dom+'_'+str(yest_dt).replace(" ","_")

			try:
				curr_nc = netCDF4.Dataset(curr_fl, 'r')
				yest_nc = netCDF4.Dataset(yest_fl, 'r')
			except IOError as e:
				sys.exit(e)
			
			# Extract 3D variables.
			crainnc	= curr_nc.variables['RAINNC'][:,:,:]
			crainc	= curr_nc.variables['RAINC'][:,:,:]
			yrainnc = yest_nc.variables['RAINNC'][:,:,:]
			yrainc	= yest_nc.variables['RAINC'][:,:,:]
			t2		= curr_nc.variables['T2'][:,:,:]
			tsk		= curr_nc.variables['TSK'][:,:,:]
			u10		= curr_nc.variables['U10'][:,:,:]
			v10		= curr_nc.variables['V10'][:,:,:]

			# Extract 4D variables.
			qv		= curr_nc.variables['QVAPOR'][:,:,:,:]
			prs		= curr_nc.variables['P'][:,:,:,:]
			prsb	= curr_nc.variables['PB'][:,:,:,:]
			ta		= curr_nc.variables['T'][:,:,:,:]

			# Calculate incremental 3-hourly precip
			# Last time step from yesterday + all 8 hours today for today's incremental array
			tmp = np.zeros([ntim+2, nlat, nlon], dtype='f4')
			tmp[0] = yrainc[ntim,:,:] + yrainnc[ntim,:,:]
			tmp[1:ntim+2] = crainc[:,:,:] + crainnc[:,:,:]
			inc_precip = np.diff(tmp, axis=0)
			# In case any incremental precip is negative (it shouldn't be),clip at 0
			inc_precip = np.clip(inc_precip, a_min=0, a_max=None)
			
			# Store 3-hourly variables in output arrays
			pr_arr[dy*8:dy*8+8,:,:] = inc_precip[:,:,:]
			t2_arr[dy*8:dy*8+8,:,:]	= curr_nc.variables['T2'][:,:,:]
			tsk_arr[dy*8:dy*8+8,:,:] = curr_nc.variables['TSK'][:,:,:]
			
			# Unrotate (i.e. earth relative) U10 and V10
			u10_unrot = u10*cosa - v10*sina
			v10_unrot = v10*cosa + u10*sina
			u10_arr[dy*8:dy*8+8,:,:] = u10_unrot[:,:,:]
			v10_arr[dy*8:dy*8+8,:,:] = v10_unrot[:,:,:]

			# Close today and yesterday's files, get ready for next time step.
			curr_nc.close()
			yest_nc.close()

		# Output monthly files of 3-hourly variables.
		if not os.path.exists(out_path+str(dt.year)):
			os.makedirs(out_path+str(dt.year))
		
		if not mod:
			fileout = out_path+str(dt.year)+'/'+'wrfpost_'+dom+'_'+str(dt.year)+str(dt.month).rjust(2,'0')+'.nc'
		else:
			fileout = out_path+str(dt.year)+'/'+'wrfpost_'+mod+'_'+dom+'_'+str(dt.year)+str(dt.month).rjust(2,'0')+'.nc'

		if os.path.exists(fileout):
			os.remove(fileout)
		
		ncfile_out = netCDF4.Dataset(fileout, 'w')
		# Create output dimensions. 
		ncfile_out.createDimension('time', nhrs)
		ncfile_out.createDimension('latitude', nlat)
		ncfile_out.createDimension('longitude', nlon)

		# Create output variables. 
		lat_out = ncfile_out.createVariable('latitude', 'f4', ('latitude', 'longitude',))
		lon_out = ncfile_out.createVariable('longitude', 'f4', ('latitude', 'longitude',))
		times_out = ncfile_out.createVariable('timestamps', 'int', ('time',))
		pr_out	= ncfile_out.createVariable('PRECIP', 'f4', ('time','latitude','longitude',))
		t2_out	= ncfile_out.createVariable('T2', 'f4', ('time','latitude','longitude',))
		tsk_out	= ncfile_out.createVariable('TSK', 'f4', ('time','latitude','longitude',))
		u10_out	= ncfile_out.createVariable('U10', 'f4', ('time','latitude','longitude',))
		v10_out	= ncfile_out.createVariable('V10', 'f4', ('time','latitude','longitude',))

		setattr(times_out, 'description', 'Times as YYYYMMDDHH')

		setattr(lat_out, 'units', 'degrees_north')
		setattr(lon_out, 'units', 'degrees_east')

		setattr(pr_out,'description','incremental total precipitation')
		setattr(pr_out,'units','mm')

		setattr(t2_out,'description','temperature at 2m')
		setattr(t2_out,'units','K')

		setattr(tsk_out,'description','surface skin temperature')
		setattr(tsk_out,'units','K')

		setattr(u10_out,'description','zonal wind speed at 10m')
		setattr(u10_out,'units','m s-1')

		setattr(v10_out,'description','meridional wind speed at 10m')
		setattr(v10_out,'units','m s-1')

		# Populate output variables with data.
		times_out[:]	= times_arr[:]
		lat_out[:,:]	= xlat[:,:]
		lon_out[:,:]	= xlon[:,:]
		pr_out[:,:,:]	= pr_arr[:,:,:]
		t2_out[:,:,:]	= t2_arr[:,:,:]
		tsk_out[:,:,:]	= tsk_arr[:,:,:]
		u10_out[:,:,:]	= u10_arr[:,:,:]
		v10_out[:,:,:]	= v10_arr[:,:,:]

		ncfile_out.close()
		
		# Make 'Time' the record dimension
		# >>> ncks --mk_rec_dmn Time fileout -O fileout
		subprocess.call(['ncks', '--mk_rec_dmn', 'Time', fileout, '-O', fileout])
		print('Finished processing {0}'.format(fileout))

if __name__=='__main__':
	# BASELINE
	strt_dt = datetime.datetime(1980,1,1)
	end_dt = datetime.datetime(2014,12,1)
	#for dom in ['d01', 'd02', 'd03']:
	for dom in ['d02']:
		ret = monthly_file_maker(strt_dt, end_dt, dom)
	
	# FUTURE
	#strt_dt = datetime.datetime(2101,4,1)
	#end_dt = datetime.datetime(2101,4,1) 
	#models = ['CNRM-CM5', 'GFDL-CM3', 'inmcm4', 'IPSL-CM5A-LR', 'MPI-ESM-LR']
	#models = ['CNRM-CM5']
	#for model in models:
	#	for dom in ['d01', 'd02', 'd03']:
	#		ret = monthly_file_maker(strt_dt, end_dt, dom, model)	
