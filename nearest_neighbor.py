"""
Carve out HUC-10 Sierra Nevada watersheds from the D03 grid

Neil Berg, March 2018
"""

import os
import sys
import netCDF4
import numpy as np
import pandas as pd
from scipy.spatial import distance

def nearest_neighbor(inv_fl, pt_lat, pt_lon):
	""" Locate nearest WRF grid cell given a point lat/lon """

	ncfile = netCDF4.Dataset(inv_fl, 'r')
	xlat = ncfile.variables['XLAT'][0,:,:]
	xlon = ncfile.variables['XLONG'][0,:,:]
	nlat, nlon = np.shape(xlat)
	ncfile.close()

	# set up 2D lat/lon coordinate array 
	npts = nlat*nlon
	wrf_coords = np.zeros([2,npts])
	wrf_indices = np.zeros([2,npts])
	idx = 0
	for i in range(nlat):
		for j in range(nlon):
			wrf_coords[0,idx] = xlat[i,j]
			wrf_coords[1,idx] = xlon[i,j]
			wrf_indices[0,idx] = i
			wrf_indices[1,idx] = j
			idx += 1

	# Compute distances from station coordinate to all WRF land grid cells
	obs_pt = np.array([pt_lat,pt_lon]).reshape((2,1))
	dist = distance.cdist(wrf_coords.T,obs_pt.T)

	near_idx = dist.argmin()
	near_lat, near_lon = wrf_coords[0,near_idx], wrf_coords[1,near_idx]
	near_lat_idx, near_lon_idx = int(wrf_indices[0,near_idx]), int(wrf_indices[1,near_idx])

	return near_lat_idx, near_lon_idx

if __name__=='__main__':
	inv_fl = '/Users/nberg/data_transfers/sbc_placeworks/invariant/invariant_d02.nc'
	pt_lat = 37.27
	pt_lon = -119.34
	near_lat_idx, near_lon_idx = nearest_neighbor(inv_fl, pt_lat, pt_lon)
	print(near_lat_idx, near_lon_idx)
