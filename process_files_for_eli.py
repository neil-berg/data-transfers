"""
Process daily T2/TSK, precip, Q2, and U/V
in the baseline and future time periods (all GCMs). 

Neil Berg, July 2018
"""

import os
import sys
import subprocess

dest_dir = '/work/nberg/data_transfers/sbc_placeworks/output/'

def baseline(src_dir, strt_yr, end_yr):
	""" Process baseline data """

	for yr in range(strt_yr, end_yr+1):
		sub_dir = src_dir+'{0}/'.format(yr)
		for mo in range(1,13):
			src_fl = sub_dir+'wrfpost_d02_{0}{1:0>2d}.nc'.format(yr, mo)
			dest_fl = dest_dir+'wrfpost_d02_{0}{1:0>2d}.nc'.format(yr,mo)
			# ncks -v SFROFF,SNOW src_fl -O dest_fl
			subprocess.call(['ncks', '-v', 'SFROFF,SNOW', src_fl, '-O', dest_fl])
			print(dest_fl)

def future(src_dir, gcms, strt_yr, end_yr):
	""" Process future data """
	
	for gcm in gcms:
		for yr in range(strt_yr, end_yr+1):
			sub_dir = src_dir+'{0}/{1}/'.format(gcm, yr)
			for mo in range(1,13):
				src_fl = sub_dir+'wrfpost_{0}_d02_{1}{2:0>2d}.nc'.format(gcm, yr, mo)
				dest_fl = dest_dir+'wrfpost_d02_{0}_{1}{2:0>2d}.nc'.format(gcm, yr,mo)
				subprocess.call(['ncks', '-v', 'SFROFF,SNOW', src_fl, '-O', dest_fl])
				print(dest_fl)

if __name__=='__main__':

	#base_dir = '/data/modelOutput/WRF/AnnenbergPost/baseline/'
	#base_strt_yr = 1991
	#base_end_yr = 2000
	#ret = baseline(base_dir, base_strt_yr, base_end_yr)
	
	future_dir = '/data/modelOutput/WRF/AnnenbergPost/future/'
	gcms = ['CNRM-CM5', 'GFDL-CM3', 'inmcm4', 'IPSL-CM5A-LR', 'MPI-ESM-LR']
	future_strt_yr = 2091
	future_end_yr = 2100
	ret = future(future_dir, gcms, future_strt_yr, future_end_yr)

