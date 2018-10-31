"""
Aggregate daily RH files into monthly files and then append to the post-processed WRF files. 

Neil Berg, February 2018
"""

import os
import sys
import datetime
import calendar
import subprocess
import numpy as np

from dateutil.rrule import rrule, MONTHLY

strt_dt = datetime.datetime(2014,1,1)
end_dt = datetime.datetime(2014,12,1)

src_dir = '/data4/nberg/data_transfers/stijn_irvine/2014/rh_files/'
dest_dir = '/data4/nberg/data_transfers/stijn_irvine/2014/'

for dt in rrule(MONTHLY, dtstart=strt_dt, until=end_dt):
	daily_fls = []
	ndays = calendar.monthrange(int(dt.year),int(dt.month))[1]
	for dy in xrange(0,ndays,1):
		curr_dt = datetime.datetime(int(dt.year),int(dt.month),int(dy+1))
		curr_fl = src_dir+'RH_wrfout_d02_'+str(curr_dt).replace(" ","_")+'.nc'
		# Make time the record dimension
		subprocess.call(['ncks', '--mk_rec_dmn', 'Time', curr_fl, '-O', curr_fl])
		daily_fls.append(curr_fl)

	monthly_fl = dest_dir+'RH_wrfout_d02_'+datetime.datetime.strftime(dt, '%Y%m')+'.nc'
	subprocess.call(['ncrcat'] + daily_fls + ['-O', monthly_fl])
	subprocess.call(['ncrename', '-d', 'south_north,latitude', 
								'-d', 'west_east,longitude', 
								'-d', 'Time,time', monthly_fl, '-O', monthly_fl]) 
	subprocess.call(['ncks', '-4', monthly_fl, '-O', monthly_fl])	
	wrf_fl = dest_dir+'wrfpost_d02_'+datetime.datetime.strftime(dt, '%Y%m')+'.nc'
	# Add RH file to the wrfpost file for this month
	subprocess.call(['ncks', '-A', monthly_fl, wrf_fl])
	print(monthly_fl)
