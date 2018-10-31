; Process relative humidity in WRF ouput
; Neil Berg, February 2018

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin

year = "2014"

src_dir = "/annenberg_raw2/output/Baseline_2007_d1d2/"
files = systemfunc("ls " + src_dir + "wrfout_d02_"+ year + "*")  ; return all netCDF files
nfiles = dimsizes(files)

dest_dir = "/data4/nberg/data_transfers/stijn_irvine/2014/rh_files/"

; Loop through each file, compute RH and output in new netCDF file
do nf=0,nfiles-1

	ncfile = addfile(files(nf), "r")
	qv = ncfile->QVAPOR
	P  = ncfile->P	; perturbation
	Pb = ncfile->PB ; base state pressure
	P  = P + Pb     ; total pressure
	theta = ncfile->T    ; perturbation potential temperature (theta+t0)
	theta = theta + 300.
	TK = wrf_tk (P, theta)
	RH = wrf_rh (qv, P, TK)

	RH_sfc = RH(:,0,:,:)

	split = str_split(files(nf), "/")
	suffix = split(3)
	dest_file = dest_dir + "RH_" + suffix + ".nc"
	
	system("/bin/rm -f " + dest_file)
	ncfile_out = addfile(dest_file, "c")

	ncfile_out->RH = RH_sfc
	RH@units = "%"

	print(dest_file)

end do
end
