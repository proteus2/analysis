
yyyy = 1993
mm = "05"
lata = 10.
p1 = 100.
klev = 35

;file_ref = "res/uanuj_ra.cgwp_offline."+yyyy+"."+mm+".nc"
file_ref = "uanuj_ra.cgwp_offline."+yyyy+"."+mm+"_rescale_5.nc"

f1 = addfile(file_ref,"r")
lat = f1->lat
z = f1->z_f(klev) * 1.e-3
wgt = where(abs(lat) .le. lata, 1., 0.)
wgt = where(abs(lat) .eq. lata, 0.5, wgt)
wgt = 1.
print(wgt((/0,dimsizes(wgt)-1/)))

z100 = 7.0*log(1000./p1)

mfe = f1->mflx_east(klev,:)
mfw = f1->mflx_west(klev,:)

if (dimsizes(klev) .eq. 2) then
  mf = 0.5*dim_avg((mfe+mfw)*conform_dims(dimsizes(mfe),wgt,1))
else
  mf = 0.5*dim_avg((mfe+mfw)*wgt)
end if
mf = mf*1.e3

print(z+"   "+mf+"   "+sqrt(10./mf))

