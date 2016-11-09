load "$UTIL/ncl/header.ncl"

begin

;--------------------------------------- parameters
 load "./ftn_extr-eqwave.ncl"

 expname = "era-int_f"
 dirname = "ERA-I-nr"
 y0 = 1979

 lat0 = 25.
 p0 = (/96.,80.,66.,54.,44.,36.,29.,23.,19.,15.,12.,10./)
 nt = 1 + 419

;--------------------------------------- read data
   fnam1 = "$DATD/tem/"+dirname(0)+"/epf/"+expname(0)+".epfsa_yz_recon0.nc"
   fnam2 = "$DATD/tem/"+dirname(0)+"/epf/"+expname(0)+".epf_yz_recon0.nc"
print(fnam1+"")

   rw_i = read_eqwave_rc0(  \
            (/"Kelvin","MRG-W","MRG-E","IG-W","IG-E","Rossby","Tide"/),  \
            (/"f_y","f_z","epd"/),(/fnam1,fnam2/),lat0,p0 )

   rw_i(0:1,:,:,:,:) = rw_i(0:1,:,:,:,:) / 6371.e3 * 1.e3

   dim = dimsizes(rw_i)
   dim(2) = nt

   rw = new(dim,float)
   rw(:,:,1:,:,:) = rw_i(:,:,:nt-2,:,:)
   rw@_FillValue = 1.e20
   rw(:,:,0,:,:) = rw@_FillValue

   t = y0 + ispan(0,nt-1,1)/12. + (1./24.)

   rw!2 = "t"
   rw&t = t

printVarSummary(rw)

   fout = addfile("$DATD/eq_wave/eq_vari/Tropical_Wave_EPflux.nc","c")
   fout->epf_y = rw(0,:,:,:,:)
   fout->epf_z = rw(1,:,:,:,:)
   fout->epfd  = rw(2,:,:,:,:)

   fout->epf_y@units = "mPa"
   fout->epf_z@units = "mPa"
   fout->epfd@units  = "m/s/day"

   globatt = True
   do iw=0, dimsizes(rw&wave)-1
     globatt@$("wave_"+iw)$ = rw&wave(iw)
   end do

   fileattdef(fout, globatt)

end

