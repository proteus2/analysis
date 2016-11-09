
load "./param_ra.ncl"
load "./ftn_ra_info.ncl"

;-------------------------------------------------------------------------------

undef("read_fspec_1")
function read_fspec_1(ra:string, var:string,  \
             yyyymm1:integer, yyyymm2:integer, p:numeric, lat_b:numeric,  \
             kmax:integer, n_ome:integer, n_sym_asym:integer, f_tail:string)
local vsp, ra, var, yyyymm1, yyyymm2, p, lat_b, kmax, n_ome, n_sym_asym, f_tail,  \
      nsa, y1, y2, f_ko, s_sa, i_s, ra_names, dirname, expname, fvname, vname,  \
      l_precal, fdir, fname0, f0, vname0f, stmp, p0, lat0, dlat, lat, k0i, nk0i, itmp, wgty, iy1, iy2, wgty, no0, vspec,  \
      fname, f1, tmps, tmpa, n, iy, im1, im2, im, tmp1r, tmp1i, tmp2r, tmp2i,  \
      tmp, nt, no, nk
begin

 nsa = 1

 y1 = yyyymm1/100
 y2 = yyyymm2/100
 f_ko = "_ko"
 if ( n_sym_asym .ne. 0 ) then
   f_ko = "_wave_ko"
   s_sa = (/"s","a"/)
   i_s = 1
   if ( n_sym_asym .eq. 3 ) then
     nsa = 2
   end if
 end if
 f_ko = f_ko+f_tail

 ra_names = get_ra_names(ra)
 dirname = ra_names(ra_names@dirname)
 expname = ra_names(ra_names@file_h )
 delete(ra_names)

 ra_names = get_ra_varnames(ra,var)
 fvname = ra_names(ra_names@fvarname)
 vname  = ra_names(ra_names@varname )
 if ( str_upper(ra_names(1)) .eq. "V" ) then
   s_sa = (/"a","s"/)
   i_s = -1
 end if
 delete(ra_names)

 l_precal = False
 fdir   = "$DATD/fspec/"+dirname+"/fcoef_ko"
 fname0 = fdir+"/"+y1+"/"+expname+".fft_"+fvname+f_ko+"."+y1+".02.nc"
; if ( yyyymm1 .eq. 197902 .and. yyyymm2 .eq. 200612 ) then
 if ( yyyymm1 .eq. 198101 .and. yyyymm2 .eq. 201012 ) then
   l_precal = True
   fdir   = "$DATD/fspec/"+dirname+"/fcoef_ko/avg_81-10"
   fname0 = fdir+"/"+expname+".fc2_"+fvname+f_ko+".nc"
 end if

 f0 = addfile(fname0,"r")

 if ( l_precal ) then
   vname0f = "fc2_"+vname
 else
   vname0f = "fcr_"+vname
 end if
 if ( n_sym_asym .ne. 0 ) then
   vname0f = vname0f+"_s"
 end if
 stmp = f0->$vname0f$!0
 p0 = f0->$stmp$
 if ( abs(p0({p}) / tofloat(p) - 1.) .gt. 0.02 ) then  ; ~ 130 m
   vsp = p0({p})
   vsp@status = "p_mismatch"
   itmp = ind_nearest_coord(p,p0,0)
   if ( p0(itmp) .gt. p ) then
     vsp@p_lower = p0(itmp)
     vsp@p_upper = min(p0(max((/itmp-1,0/)):min((/itmp+1,dimsizes(p0)-1/))))
   else
     vsp@p_upper = p0(itmp)
     vsp@p_lower = max(p0(max((/itmp-1,0/)):min((/itmp+1,dimsizes(p0)-1/))))
   end if
   print( "["+expname+"]" )
   print( "   p asked: "+p+"  /  read: "+p0(itmp)+"  mismatch, interpolate." )
   return(vsp)
 end if
 print( "["+expname+"]" )
 print( "   p asked: "+p+"  /  read: "+f0->$stmp$({p}) )
 stmp = f0->$vname0f$!1
 lat0 = f0->$stmp$
 iy1 = ind( lat0 .eq. min( lat0({-lat_b+0.001:lat_b-0.001}) ) )
 iy2 = ind( lat0 .eq. max( lat0({-lat_b+0.001:lat_b-0.001}) ) )
 if (iy1 .lt. iy2) then
   iy1 = iy1 - 1
   iy2 = iy2 + 1
 else
   iy1 = iy1 + 1
   iy2 = iy2 - 1
 end if
 iy1 = min( (/ max( (/iy1,0/) ), dimsizes(lat0)-1 /) )
 iy2 = min( (/ max( (/iy2,0/) ), dimsizes(lat0)-1 /) )
 lat = f0->$stmp$(iy1:iy2)
 lat = where( abs(lat) .lt. 1.e-3, 0., lat )  ; MERRA2 has lat = -1.79751e-13 !
 stmp = f0->$vname0f$!3
 k0i = f0->$stmp$
 nk0i = dimsizes(k0i)

 itmp = dimsizes(lat)
 wgty = lat*0. + 1.
 wgty(1:itmp-2) = 0.5*(lat(2:itmp-1) - lat(0:itmp-3))  ; grid size
 wgty(0) = lat(1) - lat(0)
 wgty(itmp-1) = lat(itmp-1) - lat(itmp-2)
 wgty = where( lat - 0.5*wgty .lt. -lat_b , (lat+0.5*wgty) + lat_b , wgty )
 wgty = where( lat + 0.5*wgty .gt. lat_b  , lat_b - (lat-0.5*wgty) , wgty )
 wgty = where( wgty .lt. 0. , 0. , wgty )  ; for useless grids
 if ( n_sym_asym .ne. 0 .and. lat(0)*lat(itmp-1) .ge. 0. ) then
   itmp = minind(abs(lat))
   wgty(itmp) = wgty(itmp)*0.5 + abs(lat(itmp))
 end if
 wgty = wgty/max(wgty)
 
 no0 = dimsizes( f0->ome_fr )

 vspec = new((/nsa,no0,kmax+1/),float)


if ( l_precal ) then

 fname = fdir+"/"+expname+".fc2_"+fvname+f_ko+".nc"
 f1 = addfile(fname,"r")
 if ( n_sym_asym .eq. 0 ) then
   vspec = dim_avg_wgt_n(  \
       f1->$("fc2_"+vname)$({p},iy1:iy2,:,:kmax), wgty, 0, 0 )*2.
   vspec(:,:,0) = vspec(:,:,0)*0.5
 else
   if ( n_sym_asym .eq. 1 .or. n_sym_asym .eq. 3 ) then
     vspec(0,:,:) = dim_avg_wgt_n(  \
         f1->$("fc2_"+vname+"_"+s_sa(0))$({p},iy1:iy2,:,:kmax), wgty, 0, 0 )*2.
     vspec(0,:,0) = vspec(0,:,0)*0.5
   end if
   if ( n_sym_asym .eq. 2 .or. n_sym_asym .eq. 3 ) then
     vspec(n_sym_asym-2,:,:) = dim_avg_wgt_n(  \
         f1->$("fc2_"+vname+"_"+s_sa(1))$({p},iy1:iy2,:,:kmax), wgty, 0, 0 )*2.
     vspec(n_sym_asym-2,:,0) = vspec(n_sym_asym-2,:,0)*0.5
   end if
 end if

else

 vspec = 0.
 n = 0
 do iy=y1, y2
   im1 = 1
   im2 = 12
   if ( iy .eq. y1 ) then
     im1 = yyyymm1 - y1*100
   end if
   if ( iy .eq. y2 ) then
     im2 = yyyymm2 - y2*100
   end if
   do im=im1, im2
     fname = fdir+"/"+iy+"/"+expname+".fft_"+fvname+f_ko+"."+ \
             iy+"."+sprinti("%2.2i",im)+".nc"
     f1 = addfile(fname,"r")

     tmp1r = f1->$("fcr_"+vname)$({p},iy1:iy2,:,:kmax)
     tmp1i = f1->$("fci_"+vname)$({p},iy1:iy2,:,:kmax)
     if ( n_sym_asym .eq. 0 ) then
       vspec(0,:,:) = vspec(0,:,:) + dim_avg_wgt_n( tmp1r^2 + tmp1i^2, wgty, 0, 0 )*2.
       vspec(0,:,0) = vspec(0,:,0) - dim_avg_wgt_n( tmp1r(:,:,0)^2 + tmp1i(:,:,0)^2, wgty, 0, 0 )
     else
       tmp2r = tmp1r(iy2-iy1:0,:,:)
       tmp2i = tmp1i(iy2-iy1:0,:,:)*tofloat(i_s)
       if ( n_sym_asym .eq. 1 .or. n_sym_asym .eq. 3 ) then
         tmp = (0.5*(tmp1r + tmp2r))^2 + (0.5*(tmp1i + tmp2i))^2
         vspec(0,:,:) = vspec(0,:,:) + dim_avg_wgt_n(tmp,wgty,0,0)*2.
         vspec(0,:,0) = vspec(0,:,0) - dim_avg_wgt_n(tmp(:,:,0),wgty,0,0)
       end if
       if ( n_sym_asym .eq. 2 .or.  n_sym_asym .eq. 3 ) then
         tmp = (0.5*(tmp1r - tmp2r))^2 + (0.5*(tmp1i - tmp2i))^2
         vspec(n_sym_asym-2,:,:) = vspec(n_sym_asym-2,:,:) + dim_avg_wgt_n(tmp,wgty,0,0)*2.
         vspec(n_sym_asym-2,:,0) = vspec(n_sym_asym-2,:,0) - dim_avg_wgt_n(tmp(:,:,0),wgty,0,0)
       end if
       delete( [/tmp2r,tmp2i,tmp/] )
     end if
     delete( [/tmp1r,tmp1i/] )
     n = n + 1
   end do
 end do
 nt = n
 vspec = vspec/tofloat(nt)

end if

 no = n_ome
 nk = kmax + 1

 vsp = new((/nsa,no/2+1,nk*2-1/),float)
 vsp@_FillValue = 1.e20
 vsp@i_k0       = nk-1
 vsp = 0.

 vsp(:,:,0:nk-2) = vspec(:,0:no/2,nk-1:1)
 vsp(:,0,nk:) = vspec(:,0,1:)  ; duplicated for visualization
 vsp(:,1:no/2,nk:) = vspec(:,no0-1:no0-no/2,1:)

 ; to conserve the sum  ----------------
 ; comment the below for visualization
; vsp(:,0,:) = vsp(:,0,:)*0.5
; if (no .eq. no0) then
;   vsp(:,no/2,:) = vsp(:,no/2,:)*0.5
; end if
 ;--------------------------------------

 ; k = 0
 vsp(:,:,nk-1) = vspec(:,0:no/2,0) + vspec(:,no0-1:no0-no/2-1,0)

 vsp(:,0,nk-1) = vsp@_FillValue

 vsp@status = "0"

 delete(vspec)

return(vsp)
end

;-------------------------------------------------------------------------------

undef("read_fspec_5d")
function read_fspec_5d(ra:string, var:string, yyyymm1:integer, yyyymm2:integer,  \
                       p:numeric, lat_b:numeric, kmax:integer, n_ome:integer,  \
                       n_sym_asym:integer, f_tail:string)
local vsp5d, ra, var, yyyymm1, yyyymm2, p, lat_b, kmax, n_ome, n_sym_asym, f_tail,  \
      nr, nv, iv, vsp_1, p1, p2, vsp_l, vsp_u, zn, zn1, zn2, dim3d
begin

 nr = dimsizes(ra)
 nv = dimsizes(var)

 do ir=0, nr-1
 do iv=0, nv-1
   vsp_1 = read_fspec_1(ra(ir),var(iv),yyyymm1,yyyymm2,p,lat_b,kmax,n_ome,  \
                        n_sym_asym,f_tail)
   if ( vsp_1@status .eq. "p_mismatch" ) then
     p1 = vsp_1@p_lower
     p2 = vsp_1@p_upper
     delete(vsp_1)
     vsp_l = read_fspec_1(ra(ir),var(iv),yyyymm1,yyyymm2,p1,lat_b,kmax,n_ome,  \
                          n_sym_asym,f_tail)
     vsp_u = read_fspec_1(ra(ir),var(iv),yyyymm1,yyyymm2,p2,lat_b,kmax,n_ome,  \
                          n_sym_asym,f_tail)
     zn  = log(1000./tofloat(p ))
     zn1 = log(1000./tofloat(p1))
     zn2 = log(1000./tofloat(p2))
     vsp_1 = ((zn2-zn)*vsp_l + (zn-zn1)*vsp_u)/(zn2-zn1)
     vsp_1@_FillValue = vsp_l@_FillValue
     vsp_1@i_k0       = vsp_l@i_k0
     delete( [/vsp_l,vsp_u/] )
   end if
   if ( ir .eq. 0 .and. iv .eq. 0 ) then
     dim3d = dimsizes(vsp_1)
     vsp5d = new((/nr,nv,dim3d(0),dim3d(1),dim3d(2)/),float)
   end if
   vsp5d(ir,iv,:,:,:) = vsp_1
   delete( vsp_1 )
 end do
 end do

return(vsp5d)
end

;-------------------------------------------------------------------------------

undef("read_fspec")
function read_fspec(ra:string, var:string, yyyymm1:integer, yyyymm2:integer,  \
                    p:numeric, lat_b:numeric, kmax:integer, n_ome:integer,  \
                    n_sym_asym:integer, f_tail:string)
local vsp, ra, var, yyyymm1, yyyymm2, p, lat_b, kmax, n_ome, n_sym_asym, f_tail,  \
      vsp5d, dim5d, dim
begin

 vsp5d = read_fspec_5d(ra,var,yyyymm1,yyyymm2,p,lat_b,kmax,n_ome,n_sym_asym,  \
                       f_tail)

 dim5d = dimsizes(vsp5d)
 dim = dim5d( ind( dim5d .ne. 1 ) )
 vsp = onedtond( ndtooned(vsp5d), dim )
 vsp@_FillValue = vsp5d@_FillValue
 vsp@i_k0       = vsp5d@i_k0
 vsp@dim_k   = dimsizes(dim)-1
 vsp@dim_ome = dimsizes(dim)-2

return(vsp)
end

;-------------------------------------------------------------------------------

undef("read_fspec_ko")
function read_fspec_ko(ra, var, yyyymm1, yyyymm2, p, lat_b, kmax, n_ome)
local vsp
begin

  vsp = read_fspec( ra, var, yyyymm1, yyyymm2, p, lat_b, kmax, n_ome, 0, "" )

return(vsp)
end

;-------------------------------------------------------------------------------

undef("read_fspec_ko_sa")
function read_fspec_ko_sa(ra, var, yyyymm1, yyyymm2, p, lat_b, kmax, n_ome,  \
                          n_sym_asym)
local vsp, f_tail
begin

  f_tail = ""
  if (kmax .gt. 20) then
    f_tail = "-p"+p
  end if
  vsp = read_fspec( ra, var, yyyymm1, yyyymm2, p, lat_b, kmax, n_ome,  \
                    n_sym_asym, f_tail )
 
return(vsp)
end

