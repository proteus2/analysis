load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.8

begin

;--------------------------------------- parameters
 load "./param0.ncl"
 load "./param1.ncl"
 load "./param_ra.ncl"
 load "./ftn_get_color_S-RIP.ncl"
 load "./ftn_read_spec.ncl"

 iwav = (/0,1/)
 k_a1   = param0_wavs_k_a1  (iwav)
 k_a2   = param0_wavs_k_a2  (iwav)
 ome_a1 = param0_wavs_ome_a1(iwav)
 ome_a2 = param0_wavs_ome_a2(iwav)
 p_max = 100
 p_min = 2

; ra = param0_ra
 ra = array_append_record(param0_ra_p,param0_ra_m,0)
 var = param0_var

 yyyymm1 = param0_yyyymm1
 yyyymm2 = param0_yyyymm2
 lat_b = param0_lat_b
 p_rng = param0_p_rng
 p_all = param_ra_plevs

 h_scale = param0_h_scale

 kmax = param1_kmax_log_draw
 omax = param1_omax_log_draw

 ndays_spec = param0_ndays_spec
 wgt_spec   = param0_wgt_spec

 l_sym_asym = False  ; True
 l_normalize = False

 ranames = get_ra_names(ra)
 varnames = get_ra_varnames_3d(ra,var)
 var0 = str_upper(varnames(0,1,:))
 
 foname = "var-rel-"+var0+"_z-y"+numeric2int(lat_b,1)+"_wav"+iwav(0)+""+iwav(1)
 wks = cr_wks("fig/"+foname,600,0)
 icol_s = get_col_srip(wks)
 colind = NhlNewColor(wks,0.6,0.6,0.6)

;--------------------------------------- read data
 vstr = ":F10:"+varnames(0,1,:)+"'"

 ne = dimsizes(ra)
 nfig = 1
 if ( l_sym_asym ) then
   nfig = 2
 end if

 no = toint( (2*ndays_spec)*omax )
 nk = kmax + 1

 nw = 3

 p_all = where( p_all .le. p_rng(0) .and. p_all .ge. p_rng(1) , p_all ,  \
                p_all@_FillValue )
 nz = new(ne,integer)
 do ie=0, ne-1
   ir = toint( ranames(ie,0) )
   nz(ie) = dimsizes( ind( .not. ismissing(p_all(ir,:)) ) )
 end do
 nz99 = max(nz)
 p99 = new((/ne,nz99/),float,1.e20)
 p99 = p99@_FillValue
 do ie=0, ne-1
   ir = toint( ranames(ie,0) )
   p99(ie,0:nz(ie)-1) = p_all( ir, ind( .not. ismissing(p_all(ir,:)) ) )
 end do
 z = 18.6 + h_scale*1.e-3*log(70./p99)

 vv = new((/nw,2,ne,nz99/),float,1.e20)
 vv = vv@_FillValue

do ie=0, ne-1
do k=0, nz(ie)-1

 p = p99(ie,k)

 vsp0 = read_fspec_ko(ra(ie),var,yyyymm1,yyyymm2,p,lat_b,kmax,no)
 vsp0 = vsp0*wgt_spec
 i0 = vsp0@i_k0

 dim3 = dimsizes(vsp0)
 vsp = new((/nfig,dim3(0),dim3(1)/),float)
 vsp = vsp0

 kwn = fspan(-kmax,kmax,kmax*2+1)
 ome = fspan(0,no/2,no/2+1) / tofloat(ndays_spec)
 dome = 1./tofloat(ndays_spec)


 if ( l_normalize .and. var0 .eq. "T" ) then
   stmp = get_ra_names(ra(ie))
   fname2 = "/data18/kyh/dat/clim/"+stmp(stmp@dirname)+"/"+stmp(stmp@file_h)+".t_yp.1979-2006.nc"
   f2 = addfile(fname2,"r")
   t0 = dim_avg(f2->$varnames(:,varnames@varname,0)$({p},{-lat_b:lat_b}))
   vsp = vsp / t0^2
   vstr(0) = vstr(0)+"/:F10:T:B:0:N:"
 end if
 if ( var0 .eq. "W" ) then
   vsp = vsp*(h_scale/(p*1.e2))^2 * 1.e6
 end if
 
;--------------------------------------- process 0
 vsp = vsp/dome  ; PS density

 vsp(:,:,i0) = vsp@_FillValue

 ; to integrate around omega = 0
; vsp(:,0,:) = vsp(:,0,:)*0.5

 vsp(:,no/2,:) = vsp@_FillValue

 vv_tot1 = dim_sum(dim_sum(vsp(:,1:,i0+1:)))*dome
 vv_tot2 = dim_sum(dim_sum(vsp(:,1:,:i0-1)))*dome
 vv_tot0 = vv_tot1 + vv_tot2 + 0.5*dim_sum(vsp(:,0,:))*dome

do iw=0, 1

 nome_a1 = ind_nearest_coord(ome_a1(iw),ome,0)
 nome_a2 = ind_nearest_coord(ome_a2(iw),ome,0)

 vv1 = dim_sum( dim_sum(vsp(:,nome_a1:nome_a2,i0+k_a1(iw):i0+k_a2(iw))) )*dome
 vv2 = dim_sum( dim_sum(vsp(:,nome_a1:nome_a2,i0-k_a2(iw):i0-k_a1(iw))) )*dome
 vv0 = 0.
 if (nome_a1 .ne. 0) then
   vv0 = 0.5*( dim_sum(vsp(:,0,i0+k_a1(iw):i0+k_a2(iw))) +  \
               dim_sum(vsp(:,0,i0-k_a2(iw):i0-k_a1(iw))) )*dome
   if (nome_a1 .gt. 1) then
     vv0 = vv0 +  \
           dim_sum( dim_sum(vsp(:,1:nome_a1-1,i0+k_a1(iw):i0+k_a2(iw))) +  \
                    dim_sum(vsp(:,1:nome_a1-1,i0-k_a2(iw):i0-k_a1(iw))) )*dome
   end if
 end if

 vv(0,iw,ie,k) = vv0/vv_tot0*100.
 vv(1,iw,ie,k) = vv1/vv_tot1*100.
 vv(2,iw,ie,k) = vv2/vv_tot2*100.

end do  ; iw

end do  ; k
end do  ; ie

;--------------------------------------- cn level
 ltck = (/1,1,1,1,1,1,1,1/)*1.0
 ldsh = (/0,0,0,0,0,1,1,1/)*16
 lcol = (/ icol_s@ei , icol_s@m , icol_s@m2 , icol_s@c , icol_s@j5 ,  \
           icol_s@ei , icol_s@m2 , icol_s@j5 /)

;--------------------------------------- ctr
 wi = 72
 he = 30
; xf = 12 + (/0,1,2,0,1,2/)*(wi+6)
; yf = 10 + (/0,0,0,1,1,1/)*(he+8)
 xf = 12 + (/0,1,0,1,0,1/)*(wi+6)
 yf = 10 + (/0,0,1,1,2,2/)*(he+8)

 xspc = 4
 wi1 = (wi-xspc)*0.5
 xfp = xf + wi1 + xspc

 mstr = ""  ; vstr+":F22:"
 xstr = vstr+":F22: "+"variance [%]"
 ystr = "p [hPa]"
; ystr = "PSD{"+vstr+":F22:}"

 ytckl = (/100,70,50,30,20,10,7,5,3,2,1/)
 ytckv = 18.6 + h_scale*1.e-3*log(70./ytckl)
 tmpx = (/300,200,100,90,80,70,60,50,40,30,20,10,9,8,7,6,5,4,3,2,1/)
 ymntckv = 18.6 + h_scale*1.e-3*log(70./tmpx)
 ymax = 18.6 + h_scale*1.e-3*log(70./p_min)
 ymin = 18.6 + h_scale*1.e-3*log(70./p_max)

 ; for plotting purpose
 do ie=0, ne-1
   do k=1, nz(ie)-1
     if ( .not. ismissing(z(ie,k)) ) then
     if ( z(ie,k) .gt. ymax ) then
       if (z(ie,k-1) .eq. ymax) then
         km = k
       else
         km = k+1
         vv(:,:,ie,k) = linint1(z(ie,:nz(ie)-1),vv(:,:,ie,:nz(ie)-1),False,ymax,0)
         z(ie,k) = ymax
       end if
       if (km .le. nz(ie)-1) then
         vv(:,:,ie,km:) = vv@_FillValue
         z(ie,km:) = z@_FillValue
         nz(ie) = km
       end if
       break
     end if
     end if
   end do
   do k=1, nz(ie)-1
     kk = nz(ie)-1-k
     if ( .not. ismissing(z(ie,kk)) ) then
     if ( z(ie,kk) .lt. ymin ) then
       if (z(ie,kk+1) .eq. ymin) then
         km = kk
       else
         km = kk-1
         vv(:,:,ie,kk) = linint1(z(ie,:nz(ie)-1),vv(:,:,ie,:nz(ie)-1),False,ymin,0)
         z(ie,kk) = ymin
       end if
       if (km .ge. 0) then
         vv(:,:,ie,:km) = vv@_FillValue
         z(ie,:km) = z@_FillValue
       end if
       break
     end if
     end if
   end do
 end do


do ip=0, nfig-1

 ; k < 0
 plot1 = xylineplot(wks,vv(2,0,:,:),z,"D",ldsh,"T",ltck*1.2,"C",lcol)
 xyaddline(plot1,vv(2,1,:,:),z,"D",ldsh,"T",ltck*0.8,"C",lcol)

 set_size(plot1,xf(ip),yf(ip),wi1,he)
 axis_range(plot1,0,100,ymin,ymax,False,False)
 axis_str(plot1,"",xstr,ystr)
 setvalues plot1
  "tmYLMode" : "Explicit"
  "tmYLValues" : ytckv
  "tmYLLabels" : ytckl
  "tmYLMinorValues" : ymntckv
  "tmXBMaxTicks" : 7
 end setvalues
 setvalues plot1
  "tiXAxisOffsetXF" : wi*0.00115
 end setvalues

 draw(plot1)


 ; k > 0
 plot2 = xylineplot(wks,vv(1,0,:,:),z,"D",ldsh,"T",ltck*1.2,"C",lcol)
 xyaddline(plot2,vv(1,1,:,:),z,"D",ldsh,"T",ltck*0.8,"C",lcol)

 set_size(plot2,xfp(ip),yf(ip),wi1,he)
 axis_range(plot2,0,100,ymin,ymax,False,False)
 axis_str(plot2,"","","")
 setvalues plot2
  "tmYLMode" : "Explicit"
  "tmYLValues" : ytckv
  "tmYLLabels" : ytckl
  "tmYLMinorValues" : ymntckv
  "tmXBMaxTicks" : 7
 end setvalues
 setvalues plot2
  "tmYLLabelsOn" : False
 end setvalues

 draw(plot2)


 dr_txt2(wks,mstr,xf(ip)+1,yf(ip)+2,22,4.5,0., 1,1)

end do

 lg = lgplot(wks,ranames(:4,1),"D",ldsh(:4),"T",ltck(:4),"C",lcol(:4),2)
 set_size_lg(lg,xf(0),yf(nfig-1)+he+6,wi*0.25,2.8*5)
 draw(lg)


 frame(wks)

end

