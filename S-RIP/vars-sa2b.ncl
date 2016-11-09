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
 load "./ftn_read_spec.ncl"
 load "./ftn_disp.ncl"

 iter_k = 23
 iter_o = 7
 l_cellfill = False  ; True

 ra_p = array_append_record( param0_ra_p, (/"JRA55C"/), 0 )
 ra_m = array_append_record( param0_ra_m, (/"JRA55C_m"/), 0 )
 ra = array_append_record( ra_p, ra_m, 0 )
; var = param0_var

 yyyymm1 = param0_yyyymm1
 yyyymm2 = param0_yyyymm2
 lat_b = param0_lat_b
 p     = param0_p

 nbv     = param0_nbv
 h_scale = param0_h_scale

 kmax_draw = 10
 omax_draw = 0.5
 kmax = kmax_draw + iter_k
 omax = 1.

 ndays_spec = param0_ndays_spec
 wgt_spec   = param0_wgt_spec

 l_normalize = False

 ranames = get_ra_names(ra)

vars = (/"T","u","v","w"/)
nv = dimsizes(vars)
do iv=0, nv-1
 var = vars(iv)

 varnames = get_ra_varnames_3d(ra,var)
 var0 = str_upper(varnames(0,1,:))
 
 foname = "vars-sa2b-p"+numeric2int(p,1)+"-y"+numeric2int(lat_b,1)
 wks = cr_wks("fig/"+foname,600,0)
 colind = NhlNewColor(wks,0.6,0.6,0.6)

;--------------------------------------- read data
 vstr = ":F10:"+varnames(0,1,:)+"'"

 ne = dimsizes(ra)
 nfig = 2

 no = toint( (2*ndays_spec)*omax )
 nk = kmax + 1

do n_sym_asym=1, 2

 vsp0 = read_fspec_ko_sa(ra,var,yyyymm1,yyyymm2,p,lat_b,kmax,no,n_sym_asym)
 vsp0 = vsp0*wgt_spec
 i0 = vsp0@i_k0

 if ( .not. isdefined("vsp") ) then
   dim3 = dimsizes(vsp0)
   vsp = new((/dim3(0),nfig,dim3(1),dim3(2)/),float)

   kwn = fspan(-kmax,kmax,kmax*2+1)
   ome = fspan(0,no/2,no/2+1) / tofloat(ndays_spec)
   dome = 1./tofloat(ndays_spec)
 end if
 vsp(:,n_sym_asym-1,:,:) = vsp0

end do

 if ( l_normalize .and. var0 .eq. "T" ) then
   stmp = get_ra_names_2d(ra)
   fname2 = "/data18/kyh/dat/clim/"+stmp(:,stmp@dirname)+"/"+stmp(:,stmp@file_h)+".t_yp.1979-2006.nc"
   do ie=0, ne-1
     f2 = addfile(fname2(ie),"r")
     t0 = dim_avg(f2->$varnames(:,varnames@varname,0)$({p},{-lat_b:lat_b}))
     vsp(ie,:,:,:) = vsp(ie,:,:,:) / t0^2
   end do
   vstr(0) = vstr(0)+"/:F10:T:B:0:N:"
 end if
 if ( var0 .eq. "W" ) then
   vsp = vsp*(h_scale/(p*1.e2))^2 * 1.e6
 end if
 
;--------------------------------------- process 0
 vsp = vsp/dome  ; PS density

 vsp(:,:,:,i0) = 0.5*(vsp(:,:,:,i0-1) + vsp(:,:,:,i0+1))
 vspb0 = vsp

 do i=1, iter_k
   vspb = wgt_runave(vspb0,(/0.25,0.5,0.25/),0)
   vspb0 = vspb
 end do
 do i=1, iter_o
   vspb = wgt_runave_n(vspb0,(/0.25,0.5,0.25/),0,2)
   vspb(:,:,0,:) = 0.5*vspb0(:,:,0,:) +  \
                   0.25*(vspb0(:,:,1,:) + vspb0(:,:,1,::-1))
   vspb0 = vspb
 end do
; vspt = vsp/vspb

; vsp(:,:,:,i0) = vsp@_FillValue  ; k = 0

print(max(vsp))

 vspm = dim_avg_n(vsp,0)

; vsp(:,:,0,: ) = vsp@_FillValue
; vsp(:,:,0,i0) = vsp@_FillValue

 kwn2 = fspan(-nk+1,nk-1,nk*2-1)
 ome2 = fspan(0,no/2,no/2+1) / 90.
 i0_2 = (nk*2-1)/2

 dispc = disp_eq(kwn2,ome2,False,nbv)
; dispc@_FillValue = 999.
 dispc = where(dispc .le. 0., 1.e-20, dispc)
 dispc = nbv^2/9.801/( (2.*3.141592/(dispc*1.e3))^2 + 0.25/h_scale^2 )

 vsp_k = (/ vsp(:,0,:,:) /)
 vsp_k!1 = "ome"
 vsp_k&ome = ome
 vsp_k(:,:,:i0) = 0.
 vsp_k(:,:,i0+11:) = 0.
 vsp_k(:,{:0.04999},:) = 0.
 vsp_k(:,{0.5001:},:) = 0.
 do j=0, dim3(0)-1
   vsp_k(j,:,:) = where(dispc(0,0,:,:) .ge. 8. .and. dispc(0,0,:,:) .le. 240., vsp_k(j,:,:), 0.)
 end do
 vv_k = dim_sum_n(dim_sum(vsp_k),1)*dome

 vsp_m = vsp(:,1,:,:)
 vsp_m(:,:,i0+1:) = 0.
 vsp_m(:,:,:i0-11) = 0.
 do j=0, dim3(0)-1
   vsp_m(j,:,:) = where(dispc(1,0,:,:) .ge. 8. .and. dispc(1,0,:,:) .le. 70., vsp_m(j,:,:), 0.)
 end do
 vv_m = dim_sum_n(dim_sum(vsp_m),1)*dome

 vsp_j = (/ vsp(:,0,:,:) /)
 vsp_j!1 = "ome"
 vsp_j&ome = ome
 vsp_j(:,:,:i0) = 0.
 vsp_j(:,:,i0+6:) = 0.
 vsp_j(:,{0.05001:},:) = 0.
 vsp_j(:,0,:) = 0.  ; check
 vv_j = dim_sum_n(dim_sum(vsp_j),1)*dome

 vspb_k = (/ vspb(:,0,:,:) /)
 vspb_k!1 = "ome"
 vspb_k&ome = ome
 vspb_k(:,:,:i0) = 0.
 vspb_k(:,:,i0+11:) = 0.
 vspb_k(:,{:0.04999},:) = 0.
 vspb_k(:,{0.5001:},:) = 0.
 do j=0, dim3(0)-1
   vspb_k(j,:,:) = where(dispc(0,0,:,:) .ge. 8. .and. dispc(0,0,:,:) .le. 240., vspb_k(j,:,:), 0.)
 end do
 vvb_k = dim_sum_n(dim_sum(vspb_k),1)*dome

 vspb_m = vspb(:,1,:,:)
 vspb_m(:,:,i0+1:) = 0.
 vspb_m(:,:,:i0-11) = 0.
 do j=0, dim3(0)-1
   vspb_m(j,:,:) = where(dispc(1,0,:,:) .ge. 8. .and. dispc(1,0,:,:) .le. 70., vspb_m(j,:,:), 0.)
 end do
 vvb_m = dim_sum_n(dim_sum(vspb_m),1)*dome

 vspb_j = (/ vspb(:,0,:,:) /)
 vspb_j!1 = "ome"
 vspb_j&ome = ome
 vspb_j(:,:,:i0) = 0.
 vspb_j(:,:,i0+6:) = 0.
 vspb_j(:,{0.05001:},:) = 0.
 vspb_j(:,0,:) = 0.  ; check
 vvb_j = dim_sum_n(dim_sum(vspb_j),1)*dome

 vva_k = (vv_k - vvb_k)
 vva_m = (vv_m - vvb_m)
 vva_j = (vv_j - vvb_j)

 vva_k = vva_k/abs(avg(vva_k(6:8)))
 vva_m = vva_m/abs(avg(vva_m(6:8)))
 vva_j = vva_j/abs(avg(vva_j(6:8)))

 ; for plotting

 if ( .not. isdefined("vva") ) then
   vva = new((/3,4,ne/),float)
 end if
 vva(:,iv,:) = (/vva_k,vva_m,vva_j/)

end do

;--------------------------------------- cn level
 ldsh = (/0  ,16 ,14 /)
 ltck = (/0.5,0.6,0.6/)
 lcol = (/2,4,6,3/)
 imrk = (/1,1,1,1,1/)*2

;--------------------------------------- ctr
 wi = 60
 he = 15
 xf = 15 + (/0,0,0/)*(wi+12)
 yf = 10 + (/0,1,2/)*(he+12)

 txtstr = ranames(:,1)

 mstr = (/"Kelvin waves (T, u, w)","MRG waves (T, u, v, w)","MJO (T, u, v, w)"/)

do iw=0, 1 ;2

  lcol2 = new(nv*ne,integer)
  do ie=0, ne-1
    lcol2(ie*nv:(ie+1)*nv-1) = lcol(:)
  end do

 dumpl = new(ne*nv,graphic)
 respl = True
 respl@gsLineThicknessF = 1.0

 ii=0
 do ie=0, ne-1
 do iv=0, nv-1
   if (ii .eq. 0) then
     plot = xylineplot(wks,ie+(iv-nv/2)*0.12,(/vva(iw,iv,ie),0/),"D",0,"T",1,"C",-1)
   else
     xyaddline(plot,ie+(iv-nv/2)*0.12,(/vva(iw,iv,ie),0/),"D",0,"T",1,"C",-1)
   end if
   xpl = ie+(iv-nv/2)*0.12 + (/0,0,1,1/)*0.10
   if ( .not. ( iw .eq. 0 .and. iv .eq. 2 ) ) then
     ypl = (/0,1,1,0/)*vva(iw,iv,ie)
     respl@gsLineColor = lcol(iv)
     if (iw .eq. 0) then
       if (iv .lt. 2) then
         xpl = xpl + 0.06
       else
         xpl = xpl - 0.06
       end if
     end if
     dumpl(ii) = gsn_add_polyline(wks,plot,xpl,ypl,respl)
   end if 
   ii = ii + 1
 end do
 end do

;  xyaddmarker(plot,imrk,"S",0.8,"T",5.0,"C",lcol2)

  do ie=0, ne-2
    xyaddline(plot,ie+0.5,(/-999,999/),"D",0,"T",0.5,"C",colind)
  end do
  xyaddline(plot,ispan(-1,ne,1),0,"D",0,"T",0.5,"C",1)
  xyaddline(plot,(/-999,999/),1,"D",0,"T",0.5,"C",1)

 set_size(plot,xf(iw),yf(iw),wi,he)
 axis_range(plot,-0.5,ne-0.5,0.4,1.2,False,False)
 axis_str(plot,"","","")
 setvalues plot
  "tmXBMode"   : "Explicit"
  "tmXBValues" : ispan(0,ne-1,1)
  "tmXBLabels" : txtstr
  "tmXBLabelAngleF" : 30
  "tmXBMinorPerMajor" : 0
  "tmXTOn" : False
;  "tmXBMajorLengthF" : -0.002
 end setvalues

 draw(plot)

 dr_txt2(wks,mstr(iw),xf(iw)+1,yf(iw)-2,22,4.0,0., 1,1)

end do

 frame(wks)

end

