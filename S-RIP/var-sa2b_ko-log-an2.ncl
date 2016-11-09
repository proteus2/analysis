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
 l_ano = False

; ra = param0_ra_m
 ra = array_append_record( param0_ra_m, (/"JRA55C_m"/), 0 )
 var = param0_var

 yyyymm1 = param0_yyyymm1
 yyyymm2 = param0_yyyymm2
 lat_b = param0_lat_b
 p     = param0_p

 nbv     = param0_nbv
 h_scale = param0_h_scale

 kmax_draw = param1_kmax_draw
 omax_draw = param1_omax_draw
 kmax = kmax_draw + iter_k
 omax = 1.

 ndays_spec = param0_ndays_spec
 wgt_spec   = param0_wgt_spec

; lzs = param1_lz_draw
lzs = (/ (/ (/2.,5.,10./),    \
                       (/2.5,5.0,10./),    \
                       (/2.,10.,60./) /), \
                    (/ (/2.,4.,8./),    \
                       (/2.5,5.0,10./),    \
                       (/2.,10.,60./) /) /)  ; [km]

 l_normalize = False

 ranames = get_ra_names(ra)
 varnames = get_ra_varnames_3d(ra,var)
 var0 = str_upper(varnames(0,1,:))
 
 foname = "var-"+var0+"-sa2b_ko-log-p"+numeric2int(p,1)+"-y"+numeric2int(lat_b,1)+"-an_m"
 wks = cr_wks("fig/"+foname,600,0)
; gsn_define_colormap(wks,"prcp_3")
 gsn_define_colormap(wks,"StepSeq25")
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
 vspt = vsp/vspb

; vsp(:,:,:,i0) = vsp@_FillValue  ; k = 0

; convert to energe-preserving form
 vsp = vsp*conform_dims(dimsizes(vsp),abs(ome)*log(10.),2)
 vsp = vsp*conform_dims(dimsizes(vsp),abs(kwn)*log(10.),3)

print(max(vsp))

 vsp(:,:,0,: ) = 1.e-20
 vsp(:,:,:,i0) = 1.e-20
 vsp(:,:,dimsizes(ome)-1,:) = 1.e-20

 vspm = dim_avg_n(vsp,0)
; vspa = vsp - conform_dims(dimsizes(vsp),vspm,(/1,2,3/))
 vspa = (vsp/conform_dims(dimsizes(vsp),vspm,(/1,2,3/)) - 1.)*100.
 vspa = runave_n(vspa,3,0,2)
; vspa(:,:,:,i0) = 0.5*(vspa(:,:,:,i0-1) + vspa(:,:,:,i0+1))

 vsp(:,:,0,: ) = vsp@_FillValue
 vsp(:,:,0,i0) = vsp@_FillValue


 ; for plotting
 dim4 = dimsizes(vsp)
 dim4(0) = dim4(0)+1
 vspp = new(dim4,float)
 vspp(:ne-1,:,:,:) = vsp
 vspp(ne,:,:,:) = vspm

 nk10 = nk*2*10
 no10 = no/2*10

 kwn2 = fspan(-nk,nk,nk10+1)
 ome2 = fspan(0,no/2,no10+1) / 90.
 i0_2 = (nk10+1)/2

 dispc = disp_eq(kwn2,ome2,False,nbv)

;--------------------------------------- cn level
; cnfil = (/5,6,7,10,12,14,15,16,17,19,21/)
; cnfil = (/26,24,21,19,17,16,14,12,2,5,11/)
 cnfil = (/26,24,21,19,17,16,14,12,2,5,11,7,10/)


 fact = 0
 if ( var0 .eq. "T" ) then
;   cnlev = (/0.02,0.05,0.1,0.2,0.5,1,2,5,10,20/)*0.1
   cnlev = (/ ispan(1,12,1)*0.15,  \
              ispan(1,12,1)*0.05 /)
   lbtit = "K:S:2:N:"
   if ( l_normalize ) then
     fact = -6
     cnlev = (/0.02,0.05,0.1,0.2,0.5,1,2,5,10,20/)
     lbtit = ""
   end if
 end if
 if ( var0 .eq. "U" ) then
;   cnlev = (/0.02,0.05,0.1,0.2,0.5,1,2,5,10,20/)
   cnlev = (/ ispan(1,12,1)*2.5,  \
              ispan(1,12,1)*1.0 /)
   lbtit = "m:S:2:N: s:S:-2:N:"
 end if
 if ( var0 .eq. "V" ) then
;   cnlev = (/0.02,0.05,0.1,0.2,0.5,1,2,5,10,20/)
   cnlev = (/ ispan(1,12,1)*2.5,  \
              ispan(1,12,1)*2.5 /)
   lbtit = "m:S:2:N: s:S:-2:N:"
 end if
 if ( var0 .eq. "W" ) then
;   cnlev = (/0.1,0.2,0.5,1,2,5,10,20,50,100/)
   cnlev = (/ ispan(1,12,1)*1.0,  \
              ispan(1,12,1)*1.0 /)
   lbtit = "mm:S:2:N: s:S:-2:N:"
 end if
 lbt_hd = ""
 if (fact .ne. 0) then
   lbt_hd = "10:S:"+fact+":N: "
 end if
 lbtit = "[ "+lbt_hd+lbtit+" ]"

 ldsh = (/0  ,16 ,14 /)
 ltck = (/0.5,0.6,0.6/)
 lcol = 1

;--------------------------------------- ctr
 wi = 32
 he = 16  ; 32
; xf = 12 + (/0,1,2,0,1,2/)*(wi+6)
; yf = 10 + (/0,0,0,1,1,1/)*(he+8)
; xf = 12 + (/0,1,0,1,0,1/)*(wi+6)
; yf = 10 + (/0,0,1,1,2,2/)*(he+7)
 yf = 10 + (/0,1,2,2,2,2/)*(he+7)

 xspc = 2
 wi1 = (wi-xspc)*0.5
; xfp = xf + wi1 + xspc

; xf = 12 + (/0,1,0,1,0,1/)*(wi1+3)
 xf = 12 + (/0,0,0,1,0,1/)*(wi1+3)

 xfa = xf + wi+7

 mstr = new((/nfig,ne+1/),string)
 xstr = new(ne+1,string)
 ystr = new(ne+1,string)

 mstr = ""
 mstr(0,:ne-1) = ranames(:,1)
 mstr(0,ne) = "mean"
 xstr = ""
 ystr = ""
; xstr = where(yf(:ne) .eq. yf(ne) .and. xf(:ne) .eq. xf(0), "zonal wavenumber", xstr)
; ystr = where(xf(:ne) .eq. xf(0), "frequency [cyc day:S:-1:N:]", ystr)
 xstr(:ne-1) = where(yf(:ne-1) .eq. yf(ne-1) .and. xf(:ne-1) .eq. xf(0), "zonal wavenumber", xstr(:ne-1))
 ystr(:ne-1) = where(xf(:ne-1) .eq. xf(0), "frequency [cyc day:S:-1:N:]", ystr(:ne-1))

; xtckv = (/1,10,100/)
 xtckv = (/1,2,4,8,16/)
; xmntckv = (/1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400/)
 xmntckv = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/)
 ytckv = (/0.01,0.1,1/)
 ymntckv = (/1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400/)*0.01


;do ie=0, ne
do ie=0, ne-1

 ; k < 0
 ip = 1
 ctr1 = cnshadeplot(wks,abs(kwn(:i0-1)),ome(1:),vspp(ie,ip,1:,:i0-1)/(10.^fact),cnlev(ip,:),cnfil,"")
 cnaddline(ctr1,wks,abs(kwn(:i0-1)),ome(1:),vspt(ie,ip,1:,:i0-1),(/1.1,1.5,2.0/),0,"D",0,"T",0.7,"C",0)
if ( l_ano ) then
 cnaddline(ctr1,wks,abs(kwn(:i0-1)),ome(1:),vspa(ie,ip,1:,:i0-1),ispan(-200,200,10),1,"D",(/0,0,1/),"T",0.7,"C",(/1,-1,1/))
end if
;if ( n_sym_asym .eq. 0 .or. n_sym_asym .eq. 2 ) then
 do iw=0, 0, 2
   cnaddline(ctr1,wks,abs(kwn2(:i0_2-1)),ome2(1:),dispc(1,iw,1:,:i0_2-1),lzs(1,iw,:),0,"D",ldsh(2),"T",ltck(2),"C",lcol)
 end do
;end if
; do iw=0, 2, 2
;   if ( ip .eq. 0 .and. iw .eq. 0 ) then
;     continue
;   end if
;   cnaddline(ctr1,wks,abs(kwn2(:i0_2-1)),ome2(1:),dispc(ip,iw,1:,:i0_2-1),lzs(ip,iw,:),0,"D",ldsh(iw),"T",ltck(iw),"C",lcol)
; end do
; do iw=2, 2
;   cnaddline(ctr1,wks,abs(kwn2(:i0_2-1)),ome2(1:),dispc(0,iw,1:,:i0_2-1),lzs(0,iw,:),0,"D",ldsh(0),"T",ltck(0),"C",lcol)
; end do

if ( l_cellfill ) then
 setvalues ctr1
  "cnFillMode" : "CellFill"
 end setvalues
end if

 set_size(ctr1,xfa(ie),yf(ie),wi1,he)
 axis_range(ctr1,kmax_draw,1,0.01,omax_draw+0.001,True,True)
 axis_str(ctr1,"",xstr(ie),"")
 setvalues ctr1
  "tmXBMode" : "Explicit"
  "tmXBValues" : xtckv
  "tmXBLabels" : xtckv*(-1)
  "tmXBMinorValues" : xmntckv
  "tmYLMode" : "Explicit"
  "tmYLValues" : ytckv
  "tmYLLabels" : ytckv
  "tmYLMinorValues" : ymntckv
 end setvalues
 setvalues ctr1
  "tiXAxisOffsetXF" : wi*0.0012
 end setvalues
 if (xfa(ie) .ne. xfa(0)) then
 setvalues ctr1
  "tmYLLabelsOn" : False
 end setvalues
 end if

 cnfinalize(ctr1,wks)

 draw(ctr1)


 ; k > 0
 ip = 0
 ctr2 = cnshadeplot(wks,kwn(i0+1:),ome(1:),vspp(ie,ip,1:,i0+1:)/(10.^fact),cnlev(ip,:),cnfil,"")
 cnaddline(ctr2,wks,kwn(i0+1:),ome(1:),vspt(ie,ip,1:,i0+1:),(/1.1,1.5,2.0/),0,"D",0,"T",0.7,"C",0)
if ( l_ano ) then
 cnaddline(ctr2,wks,kwn(i0+1:),ome(1:),vspa(ie,ip,1:,i0+1:),ispan(-200,200,10),1,"D",(/0,0,1/),"T",0.7,"C",(/1,-1,1/))
end if
; do iw=0, 0
;   cnaddline(ctr2,wks,kwn2(i0_2+1:),ome2(1:),dispc(ip,iw,1:,i0_2+1:),lzs(ip,iw,:),0,"D",ldsh(iw),"T",ltck(iw),"C",lcol)
; end do
;if ( n_sym_asym .eq. 0 .or. n_sym_asym .eq. 1 ) then
 do iw=0, 0
   cnaddline(ctr2,wks,kwn2(i0_2+1:),ome2(1:),dispc(0,iw,1:,i0_2+1:),lzs(0,iw,:),0,"D",ldsh(0),"T",ltck(0),"C",lcol)
 end do
;end if
;if ( n_sym_asym .eq. 0 .or. n_sym_asym .eq. 2 ) then
; do iw=0, 1, 2
;   cnaddline(ctr2,wks,kwn2(i0_2+1:),ome2(1:),dispc(1,iw,1:,i0_2+1:),lzs(1,iw,:),0,"D",ldsh(2),"T",ltck(2),"C",lcol)
; end do
;end if
if ( l_cellfill ) then
 setvalues ctr2
  "cnFillMode" : "CellFill"
 end setvalues
end if

 set_size(ctr2,xf(ie),yf(ie),wi1,he)
 axis_range(ctr2,1,kmax_draw,0.01,omax_draw+0.001,True,True)
 axis_str(ctr2,"",xstr(ie),ystr(ie))
 setvalues ctr2
  "tmXBMode" : "Explicit"
  "tmXBValues" : xtckv
  "tmXBLabels" : xtckv
  "tmXBMinorValues" : xmntckv
  "tmYLMode" : "Explicit"
  "tmYLValues" : ytckv
  "tmYLLabels" : ytckv
  "tmYLMinorValues" : ymntckv
 end setvalues
 setvalues ctr2
  "tiXAxisOffsetXF" : wi*0.0012
 end setvalues
 if (xf(ie) .ne. xf(0)) then
 setvalues ctr2
  "tmYLLabelsOn" : False
 end setvalues
 end if

 cnfinalize(ctr2,wks)

 draw(ctr2)


 dr_txt2(wks,mstr(ip,ie),xf(ie)+1,yf(ie)-2,22,4.0,0., 1,1)
 dr_txt2(wks,mstr(ip,ie),xfa(ie)+1,yf(ie)-2,22,4.0,0., 1,1)

end do

 dr_txt2(wks,"(a)",xf(0)-2,yf(0)-5,25,5.0,0., 1,1)
 dr_txt2(wks,"(b)",xfa(0)-2,yf(0)-5,25,5.0,0., 1,1)

xfl = (/xf(0),xfa(0)/)
do ip=0, 1

 label = create "lab" labelBarClass wks
   "lbOrientation" : "Horizontal"
   "lbBoxCount"   : dimsizes(cnfil)
   "lbFillColors" : cnfil
   "lbMonoFillPattern" : True
   "lbAutoManage"   : False
   "lbLabelsOn"     : True
   "lbLabelStrings" : cnlev(ip,:)
   "lbLabelAutoStride"  : "True"
;   "lbLabelStride"  : 1
   "lbLabelAlignment" : "InteriorEdges"
   "lbPerimOn"        : False
   "lbTitleString"    : lbtit
   "lbTitlePosition"  : "Right"
   "lbTitleDirection" : "Across"
 end create
 set_size_lb(label,xfl(ip),yf(ne-1)+he+8,wi*1.0,3)
 draw(label)

end do

 frame(wks)

end

