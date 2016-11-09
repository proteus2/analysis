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
 load "./ftn_get_color_S-RIP.ncl"
 load "./ftn_read_spec.ncl"

 param0_i_wav = 1
 k_a1 = param0_wavs_k_a1(param0_i_wav)
 k_a2 = param0_wavs_k_a2(param0_i_wav)
 omin9 = 0.04
 omax9 = 2.

; ra = param0_ra
 ra = array_append_record(param0_ra_p,param0_ra_m,0)
 var = "T"

 yyyymm1 = param0_yyyymm1
 yyyymm2 = param0_yyyymm2
 lat_b = param0_lat_b
 p     = 10 ;param0_p

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
 
 foname = "z"
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

 vsp0 = read_fspec_ko(ra,var,yyyymm1,yyyymm2,p,lat_b,kmax,no)
 vsp0 = vsp0*wgt_spec
 i0 = vsp0@i_k0

 dim3 = dimsizes(vsp0)
 vsp = new((/dim3(0),nfig,dim3(1),dim3(2)/),float)
 vsp = vsp0

 kwn = fspan(-kmax,kmax,kmax*2+1)
 ome = fspan(0,no/2,no/2+1) / tofloat(ndays_spec)
 dome = 1./tofloat(ndays_spec)


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

 vsp(:,:,:,i0) = vsp@_FillValue  ; k = 0

 ; to integrate around omega = 0
; vsp(:,:,0,:) = vsp(:,:,0,:)*0.5

 vsp(:,:,no/2,:) = vsp@_FillValue

 vsp_o1 = dim_sum(vsp(:,:,:,i0+k_a1:i0+k_a2))
 vsp_o2 = dim_sum(vsp(:,:,:,i0-k_a2:i0-k_a1))

 vsp_o1!2 = "ome"
 vsp_o1&ome = ome
 vsp_o2!2 = "ome"
 vsp_o2&ome = ome

; convert to energe-preserving form
 vsp_o1 = vsp_o1*conform_dims(dimsizes(vsp_o1),abs(ome)*log(10.),2)
 vsp_o2 = vsp_o2*conform_dims(dimsizes(vsp_o2),abs(ome)*log(10.),2)
; vsp_k = vsp_k*conform_dims(dimsizes(vsp_k),abs(kwn)*log(10.),2)

print(max(vsp_o1))
print(max(vsp_o2))

print(abs(ome(1:))+"  "+dim_avg_n(vsp_o1(:,0,1:),0))

;--------------------------------------- cn level
 fact = 0
 if ( var0 .eq. "T" ) then
   yunt = "K:S:2:N:"
   if ( l_normalize ) then
     fact = -6
     yunt = ""
   end if
 end if
 if ( var0 .eq. "U" ) then
   yunt = "m:S:2:N: s:S:-2:N:"
 end if
 if ( var0 .eq. "V" ) then
   yunt = "m:S:2:N: s:S:-2:N:"
 end if
 if ( var0 .eq. "W" ) then
   yunt = "mm:S:2:N: s:S:-2:N:"
 end if
 lbt_hd = ""
 if (fact .ne. 0) then
   lbt_hd = "10:S:"+fact+":N: "
 end if
; yunt = "[ "+lbt_hd+yunt+" / cyc day:S:-1:N: ]"
 yunt = "["+lbt_hd+yunt+"]"

 ltck = (/1,1,1,1,1,1,1,1/)*0.8
 ldsh = (/0,0,0,0,0,14,14,14/)
 lcol = (/ icol_s@ei , icol_s@m , icol_s@m2 , icol_s@c , icol_s@j5 ,  \
           icol_s@ei , icol_s@m2 , icol_s@j5 /)

;--------------------------------------- ctr
 wi = 56
 he = 30
; xf = 12 + (/0,1,2,0,1,2/)*(wi+6)
; yf = 10 + (/0,0,0,1,1,1/)*(he+8)
 xf = 12 + (/0,1,0,1,0,1/)*(wi+6)
 yf = 10 + (/0,0,1,1,2,2/)*(he+8)

 xspc = 6
 wi1 = (wi-xspc)*0.5
 xfp = xf + wi1 + xspc

 mstr = vstr+":F22:"
 xstr = "frequency [cyc day:S:-1:N:]"
 ystr = "PSD "+yunt
; ystr = "PSD{"+vstr+":F22:} "+yunt

 xtckv = (/0.1,0.2,0.5,1,2,5,10,20,50/)*0.1
 xmntckv = (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50/)*0.1
 ytckv = (/0.01,0.1,1/)
 ymntckv = (/1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400/)*0.01


do ip=0, nfig-1

 ; k < 0

 tmp = vsp_o2(:,ip,1:)/(10.^fact)
 tmp!1 = "ome"
 tmp&ome = ome(1:)
 tmp(:,{0.8:}) = 0.
 ymax = max(tmp)

; plot1 = xylineplot(wks,ome(1:),vsp_o2(:,ip,1:)/(10.^fact),"D",ldsh,"T",ltck,"C",lcol)
 plot1 = xylineplot(wks,(/10,20/),(/0,ymax/),"D",0,"T",0.2,"C",-1)
 xyaddline(plot1,ome(1:),vsp_o2(:,ip,1:)/(10.^fact),"D",ldsh,"T",ltck,"C",lcol)

 set_size(plot1,xf(ip),yf(ip),wi1,he)
 axis_range(plot1,omin9,omax9,0,999,True,False)
; axis_range(plot1,omax9,omin9,0,999,True,False)
; axis_range(plot1,omax9,omin9,0.01,999,True,True)
 axis_str(plot1,"",xstr,ystr)
 setvalues plot1
  "tmXBMode" : "Explicit"
  "tmXBValues" : xtckv
  "tmXBLabels" : xtckv
  "tmXBMinorValues" : xmntckv
;  "tmYLMode" : "Explicit"
;  "tmYLValues" : ytckv
;  "tmYLLabels" : ytckv
;  "tmYLMinorValues" : ymntckv
 end setvalues
 setvalues plot1
  "tiXAxisOffsetXF" : wi*0.00115
 end setvalues

 draw(plot1)


 ; k > 0

 tmp = vsp_o1(:,ip,1:)/(10.^fact)
 tmp!1 = "ome"
 tmp&ome = ome(1:)
 tmp(:,{0.8:}) = 0.
 ymax = max(tmp)

; plot2 = xylineplot(wks,ome(1:),vsp_o1(:,ip,1:)/(10.^fact),"D",ldsh,"T",ltck,"C",lcol)
 plot2 = xylineplot(wks,(/10,20/),(/0,ymax/),"D",0,"T",0.2,"C",-1)
 xyaddline(plot2,ome(1:),vsp_o1(:,ip,1:)/(10.^fact),"D",ldsh,"T",ltck,"C",lcol)

 set_size(plot2,xfp(ip),yf(ip),wi1,he)
 axis_range(plot2,omin9,omax9,0,999,True,False)
; axis_range(plot2,omin9,omax9,0.01,999,True,True)
 axis_str(plot2,"","","")
 setvalues plot2
  "tmXBMode" : "Explicit"
  "tmXBValues" : xtckv
  "tmXBLabels" : xtckv
  "tmXBMinorValues" : xmntckv
;  "tmYLMode" : "Explicit"
;  "tmYLValues" : ytckv
;  "tmYLLabels" : ytckv
;  "tmYLMinorValues" : ymntckv
 end setvalues
; setvalues plot2
;  "tmYLLabelsOn" : False
; end setvalues

 draw(plot2)


 dr_txt2(wks,mstr,xf(ip)+1,yf(ip)+2,22,4.5,0., 1,1)

end do

 lg = lgplot(wks,ranames(:4,1),"D",ldsh(:4),"T",ltck(:4),"C",lcol(:4),2)
 set_size_lg(lg,xf(0),yf(nfig-1)+he+6,wi*0.25,2.8*5)
 draw(lg)


 frame(wks)

end

