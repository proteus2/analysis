load "$UTIL/ncl/header.ncl"

xspace = 83   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

;axisfontratio = 0.9

begin

;--------------------------------------- parameters
 expname = "era-int"+"_f"
 dirname = "ERA-I"+"-nr"
; varname = (/"U","V","W","T"/)+"_GDS0_ISBL"
 varname = (/"u","v","w","t"/)

; expname = "merra"
; dirname = "MERRA"
; varname = (/"u","v","omega","t"/)

 yy = "2004";"1983"
 mm = "01";"05"
 ndays = 90
 p0 = 80  ; 96, 80, 66
; nn = 120 + (/0, 119/)
 nn = 120 + (/0, 119+80/) - 40
 lat0 = 20

 nkc = 6
 perc = 1  ; [day]  ; cut-off period

 nt = ndays*4
 if (expname .eq. "merra") then  ; 8-hr data
   nt = nt*2
   nn = nn*2
   nn(1) = nn(1) + 1
 end if

 wks = cr_wks("fig/kw_xt-"+yy+mm,600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colindr = NhlNewColor(wks,1.0,0.0,0.0)
 colindb = NhlNewColor(wks,0.0,0.0,1.0)

;--------------------------------------- read data
 nv = dimsizes(varname)
 ntr = nn(1) - nn(0) + 1

 fdir  = "$DATD/fspec/"+dirname+"/fcoef_ko-wav4"
 fname = fdir+"/"+yy+"/"+expname+".fft_var4_ko."+yy+"."+mm+".nc"
print(fname+"")

; axis
 f0 = addfile(fname,"r")
 s_v0 = "fcr_"+varname(0)+"_k"
 no = dimsizes(f0->$(f0->$s_v0$!2)$) / 2
 nk = dimsizes(f0->$(f0->$s_v0$!3)$)
 lath = f0->$(f0->$s_v0$!1)$({0:lat0})
 ny = dimsizes(lath)
 ny2 = ny*2 - 1
 nx = (nk-1)*2
; lat = array_append_record( lath(ny-1:1)*(-1.), lath, 0 )

 sw = "k"   ; k  /  r_s, r_a  /  mrg  /  ig_s, ig_a
 nsw = dimsizes(sw)  ; should be 1 or 2
 cf = new((/2,nsw,nv,ny,nt,nk/),float)
 cf = 0.
 si = (/"r","i"/)
 f0 = addfile(fname,"r")
 do iv=0, nv-1
   do i=0, 1
   do j=0, nsw-1
     s_v0 = "fc"+si(i)+"_"+varname(iv)+"_"+sw(j)
     cf(i,j,iv,:,:no,:) = f0->$s_v0$({p0},{0:lat0},:no,:)
     cf(i,j,iv,:,nt-no:,:) = f0->$s_v0$({p0},{0:lat0},no:,:)
   end do
   end do
 end do
 cf(:,:,:,:,ndays/perc:nt-ndays/perc,:) = 0.
 cf(:,:,:,:,:,nkc+1:) = 0.

;--------------------------------------- process 0
 vv = new((/nsw,nv,ny,nt,nx/),float)
 vv = 0.
 cfx = new((/2,nsw,ny,nt,nx/),float)
 do iv=0, nv-1
   cfx = 0.
   do i=1, nkc
     cfx(:,:,:,:,i) = cfftb(cf(:,:,iv,:,:,i),0)
   end do
   do j=0, nsw-1
     vv(j,iv,:,:,:) = cfftb(cfx(:,j,:,:,:),1)*2. * (nx*nt)
   end do
 end do

; vv2 = new((/nv,ny2,nx+1/),float)
; if (nsw .eq. 1) then
;   vv2(:,ny-1:,:nx-1) = vv(0,:,:,n0,:)
;   if (sw(0) .eq. "k") then
;     vv2(:,0:ny-2,:nx-1) = vv(0,:,ny-1:1,n0,:)
;   end if
;   if (sw(0) .eq. "mrg") then
;     vv2(:,0:ny-2,:nx-1) = vv(0,:,ny-1:1,n0,:)*(-1.)
;   end if
; else
;   vv2(:,ny-1:,:nx-1) = vv(0,:,:,n0,:) + vv(1,:,:,n0,:)
;   vv2(:,0:ny-2,:nx-1) = vv(0,:,ny-1:1,n0,:) + vv(1,:,ny-1:1,n0,:)*(-1.)
; end if
; vv2(1,0:ny-2,:nx-1) = vv2(1,0:ny-2,:nx-1)*(-1.)  ; v
; vv2(:,:,nx) = vv2(:,:,0)  ; cyclic

 vv2 = new((/nv,ntr,nx+1/),float)
 vv2(:,:,:nx-1) = ( dim_sum_n( vv(0,:,1:,nn(0):nn(1),:), 1 )*2. +  \
                    vv(0,:,0,nn(0):nn(1),:) )/tofloat(ny2)
 vv2(1,:,:) = 0.
 vv2(:,:,nx) = vv2(:,:,0)  ; cyclic

 lon = fspan(0.,360.,nx+1)
 tmp = fspan(1.-10.,31.+10.,ntr+1)
 t = tmp(:ntr-1)

;--------------------------------------- cn level
 cnlev = fspan(-1,1,21)*5.*3.
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnlevh = fspan(-1,1,21)*3.*3.
 cnfac = 1.0

;--------------------------------------- contour
 wi = 120
 he = 150
 xf = 8 + (/0,0,0,0/)*(wi+8)
 yf = 10 + (/0,1,2,3/)*(he+5)


 mstr = ""
 xstr = ""
 ystr = ""
; xstr = "longitude"
 ystr = "t (day)"
 lbtxt = "(10:S:-3:N: m:S:2:N: s:S:-2:N:)" ;"(mPa s:S:-1:N:)"

 iv = 2  ; omega
; contour1 = cnshadeplot(wks,lon,t,vv2(iv,:,:)*1.e3,cnlev*cnfac,cnfil,lbtxt)
 mf = vv2(2,:,:)*vv2(0,:,:)*6.6e3/(p0*100.)*(-1.)
 contour1 = cnshadeplot(wks,lon,t,mf*1.e3,cnlev*cnfac,cnfil,lbtxt)
 iv = 3  ; T
; cnaddline(contour1,wks,lon,t,vv2(iv,:,:),cnlevh*cnfac,0,"D",0,"T",0.5,"C",(/colindr,-1,colindb/))
 cnaddline(contour1,wks,lon,t,vv2(iv,:,:)*abs(vv2(iv,:,:)),cnlevh*cnfac,0,"D",(/0,0,1/),"T",0.5,"C",(/1,-1,1/))
 setvalues contour1
;  "pmLabelBarOrthogonalPosF" : 0.15
;  "pmLabelBarParallelPosF"   : 0.50-0.40
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.6
  "pmLabelBarHeightF"        : 0.05
  "lbLabelFontHeightF" : 0.010
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.10
  "lbTitleFontHeightF" : 0.010
 end setvalues
 set_size(contour1,xf,yf,wi,he)
 axis_range(contour1,0,360,1-10,31+10,False,False)
 axis_str(contour1,mstr,xstr,ystr)
 setvalues contour1
  "tmXBMode"   : "Explicit"
  "tmXBValues" : ispan(0,360,60)
  "tmXBLabels" : (/"0","60E","120E","180","120W","60W","0"/)
  "tmXBMinorValues" : ispan(0,360,10)
 end setvalues

 cnfinalize(contour1,wks)

; vf = create "epf_vec" vectorFieldClass defaultapp
;  "vfDataArray" : vv2(0:1,::4,:)
;  "vfXArray"    : lon
;  "vfYArray"    : t(::4)
; end create
; ref = 2.*cnfac
; vec1 = create "vector" vectorPlotClass wks
;  "vcVectorFieldData" : vf
;  "vcRefMagnitudeF" : ref
;  "vcRefLengthF"    : 0.005
;  "vcRefAnnoOn"     : False
;;  "vcLineArrowThicknessF" :  thc(nnn)
;;  "vcLineArrowColor"      :  col(nnn)
;;  "vcRefAnnoString1" : ref*tmpc(nnn)
;  "vcRefAnnoString2On" : False
;;  "vcRefAnnoParallelPosF" :  1. - nnn*0.1
; end create
; set_size(vec1,xf,yf,wi,he)
; axis_range(vec1,0,360,1,31,False,False)
; axis_str(vec1,"","","")
;
; overlay(contour1,vec1)
 draw(contour1)
; NhlRemoveOverlay(contour1,vec1,False)

 frame(wks)

end

