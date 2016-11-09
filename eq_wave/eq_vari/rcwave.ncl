load "$UTIL/ncl/header.ncl"

xspace = 83   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

;axisfontratio = 0.9

begin

;--------------------------------------- parameters
; expname = "era-int"
; dirname = "ERA-I"
; varname = (/"U","V","W","T"/)+"_GDS0_ISBL"

 expname = "merra"
 dirname = "MERRA"
 varname = (/"u","v","omega","t"/)

 yy = "1979";"1983"
 mm = "12";"05"
 ndays = 90
 p0 = 70
 n0 = 120 + 90  ; 90, 0 - 120
 lat0 = 25
 ymax = 20

 nkc = 6
 perc = 1  ; [day]  ; cut-off period

 nt = ndays*4
 if (expname .eq. "merra") then  ; 8-hr data
   nt = nt*2
   n0 = n0*2
 end if

 wks = cr_wks("fig/rcwave",600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colindr = NhlNewColor(wks,1.0,0.0,0.0)
 colindb = NhlNewColor(wks,0.0,0.0,1.0)

;--------------------------------------- read data
 nv = dimsizes(varname)
 nw = 4

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
 lat = array_append_record( lath(ny-1:1)*(-1.), lath, 0 )

 cf = new((/2,nw+2,nv,ny,nt,nk/),float)
 cf = 0.
 si = (/"r","i"/)
 sw = (/"k","r_s","r_a","mrg","ig_s","ig_a"/)
 f0 = addfile(fname,"r")
 do iv=0, nv-1
   do i=0, 1
   do j=0, nw+1
     s_v0 = "fc"+si(i)+"_"+varname(iv)+"_"+sw(j)
     cf(i,j,iv,:,:no,:) = f0->$s_v0$({p0},{0:lat0},:no,:)
     cf(i,j,iv,:,nt-no:,:) = f0->$s_v0$({p0},{0:lat0},no:,:)
   end do
   end do
 end do
 cf(:,:,:,:,ndays/perc:nt-ndays/perc,:) = 0.
 cf(:,:,:,:,:,nkc+1:) = 0.

;--------------------------------------- process 0
 vv = new((/nw+2,nv,ny,nx/),float)
 vv = 0.
 cfx = new((/2,nw+2,ny,nt,nx/),float)
 do iv=0, nv-1
   cfx = 0.
   do i=1, nkc
     cfx(:,:,:,:,i) = cfftb(cf(:,:,iv,:,:,i),0)
   end do
   do j=0, nw+1
     vv(j,iv,:,:) = cfftb(cfx(:,j,:,n0,:),1)*2. * (nx*nt)
   end do
 end do

 vv2 = new((/nw,nv,ny2,nx+1/),float)
 vv2(0,:,ny-1:,:nx-1) = vv(0,:,:,:)
 vv2(1,:,ny-1:,:nx-1) = vv(1,:,:,:) + vv(2,:,:,:)
 vv2(2,:,ny-1:,:nx-1) = vv(3,:,:,:)
 vv2(3,:,ny-1:,:nx-1) = vv(4,:,:,:) + vv(5,:,:,:)
 vv2(0,:,0:ny-2,:nx-1) = vv(0,:,ny-1:1,:)
 vv2(1,:,0:ny-2,:nx-1) = vv(1,:,ny-1:1,:) + vv(2,:,ny-1:1,:)*(-1.)
 vv2(2,:,0:ny-2,:nx-1) = vv(3,:,ny-1:1,:)*(-1.)
 vv2(3,:,0:ny-2,:nx-1) = vv(4,:,ny-1:1,:) + vv(5,:,ny-1:1,:)*(-1.)
 vv2(:,1,0:ny-2,:nx-1) = vv2(:,1,0:ny-2,:nx-1)*(-1.)  ; v
 vv2(:,:,:,nx) = vv2(:,:,:,0)  ; cyclic

 lon = fspan(0.,360.,nx+1)

; tmp = vv2(:,2,:,:)
; do j=1, ny2-2
;   vv2(:,2,j,:) = (2.*tmp(:,j,:) + tmp(:,j+1,:) + tmp(:,j-1,:))/4.
; end do

;--------------------------------------- cn level
 cnlev = fspan(-1,1,21)*5.
; cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)
 cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnlevh = fspan(-1,1,21)*5.
 cnfac = (/1.0,1.0,1.0,1.0/)

;--------------------------------------- contour
 wi = 120
 he = 20
 xf = 8 + (/0,0,0,0/)*(wi+8)
 yf = 10 + (/0,1,2,3/)*(he+5)

 mstr  = new(nw,string)
 xstr  = new(nw,string)
 ystr  = new(nw,string)
 lbopt = new(nw,string)

 mstr = "" ;(/"Kelvin","Rossby","MRG","IG"/)
 xstr(:) = ""
 ystr(:) = ""
; xstr(3) = "longitude"
; ystr(:) = "latitude"
 lbtxt = (/"","","","(mPa s:S:-1:N:)"/)

do il=0, nw-1

 iv = 2  ; omega
 contour1 = cnshadeplot(wks,lon,lat,vv2(il,iv,:,:)*1.e3,cnlev*cnfac(il),cnfil,lbtxt(il))
 iv = 3  ; T
 cnaddline(contour1,wks,lon,lat,vv2(il,iv,:,:),cnlevh*cnfac(il),0,"D",0,"T",0.5,"C",(/colindr,-1,colindb/))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.15
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
 set_size(contour1,xf(il),yf(il),wi,he)
 axis_range(contour1,0,360,-ymax,ymax,False,False)
 axis_str(contour1,mstr(il),xstr(il),ystr(il))
 setvalues contour1
  "tmXBMode"   : "Explicit"
  "tmXBValues" : ispan(0,360,60)
  "tmXBLabels" : (/"0","60E","120E","180","120W","60W","0"/)
  "tmXBMinorValues" : ispan(0,360,10)
  "tmYLMode"   : "Explicit"
  "tmYLValues" : ispan(-30,30,10)
  "tmYLLabels" : (/"30S","20S","10S","EQ","10N","20N","30N"/)
  "tmYLMinorValues" : ispan(-30,30,5)
 end setvalues

 cnfinalize(contour1,wks)

 vf = create "epf_vec" vectorFieldClass defaultapp
  "vfDataArray" : vv2(il,0:1,1::4,:)
  "vfXArray"    : lon
  "vfYArray"    : lat(1::4)
 end create
 ref = 2.*cnfac
 vec1 = create "vector" vectorPlotClass wks
  "vcVectorFieldData" : vf
  "vcRefMagnitudeF" : ref(il)
  "vcRefLengthF"    : 0.01
  "vcRefAnnoOn"     : False
;  "vcLineArrowThicknessF" :  thc(nn)
;  "vcLineArrowColor"      :  col(nn)
;  "vcRefAnnoString1" : ref*tmpc(nn)
  "vcRefAnnoString2On" : False
;  "vcRefAnnoParallelPosF" :  1. - nn*0.1
 end create
 set_size(vec1,xf(il),yf(il),wi,he)
 axis_range(vec1,0,360,-ymax,ymax,False,False)
 axis_str(vec1,"","","")

 overlay(contour1,vec1)
 draw(contour1)
 NhlRemoveOverlay(contour1,vec1,False)

end do

 frame(wks)

end

