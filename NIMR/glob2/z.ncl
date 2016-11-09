load "$UTIL/ncl/header.ncl"
xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"
axisfontratio = 1.2

begin

;--------------------------------------- parameters
 jobname = (/"tc","tk"/)
 expname = (/"AOL60CGW","AOL60CTL"/)
 ppid    = "e"
 varname = (/"epd","du_gw"/)
 vvarname = (/"f_y","f_z"/)

 y0 = 1951
 y1 = 1960
 y2 = 1999
 mms = "djf"

 wks = cr_wks("fig/clim-wf_yz-"+mms+"-chg-diff",600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 nv = dimsizes(varname)+1
 nyr = y2 - y1 + 1

 mm3s = (/"djf","mam","jja","son"/)
 mm3i = (/ (/1,2,12/), (/3,4,5/), (/6,7,8/), (/9,10,11/) /)
 do is=0, dimsizes(mm3s)-1
   if ( mms .eq. mm3s(is) ) then
     mm = mm3i(is,:)
   end if
 end do
 nm = dimsizes(mm)

 fnamej = "../dat/"+expname+"/xm.uan"+jobname+"a.p"+ppid+".nc"
 fnamet = "../dat/"+expname+"/tem/uan"+jobname+".tem-wc_yzt.nc"
 fj = addfile(fnamej(0),"r")
 fjx = addfile(fnamej(1),"r")
 ft = addfile(fnamet(0),"r")
 ftx = addfile(fnamet(1),"r")
 n0 = (y1 - y0)*12

 nf = n0 + 98*12

 ; get y, z
 tmp0 = ft->$varname(0)$(nf:nf+nyr*12-1,:,:) - ft->$varname(0)$(n0:n0+nyr*12-1,:,:)
 tmp1 = fj->$varname(1)$(nf:nf+nyr*12-1,:,:)*72. - fj->$varname(1)$(n0:n0+nyr*12-1,:,:)*72.
 tmpv0 = ft->$vvarname(0)$(nf:nf+nyr*12-1,:,:) - ft->$vvarname(0)$(n0:n0+nyr*12-1,:,:)
 tmpv1 = ft->$vvarname(1)$(nf:nf+nyr*12-1,:,:) - ft->$vvarname(1)$(n0:n0+nyr*12-1,:,:)
 if ( isfilevar(fj,"dudt_gwconv") ) then
   tmp1 = tmp1 + fj->dudt_gwconv(nf:nf+nyr*12-1,:,:)*86400. - fj->dudt_gwconv(n0:n0+nyr*12-1,:,:)*86400.
 end if

 tmp = ft->$varname(0)$(n0:n0+nyr*12-1,:,:)
 lat = tmp&$tmp!2$
 z = tmp&$tmp!1$ / 1.e3
 ny = dimsizes(lat)
 nz = dimsizes(z)

 tmp0x = ftx->$varname(0)$(nf:nf+nyr*12-1,:,:) - ftx->$varname(0)$(n0:n0+nyr*12-1,:,:)
 tmp1x = fjx->$varname(1)$(nf:nf+nyr*12-1,:,:)*72. - fjx->$varname(1)$(n0:n0+nyr*12-1,:,:)*72.
 tmpv0x = ftx->$vvarname(0)$(nf:nf+nyr*12-1,:,:) - ftx->$vvarname(0)$(n0:n0+nyr*12-1,:,:)
 tmpv1x = ftx->$vvarname(1)$(nf:nf+nyr*12-1,:,:) - ftx->$vvarname(1)$(n0:n0+nyr*12-1,:,:)
 if ( isfilevar(fjx,"dudt_gwconv") ) then
   tmp1x = tmp1x + fjx->dudt_gwconv(nf:nf+nyr*12-1,:,:)*86400. - fjx->dudt_gwconv(n0:n0+nyr*12-1,:,:)*86400.
 end if

 vara = new((/nv,nyr,nz,ny,nm/),float)
 varax = vara
 veca = new((/2,nyr,nz,ny,nm/),float)
 vecax = veca
 do im=0, nm-1
   vara(0,:,:,:,im) = (/ tmp0(ispan(mm(im)-1,nyr*12-1,12),:,:) /)
   vara(1,:,:,:,im) = tmp1(ispan(mm(im)-1,nyr*12-1,12),:,:)
   vara(2,:,:,:,im) = vara(0,:,:,:,im) + vara(1,:,:,:,im)
   veca(:,:,:,:,im) = (/ tmpv0(ispan(mm(im)-1,nyr*12-1,12),:,:), \
                         tmpv1(ispan(mm(im)-1,nyr*12-1,12),:,:) /)/6371.*1.e3
 end do
 do im=0, nm-1
   varax(0,:,:,:,im) = (/ tmp0x(ispan(mm(im)-1,nyr*12-1,12),:,:) /)
   varax(1,:,:,:,im) = tmp1x(ispan(mm(im)-1,nyr*12-1,12),:,:)
   varax(2,:,:,:,im) = varax(0,:,:,:,im) + varax(1,:,:,:,im)
   vecax(:,:,:,:,im) = (/ tmpv0x(ispan(mm(im)-1,nyr*12-1,12),:,:), \
                         tmpv1x(ispan(mm(im)-1,nyr*12-1,12),:,:) /)/6371.*1.e3
 end do

 delete(tmp0)
 delete(tmp1)
 delete(tmpv0)
 delete(tmpv1)
 delete(tmp0x)
 delete(tmp1x)
 delete(tmpv0x)
 delete(tmpv1x)

;--------------------------------------- process 0
 vary = dim_avg(vara)
 vecy = dim_avg(veca)
 varyx = dim_avg(varax)
 vecyx = dim_avg(vecax)

 vavg = dim_avg_n(vary,1)
 vec = dim_avg_n(vecy,1)
 vavgx = dim_avg_n(varyx,1)
 vecx = dim_avg_n(vecyx,1)

; vstd = dim_stddev_n(vary,1)
; vstdx = dim_stddev_n(varyx,1)
;
; cl = 100.*(1. - ttest(vavg,vstd^2,nyr, \
;                       vavgx,vstdx^2,nyr,False,False) )

 vavg = vavg - vavgx
 vec  = vec  - vecx

 do k=19, nz-1
   vec(:,k,:) = vec(:,k,:)*exp(z(k)/7.)
 end do
 vec(0,:,:) = vec(0,:,:)/111.

 vec@_FillValue = 1.e20

;--------------------------------------- cn level
 cnlev = ispan(-10,10,1)*0.04

; cnlev2 = array_append_record( ispan(-30,-10,1), ispan(10,30,1), 0 )*0.2
 cnlev2 = cnlev*10.

 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)
 
;--------------------------------------- contour
 wi = 32
 he0 = 40
 xf = 10 + (/0,1,0,1/)*40
 yf = 10 + (/0,0,1,1/)*58
 ymin = (/8,8,8,8/)
 he = he0*((50-ymin)/50.)

 ystr  = new(nv,string)

 mstr = (/"EPD, EP flux","GWD","EPD + GWD"/)
 xstr = "latitude (deg)"
 ystr(:) = ""
 ystr((/0,2/)) = "hybrid height (km)"
 lbtit = "(m s:S:-1:N: day:S:-1:N:)"

do iv=0, nv-1

 contour1 = cnshadeplot(wks,lat,z,vavg(iv,:,:),cnlev,cnfil,"")
; cnaddline(contour1,wks,lat,z,cl(iv,:,:),(/90,95/),0,"D",0,"T",1.2,"C",1)
 cnaddline(contour1,wks,lat,z,vavg(iv,:,:),cnlev2,0,"D",0,"T",0.5,"C",(/1,-1,1/))

 set_size(contour1,xf(iv),yf(iv),wi,he(iv))
 axis_range(contour1,-90,90,ymin(iv),50,False,False)
 axis_str(contour1,mstr(iv),xstr,ystr(iv))

 cnfinalize(contour1,wks)

 if (iv .eq. nv-1) then
 label = create "lab" labelBarClass wks
  "lbOrientation" : "Vertical"
  "lbBoxCount"   : dimsizes(cnfil)
  "lbFillColors" : cnfil
  "lbMonoFillPattern" : True
  "lbLabelsOn" : True
  "lbLabelStrings" : cnlev
  "lbLabelStride" : 1
  "lbAutoManage" : False
  "lbLabelAlignment" : "InteriorEdges"
  "lbPerimOn" : False
  "lbTitleString" : ""
  "lbTitlePosition" : "Bottom"
 end create
 set_size_lb(label,xf(1)+2,yf(nv-1)-2,6,he(nv-1)*1.3)
 draw(label)
 dr_txt(wks,lbtit,xf(1)+2+2.5,yf(nv-1)-2+(he(nv-1)*1.3)*0.9,22,5,0.)
 end if


 if (iv .eq. 0) then

  vf = create "epf_vec" vectorFieldClass defaultapp
   "vfDataArray" : vec(:,22::2,::8)
   "vfXArray"    : lat(::8)
   "vfYArray"    : z(22::2)
   "vfMissingUValueV" : 1.e20
   "vfMissingVValueV" : 1.e20
  end create
 
  vc1 = create "vector" vectorPlotClass wks
   "vcVectorFieldData" : vf
   "vcRefMagnitudeF" : 1.e6 * 0.1
   "vcRefLengthF"    : 0.05
   "vcRefAnnoOn"     : False
 ;  "vcLineArrowThicknessF" :  thc(nn)
 ;  "vcLineArrowColor"      :  col(nn)
 ;  "vcRefAnnoString1" : ref*tmpc(nn)
   "vcRefAnnoString2On" : False
 ;  "vcRefAnnoParallelPosF" :  1. - nn*0.1
   "vcGlyphStyle" : "FillArrow"
   "vcFillArrowEdgeThicknessF" : 0.1
   "vcFillArrowWidthF" : 0.1
   "vcPositionMode" : "ArrowTail"
  end create

  set_size(vc1,xf(iv),yf(iv),wi,he(iv))
  axis_range(vc1,-90,90,ymin(iv),50,False,False)
  axis_str(vc1,"","","")
 
  overlay(contour1,vc1)

 end if

 draw(contour1)

 if (iv .eq. 0) then
   NhlRemoveOverlay(contour1,vc1,False)
 end if

end do  ; iv


 frame(wks)

end

