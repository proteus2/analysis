load "$UTIL/ncl/header.ncl"
xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"
axisfontratio = 1.2

begin

;--------------------------------------- parameters
 jobname = "tc"
 expname = "AOL60CGW"
 ppid    = "j"
 varname = (/"theta","u","w_res"/)
 vvarname = (/"v_res","w_res"/)

 y0 = 1951
 y1 = 1960;+98
 y2 = 1999;+98
 mms = "jja"

 wks = cr_wks("fig/clim-circ_yz-"+mms+"_"+y1+"-"+y2+"-"+jobname,600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 nv = dimsizes(varname)
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
 fj = addfile(fnamej,"r")
 ft = addfile(fnamet,"r")
 n0 = (y1 - y0)*12

 ; get y, z
 tmp0 = fj->$varname(0)$(n0:n0+nyr*12-1,:,:)
 tmp1 = fj->$varname(1)$(n0:n0+nyr*12-1,:,:)
 tmp2 = ft->$varname(2)$(n0:n0+nyr*12-1,:,:)
 tmpp = fj->p(n0:n0+nyr*12-1,:,:)
 tmpv0 = ft->$vvarname(0)$(n0:n0+nyr*12-1,:,:)
 tmpv1 = ft->$vvarname(1)$(n0:n0+nyr*12-1,:,:)
 lat = tmp0&$tmp0!2$
 z = tmp1&$tmp1!1$ / 1.e3
 ny = dimsizes(lat)
 nz = dimsizes(z)

 tmpt = tmp0
 do k=1, nz-1
   tmpt(:,k,:) = exp(0.5*(log(tmp0(:,k-1,:)*tmp0(:,k,:))))
 end do
 tmp0 = tmpt * (tmpp/1.e5)^0.286

 vara = new((/nv,nyr,nz,ny,nm/),float)
 veca = new((/2,nyr,nz,ny,nm/),float)
 do im=0, nm-1
   vara(0,:,:,:,im) = (/ tmp0(ispan(mm(im)-1,nyr*12-1,12),:,:) /) - 273.15
   vara(1,:,:,:,im) = tmp1(ispan(mm(im)-1,nyr*12-1,12),:,:)
   vara(2,:,:,:,im) = tmp2(ispan(mm(im)-1,nyr*12-1,12),:,:)*1.e3
   veca(:,:,:,:,im) = (/ tmpv0(ispan(mm(im)-1,nyr*12-1,12),:,:), \
                         tmpv1(ispan(mm(im)-1,nyr*12-1,12),:,:) /)*1.e3
 end do

 delete(tmp0)
 delete(tmp1)
 delete(tmp2)
 delete(tmpp)
 delete(tmpv0)
 delete(tmpv1)

;--------------------------------------- process 0
 vary = dim_avg(vara)
 vecy = dim_avg(veca)

 vavg = dim_avg_n(vary,1)
 vec = dim_avg_n(vecy,1)

 do k=19, nz-1
   vec(:,k,:) = vec(:,k,:)*exp(-z(k)/(7.*2.))
 end do
 vec(0,:,:) = vec(0,:,:)/111.

 vec@_FillValue = 1.e20

;--------------------------------------- cn level
 cnlev = new((/nv,21/),float)
 cnlev(0,:) = ispan(-10,10,1)*4-20.
 cnlev(1,:) = ispan(-10,10,1)*5.
 cnlev(2,:) = ispan(-10,10,1)*0.4

 cnlev2 = new((/nv,42/),float)
; cnlev2(0,:) = ispan(-20,20,1)*4-20.
 cnlev2(0,:) = array_append_record( ispan(-30,-10,1), ispan(10,30,1), 0 )*4.-20.
 cnlev2(1,:) = array_append_record( ispan(-30,-10,1), ispan(10,30,1), 0 )*5.
 cnlev2(2,:) = array_append_record( ispan(-30,-10,1), ispan(10,30,1), 0 )*0.4

 cnfil = new((/nv,22/),integer)
 cnfil(0,:) = (/17,18,19,20,21,22,23,24,25,26,0,11,10,9,8,7,6,5,4,3,2,1/)
 cnfil(1,:) = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)
 cnfil(2,:) = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)
 
;--------------------------------------- contour
 wi = 32
 he0 = 40
 xf = 10 + (/0,1,0,1/)*40
 yf = 10 + (/0,0,1,1/)*58
 ymin = (/0,0,8,8/)
 he = he0*((50-ymin)/50.)

 ystr  = new(nv,string)

 mstr = (/"T","u","w*, (v*,w*)"/)
 xstr = "latitude (deg)"
 ystr(:) = ""
 ystr((/0,2/)) = "hybrid height (km)"
 lbtit = (/"T (:S:o:N:C)","u (m s:S:-1:N:)","w* (mm s:S:-1:N:)"/)

do iv=0, nv-1

 contour1 = cnshadeplot(wks,lat,z,vavg(iv,:,:),cnlev(iv,:),cnfil(iv,:),"")
; cnaddline(contour1,wks,lat,z,varx(iv,:,:),cnlevx,0,"D",0,"T",1.2,"C",1)
 cnaddline(contour1,wks,lat,z,vavg(iv,:,:),cnlev2(iv,:),0,"D",0,"T",0.5,"C",0)

 set_size(contour1,xf(iv),yf(iv),wi,he(iv))
 axis_range(contour1,-90,90,ymin(iv),50,False,False)
 axis_str(contour1,mstr(iv),xstr,ystr(iv))

 cnfinalize(contour1,wks)

 label = create "lab" labelBarClass wks
  "lbOrientation" : "Vertical"
  "lbBoxCount"   : dimsizes(cnfil(iv,:))
  "lbFillColors" : cnfil(iv,:)
  "lbMonoFillPattern" : True
  "lbLabelsOn" : True
  "lbLabelStrings" : cnlev(iv,:)
  "lbLabelStride" : 1
  "lbAutoManage" : False
  "lbLabelAlignment" : "InteriorEdges"
  "lbPerimOn" : False
  "lbTitleString" : "";lbtit(iv)
  "lbTitlePosition" : "Bottom"
 end create
 set_size_lb(label,xf(1)+2+iv*12,yf(nv-1)-2,6,he(nv-1)*1.3)
 draw(label)
 dr_txt(wks,lbtit(iv),xf(1)+2+iv*12+2.5,yf(nv-1)-2+(he(nv-1)*1.3)*0.9,22,5,0.)


 if (iv .eq. 2) then

  vf = create "epf_vec" vectorFieldClass defaultapp
   "vfDataArray" : vec(:,24+4:52:2,::6)
   "vfXArray"    : lat(::6)
   "vfYArray"    : z(24+4:52:2)
   "vfMissingUValueV" : 1.e20
   "vfMissingVValueV" : 1.e20
  end create
 
  vc1 = create "vector" vectorPlotClass wks
   "vcVectorFieldData" : vf
   "vcRefMagnitudeF" : 0.5
   "vcRefLengthF"    : 0.05
   "vcRefAnnoOn"     : False
 ;  "vcLineArrowThicknessF" :  thc(nn)
 ;  "vcLineArrowColor"      :  col(nn)
 ;  "vcRefAnnoString1" : ref*tmpc(nn)
   "vcRefAnnoString2On" : False
 ;  "vcRefAnnoParallelPosF" :  1. - nn*0.1
   "vcGlyphStyle" : "CurlyVector"
   "vcFillArrowEdgeThicknessF" : 0.1
   "vcFillArrowWidthF" : 0.1
   "vcPositionMode" : "ArrowTail"
   "vcLineArrowHeadMaxSizeF" : 0.05
   "vcLineArrowHeadMinSizeF" : 0.001
  end create

  set_size(vc1,xf(iv),yf(iv),wi,he(iv))
  axis_range(vc1,-90,90,ymin(iv),50,False,False)
  axis_str(vc1,"","","")
 
  overlay(contour1,vc1)

 end if

 draw(contour1)

 if (iv .eq. 2) then
   NhlRemoveOverlay(contour1,vc1,False)
 end if

end do  ; iv


 frame(wks)

end

