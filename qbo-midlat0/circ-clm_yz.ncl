load "$UTIL/ncl/header.ncl"
xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"
axisfontratio = 1.2
begin

;--------------------------------------- parameters
 casecode = "ra"
; casecode = "hg"

 if (casecode .eq. "ra") then
   load "./param-ra-0.ncl"
 else
   load "./param-hg-0.ncl"
 end if

 cexp    = param0_cexp  ; (/"CGW","CTL"/)
 expname = param0_expname  ;(/"tc" ,"tk" /)

 varname = (/"v_res","w_res",param0_u0varname/)

 y0 = param0_y0
 y1 = param0_y1
 y2 = param0_y2

 lat0 = param0_lat0  ; 10.
 p1 = param0_p1  ; 70
 p2 = param0_p2  ; 30

 nma = param0_nma

 if (casecode .eq. "ra") then
   cexp_v  = param0_cexp_v
   expname_v = param0_expname_v
   fname_v = param0_fname_uv
;   y9 = param0_y9
;   pp = param0_pp
;   ppi0 = param0_ppi0
 end if
 
 vnamef = "circ-clm"
 if (casecode .eq. "hg") then
   foutname = vnamef+"_yz-"+nma+"m_hg_p"+p1+"-"+p2+"_"+y1+"-"+y2
 else
   if (expname .eq. "era-int_f") then
     foutname = vnamef+"_yz-"+nma+"m_"+cexp+"_p"+pp(p1-ppi0)+"-"+pp(p2-ppi0)+"_"+y1+"-"+y2
   else
     foutname = vnamef+"_yz-"+nma+"m_"+cexp+"_p"+p1+"-"+p2+"_"+y1+"-"+y2
   end if
 end if

;--------------------------------------- read data
 nvi = dimsizes(varname)
 nv = nvi
 nt = (y2-y1+1)*12
 it0 = (y1-y0)*12
 it9 = (y2+1-y0)*12 - 1

 iv_vr = 0
 iv_wr = 1
 iv_u = 2

 if (casecode .eq. "ra") then
   fname0 = "$DATD/tem/"+cexp_v+"/2000/"+expname_v+".tem_ypt.2000.01.nc"
   fname1 = "$DATD/clim/"+cexp_v+"/monthly/2000/"+expname_v+".u_yp.200001.nc"
 end if
 if (casecode .eq. "hg") then
;   fname0 = "/hippo0/HG2CMIP/L60"+cexp+"/pj/uan"+expname+"a.pj_"+y1+"01.nc"
 end if
 ff0 = addfile(fname0,"r")
 lat = ff0->lat({-30:85})
 p   = ff0->p({500:10})
 z   = 7.0*log(1000./p)
 ny = dimsizes(lat)
 nz = dimsizes(p)

 ff1 = addfile(fname1,"r")
 axisname = ff1->$varname(iv_u)$!2
 lat2 = ff1->$axisname$({-30:85})
 axisname = ff1->$varname(iv_u)$!1
 p2a  = ff1->$axisname$({500:10})

 if ( any(lat .ne. lat2) .or. any(p .ne. p2a) ) then
   print("Grids are different.")
   exit
 end if

 var = new((/nv,nt,nz,ny/),float,1.e32)
 l_firstmissing = False

 fname = new(nt,string)
 fnameu = new(nt,string)

 do iv=0, nvi-1

   n = 0
   do yy=y1, y2
   do m=1, 12
     mm = sprinti("%2.2i",m)
     yyyymm = yy+sprinti("%2.2i",m)
     if (casecode .eq. "ra") then
       fname(n) = "$DATD/tem/"+cexp_v+"/"+yy+"/"+expname_v+".tem_ypt."+yy+"."+mm+".nc"
       if (iv .eq. iv_u) then
         fname(n) = "$DATD/clim/"+cexp_v+"/monthly/"+yy+"/"+expname_v+".u_yp."+yyyymm+".nc"
       end if
     end if
     if (casecode .eq. "hg") then
       fname(n) = "/hippo0/HG2CMIP/L60"+cexp+"/pj/uan"+expname+"a.pj_"+yyyymm+".nc"
     end if
     if ( yyyymm .eq. "197901" .and. .not. isfilepresent(fname(n)) ) then
       l_firstmissing = True
     end if
     n = n + 1
   end do
   end do

   if ( l_firstmissing ) then
     fname(0) = fname(1)
   end if

   ff = addfiles(fname,"r")
   if (casecode .eq. "ra") then
     if (iv .eq. iv_u) then
       ListSetType(ff, "cat")
     else
       ListSetType(ff, "join")
     end if
   end if
   if (casecode .eq. "hg") then
     ListSetType(ff, "cat")
   end if
   var(iv,:,:,:) = (/ ff[:]->$varname(iv)$(:,{500:10},{-30:85}) /)
   delete([/ff/])

 end do

 if ( l_firstmissing ) then
   var(:,0,:,:) = var@_FillValue
 end if

;--------------------------------------- process 0
 if (nma .eq. 3) then
   ; Dec. to the first
   tmpv = var
   var(:,0,:,:) = tmpv(:,nt-1,:,:)
   var(:,1:,:,:) = tmpv(:,:nt-2,:,:)
   delete(tmpv)
 end if

 varc = new((/nv,12,nz,ny/),float,1.e20)
 do im=0, 11
   varc(:,im,:,:) = dim_avg_n(var(:,im::12,:,:), 1)
 end do

 var8m = new((/12/nma,nv,nz,ny/),float)
 do im=0, 12/nma-1
   imon = im*nma
   var8m(im,:,:,:) = dim_avg_n(varc(:,imon:imon+nma-1,:,:),1)
 end do

;--------------------------------------- cn level
 cnlev = ispan(-10,10,1)*5

; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 cnlab = new(12,string)
 cnlab(:) = ""
 cnlab(0) = "U [m s:S:-1:N:]"

 mstr = cexp
;p xstr = (/"","lat [deg]","","lat [deg]"/)
;p ystr = (/"z [km]","z [km]","",""/)
 xstr = (/"lat [deg]","lat [deg]","lat [deg]","lat [deg]"/)
 ystr = (/"z [km]","","",""/)

 wi = 50
 he = 30
 xf = 15
 yf = 10


var8m(:,iv_wr,:,:) = var8m(:,iv_wr,:,:)*111.

do im=0, 12/nma-1

 wks = cr_wks("fig/"+nma+"m"+(im+1)+"/"+foutname,600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.3,0.3,0.3)

 vw = create "f1" vectorFieldClass defaultapp
  "vfDataArray" : var8m(im,(/iv_vr,iv_wr/),:,::4)
  "vfXArray" : lat(::4)
  "vfYArray" : z
;  "vfMissingUValueV" : 1.e15
;  "vfMissingVValueV" : 1.e15
 end create
 vec1 = create "vector" vectorPlotClass wks
  "vcVectorFieldData" : vw
  "vcRefMagnitudeF" : 1.0
  "vcRefLengthF"    : 0.05
  "vcRefAnnoOn"     : True
  "vcGlyphStyle"    : "FillArrow"
  "vcFillArrowHeadYF" : 0.18
  "vcFillArrowHeadXF" : 0.5
  "vcMinMagnitudeF" : 0.05

;  "vcRefAnnoArrowAngleF" : 90.
;  "vcRefAnnoFontHeightF" : 0.024
  "vcRefAnnoString1On"   : False
;  "vcRefAnnoString1"     : "1 mPa"
;  "vcRefAnnoString2On"   : False
  "vcRefAnnoString2"     : "1 m s:S:-1:N:"
;  "vcRefAnnoPerimSpaceF" : 0.6
;  "vcRefAnnoOrthogonalPosF" : 0.15
;  "vcRefAnnoParallelPosF"   :  1. - nn*0.1

 end create
 set_size(vec1,xf(0),yf(0),wi,he)
 axis_range(vec1,-25,80,-999,30,False,False)

 contour1 = cnshadeplot(wks,lat,z,var8m(im,iv_u,:,:),  \
            cnlev,cnfil,cnlab(0))
 setvalues contour1
  "cnMissingValFillColor" : colind
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.03
 end setvalues
 set_size(contour1,xf(0),yf(0),wi,he)
 axis_range(contour1,-25,80,-999,30,False,False)

 cnfinalize(contour1,wks)

 overlay(contour1,vec1)

 draw(contour1)
; cnfinalize(vec1,wks)
; draw(vec1)

 frame(wks)
 delete(wks)

end do

end

