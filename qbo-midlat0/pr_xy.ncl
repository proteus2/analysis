load "$UTIL/ncl/header.ncl"
xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"
axisfontratio = 1.2
begin

 load "./ftn_qbo_phase.ncl"

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

 varname = param0_varname_pr

 y0 = param0_y0
 y1 = param0_y1
 y2 = param0_y2

 lat0 = param0_lat0  ; 10.
 p1 = param0_p1  ; 70
 p2 = param0_p2  ; 30

 nma = param0_nma
 nph = param0_nph

 if (casecode .eq. "ra") then
   cexp_v  = param0_cexp_v
   expname_v = param0_expname_v
   fname_v = param0_fname_pr
   y9 = param0_y9
   pp = param0_pp
   ppi0 = param0_ppi0
 end if
 
 vnamef = "pr"
 if (casecode .eq. "hg") then
   foutname = vnamef+"_xy_"+nma+"m"+nph+"p_hg_p"+p1+"-"+p2+"_"+y1+"-"+y2
 else
   if (expname(0) .eq. "era-int_f") then
     foutname = vnamef+"_xy_"+nma+"m"+nph+"p_"+cexp(0)+"_p"+pp(p1-ppi0)+"-"+pp(p2-ppi0)+"_"+y1+"-"+y2
   else
     foutname = vnamef+"_xy_"+nma+"m"+nph+"p_"+cexp(0)+"_p"+p1+"-"+p2+"_"+y1+"-"+y2
   end if
 end if

;--------------------------------------- read data
 nt = (y2-y1+1)*12
 it0 = (y1-y0)*12
 it9 = (y2+1-y0)*12 - 1

 if (cexp .eq. "ERA-I" .or. cexp .eq. "ERA-I-nr") then
   fname = "/data11/data-arch/ERA-I/2000/01/era-int."+fname_v+".fcst.00_12.sfc.200001.nc"
 end if
 if (casecode .eq. "hg") then
   fname = "/hippo0/HG2CMIP/L60"+cexp(0)+"/pj/uan"+expname(0)+"a.pj_"+y1+"01.nc"
 end if
 ff = addfile(fname,"r")
 axisname = ff->$varname$!2
 lon = ff->$axisname$
 axisname = ff->$varname$!1
 lat = ff->$axisname$
 nx = dimsizes(lon)
 ny = dimsizes(lat)

 var = new((/nt,ny,nx/),float,1.e20)

 fname00 = new(nt,string)
 fname12 = new(nt,string)

 n = 0
 do yy=y1, y2
 do m=1, 12
   mm = sprinti("%2.2i",m)
   yyyymm = yy+sprinti("%2.2i",m)
   if (cexp .eq. "ERA-I" .or. cexp .eq. "ERA-I-nr") then
     fname_dir = "/data11/data-arch/"+cexp_v+"/monthly/"+yy+"/"+mm
     fname00(n) = fname_dir+"/"+  \
                  expname_v+"."+fname_v+".fcst.00_12.sfc."+yyyymm+".mm.nc"
     fname12(n) = fname_dir+"/"+  \
                  expname_v+"."+fname_v+".fcst.12_12.sfc."+yyyymm+".mm.nc"
   end if
   if (casecode .eq. "hg") then
     fname00(n) = "/hippo0/HG2CMIP/L60"+cexp+"/pj/uan"+expname+"a.pj_"+yyyymm+".nc"
   end if
   n = n + 1
 end do
 end do

 if (cexp .eq. "ERA-I" .or. cexp .eq. "ERA-I-nr") then
   ff00 = addfiles(fname00,"r")
   ListSetType(ff00, "join")
   var = ff00[:]->$varname$
   ff12 = addfiles(fname12,"r")
   ListSetType(ff12, "join")
   var = var + ff12[:]->$varname$
   var = var*1.e3 * 30.  ; [mm/mon]
 end if
 if (casecode .eq. "hg") then
   ff = addfiles(fname00,"r")
   ListSetType(ff, "cat")
   var = ff[:]->$varname$(:,:,:)
 end if
 
;--------------------------------------- process 0
 varc = new((/12,ny,nx/),float,1.e20)
 do im=0, 11
   varc(im,:,:) = dim_avg_n(var(im::12,:,:), 0)
 end do

 vara = var
 do im=0, 11
   vara(im::12,:,:) = var(im::12,:,:) -  \
       conform_dims(dimsizes(var(im::12,:,:)),varc(im,:,:),(/1,2/))
 end do

 nt_phs4 = new((/12/nma,nph/),integer)

 season = (/"DJF","MAM","JJA","SON"/)

 var8m = new((/12/nma,nph,ny,nx/),float)
 std8m = var8m
 prob = var8m
 do im=0, 12/nma-1

   it_phs = qbo_phase_ind(season(im),nph,param0_fname_u0,param0_u0varname,  \
                          p1,p2,lat0,(/param0_umin_p1,param0_umax_p1/),  \
                                     (/param0_umin_p2,param0_umax_p2/),  \
                          y0,y1,y2)
   nt_phs = it_phs@nt_phase

   do ip=0, nph-1
     its = it_phs(ip,:nt_phs(ip)-1)
     tmp = vara(its,:,:)
     var8m(im,ip,:,:) = dim_avg_n(tmp,0)
     std8m(im,ip,:,:) = sqrt( ( dim_sum_n(tmp^2,0) -  \
         nt_phs(ip)*var8m(im,ip,:,:)^2 )/tofloat(nt_phs(ip)-1) )
     prob(im,ip,:,:) = (var8m(im,ip,:,:) - 0.)*  \
                       sqrt(nt_phs(ip))/std8m(im,ip,:,:)
     prob(im,ip,:,:) = student_t(prob(im,ip,:,:),nt_phs(ip)-1)
     delete([/tmp,its/])
   end do
   nt_phs4(im,:) = nt_phs
   delete([/it_phs,nt_phs/])

 end do

;--------------------------------------- cn level
 cnlev = fspan(-10,10,21)*15
 if (nph .eq. 4) then
   cnlev = fspan(-10,10,21)*8
 end if
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)

;--------------------------------------- contour
 cnlab = new(12,string)
 cnlab(:) = ""
 cnlab(nph-1) = "Prec. [mm mon:S:-1:N:]"

 mstr = cexp

 wi = 60
 he = 30
 xf = new(nph,float)
 yf = new(nph,float)
 xf = 15
 yf = 10
 xf(nph/2:) = xf(nph/2:) + (wi+8)
 yf(:nph/2-1) = yf(:nph/2-1) + ispan(0,nph/2-1,1)*(he+10)
 yf(nph/2:  ) = yf(nph/2:  ) + ispan(0,nph/2-1,1)*(he+10)


do im=0, 12/nma-1

 wks = cr_wks("fig/"+nma+"m"+(im+1)+"/"+foutname,600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.3,0.3,0.3)

do ip=0, nph-1

 map = create "mp" mapPlotClass wks
;  "mpProjection"     : "Mollweide" ;"Robinson"
  "mpCenterLonF"     : 125.
  "mpOutlineOn"      : True
  "mpGridAndLimbOn"  : False
  "mpPerimOn"        : True
  "mpPerimDrawOrder" : "PostDraw"
  "mpGeophysicalLineThicknessF" : 0.8

;  "mpShapeMode"     : "FixedAspectFitBB"
  "mpShapeMode"     : "FreeAspect"
  "mpGeophysicalLineColor" : colind ;"black"  ;"Orange4"
  "mpDataBaseVersion" : "Ncarg4_0"

  "mpLimitMode"    : "LatLon"
  "mpMinLonF"      : 0. ;85.
  "mpMaxLonF"      : 360. ;165.
  "mpMinLatF"      : -25. ;0.
  "mpMaxLatF"      : 80. ;50.
 end create
 set_size_mp(map,xf(ip),yf(ip),wi,he)

 contour1 = cnshadeplot(wks,lon,lat,var8m(im,ip,:,:),cnlev,cnfil,cnlab(ip))
 cnaddline(contour1,wks,lon,lat,prob(im,ip,:,:)*100.,5.,0,"D",0,"T",1.0,"C",1)
 setvalues contour1
  "cnMissingValFillColor" : colind
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.06
 end setvalues
 set_size(contour1,xf(ip),yf(ip),wi,he)
 axis_str(contour1,"phase "+(ip+1),"","")

 overlay(map,contour1)
; overlay(map,vec1)
; NhlAddOverlay(map,contour1,map)

 draw(map)

 dr_txt(wks,nt_phs4(im,ip)+"",xf(ip)+wi-3,yf(ip)+he-2.5,22,8,0.)

end do

 frame(wks)
 delete(wks)

end do

end

