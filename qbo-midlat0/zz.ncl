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
;   load "./param-ra-0.ncl"
   load "./param-ra2-0.ncl"
 else
   load "./param-hg-0.ncl"
 end if

 cexp    = param0_cexp  ; (/"CGW","CTL"/)
 expname = param0_expname  ;(/"tc" ,"tk" /)

 varname = param0_varname_z
 p_v = (/70,100,200,500,850/)

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
   fname_v = param0_fname_z
;   y9 = param0_y9
;   pp = param0_pp
;   ppi0 = param0_ppi0
 end if
 
 vnamef = "z"
 if (casecode .eq. "hg") then
   foutname = vnamef+"_xy-5p_"+nma+"m"+nph+"p_"+cexp+"_p"+p1+"-"+p2+"_"+y1+"-"+y2
 else
   if (expname .eq. "era-int_f") then
     foutname = vnamef+"_xy-5p_"+nma+"m"+nph+"p_"+cexp+"_p"+pp(p1-ppi0)+"-"+pp(p2-ppi0)+"_"+y1+"-"+y2
   else
     foutname = vnamef+"_xy-5p_"+nma+"m"+nph+"p_"+cexp+"_p"+p1+"-"+p2+"_"+y1+"-"+y2
   end if
 end if

;--------------------------------------- read data
 nt = (y2-y1+1)*12
 nz = dimsizes(p_v)
 it0 = (y1-y0)*12
 it9 = (y2+1-y0)*12 - 1

 if (cexp .eq. "ERA-I" .or. cexp .eq. "ERA-I-nr") then
   fname0 = "/data11/data-arch/ERA-I/2000/01/era-int."+fname_v+".anal.00.pl.200001.nc"
 else if (cexp .eq. "JRA55") then
   fname0 = "/data11/data-arch/JRA55/2000/01/jra55."+fname_v+".anal.pl.200001.nc"
 end if
 end if
 if (casecode .eq. "hg") then
   fname0 = "/hippo0/HG2AMIP/L60"+cexp+"/pj/uan"+expname+"a.pj_"+y1+"01.nc"
 end if
 ff0 = addfile(fname0,"r")
 axisname = ff0->$varname$!3
 lon = ff0->$axisname$
 axisname = ff0->$varname$!2
 lat = ff0->$axisname$({-30:85})
 nx = dimsizes(lon)
 ny = dimsizes(lat)

 var = new((/nt,nz,ny,nx/),float,1.e20)
 if (casecode .eq. "hg") then
   var@_FillValue = 0.
 end if
 l_firstmissing = False

 fname = new(nt,string)

 n = 0
 do yy=y1, y2
 do m=1, 12
   mm = sprinti("%2.2i",m)
   yyyymm = yy+sprinti("%2.2i",m)
   if ( cexp .eq. "ERA-I" .or. cexp .eq. "ERA-I-nr" .or.  \
        cexp .eq. "JRA55" ) then
     fname(n) = "/data11/data-arch/"+cexp_v+"/monthly/"+yy+"/"+mm+"/"+  \
                expname_v+"."+fname_v+".anal.pl."+yyyymm+".mm.nc"
   end if
   if (casecode .eq. "hg") then
     fname(n) = "/hippo0/HG2AMIP/L60"+cexp+"/pj/uan"+expname+"a.pj_"+yyyymm+".nc"
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
 if ( cexp .eq. "ERA-I" .or. cexp .eq. "ERA-I-nr" .or.  \
      cexp .eq. "JRA55" ) then
   ListSetType(ff, "join")
 end if
 if (casecode .eq. "hg") then
   ListSetType(ff, "cat")
 end if
 var = (/ ff[:]->$varname$(:,{p_v},{-30:85},:) /)
 delete(ff)

 if (cexp .eq. "ERA-I" .or. cexp .eq. "ERA-I-nr") then
   var = var/9.81   ; ERA provides geopotential, not height
 end if

 if ( l_firstmissing ) then
   var(0,:,:,:) = var@_FillValue
 end if

;--------------------------------------- process 0
 varc = new((/12,nz,ny,nx/),float,1.e20)
 do im=0, 11
   varc(im,:,:,:) = dim_avg_n(var(im::12,:,:,:), 0)
 end do

 vara = var
 do im=0, 11
   vara(im::12,:,:,:) = var(im::12,:,:,:) -  \
       conform_dims(dimsizes(var(im::12,:,:,:)),varc(im,:,:,:),(/1,2,3/))
 end do

 nt_phs4 = new((/12/nma,nph/),integer)

 season = (/"DJF","MAM","JJA","SON"/)

 var8m = new((/12/nma,nph,nz,ny,nx/),float)
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
     tmp = vara(its,:,:,:)
     var8m(im,ip,:,:,:) = dim_avg_n(tmp,0)
     std8m(im,ip,:,:,:) = sqrt( ( dim_sum_n(tmp^2,0) -  \
         nt_phs(ip)*var8m(im,ip,:,:,:)^2 )/tofloat(nt_phs(ip)-1) )
     prob(im,ip,:,:,:) = (var8m(im,ip,:,:,:) - 0.)*  \
                       sqrt(nt_phs(ip))/std8m(im,ip,:,:,:)
     prob(im,ip,:,:,:) = student_t(prob(im,ip,:,:,:),nt_phs(ip)-1)
     delete([/tmp,its/])
   end do
   nt_phs4(im,:) = nt_phs
   delete([/it_phs,nt_phs/])

 end do

;--------------------------------------- cn level
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 cnlab = new(12,string)
 cnlab(:) = ""
 cnlab(nph-1) = "Z [m]"

 mstr = cexp

 wi = 50
 he = 25
 xf = new(nph,float)
;p yf = new(nph,float)
 yf = new(nz,float)
 xf = 15
 yf = 10
;p xf(nph/2:) = xf(nph/2:) + (wi+8)
;p yf(:nph/2-1) = yf(:nph/2-1) + ispan(0,nph/2-1,1)*(he+10)
;p yf(nph/2:  ) = yf(nph/2:  ) + ispan(0,nph/2-1,1)*(he+10)
 xf = 15 + ispan(0,nph-1,1)*(wi+4)
 yf = 10 + ispan(0,nz-1,1)*(he+8)


do im=0, 12/nma-1

; wks = cr_wks("fig/"+nma+"m"+(im+1)+"/"+foutname,600,0)
 wks = cr_wks("zz-"+im,600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.3,0.3,0.3)

do ip=0, nph-1
do k=0, nz-1

 cnlev = fspan(-10,10,21)*8
 if (nph .eq. 4) then
   cnlev = fspan(-10,10,21)*4
 end if
 if (p_v(k) .eq. 500) then
   cnlev = fspan(-10,10,21)*6
   if (nph .eq. 4) then
     cnlev = fspan(-10,10,21)*3
   end if
 end if
 if (p_v(k) .gt. 500) then
   cnlev = fspan(-10,10,21)*4
   if (nph .eq. 4) then
     cnlev = fspan(-10,10,21)*2
   end if
 end if

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
 set_size_mp(map,xf(ip),yf(k),wi,he)

 contour1 = cnshadeplot(wks,lon,lat,var8m(im,ip,k,:,:),cnlev,cnfil,cnlab(ip))
 cnaddline(contour1,wks,lon,lat,prob(im,ip,k,:,:)*100.,5.,0,"D",0,"T",1.0,"C",1)
 setvalues contour1
  "cnMissingValFillColor" : colind
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.06
 end setvalues
 set_size(contour1,xf(ip),yf(k),wi,he)
if (k .eq. 0) then
 axis_str(contour1,"phase "+(ip+1),"","")
end if

 overlay(map,contour1)
; overlay(map,vec1)
; NhlAddOverlay(map,contour1,map)

 draw(map)

if (k .eq. 0) then
 dr_txt(wks,nt_phs4(im,ip)+"",xf(ip)+wi-3,yf(k)+he-2.5,22,8,0.)
end if

end do
end do

 frame(wks)
 delete(wks)

end do

end

