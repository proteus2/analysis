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

; varname = (/"f_x","f_y",param0_varname_t/)
 varname = (/"f_x","f_y","f_z"/)
 p_v = (/70,100,200,500,850/)
 if (casecode .eq. "hg") then
 if ( param0_lowlev ) then
   varname(2) = param0_varname_ll
;   delete(p_v)
;   p_v = param0_z_v_ll
 end if
 end if

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
;   fname_v = (/"tem3d_qg_st-p85gp_xyp","tem3d_qg_st-p85gp_xyp",param0_fname_t/)
   fname_v = (/"tem3d_qg_st-p85gp_xyp","tem3d_qg_st-p85gp_xyp","tem3d_qg_st-p85gp_xyp"/)
;   y9 = param0_y9
;   pp = param0_pp
;   ppi0 = param0_ppi0
 end if
 
 vnamef = "waf9"
 if (casecode .eq. "hg") then
   foutname = vnamef+"_xy-5p_"+nma+"m"+nph+"p_hg_p"+p1+"-"+p2+"_"+y1+"-"+y2
   if ( param0_lowlev ) then
     foutname = vnamef+"_xy-z"+(p_v/1.e3)+"_"+nma+"m"+nph+"p_hg_p"+p1+"-"+p2+"_"+y1+"-"+y2
   end if
 else
   if (expname .eq. "era-int_f") then
     foutname = vnamef+"_xy-5p_"+nma+"m"+nph+"p_"+cexp(0)+"_p"+pp(p1-ppi0)+"-"+pp(p2-ppi0)+"_"+y1+"-"+y2
   else
     foutname = vnamef+"_xy-5p_"+nma+"m"+nph+"p_"+cexp(0)+"_p"+p1+"-"+p2+"_"+y1+"-"+y2
   end if
 end if

;--------------------------------------- read data
 nvi = dimsizes(varname)
 nv = nvi
 nt = (y2-y1+1)*12
 nz = dimsizes(p_v)
 it0 = (y1-y0)*12
 it9 = (y2+1-y0)*12 - 1

 iv_fx = 0
 iv_fy = 1
 iv_fz = 2

 if (cexp .eq. "ERA-I" .or. cexp .eq. "ERA-I-nr") then
   fname0 = "$DATD/tem3d/ERA-I/2000/01/"+expname+"."+fname_v(iv_fx)+".2000.01.nc"
;   fname1 = "/data11/data-arch/ERA-I/monthly/2000/01/era-int."+fname_v(iv_t)+".anal.pl.200001.mm.nc"
 end if
 if (casecode .eq. "hg") then
   fname0 = "/hippo0/HG2CMIP/L60"+cexp(0)+"/pj/uan"+expname+"a.pj_"+y1+"01.nc"
 end if
 ff0 = addfile(fname0,"r")
 axisname = ff0->$varname(iv_fx)$!2
 lon = ff0->$axisname$
 axisname = ff0->$varname(iv_fx)$!1
 lat = ff0->$axisname$({10:85})
 nx = dimsizes(lon)
 ny = dimsizes(lat)

; ff1 = addfile(fname1,"r")
; axisname = ff1->$varname(iv_fx)$!2
; lon2 = ff1->$axisname$
; axisname = ff1->$varname(iv_fx)$!1
; lat2 = ff1->$axisname$({10:85})
;
; if ( any(lat .ne. lat2) .or. any(lon .ne. lon2) ) then
;   print("Grids are different.")
;   exit
; end if

 var = new((/nv,nt,nz,ny,nx/),float,1.e20)
; l_firstmissing = False
 l_firstmissing = True

 fname = new(nt,string)
 do iv=0, nvi-1

   n = 0
   do yy=y1, y2
   do m=1, 12
     mm = sprinti("%2.2i",m)
     yyyymm = yy+sprinti("%2.2i",m)
     if (cexp .eq. "ERA-I" .or. cexp .eq. "ERA-I-nr") then
       fname(n) = "$DATD/tem3d/ERA-I/"+yy+"/"+mm+"/"+  \
                  expname+"."+fname_v(iv)+"."+yy+"."+mm+".nc"
;       if (iv .eq. iv_t) then
;         fname(n) = "/data11/data-arch/"+cexp_v+"/monthly/"+yy+"/"+mm+"/"+  \
;                    expname_v+"."+fname_v(iv)+".anal.pl."+yyyymm+".mm.nc"
;       end if
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
   if (cexp .eq. "ERA-I" .or. cexp .eq. "ERA-I-nr") then
     ListSetType(ff, "join")
   end if
   if (casecode .eq. "hg") then
     ListSetType(ff, "cat")
   end if
   var(iv,:,:,:,:) = (/ ff[:]->$varname(iv)$(:,{p_v},{10:85},:) /)
   delete(ff)

 end do
 var(iv_fx,:,:,:,:) = var(iv_fx,:,:,:,:)*1.e3
 var(iv_fy,:,:,:,:) = var(iv_fy,:,:,:,:)*1.e3
 var(iv_fz,:,:,:,:) = var(iv_fz,:,:,:,:)*1.e3

 if ( l_firstmissing ) then
   var(:,0,:,:,:) = var@_FillValue
 end if

;--------------------------------------- process 0
 varc = new((/nv,12,nz,ny,nx/),float,1.e20)
 do im=0, 11
   varc(:,im,:,:,:) = dim_avg_n(var(:,im::12,:,:,:), 1)
 end do

 vara = var
; do im=0, 11
;   vara(:,im::12,:,:,:) = var(:,im::12,:,:,:) -  \
;       conform_dims(dimsizes(var(:,im::12,:,:,:)),varc(:,im,:,:,:),(/0,2,3,4/))
; end do

 nt_phs4 = new((/12/nma,nph/),integer)

 season = (/"DJF","MAM","JJA","SON"/)

 var8m = new((/12/nma,nph,nv,nz,ny,nx/),float)
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
     tmp = vara(:,its,:,:,:)
     var8m(im,ip,:,:,:,:) = dim_avg_n(tmp,1)
     std8m(im,ip,:,:,:,:) = sqrt( ( dim_sum_n(tmp^2,1) -  \
         nt_phs(ip)*var8m(im,ip,:,:,:,:)^2 )/tofloat(nt_phs(ip)-1) )
     prob(im,ip,:,:,:,:) = (var8m(im,ip,:,:,:,:) - 0.)*  \
                       sqrt(nt_phs(ip))/std8m(im,ip,:,:,:,:)
     prob(im,ip,:,:,:,:) = student_t(prob(im,ip,:,:,:,:),nt_phs(ip)-1)
     delete([/tmp,its/])
   end do
   nt_phs4(im,:) = nt_phs
   delete([/it_phs,nt_phs/])

 end do

; fcoef = ezfftf(var8m)
; fcoef(:,:,:,:,:,:,10:) = 0.
; var8m = ezfftb(fcoef, fcoef@xbar)
; var8m(:,:,iv_t,:,:,:) = wgt_runave_n_Wrap(var8m(:,:,iv_t,:,:,:),(/1.,2.,1./),0,3)

;--------------------------------------- cn level
 cnlev = fspan(-10,10,21)*0.1
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 cnlab = new((/12,nz/),string)
 cnlab = ""
; cnlab(nph-1) = "M [10:S:2:N: mPa]"
 cnlab(0,nz-1) = "M [10:S:3:N: mPa]"

 mstr = cexp

 wi = 70 ;50
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

 wks = cr_wks("fig/"+nma+"m"+(im+1)+"/"+foutname,600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.3,0.3,0.3)

do ip=0, 2 ;nph-1
do k=0, nz-1

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
  "mpMinLatF"      : 15. ;0.
  "mpMaxLatF"      : 80. ;50.
 end create
 set_size_mp(map,xf(ip),yf(k),wi,he)

 wind = create "f1" vectorFieldClass defaultapp
  "vfDataArray" : var8m(im,ip,(/iv_fx,iv_fy/),k,::4,::8) * 1.e-5 * (1.e3/p_v(k))
  "vfXArray" : lon(::8)
  "vfYArray" : lat(::4)
  "vfMissingUValueV" : 1.e15
  "vfMissingVValueV" : 1.e15
 end create
 vec1 = create "vector" vectorPlotClass wks
  "vcVectorFieldData" : wind
  "vcRefMagnitudeF" : 2.
  "vcRefLengthF"    : 0.08
  "vcRefAnnoOn"     : False
  "vcGlyphStyle"    : "FillArrow"
  "vcFillArrowHeadYF" : 0.18
  "vcFillArrowHeadXF" : 0.5
;  "vcMinMagnitudeF" : 0.25
 end create
 set_size(vec1,xf(ip),yf(k),wi,he)

; contour1 = cnshadeplot(wks,lon,lat,var8m(im,ip,iv_t,k,:,:),cnlev,cnfil,cnlab(ip))
 contour1 = cnshadeplot(wks,lon,lat,var8m(im,ip,iv_fz,k,:,:)*1.e-3*(1.e3/p_v(k)),cnlev,cnfil,cnlab(ip,k))
; cnaddline(contour1,wks,lon,lat,prob(im,ip,iv_fx,k,:,:)*100.,5.,0,"D",0,"T",1.0,"C",1)
; cnaddline(contour1,wks,lon,lat,prob(im,ip,iv_fy,k,:,:)*100.,5.,0,"D",1,"T",1.0,"C",1)
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

; overlay(map,contour1)
 overlay(map,vec1)
 NhlAddOverlay(map,contour1,map)

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

