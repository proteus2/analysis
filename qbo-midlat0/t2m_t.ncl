load "$UTIL/ncl/header.ncl"
xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"
axisfontratio = 1.2
begin

 load "./ftn_qbo_phase.ncl"

;--------------------------------------- parameters
 season = "DJF"  ; "JJA" ; "DJF"
; lons = (/57,80/)
; lats = (/45,65/)
 lons = (/55,85/)
 lats = (/45,70/)

 casecode = "ra"
; casecode = "hg"

 if (casecode .eq. "ra") then
   load "./param-ra-0.ncl"
 else
   load "./param-hg-0.ncl"
 end if

 cexp    = param0_cexp  ; (/"CGW","CTL"/)
 expname = param0_expname  ;(/"tc" ,"tk" /)

 varname = param0_varname_t2m

 y0 = param0_y0
 y1 = param0_y1
 y2 = param0_y2

 lat0 = param0_lat0  ; 10.
 p1 = param0_p1  ; 70
 p2 = param0_p2  ; 30

 nph = param0_nph

 if (casecode .eq. "ra") then
   cexp_v  = param0_cexp_v
   expname_v = param0_expname_v
   fname_v = param0_fname_t
   y9 = param0_y9
   pp = param0_pp
   ppi0 = param0_ppi0
 end if
 
 vnamef = "t2m"
 if (casecode .eq. "hg") then
   foutname = vnamef+"_t-"+season+"."+nph+"p_hg_p"+p1+"-"+p2+"_"+y1+"-"+y2
   if ( param0_lowlev ) then
;     foutname = vnamef+"_t-z"+(p_v/1.e3)+"_"+season+"."+nph+"p_hg_p"+p1+"-"+p2+"_"+y1+"-"+y2
   end if
 else
   if (expname(0) .eq. "era-int_f") then
     foutname = vnamef+"_t-"+season+"."+nph+"p_"+cexp(0)+"_p"+pp(p1-ppi0)+"-"+pp(p2-ppi0)+"_"+y1+"-"+y2
   else
     foutname = vnamef+"_t-"+season+"."+nph+"p_"+cexp(0)+"_p"+p1+"-"+p2+"_"+y1+"-"+y2
   end if
 end if

;--------------------------------------- read data
 it_phs = qbo_phase_ind(season,nph,param0_fname_u0,param0_u0varname,  \
                        p1,p2,lat0,(/param0_umin_p1,param0_umax_p1/),  \
                                   (/param0_umin_p2,param0_umax_p2/),  \
                        y0,y1,y2)
 nt_phs = it_phs@nt_phase
 phs    = it_phs@phases
 deg    = where(it_phs@degree .ge. 360-45, it_phs@degree - 360., it_phs@degree)
 it_all = it_phs@i_allphases
 nts = dimsizes(it_all)

 ne = dimsizes(expname)
 nt = (y2-y1+1)*12
 it0 = (y1-y0)*12
 it9 = (y2+1-y0)*12 - 1

 if (cexp .eq. "ERA-I" .or. cexp .eq. "ERA-I-nr") then
   fname = "$DATD/clim/"+cexp_v+"/monthly/t2_2.5deg/t2.mon.2000.nc"
 end if
 if (casecode .eq. "hg") then
   fname = "/hippo0/HG2CMIP/L60"+cexp(0)+"/pj/uan"+expname(0)+"a.pj_"+y1+"01.nc"
 end if
 ff = addfile(fname,"r")
 axisname = ff->$varname$!2
 lon = ff->$axisname$ ;({lons(0):lons(1)})
 axisname = ff->$varname$!1
 lat = ff->$axisname$ ;({lats(0):lats(1)})
 nx = dimsizes(lon)
 ny = dimsizes(lat)

 var = new((/ne,nt,ny,nx/),float,1.e20)
 fname_dir = "$DATD/clim/"+cexp_v+"/monthly/t2_2.5deg"
 fnames = fname_dir+"/t2.mon."+ispan(y1,y2,1)+".nc"
 ffs = addfiles(fnames,"r")
 ListSetType(ffs, "cat")
 ie = 0
 var(ie,:,:,:) = ffs[:]->$varname$ ;(:,{lats(0):lats(1)},{lons(0):lons(1)})

 coslat = cos(lat({lats(0):lats(1)})*3.141592/180.)
 varl = dim_avg_wgt_Wrap(dim_avg_Wrap(var(:,:,{lats(0):lats(1)},{lons(0):lons(1)})),coslat,0)
 
;--------------------------------------- process 0
 varlc = new((/ne,12/),float,1.e20)
 varc = new((/ne,12,ny,nx/),float,1.e20)
 do im=0, 11
   varlc(:,im) = dim_avg(varl(:,im::12))
   varc(:,im,:,:) = dim_avg_n(var(:,im::12,:,:), 1)
 end do

 varla = varl
 vara = var
 do im=0, 11
   varla(:,im::12) = varl(:,im::12) -  \
       conform_dims(dimsizes(varl(:,im::12)),varlc(:,im),0)
   vara(:,im::12,:,:) = var(:,im::12,:,:) -  \
       conform_dims(dimsizes(var(:,im::12,:,:)),varc(:,im,:,:),(/0,2,3/))
 end do

 var8m = new((/ne,nph,ny,nx/),float)
 std8m = var8m
 prob = var8m
 ie = 0
 varla_s = varla(ie,it_all)
 do ip=0, nph-1
   its = it_phs(ip,:nt_phs(ip)-1)
   tmp = vara(ie,its,:,:)
   var8m(ie,ip,:,:) = dim_avg_n(tmp,0)
   std8m(ie,ip,:,:) = sqrt( ( dim_sum_n(tmp^2,0) -  \
       nt_phs(ip)*var8m(ie,ip,:,:)^2 )/tofloat(nt_phs(ip)-1) )
   prob(ie,ip,:,:) = (var8m(ie,ip,:,:) - 0.)*  \
                     sqrt(nt_phs(ip))/std8m(ie,ip,:,:)
   prob(ie,ip,:,:) = student_t(prob(ie,ip,:,:),nt_phs(ip)-1)
   delete([/tmp,its/])
 end do

;--------------------------------------- cn level
 cnlev = fspan(-10,10,21)*0.5
 if (nph .eq. 4) then
   cnlev = fspan(-10,10,21)*0.2
 end if
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- contour
 cnlab = new(12,string)
 cnlab = ""
 cnlab(nph-1) = "T [K]"

 mstr = cexp

 wi = 50
 he = 25
 xf = new(nph,float)
;p yf = new(nph,float)
 xf = 15
 yf = 10
 wi2 = 100
 yf2 = yf + he + 10
;p xf(nph/2:) = xf(nph/2:) + (wi+8)
;p yf(:nph/2-1) = yf(:nph/2-1) + ispan(0,nph/2-1,1)*(he+10)
;p yf(nph/2:  ) = yf(nph/2:  ) + ispan(0,nph/2-1,1)*(he+10)
 xf = 15 + ispan(0,nph-1,1)*(wi+4)


 wks = cr_wks("fig/"+season+"/"+foutname,600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.3,0.3,0.3)
 colgrn = NhlNewColor(wks,0.0,1.0,0.0)
 colred = NhlNewColor(wks,1.0,0.0,0.0)
 colblu = NhlNewColor(wks,0.0,0.0,1.0)
 colcya = NhlNewColor(wks,0.0,1.0,1.0)

do ip=0, 2 ;0, nph-1

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
  "mpMinLatF"      : 20. ;0.
  "mpMaxLatF"      : 85. ;50.
 end create
 set_size_mp(map,xf(ip),yf,wi,he)

 contour1 = cnshadeplot(wks,lon,lat,var8m(0,ip,:,:),cnlev,cnfil,cnlab(ip))
 cnaddline(contour1,wks,lon,lat,prob(0,ip,:,:)*100.,5.,0,"D",0,"T",1.0,"C",1)
 setvalues contour1
  "cnMissingValFillColor" : colind
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.06
 end setvalues
 set_size(contour1,xf(ip),yf,wi,he)
 axis_str(contour1,"phase "+(ip+1),"","")

 xpg = (/lons(0),lons(1),lons(1),lons(0),lons(0)/)
 ypg = (/lats(0),lats(0),lats(1),lats(1),lats(0)/)
 resPG = True
 resPG@gsLineColor = colgrn
 dumbx = gsn_add_polyline(wks,contour1,xpg,ypg,resPG)

 overlay(map,contour1)
; overlay(map,vec1)
; NhlAddOverlay(map,contour1,map)

 draw(map)

end do

 t = ispan(0,nts-1,1)

 col = (/0,colred,colgrn,colblu,colcya/)
 col(4) = colind

 plot = xylineplot(wks,deg,varla_s,"D",0,"T",1.0,"C",-1)
 do ip=0, nph-1
   indphs = ind(phs .eq. ip+1)
   xyaddline(plot,deg(indphs),varla_s(indphs),"D",0,"T",1.0,"C",-1)
   delete(indphs)
 end do
 xyaddmarker(plot,16,"S",0.8,"T",1.0,"C",col)
 do ip=0, nph-1
   indphs = ind(phs .eq. ip+1)
   xyaddline(plot,deg(indphs),avg(varla_s(indphs)),"D",0,"T",1.0,"C",col(ip+1))
   delete(indphs)
 end do
 set_size(plot,xf(0),yf2,wi2,he)
 axis_range(plot,-45,360-45,-999,999,False,False)
 axis_str(plot,"","phase [degree]","T anomaly [K]")

; dr_txt(wks,p(ip)+" hPa",xf(0)+8,yf(ip)+3,22,8,0.)

 draw(plot)

 frame(wks)

end

