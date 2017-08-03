load "$UTIL/ncl/header.ncl"
xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"
axisfontratio = 1.2
begin

;--------------------------------------- parameters
 load "./param-qbodata-obs.ncl"
 load "./param-1-ra.ncl"
 load "./ftn_qbo_phase.ncl"

 season = par1_season
 y1 = par1_y1
 y2 = par1_y2

 varname = par1_varname_z
 p_v = (/70,100,200,500,850/)

 cexp    = par1_expgroup_v
 expname = par1_expname_v
 fname_v = par1_fname_z

 p1 = par0_p_q(0)
 p2 = par0_p_q(1:)
 if (par0_case_q .ne. "obs") then
   lat0 = par0_latbnd_q
 end if

 vnamef = "z"
 foutname = vnamef+"_xy-5p_"+season+"_"+cexp+"_p"+p1+"-"+p2+"_"+y1+"-"+y2

;-------------------------------------------------
 nma = par1_nma

;--------------------------------------- read data
 nt = (y2-y1+1)*nma
 nz = dimsizes(p_v)

 if (cexp .eq. "ERA-I") then
   fname0 = "/data11/data-arch/"+cexp+"/2000/01/"+expname+"."+fname_v+".anal.00.pl.200001.nc"
 end if
 if (par1_case_v .eq. "hg") then
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
 if (par1_case_v .eq. "hg") then
   var@_FillValue = 0.
 end if
 l_firstmissing = False

 fname = new(nt,string)

 n = 0
 do yy=y1, y2
 do im=0, nma-1
   m = par1_imon0(par1_iseason,im)+1
   mm = sprinti("%2.2i",m)
   yyyymm = yy+sprinti("%2.2i",m)
   if (cexp .eq. "ERA-I") then
     fname(n) = "/data11/data-arch/"+cexp+"/monthly/"+yy+"/"+mm+"/"+  \
                expname+"."+fname_v+".anal.pl."+yyyymm+".mm.nc"
   end if
   if (par1_case_v .eq. "hg") then
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
 if (cexp .eq. "ERA-I") then
   ListSetType(ff, "join")
 end if
 if (par1_case_v .eq. "hg") then
   ListSetType(ff, "cat")
 end if
 var = (/ ff[:]->$varname$(:,{p_v},{-30:85},:) /)
 delete(ff)

 if (cexp .eq. "ERA-I") then
   var = var/9.81   ; ERA provides geopotential, not height
 end if

 if ( l_firstmissing ) then
   var(0,:,:,:) = var@_FillValue
 end if

;--------------------------------------- process 0
 varc = new((/nma,nz,ny,nx/),float,1.e20)
 do im=0, nma-1
   varc(im,:,:,:) = dim_avg_n(var(im::nma,:,:,:), 0)
 end do

 vara = var
 do im=0, nma-1
   vara(im::nma,:,:,:) = var(im::nma,:,:,:) -  \
       conform_dims(dimsizes(var(im::nma,:,:,:)),varc(im,:,:,:),(/1,2,3/))
 end do

 it_phs = qbo_phase_ind(season,par0_fname_q,par0_varname_q,  \
                        p1,p2,10,par0_y0_q,y1,y2)
;                        p1,p2,10,par0_y0_q,par0_y1_q,par0_y2_q)
 nph = it_phs@n_phase
 nt_phs = it_phs@nt_phase
 phs = it_phs@phases

 var8m = new((/nph,nz,ny,nx/),float)
 std8m = var8m
 prob = var8m
 do ip=0, nph-1
   its = it_phs(ip,:nt_phs(ip)-1)
   tmp = vara(its,:,:,:)
   var8m(ip,:,:,:) = dim_avg_n(tmp,0)
   std8m(ip,:,:,:) = sqrt( ( dim_sum_n(tmp^2,0) -  \
       nt_phs(ip)*var8m(ip,:,:,:)^2 )/tofloat(nt_phs(ip)-1) )
   prob(ip,:,:,:) = (var8m(ip,:,:,:) - 0.)*  \
                     sqrt(nt_phs(ip))/std8m(ip,:,:,:)
   prob(ip,:,:,:) = student_t(prob(ip,:,:,:),nt_phs(ip)-1)
   delete([/tmp,its/])
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


 wks = cr_wks("fig/"+season+"/"+foutname,600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.3,0.3,0.3)

do ip=0, nph-1
do k=0, nz-1

 cnlev = fspan(-10,10,21)*8
 if (nph .eq. 3) then
   cnlev = fspan(-10,10,21)*4
 end if
 if (p_v(k) .eq. 500) then
   cnlev = fspan(-10,10,21)*6
   if (nph .eq. 3) then
     cnlev = fspan(-10,10,21)*3
   end if
 end if
 if (p_v(k) .gt. 500) then
   cnlev = fspan(-10,10,21)*4
   if (nph .eq. 3) then
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

 contour1 = cnshadeplot(wks,lon,lat,var8m(ip,k,:,:),cnlev,cnfil,cnlab(ip))
 cnaddline(contour1,wks,lon,lat,prob(ip,k,:,:)*100.,5.,0,"D",0,"T",1.0,"C",1)
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
 dr_txt(wks,nt_phs(ip)+"",xf(ip)+wi-3,yf(k)+he-2.5,22,8,0.)
end if

end do
end do

 frame(wks)
 delete(wks)

end

