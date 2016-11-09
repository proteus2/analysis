load "header.ncl"

xspace = 170
yspace = 240
load "figsize.ncl"
load "attribute.ncl"

axisfontratio = 1.0

begin

;--------------------------------------- parameters
 year1 = stringtoint(getenv("year1"))
 year2 = stringtoint(getenv("year2"))
 yearf = stringtoint(getenv("yearf"))

;--------------------------------------- data reading
 nt = year2 - year1 + 1
 time = fspan(year1,year2,nt)

 easmi = new(nt,float)
 easmi = 0.

 n = 0
 do iy=year1, year2
   do im=6, 8
     file_u = "../dat/"+iy+"/merra.u_m.assm.pl"+850+"."+ \
              iy+sprinti("%0.2i",im)+".nc"
     file_v = "../dat/"+iy+"/merra.v_m.assm.pl"+850+"."+ \
              iy+sprinti("%0.2i",im)+".nc"
     if ( .not. (isfilepresent(file_u) .and. isfilepresent(file_v)) ) then
       print("Some files were not downloaded.")
       status_exit(20)
     end if
     fu = addfile(file_u,"r")
     fv = addfile(file_v,"r")
     if ( .not. isdefined("um") ) then
       um = tofloat(fu->u(0,0,:,:))
       vm = tofloat(fv->v(0,0,:,:))
       um = 0.
       vm = 0.
       lon = tofloat(fu->u&longitude)
       lat = tofloat(fu->u&latitude )
     end if
     ; calculate EASMI
     u = fu->u(0,0,:,:)
     v = fv->v(0,0,:,:)
     easmi(n) = easmi(n) + ( avg(u({22.5:32.5},{110:140})) - \
                             avg(u({ 5.0:15.0},{ 90:130})) )/3.
     um = um + u/(tofloat(nt*3))
     vm = vm + v/(tofloat(nt*3))
   end do
   n = n + 1
 end do
print(avg(easmi))

 easmi = easmi - avg(easmi)

 do im=6, 8
   file_u = "../dat/"+yearf+"/merra.u_m.assm.pl"+850+"."+ \
            yearf+sprinti("%0.2i",im)+".nc"
   file_v = "../dat/"+yearf+"/merra.v_m.assm.pl"+850+"."+ \
            yearf+sprinti("%0.2i",im)+".nc"
   if ( .not. (isfilepresent(file_u) .and. isfilepresent(file_v)) ) then
     print("Some files were not downloaded.")
     status_exit(20)
   end if
   fu = addfile(file_u,"r")
   fv = addfile(file_v,"r")
   if ( .not. isdefined("um1") ) then
     um1 = tofloat(fu->u(0,0,:,:))
     vm1 = tofloat(fv->v(0,0,:,:))
     um1 = 0.
     vm1 = 0.
   end if
   u = fu->u(0,0,:,:)
   v = fv->v(0,0,:,:)
   um1 = um1 + u/3.
   vm1 = vm1 + v/3.
 end do

;--------------------------------------- process 0
 um1 = um1 - um
 vm1 = vm1 - vm

 wnd = (/um(::5,::5),vm(::5,::5)/)
 wind = create "f1" vectorFieldClass defaultapp
  "vfDataArray" : wnd
  "vfXArray" : lon(::5)
  "vfYArray" : lat(::5)
  "vfMissingUValueV" : 1.e15
  "vfMissingVValueV" : 1.e15
 end create

 wnd1 = (/um1(::5,::5),vm1(::5,::5)/)
 wind1 = create "f1" vectorFieldClass defaultapp
  "vfDataArray" : wnd1
  "vfXArray" : lon(::5)
  "vfYArray" : lat(::5)
  "vfMissingUValueV" : 1.e15
  "vfMissingVValueV" : 1.e15
 end create

;--------------------------------------- cn level
 cnlev = fspan(-10,15,26)*1.0
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2,16,15,14,13,12/)

;--------------------------------------- drawing
 wi = (/nt*2,36/)
 he = (/20,30/)
 xf = (/10,10,50,90/)
 yf = (/10,42/)

 outfile = getenv("out1")
 wks = cr_wks(outfile,600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.6,0.6,0.6)
 plot = xylineplot(wks,time,easmi,"D",0,"T",1.2,"C",1)
 xyaddmarker(plot,16,"S",0.5,"T",1.0,"C",1)
 xyaddline(plot,time,0.,"D",0,"T",0.6,"C",colind)
 set_size(plot,xf(0),yf(0),wi(0),he(0))
 axis_range(plot,year1-0.5,year2+0.5,-999,999,False,False)
 setvalues plot
  "tmXBMode" : "Manual"
  "tmXBTickStartF"    : 1975
  "tmXBTickSpacingF"  : 5
  "tmXBMinorPerMajor" : 4
 end setvalues
 axis_str(plot,"East Asian Summer Monsoon Index","year","(m s:S:-1:N:)")
 draw(plot)

 map = create "mp" mapPlotClass wks
  "mpProjection"     : "Mollweide" ;"Robinson"
  "mpCenterLonF"     : 125.
  "mpOutlineOn"      : True
  "mpGridAndLimbOn"  : False
  "mpPerimOn"        : True
  "mpPerimDrawOrder" : "PostDraw"

  "mpShapeMode"     : "FixedAspectFitBB"
  "mpGeophysicalLineColor" : "Orange4"
  "mpDataBaseVersion" : "Ncarg4_0"

  "mpLimitMode"    : "LatLon"
  "mpMinLonF"      : 85.;90.
  "mpMaxLonF"      : 165.;140.
  "mpMinLatF"      : 0.;5.
  "mpMaxLatF"      : 50.;45.
 end create
 set_size_mp(map,xf(1),yf(1),wi(1),he(1))

 vec1 = create "vector" vectorPlotClass wks
  "vcVectorFieldData" : wind
  "vcRefMagnitudeF" : 10.
  "vcRefLengthF"    : 0.03
  "vcRefAnnoOn"     : False
 end create
 set_size(vec1,xf(1),yf(1),wi(1),he(1))
 axis_str(vec1,"MERRA","","")

 contour1 = cnshadeplot(wks,lon,lat,um,cnlev,cnfil,"u (m s:S:-1:N:)")
 setvalues contour1
  "cnMissingValFillColor" : colind
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.06
 end setvalues
 set_size(contour1,xf(1),yf(1),wi(1),he(1))

 resp = True
 resp@gsLineThicknessF = 1.2
 dum1 = gsn_add_polyline(wks,contour1,(/90,90,130,130,90/),(/5,15,15,5,5/),resp)
 dum2 = gsn_add_polyline(wks,contour1,(/110,110,140,140,110/),(/22.5,32.5,32.5,22.5,22.5/),resp)

 overlay(map,vec1)
 NhlAddOverlay(map,contour1,map)

 draw(map)

 NhlRemoveOverlay(map,contour1,False)
 NhlRemoveOverlay(map,vec1    ,False)

 set_size_mp(map,xf(2),yf(1),wi(1),he(1))

 vec2 = create "vector" vectorPlotClass wks
  "vcVectorFieldData" : wind1
  "vcRefMagnitudeF" : 10.
  "vcRefLengthF"    : 0.03
  "vcRefAnnoOn"     : False
 end create
 set_size(vec2,xf(2),yf(1),wi(1),he(1))

 axis_str(vec2,"MERRA "+yearf,"","")

 contour2 = cnshadeplot(wks,lon,lat,um1,cnlev,cnfil,"u (m s:S:-1:N:)")
 setvalues contour2
  "cnMissingValFillColor" : colind
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.06
 end setvalues
 set_size(contour2,xf(2),yf(1),wi(1),he(1))

 dum1 = gsn_add_polyline(wks,contour2,(/90,90,130,130,90/),(/5,15,15,5,5/),resp)
 dum2 = gsn_add_polyline(wks,contour2,(/110,110,140,140,110/),(/22.5,32.5,32.5,22.5,22.5/),resp)

 overlay(map,vec2)
 NhlAddOverlay(map,contour2,map)

 draw(map)

 frame(wks)

 status_exit(10)

end

