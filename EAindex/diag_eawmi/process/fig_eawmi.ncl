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

 eawmi = new(nt,float)
 eawmi = 0.

 mm = (/12,1,2/)

 n = 0
 do iy=year1, year2
   iy2 = iy
   do im=0, 2
     if (mm(im) .eq. 1) then
       iy2 = iy + 1
     end if
     file_u = "../dat/"+iy2+"/merra.u_m.assm.pl"+300+"."+ \
              iy2+sprinti("%0.2i",mm(im))+".nc"
     file_v = "../dat/"+iy2+"/merra.v_m.assm.pl"+300+"."+ \
              iy2+sprinti("%0.2i",mm(im))+".nc"
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
     ; calculate EAWMI
     u = fu->u(0,0,:,:)
     v = fv->v(0,0,:,:)
     eawmi(n) = eawmi(n) + ( avg(u({27.5:37.5},{110:170})) - \
                             avg(u({50.0:60.0},{ 80:140})) )/3.
     um = um + u/(tofloat(nt*3))
     vm = vm + v/(tofloat(nt*3))
   end do
   n = n + 1
 end do

 eawmi = eawmi - avg(eawmi)

 iy2 = yearf
 do im=0, 2
   if (mm(im) .eq. 1) then
     iy2 = yearf + 1
   end if
   file_u = "../dat/"+iy2+"/merra.u_m.assm.pl"+300+"."+ \
            iy2+sprinti("%0.2i",mm(im))+".nc"
   file_v = "../dat/"+iy2+"/merra.v_m.assm.pl"+300+"."+ \
            iy2+sprinti("%0.2i",mm(im))+".nc"
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

 wnd = (/um(::8,::8),vm(::8,::8)/)
 wind = create "f1" vectorFieldClass defaultapp
  "vfDataArray" : wnd
  "vfXArray" : lon(::8)
  "vfYArray" : lat(::8)
  "vfMissingUValueV" : 1.e15
  "vfMissingVValueV" : 1.e15
 end create

 wnd1 = (/um1(::8,::8),vm1(::8,::8)/)
 wind1 = create "f1" vectorFieldClass defaultapp
  "vfDataArray" : wnd1
  "vfXArray" : lon(::8)
  "vfYArray" : lat(::8)
  "vfMissingUValueV" : 1.e15
  "vfMissingVValueV" : 1.e15
 end create

;--------------------------------------- cn level
 cnlev = fspan(4,72,18)
 cnfil = (/0,11,10,9,6,5,4,3,2,26,25,24,23,22,21,20,19,18,17/)

 cnleva = fspan(-10,10,21)
 cnfila = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

;--------------------------------------- drawing
 wi = (/nt*2,36/)
 he = (/20,30/)
 xf = (/10,10,50,90/)
 yf = (/10,42/)

 outfile = getenv("out1")
 wks = cr_wks(outfile,600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.6,0.6,0.6)
 plot = xylineplot(wks,time,eawmi,"D",0,"T",1.2,"C",1)
 xyaddmarker(plot,16,"S",0.5,"T",1.0,"C",1)
 xyaddline(plot,time,0.,"D",0,"T",0.6,"C",colind)
 set_size(plot,xf(0),yf(0),wi(0),he(0))
 axis_range(plot,year1-0.5,year2+0.5,-999,999,False,False)
 setvalues plot
  "tmXBMode"   : "Explicit"
  "tmXBValues" : ispan(1975,2100,5)
  "tmXBLabels" : sprinti("%2.2i",mod(ispan(1975,2100,5),100))+"/"+ \
                 sprinti("%2.2i",mod(ispan(1976,2101,5),100))
  "tmXBMinorValues" : ispan(1975,2100,1)
 end setvalues
 axis_str(plot,"East Asian Winter Monsoon Index","year","(m s:S:-1:N:)")
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
  "mpMinLonF"      : 85.
  "mpMaxLonF"      : 165.
  "mpMinLatF"      : 10.
  "mpMaxLatF"      : 60.
 end create
 set_size_mp(map,xf(1),yf(1),wi(1),he(1))

 vec1 = create "vector" vectorPlotClass wks
  "vcVectorFieldData" : wind
  "vcRefMagnitudeF" : 30.
  "vcRefLengthF"    : 0.02
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
 dum1 = gsn_add_polyline(wks,contour1,(/80,80,140,140,80/),(/50,60,60,50,50/),resp)
 dum2 = gsn_add_polyline(wks,contour1,(/110,110,170,170,110/),(/27.5,37.5,37.5,27.5,27.5/),resp)

 overlay(map,vec1)
 NhlAddOverlay(map,contour1,map)

 draw(map)

 NhlRemoveOverlay(map,contour1,False)
 NhlRemoveOverlay(map,vec1    ,False)

 set_size_mp(map,xf(2),yf(1),wi(1),he(1))

 vec2 = create "vector" vectorPlotClass wks
  "vcVectorFieldData" : wind1
  "vcRefMagnitudeF" : 5.
  "vcRefLengthF"    : 0.02
  "vcRefAnnoOn"     : False
 end create
 set_size(vec2,xf(2),yf(1),wi(1),he(1))
 axis_str(vec2,"MERRA "+yearf+"/"+(yearf+1),"","")

 contour2 = cnshadeplot(wks,lon,lat,um1,cnleva,cnfila,"u (m s:S:-1:N:)")
 setvalues contour2
  "cnMissingValFillColor" : colind
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.06
 end setvalues
 set_size(contour2,xf(2),yf(1),wi(1),he(1))

 dum1 = gsn_add_polyline(wks,contour2,(/80,80,140,140,80/),(/50,60,60,50,50/),resp)
 dum2 = gsn_add_polyline(wks,contour2,(/110,110,170,170,110/),(/27.5,37.5,37.5,27.5,27.5/),resp)

 overlay(map,vec2)
 NhlAddOverlay(map,contour2,map)

 draw(map)

 frame(wks)

 status_exit(10)

end

