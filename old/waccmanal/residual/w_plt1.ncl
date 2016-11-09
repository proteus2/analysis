load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"

begin

 setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize" : 55000000
 end setvalues


do i11 = 0, 2
do i33 = 1, 12


 ;  0-CTL  ;  1-SGWDC  ;  2-LSGWDC
 idexp = i11;0
 cexp = (/"CTL","SGWDC","LSGWDC"/)

 ;  0-EPD  ;  1-BGWD  ;  2-GWDO  ;  3-GWDC
 cvar = (/"EPD","BGWD","GWDO","GWDC"/)
 idash = (/0,1,2,3/)

 mon = i33;1
 cmon = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

 ;  7- 13 km  ;  8- 15 km  ;  9- 17 km
 iheight = 8
 cheight = (/"13","15","17"/)

 ymax = (/0.4,0.4/)

 mainstr = (/cheight(iheight-7)+" km, "+cmon(mon-1)+", "+cexp(idexp),""/)

 ofname = "w"+cheight(iheight-7)+"-"+sprinti("%2.2i",mon)+"-"+cexp(idexp)

 xstr = (/"lat (deg)"/)
 ystr = (/"W* (mm s:S:-1:N:)"/)

;--------------------------------------- read base data

 fdir = "res/"

 if (idexp .eq. 0) then
   fnam1 = "control1.nc"
 end if
 if (idexp .eq. 1) then
   fnam1 = "sgwdclm2.nc"
 end if
 if (idexp .eq. 2) then
   fnam1 = "lsgwdclm.nc"
 end if

 missv = 1.e20

 f1 = addfile(fdir+fnam1,"r")

 lat= f1->lat
 z  = f1->z

 ny = dimsizes(lat)

 var1v = new((/4,ny/),float)
 var2v = new((/4,ny/),float)
 var1w = new((/4,ny/),float)
 var2w = new((/4,ny/),float)

 var1v(0,:) = f1->vsepd(mon-1,iheight,:)
 var1w(0,:) = f1->wsepd(mon-1,iheight,:)*1.e3
 var1v(1,:) = f1->vsgwd(mon-1,iheight,:)
 var1w(1,:) = f1->wsgwd(mon-1,iheight,:)*1.e3
 var1v(2,:) = f1->vsgwo(mon-1,iheight,:)
 var1w(2,:) = f1->wsgwo(mon-1,iheight,:)*1.e3
 if (idexp .ne. 0) then
   var1v(3,:) = f1->vsgwc(mon-1,iheight,:)
   var1w(3,:) = f1->wsgwc(mon-1,iheight,:)*1.e3
 else
   var1v(3,:) = f1->vsgwd(mon-1,iheight,:)*0.
   var1w(3,:) = f1->wsgwd(mon-1,iheight,:)*0.
 end if


 vardraw = new((/2,4,ny/),float,missv)
 vardraw(0,:,:) = (/var1w/)
; vardraw(1,:,:) = (/var1w/)

print(min(vardraw(0,:,:)))

;---------------------------------------------------

 wks = create "ps" psWorkstationClass defaultapp
     "wkPSFileName"   : ofname+".eps"
     "wkPSFormat"     : "eps"
     "wkPSResolution" : 600
     "wkColorMap"        : "posneg_2"
     "wkBackgroundColor" : (/1.,1.,1./)
 end create
 colind = NhlNewColor(wks,0.5,0.5,0.5)

;--------------------------------------- cn level

 ncn = 41
 cnlev = new(ncn,float)

 cnfil = new(ncn+1,integer,0)
 cnfil = (/ 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, \
               7, 7, 7, 7, 9, 9, 9, 0, 0, 0, \
               0, 0, 0,14,14,14,16,16,16,16, \
              16,16,16,16,16,16,16,16,16,16,16/)

;--------------------------------------- field

  pair = new((/2,4/),graphic)

  do n=0, 1
  do j=0, 3
  pair(n,j) = create "w" coordArraysClass defaultapp
    "caXArray" : lat
    "caYArray" : vardraw(n,j,:)
    "caYMissingV" : missv
  end create
  end do
  end do

  pair0 = create "w" coordArraysClass defaultapp
    "caXArray" : lat
    "caYArray" : vardraw(0,0,:)*0.
    "caYMissingV" : missv
  end create

;--------------------------------------- contour
; wprt1  
;---------------------------------------

  if (idexp .eq. 0) then
    nline = 3
  else
    nline = 4
  end if

  do n=0, 0
  do j=0, nline-1

  plot =  create "lplot" xyPlotClass wks
   "xyCoordData" : pair(n,j)

   "trXMaxF" : 60
   "trXMinF" : -60
   "trYMaxF" : ymax(n)
   "trYMinF" : -ymax(n)

   "tmXBMode" : "Manual"
   "tmXBTickSpacingF"  : 30
   "tmXBMinorPerMajor" : 2
   "tmXBMajorLengthF"  : -0.008
   "tmXBMinorLengthF"  : -0.005
   "tmXBLabelFont"        : 21
   "tmXBLabelFontHeightF" : 0.014

   "tmYLMode" : "Manual"
   "tmYLTickSpacingF"  : 0.1
   "tmYLMinorPerMajor" : 4
   "tmYLMajorLengthF"  : -0.008
   "tmYLMinorLengthF"  : -0.005
   "tmYLLabelFont"        : 21
   "tmYLLabelFontHeightF" : 0.014

   "tiMainString"       : mainstr(n)
   "tiMainFont"         : 22
   "tiMainFontHeightF"  : 0.018

   "tiXAxisString"      : xstr(n)
   "tiXAxisFont"        : 21
   "tiXAxisFontHeightF" : 0.017

   "tiYAxisString"      : ystr(n)
   "tiYAxisFont"        : 21
   "tiYAxisFontHeightF" : 0.017

   "vpXF" : 0.2
   "vpYF" : 0.9 - n*0.25
   "vpHeightF" : 0.5
   "vpWidthF"  : 0.5
  end create

  if (j .eq. 0) then
    setvalues plot
     "xyCoordData" : pair0
    end setvalues
    getvalues plot
     "xyCoordDataSpec" : spec
    end getvalues
    setvalues spec
     "xyLineColor" : colind
    end setvalues
    draw(plot)
    setvalues plot
     "xyCoordData" : pair(n,j)
    end setvalues
  end if

  getvalues plot
   "xyCoordDataSpec" : spec
  end getvalues
  setvalues spec
   "xyMonoDashPattern" : True
   "xyDashPattern"     : idash(j)
   "xyMonoLineThickness" : True
   "xyLineThicknessF"    : 1.5
  end setvalues
  draw(plot)


  end do
  end do


  lg1 = create "lg1" legendClass wks
    "lgAutoManage" : False
    "lgItemCount"        : nline
    "lgLabelStrings"     : "   "+(/"EPD","BGWD","GWDO","GWDC"/)

    "lgMonoDashIndex"      : False
    "lgMonoLineThickness"  : True
    "lgMonoLineColor"      : True
    "lgLineThicknessF"     : 1.5
    "lgDashIndexes"        : idash
    "lgLabelFont"        : 21
    "lgLabelFontHeightF" : 0.012
    "lgPerimOn"          : False

    "vpXF"        : 0.21
    "vpYF"        : 0.90
    "vpWidthF"    : 0.10
    "vpHeightF"   : 0.08
  end create
  draw(lg1)



   frame(wks)


end do
end do


 end

