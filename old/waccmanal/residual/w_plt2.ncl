load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"

begin

 setvalues NhlGetWorkspaceObjectId()
  "wsMaximumSize" : 55000000
 end setvalues


do i22 = 0, 3
do i33 = 1, 12


 idash = (/1,0/)

 ;  0-EPD  ;  1-BGWD  ;  2-GWDO  ;  3-GWDC
 idvar = i22
 cvar = (/"EPD","BGWD","GWDO","GWDC"/)

 mon = i33;1
 cmon = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

 ;  7- 13 km  ;  8- 15 km  ;  9- 17 km
 iheight = 8
 cheight = (/"13","15","17"/)

 ymax = (/0.3,0.01,0.03,0.012/)

 mainstr = (/cheight(iheight-7)+" km, "+cvar(idvar)+", "+cmon(mon-1)/)

 ofname = "w"+cheight(iheight-7)+"_"+cvar(idvar)+"-"+sprinti("%2.2i",mon)

 xstr = (/"lat (deg)"/)
 ystr = (/"W* (mm s:S:-1:N:)"/)

;--------------------------------------- read base data

 fdir = "res/"

 fnam1 = "control1.nc"
; fnam2 = "sgwdclm2.nc"
 fnam2 = "lsgwdclm.nc"

 missv = 1.e20

 f1 = addfile(fdir+fnam1,"r")
 f2 = addfile(fdir+fnam2,"r")

 lat= f1->lat
 z  = f1->z

 ny = dimsizes(lat)

 var1v = new((/2,ny/),float)
 var1w = new((/2,ny/),float)

 if (idvar .eq. 0) then
   var1v(0,:) = f1->vsepd(mon-1,iheight,:)
   var1w(0,:) = f1->wsepd(mon-1,iheight,:)*1.e3
   var1v(1,:) = f2->vsepd(mon-1,iheight,:)
   var1w(1,:) = f2->wsepd(mon-1,iheight,:)*1.e3
 end if
 if (idvar .eq. 1) then
   var1v(0,:) = f1->vsgwd(mon-1,iheight,:)
   var1w(0,:) = f1->wsgwd(mon-1,iheight,:)*1.e3
   var1v(1,:) = f2->vsgwd(mon-1,iheight,:)
   var1w(1,:) = f2->wsgwd(mon-1,iheight,:)*1.e3
 end if
 if (idvar .eq. 2) then
   var1v(0,:) = f1->vsgwo(mon-1,iheight,:)
   var1w(0,:) = f1->wsgwo(mon-1,iheight,:)*1.e3
   var1v(1,:) = f2->vsgwo(mon-1,iheight,:)
   var1w(1,:) = f2->wsgwo(mon-1,iheight,:)*1.e3
 end if
 if (idvar .eq. 3) then
   var1v(1,:) = f2->vsgwc(mon-1,iheight,:)
   var1w(1,:) = f2->wsgwc(mon-1,iheight,:)*1.e3
   var1v(0,:) = var1v(1,:)*0.
   var1w(0,:) = var1w(1,:)*0.
 end if


 vardraw = new((/2,2,ny/),float,missv)
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

  pair = new((/2,2/),graphic)

  do n=0, 0
  do j=0, 1
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

  if (idvar .eq. 3) then
    nline = 1
  else
    nline = 2
  end if

  do n=0, 0
  do j=2-nline, 1

  plot =  create "lplot" xyPlotClass wks
   "xyCoordData" : pair(n,j)

   "trXMaxF" : 60
   "trXMinF" : -60
   "trYMaxF" : ymax(idvar)
   "trYMinF" : -ymax(idvar)

   "tmXBMode" : "Manual"
   "tmXBTickSpacingF"  : 30
   "tmXBMinorPerMajor" : 2
   "tmXBMajorLengthF"  : -0.008
   "tmXBMinorLengthF"  : -0.005
   "tmXBLabelFont"        : 21
   "tmXBLabelFontHeightF" : 0.014

   "tmYLMode" : "Automatic" ;"Manual"
;   "tmYLTickSpacingF"  : 0.25
;   "tmYLMinorPerMajor" : 4
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
   "vpHeightF" : 0.2
   "vpWidthF"  : 0.4
  end create

  if (j .eq. 2-nline) then
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


  lbstr = (/"CTL","LSGWDC"/)

  lg1 = create "lg1" legendClass wks
    "lgAutoManage" : False
    "lgItemCount"        : nline
    "lgLabelStrings"     : "   "+lbstr(2-nline:)

    "lgMonoDashIndex"      : False
    "lgMonoLineThickness"  : True
    "lgMonoLineColor"      : True
    "lgLineThicknessF"     : 1.5
    "lgDashIndexes"        : idash(2-nline:)
    "lgLabelFont"        : 21
    "lgLabelFontHeightF" : 0.012
    "lgPerimOn"          : False

    "vpXF"        : 0.21
    "vpYF"        : 0.90
    "vpWidthF"    : 0.10
    "vpHeightF"   : 0.05
  end create
  draw(lg1)



   frame(wks)


end do
end do


 end

