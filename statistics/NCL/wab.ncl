begin

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;--------------------------------------------

wks = create "1sttry" psWorkstationClass defaultapp
      "wkPSFileName" : "wab.ps"
      "wkColorMap"   : "hotres"
end create

;---------------------------------------

 file1 = "res/wabpic.nc"
 f1    = addfile(file1,"r")
 ps1   = f1->PS
 xout  = f1->Xout
 peri  = f1->FScale

 file2 = "res/wabpicsig.nc"
 f2    = addfile(file2,"r")
 sig   = f2->Sig
 ps_theor = f2->PS_theor
 sig_glob = f2->Sig_glob

 file3 = "res/sa1b2pic.nc"
 f3    = addfile(file3,"r")
 fps1  = f3->PS_Nino*(2.*3.141593)/12.
 ffreq = f3->freq


 psd = new((/2,dimsizes(peri),dimsizes(xout)-1/),float)
 x   = new(dimsizes(xout)-1,float)
 do i=0, dimsizes(xout)-2
   x(i) = xout(i)
   psd(0,:,i) = ps1(i,:)
   psd(1,:,i) = psd(0,:,i)/sig(:)
 end do

 ps_glob = new((/3,dimsizes(peri)/),float)
 ps_glob(0,:) = ps1(dimsizes(xout)-1,:)
 ps_glob(1,:) = ps_theor(:)
 ps_glob(2,:) = sig_glob(:)

 fperi = ffreq
 fperi(0) = 999.
 fperi(1:) = 1./ffreq(1:)/12.

;---------------------------------- cn level

 cnlev = new(6,float)
 cnfil = new(7,integer)
 do i=0,5
  cnlev(i) = 2.^(i-5 -1)
 end do
 cnfil = (/1,4,9,16,36,64,90/)*2


;------------------------------------------

 field = new(2,graphic)
 do i=0, 1
 field(i) = create "w" scalarFieldClass defaultapp
   "sfDataArray" : psd(i,:,:)
   "sfYArray" : peri
   "sfXArray" : x
 end create
 end do

 yf = (/0.95,0.95,0.31/)
 do i=0, 1
  cont1 = create "contour" contourPlotClass wks
   "cnScalarFieldData" : field(i)
   "cnFillOn"        :  True           ; color fill
   "cnLevelSelectionMode" : "ExplicitLevels"
   "cnLevels" : cnlev
   "cnFillColors" : cnfil
   "cnLinesOn" : False
   "cnInfoLabelOn" : False
   "cnLineLabelPlacementMode" : "Constant"
   "cnLineLabelFont" : 21
   "cnLineLabelFontHeightF" : 0.011
   "cnLineLabelsOn" : False

   "tmXBMode" : "Manual"
   "tmXBTickSpacingF" : 20
   "tmXBMinorPerMajor" : 3
   "tmXBMajorLengthF" : -0.007
   "tmXBMinorLengthF" : -0.004
   "tmXBLabelFont" : 21
   "tmXBLabelFontHeightF" : 0.011
   "tmXBFormat" : "@*+^sg"

   "tmYLMode" : "Explicit"
   "tmYLValues" : 2.^(/0,1,2,3,4,5,6,7,8,9/)
   "tmYLLabels" : 2.^(/0,1,2,3,4,5,6,7,8,9/)
   "tmYLMinorLengthF" : -0.004
   "tmYLMajorLengthF" : -0.007
   "tmYLLabelFont" : 21
   "tmYLLabelFontHeightF" : 0.011
   "tmYLFormat" : "@*+^sg"

   "vpXF" : 0.1
   "vpYF" : yf(i)
   "vpWidthF" : 0.65
   "vpHeightF" : 0.30

   "trXMinF" : 0
   "trXMaxF" : max(x)
   "trYMinF" : 1
   "trYMaxF" : 32
   "trYReverse" : "True"

   "tiMainString" : "Wavelet power spectrum of Nino3.4 index"
   "tiMainFont"  : 22
   "tiMainFontHeightF" : 0.017

   "tiXAxisString" : "time (yr)"
   "tiXAxisFont"  : 21
   "tiXAxisFontHeightF" : 0.015

   "tiYAxisString" : "period (yr)"
   "tiYAxisFont"  : 21
   "tiYAxisFontHeightF" : 0.015
  end create
  if (i .eq. 1) then
    setvalues cont1
      "cnFillOn"     : False
      "cnLevels"     : 1.
      "cnLinesOn"    : True
    end setvalues
  end if
  draw(cont1)
  delete(cont1)
  end do


  data = new(4,graphic)
  do i=0, 2
  data(i) = create "00" coordArraysClass defaultapp
    "caXArray" : peri
    "caYArray" : ps_glob(i,:)
  end create
  end do
  data(3) = create "00" coordArraysClass defaultapp
    "caXArray" : fperi
    "caYArray" : fps1
  end create

  lthc = (/2,1,1,1/)
  ldsh = (/0,1,1,0/)
  do i=0, 3
  plot0 = create "plot0" xyPlotClass wks
    "xyCoordData" : data(i)

    "trYMaxF" : 0.5
    "trYMinF" : 0
    "trXMaxF" : 32
    "trXMinF" : 1
    "trXLog"  : True
    "trXReverse" : True

    "tmXBMode" : "Explicit"
    "tmXBValues" : 2.^(/0,1,2,3,4,5,6,7,8,9/)
    "tmXBLabels" : 2.^(/0,1,2,3,4,5,6,7,8,9/)
    "tmXBLabelFont" : 21
    "tmXBLabelFontHeightF" : 0.011
    "tmXBMajorLengthF" :  -0.008
    "tmXBMinorLengthF" :  -0.004
    "tmXBFormat" : "@*+^sg"

    "tmYLMinorLengthF" : -0.004
    "tmYLMajorLengthF" : -0.008
    "tmYLMode" : "MANUAL"
    "tmYLTickSpacingF" : 0.05
    "tmYLMinorPerMajor" : 4
    "tmYLLabelFont" : 21
    "tmYLLabelFontHeightF" : 0.011
;    "tmYLFormat" : "@*+^sg"

    "tiMainFont" : 22
    "tiMainFontHeightF" : 0.020
    "tiMainString" : "Global avg."
    "tiXAxisFont" : 21
    "tiXAxisFontHeightF" : 0.017
    "tiXAxisString" : "period (yr)"
    "tiYAxisFont" : 21
    "tiYAxisFontHeightF" : 0.017
    "tiYAxisString" : "PSD"

    "vpYF" : 0.50
    "vpXF" : 0.1
    "vpHeightF" : 0.30
    "vpWidthF" : 0.5
  end create
  getvalues plot0
    "xyCoordDataSpec" : spec
  end getvalues
  setvalues spec
    "xyDashPattern"    : ldsh(i)
    "xyLineThicknessF" : lthc(i)
  end setvalues
  draw(plot0)
  delete(plot0)
  end do




   label = create "lab" labelBarClass wks
     "lbOrientation" : "vertical"
     "lbBoxCount" : 7
     "lbFillColors" : cnfil
     "lbMonoFillPattern" : True
     "lbLabelsOn" : True
     "lbLabelStrings" : log(cnlev)/log(2.)
     "lbAutoManage" : False
     "lbLabelFont" : 21
     "lbLabelFontHeightF" : 0.010
     "lbLabelPosition" : "Bottom";,Right"
     "lbLabelAlignment" : "InteriorEdges"
     "lbLabelAngleF" : 0.
     "lbPerimOn" : False
     "lbTitleString" : "log:B:2:N: (PSD)"
     "lbTitlePosition" : "Bottom";,"Top"
     "lbTitleDirection" : "Across"
     "lbTitleFont" : 22
     "lbTitleFontHeightF" : 0.010
     "vpXF" : 0.80
     "vpYF" : 0.95
     "vpWidthF"  : 0.06
     "vpHeightF" : 0.25
   end create
   draw(label)


;-----------------------------------------------

frame(wks)

end



