begin

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;--------------------------------------------

wks = new(2,graphic)
wks(0) = create "1sttry" psWorkstationClass defaultapp
      "wkPSFileName" : "sa1a1.ps"
;      "wkColorMap"   : "ViBlGrWhYeOrRe"
end create
wks(1) = create "1sttry" psWorkstationClass defaultapp
      "wkPSFileName" : "sa1a2.ps"
;      "wkColorMap"   : "ViBlGrWhYeOrRe"
end create

;---------------------------------------

 file1 = "res/sa1a_0.nc"
 f1   = addfile(file1,"r")
 t    = f1->t
 x01  = f1->X1
 x02  = f1->X2
 fr0  = f1->freq
 ps01 = f1->PS1
 ps02 = f1->PS2

 file2 = "res/sa1a_hlag.nc"
 f2   = addfile(file2,"r")
 frh  = f2->freq
 psh1 = f2->PS1
 psh2 = f2->PS2

 file3 = "res/sa1a_llag.nc"
 f3   = addfile(file3,"r")
 frl  = f3->freq
 psl1 = f3->PS1
 psl2 = f3->PS2

 file4 = "res/sa1a_wn.nc"
 f4   = addfile(file4,"r")
 xw1  = f4->X1
 xw2  = f4->X2
 psw1 = f4->PS1
 psw2 = f4->PS2


;------------------------------------------

 data0 = new(2,graphic)
 dataw = new(2,graphic)
 data0(0)= create "00" coordArraysClass defaultapp
   "caXArray" : t
   "caYArray" : x01
 end create
 dataw(0)= create "00" coordArraysClass defaultapp
   "caXArray" : t
   "caYArray" : xw1
 end create
 data0(1)= create "00" coordArraysClass defaultapp
   "caXArray" : t
   "caYArray" : x02
 end create
 dataw(1)= create "00" coordArraysClass defaultapp
   "caXArray" : t
   "caYArray" : xw2
 end create

 do i=0, 1
 plot0 = create "plot0" xyPlotClass wks(i)
  "xyCoordData" : dataw(i)

  "trYMaxF" : 10
  "trYMinF" : -10
  "trXMaxF" : 1000
  "trXMinF" : 0

  "tmXBMode" : "MANUAL"
  "tmXBTickSpacingF" : 200.
  "tmXBMinorPerMajor" : 3
  "tmXBLabelFont" : 21
  "tmXBLabelFontHeightF" : 0.011
  "tmXBMajorLengthF" :  -0.007
  "tmXBMinorLengthF" :  -0.004

  "tmYLMinorLengthF" : -0.004
  "tmYLMajorLengthF" : -0.008
  "tmYLMode" : "MANUAL"
  "tmYLTickSpacingF" : 5
  "tmYLMinorPerMajor" : 4
  "tmYLLabelFont" : 21
  "tmYLLabelFontHeightF" : 0.011
  "tmYLFormat" : "@*+^sg"

  "tiMainFont" : 22
  "tiMainFontHeightF" : 0.020
  "tiMainString" : "Idealized Exp. "+(i+1)
  "tiXAxisFont" : 21
  "tiXAxisFontHeightF" : 0.017
  "tiXAxisString" : "t"
  "tiYAxisFont" : 21
  "tiYAxisFontHeightF" : 0.017
  "tiYAxisString" : "X"

  "vpYF" : 0.94
  "vpXF" : 0.10
  "vpHeightF" : 0.20
  "vpWidthF" : 0.80

 end create
 getvalues plot0
  "xyCoordDataSpec" : spec
 end getvalues
 setvalues spec
  "xyLineColor" : "grey"
 end setvalues
 draw(plot0)

 setvalues plot0
  "xyCoordData" : data0(i)
 end setvalues
 draw(plot0)

 delete(plot0)
 end do

;-----------------------------------------------
 data1 = new((/2,4/),graphic)
 data1(0,0)= create "00" coordArraysClass defaultapp
   "caXArray" : frh
   "caYArray" : psh1
 end create
 data1(0,1)= create "00" coordArraysClass defaultapp
   "caXArray" : fr0
   "caYArray" : ps01
 end create
 data1(0,2)= create "00" coordArraysClass defaultapp
   "caXArray" : frl
   "caYArray" : psl1
 end create
 data1(0,3)= create "00" coordArraysClass defaultapp
   "caXArray" : fr0
   "caYArray" : psw1
 end create
 data1(1,0)= create "00" coordArraysClass defaultapp
   "caXArray" : frh
   "caYArray" : psh2
 end create
 data1(1,1)= create "00" coordArraysClass defaultapp
   "caXArray" : fr0
   "caYArray" : ps02
 end create
 data1(1,2)= create "00" coordArraysClass defaultapp
   "caXArray" : frl
   "caYArray" : psl2
 end create
 data1(1,3)= create "00" coordArraysClass defaultapp
   "caXArray" : fr0
   "caYArray" : psw2
 end create

 xf = (/0.1,0.55,0.1,0.55/)
 yf = (/0.6,0.6,0.26,0.26/)
 str = (/"maxlag = nt-1","maxlag = nt/2","maxlag = nt/10","with white noise"/)
 xst = (/"","","frequency","frequency"/)
 yst = (/"PSD","","PSD",""/)
 ymax = (/800,120/)
 dy = (/200,40/)
 do i=0, 1
 do j=0, 3
 plot1 = create "plot1" xyPlotClass wks(i)
  "xyCoordData" : data1(i,j)

  "trYMaxF" : ymax(i)
  "trYMinF" : 0
  "trXMaxF" : 0.1
  "trXMinF" : 0

  "tmXBMode" : "MANUAL"
  "tmXBTickSpacingF" : 0.02
  "tmXBMinorPerMajor" : 3
  "tmXBLabelFont" : 21
  "tmXBLabelFontHeightF" : 0.011
  "tmXBMajorLengthF" :  -0.007
  "tmXBMinorLengthF" :  -0.004

  "tmYLMinorLengthF" : -0.004
  "tmYLMajorLengthF" : -0.008
  "tmYLMode" : "MANUAL"
  "tmYLTickSpacingF" : dy(i)
  "tmYLMinorPerMajor" : 3
  "tmYLLabelFont" : 21
  "tmYLLabelFontHeightF" : 0.011
  "tmYLFormat" : "@*+^sg"

  "tiMainFont" : 21
  "tiMainFontHeightF" : 0.020
  "tiMainString" : str(j)
  "tiXAxisFont"  : 21
  "tiXAxisFontHeightF" : 0.017
  "tiXAxisString" : xst(j)
  "tiYAxisFont" : 21
  "tiYAxisFontHeightF" : 0.017
  "tiYAxisString" : yst(j)

  "vpYF" : yf(j)
  "vpXF" : xf(j)
  "vpHeightF" : 0.20
  "vpWidthF" : 0.35

 end create
 draw(plot1)
 delete(plot1)
 end do
 end do

;-------------------------------------------------

frame(wks(0))
frame(wks(1))

end



