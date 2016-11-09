begin

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;--------------------------------------------

wks = new(2,graphic)
wks(0) = create "1sttry" psWorkstationClass defaultapp
      "wkPSFileName" : "sa1b2co2.ps"
;      "wkColorMap"   : "ViBlGrWhYeOrRe"
end create
wks(1) = create "1sttry" psWorkstationClass defaultapp
      "wkPSFileName" : "sa1b2pic.ps"
;      "wkColorMap"   : "ViBlGrWhYeOrRe"
end create

;---------------------------------------

 file1 = "res/sa1b2co2.nc"
 f1   = addfile(file1,"r")
 t1   = f1->time
 x11  = f1->Nino3
 x12  = f1->SSTmean
 x13  = f1->Tamean
 fr1  = f1->freq
 p11  = f1->PS_Nino
 p12  = f1->PS_SST
 p13  = f1->PS_Ta

 file2 = "res/sa1b2pic.nc"
 f2   = addfile(file2,"r")
 t2   = f2->time
 x21  = f2->Nino3
 x22  = f2->SSTmean
 x23  = f2->Tamean
 fr2  = f2->freq 
 p21  = f2->PS_Nino
 p22  = f2->PS_SST
 p23  = f2->PS_Ta

;------------------------------------------

 datax = new((/2,3/),graphic)
 datap = new((/2,3/),graphic)
 datax(0,0)= create "00" coordArraysClass defaultapp
   "caXArray" : t1
   "caYArray" : x11
 end create
 datax(0,1)= create "00" coordArraysClass defaultapp
   "caXArray" : t1
   "caYArray" : x12
 end create
 datax(0,2)= create "00" coordArraysClass defaultapp
   "caXArray" : t1
   "caYArray" : x13
 end create
 datax(1,0)= create "00" coordArraysClass defaultapp
   "caXArray" : t2 
   "caYArray" : x21 
 end create 
 datax(1,1)= create "00" coordArraysClass defaultapp
   "caXArray" : t2
   "caYArray" : x22
 end create 
 datax(1,2)= create "00" coordArraysClass defaultapp
   "caXArray" : t2
   "caYArray" : x23
 end create
 datap(0,0)= create "00" coordArraysClass defaultapp
   "caXArray" : fr1 
   "caYArray" : p11 
 end create 
 datap(0,1)= create "00" coordArraysClass defaultapp
   "caXArray" : fr1
   "caYArray" : p12
 end create 
 datap(0,2)= create "00" coordArraysClass defaultapp
   "caXArray" : fr1
   "caYArray" : p13
 end create 
 datap(1,0)= create "00" coordArraysClass defaultapp
   "caXArray" : fr2
   "caYArray" : p21
 end create  
 datap(1,1)= create "00" coordArraysClass defaultapp
   "caXArray" : fr2 
   "caYArray" : p22
 end create  
 datap(1,2)= create "00" coordArraysClass defaultapp
   "caXArray" : fr2 
   "caYArray" : p23
 end create

;-----------------------------------------------
 yf = (/0.95,0.63,0.31/)
 str = (/"Nino3.4","SST mean","Ta mean"/)
 xst = (/"","","frequency","frequency"/)
 yst = (/"PSD","","PSD",""/)
 ymax = (/800,120/)
 xmin = (/2051,1871/)
 dy = (/200,40/)
;-------------------------------------------------

 do i=0, 1
 do j=0, 2
 plot0 = create "plot0" xyPlotClass wks(i)
  "xyCoordData" : datax(i,j)

;  "trYMaxF" : 10
;  "trYMinF" : -10
;  "trXMaxF" : 1000
;  "trXMinF" : xmin(i)

  "tmXBMode" : "automatic"
  "tmXBTickSpacingF" : 200.
  "tmXBMinorPerMajor" : 3
  "tmXBLabelFont" : 21
  "tmXBLabelFontHeightF" : 0.011
  "tmXBMajorLengthF" :  -0.007
  "tmXBMinorLengthF" :  -0.004

  "tmYLMinorLengthF" : -0.004
  "tmYLMajorLengthF" : -0.008
  "tmYLMode" : "automatic"
  "tmYLTickSpacingF" : 5
  "tmYLMinorPerMajor" : 4
  "tmYLLabelFont" : 21
  "tmYLLabelFontHeightF" : 0.011
  "tmYLFormat" : "@*+^sg"

  "tiMainFont" : 22
  "tiMainFontHeightF" : 0.020
  "tiMainString" : str(j)
  "tiXAxisFont" : 21
  "tiXAxisFontHeightF" : 0.017
  "tiXAxisString" : "t"
  "tiYAxisFont" : 21
  "tiYAxisFontHeightF" : 0.017
  "tiYAxisString" : ""

  "vpYF" : yf(j)
  "vpXF" : 0.05
  "vpHeightF" : 0.20
  "vpWidthF" : 0.40
 end create
 draw(plot0)

 setvalues plot0
  "xyCoordData" : datap(i,j)
  "vpXF" : 0.58
  "tiXAxisString" : "frequency (cycle/mon)"
  "tiYAxisString" : "PS (K:S:2:N: / rad mon:S:-1:N:)"
 end setvalues
 draw(plot0)

 delete(plot0)
 end do
 end do

frame(wks(0))
frame(wks(1))

end



