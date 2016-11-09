begin

wks = create "1sttry" psWorkstationClass defaultapp
      "wkPSFileName" :   "m2p.ps"
      "wkColorMap"   : "testcmap"
end create

;---------------------------------------

f1  = addfile("/usr/users/kyh/anal/work/sep60/data/m2.nc","r")
m2d = f1->m2
z   = f1->z

f2  = addfile("/usr/users/kyh/anal/work/sep60/data/ko60.nc","r")
psd = f2->psd

f3  = addfile("/usr/users/kyh/anal/work/sep60/data/co60.nc","r")
csd = f3->co_sd

o1 = 197
k1 = 18
o2 = 198
k2 = 20
mm1 = m2d(:,o1,k1)
mm2 = m2d(:,o2,k2)
;mm1 = sqrt(sqrt(mm1^2))
;mm2 = sqrt(sqrt(mm2^2))


pp1 = psd(:,o1,k1)
pp2 = psd(:,o2,k2) 

cc1 = csd(:,o1,k1)
cc2 = csd(:,o2,k2)


;----------------------------------------

data1 = create "T" coordArraysClass defaultapp
      "caXArray"  : mm1 * 10.^8
      "caYArray"  : z
end create

data2 = create "T" coordArraysClass defaultapp
      "caXArray"  : mm2 * 10.^8
      "caYArray"  : z
end create

line0 = new(202,"float")
line0(:) = 0.
data00 = create "T" coordArraysClass defaultapp
      "caXArray"  : line0
      "caYArray"  : z
end create

data1p = create "T" coordArraysClass defaultapp
      "caXArray"  : pp1 / 10.^6
      "caYArray"  : z
end create

data2p = create "T" coordArraysClass defaultapp
      "caXArray"  : pp2 / 10.^6
      "caYArray"  : z
end create

data1c = create "T" coordArraysClass defaultapp
      "caXArray"  : cc1 / 10.^6
      "caYArray"  : z 
end create 
 
data2c = create "T" coordArraysClass defaultapp
      "caXArray"  : cc2 / 10.^6
      "caYArray"  : z
end create

;----------------------------------------

plot1 =  create "lplot" xyPlotClass wks

  "xyCoordData" : data1

  "trYMaxF" : 60.
  "trYMinF" : 0.
  "trXMaxF" : 2.5
  "trXMinF" : -0.5

  "tmXBMode" : "MANUAL"
  "tmXBTickSpacingF" : 0.5
  "tmXBMinorPerMajor" : 4
  "tmXBLabelFont" : 21
  "tmXBLabelFontHeightF" : 0.012
  "tmXBMajorLengthF" : -0.007
  "tmXBMinorLengthF" : -0.004

  "tmYLMinorLengthF" : -0.004
  "tmYLMajorLengthF" : -0.007
  "tmYLMode" : "MANUAL"
  "tmYLTickSpacingF" : 10.
  "tmYLMinorPerMajor" : 4
  "tmYLLabelFont" : 21
  "tmYLLabelFontHeightF" : 0.012


  "tiMainFont" : 25
  "tiMainFontHeightF" : 0.017
    "tiMainString" : "L:B:x:N: = 38.89 km, T = 21.18 min, c:B:p:N: = 30.6 ms:S:-1:N:"
  "tiXAxisFont" : 21
  "tiXAxisFontHeightF" : 0.014
  "tiXAxisString" : "m:S:2:N: [10:S:-8:N: (cycle/m):S:2:N:]"
  "tiYAxisFont" : 21
  "tiYAxisFontHeightF" : 0.014
  "tiYAxisString" : "z (km)"

  "vpYF" : 0.95
  "vpXF" : 0.1
  "vpHeightF" : 0.35
  "vpWidthF" : 0.3

end create
draw(plot1)
setvalues plot1
  "xyCoordData" : data00
end setvalues
getvalues plot1
  "xyCoordDataSpec" : spec0
end getvalues
setvalues spec0
  "xyMonoDashPattern" : True
  "xyDashPattern" : 1
end setvalues
draw(plot1)

setvalues plot1
  "xyCoordData" : data2
  "trXMaxF" : 3.0
  "tiMainString" : "L:B:x:N: = 35 km, T = 20 min, c:B:p:N: = 29.2 ms:S:-1:N:"
  "vpYF" : 0.45
  "vpXF" : 0.1
end setvalues
draw(plot1)

setvalues plot1
  "xyCoordData" : data00
end setvalues
getvalues plot1
  "xyCoordDataSpec" : spec0
end getvalues
setvalues spec0
  "xyMonoDashPattern" : True
  "xyDashPattern" : 1
end setvalues
draw(plot1)

;-------------------------

plot1p =  create "lplot" xyPlotClass wks

  "xyCoordData" : data1p

  "trYMaxF" : 60.
  "trYMinF" : 0.
  "trXMaxF" : 3
  "trXMinF" : -3

  "tmXBMode" : "MANUAL"
  "tmXBTickSpacingF" : 1
  "tmXBMinorPerMajor" : 4
  "tmXBLabelFont" : 21
  "tmXBLabelFontHeightF" : 0.012
  "tmXBMajorLengthF" : -0.007
  "tmXBMinorLengthF" : -0.004
  "tmXBFormat" : "@*+^sg"

  "tmYLMinorLengthF" : -0.004
  "tmYLMajorLengthF" : -0.007
  "tmYLMode" : "MANUAL"
  "tmYLTickSpacingF" : 10.
  "tmYLMinorPerMajor" : 4
  "tmYLLabelFont" : 21
  "tmYLLabelFontHeightF" : 0.012


  "tiMainFont" : 21
  "tiMainFontHeightF" : 0.014
  "tiMainFont" : 21
  "tiMainFontHeightF" : 0.014
  "tiMainString" : "CoSD [10:S:6:N: (m:S:2:N:s:S:-2:N:/cycle:S:2:N:m:S:-1:N:s:S:-1:N:)]"
  "tiXAxisFont" : 21
  "tiXAxisFontHeightF" : 0.014
  "tiXAxisString" : "PSD [10:S:6:N: (m:S:2:N:s:S:-2:N:/cycle:S:2:N:m:S:-1:N:s:S:-1:N:)]"
  "tiYAxisFont" : 21
  "tiYAxisFontHeightF" : 0.014
  "tiYAxisString" : "z (km)"

  "vpYF" : 0.95
  "vpXF" : 0.6
  "vpHeightF" : 0.35
  "vpWidthF" : 0.3

end create
draw(plot1p)

setvalues plot1p
  "xyCoordData" : data00
end setvalues
draw(plot1p)

setvalues plot1p
  "xyCoordData" : data1c
end setvalues
getvalues plot1p
  "xyCoordDataSpec" : spec
end getvalues
setvalues spec
  "xyMonoDashPattern" : True
  "xyDashPattern" : 1
  "xyMonoLineColor" : True
  "xyLineColor" : 180
end setvalues
draw(plot1p)

setvalues plot1p
  "xyCoordData" : data2p
  "trXMaxF" : 3
  "trXMinF" : -3
  "tmXBTickSpacingF" : 1
  "tmXBMinorPerMajor" : 4
  "vpYF" : 0.45
  "vpXF" : 0.6
  "tiMainString" : "CoSD [10:S:6:N: (m:S:2:N:s:S:-2:N:/cycle:S:2:N:m:S:-1:N:s:S:-1:N:)]"
end setvalues
draw(plot1p)

setvalues plot1p
  "xyCoordData" : data00
end setvalues 
draw(plot1p)

setvalues plot1p
  "xyCoordData" : data2c
end setvalues
getvalues plot1p
  "xyCoordDataSpec" : spec
end getvalues
setvalues spec
  "xyMonoDashPattern" : True
  "xyDashPattern" : 1
  "xyMonoLineColor" : True
  "xyLineColor" : 180
end setvalues 
draw(plot1p)


;--------------------------------

frame(wks)


end
