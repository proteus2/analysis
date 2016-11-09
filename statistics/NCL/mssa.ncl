begin

 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;--------------------------------------------

wks = create "1sttry" psWorkstationClass defaultapp
      "wkPSFileName" : "mssa.ps"
;      "wkColorMap"   : "hotres"
end create

;---------------------------------------

 file1 = "res/mssa31_compo.nc"
 f1    = addfile(file1,"r")
 pc1   = f1->PC
 rc1   = f1->RC
 re1   = f1->recon
 inino = f1->data
 x     = f1->x
 t     = f1->t
 t2    = f1->t2
 com   = f1->com

 file2 = "res/mssa31_eigen.nc"
 f2    = addfile(file2,"r")
 eof   = f2->ei_vec
 val   = f2->ei_val_100
 lag   = f2->lag


  xf = (/0.1, 0.55,0.1,0.1,0.1/)
  yf = (/0.97,0.97,0.75,0.53,0.31/)
  wi = (/0.35,0.35,0.35,0.35,0.35/)
  he = (/0.15,0.15,0.15,0.15,0.15/)

  xmin = (/0,0,0,0,0/)
  xmax = (/100,max(com),max(lag),100,100/)
  ymin = (/-10,    0,-0.3,-40,-2/)
  ymax = (/ 10,10, 0.3, 40, 2/)
  xtick = (/20,5, 1, 20,20/)
  xmtck = (/3, 4, 0,  3, 3/)
  ytick = (/10,    5,0.1,20,1/)
  ymtck = (/0,    4,  4,1,  1/) 

  mark = (/0,2,0,0,0/)
  ystr = (/"Recon. series","Eigenvalue","EOF","PC","RC"/)
  xstr = (/"t (yr)","component","lag (yr)","t (yr)","t (yr)"/)

  data1 = new(3,graphic)
  data1(0) = create "00" coordArraysClass defaultapp
    "caXArray" : t
    "caYArray" : re1(:,0)
  end create
  data1(1) = create "00" coordArraysClass defaultapp
    "caXArray" : com
    "caYArray" : val
  end create
  do i=0, 1
  plot0 = create "plot0" xyPlotClass wks
    "xyCoordData" : data1(i)

    "trYMaxF" : ymax(i)
    "trYMinF" : ymin(i)
    "trXMaxF" : xmax(i)
    "trXMinF" : xmin(i)

    "tmXBMode" : "manual"
    "tmXBTickSpacingF" : xtick(i)
    "tmXBMinorPerMajor" : xmtck(i)
    "tmXBLabelFont" : 21
    "tmXBLabelFontHeightF" : 0.011
    "tmXBMajorLengthF" :  -0.008
    "tmXBMinorLengthF" :  -0.004
    "tmXBFormat" : "@*+^sg"

    "tmYLMinorLengthF" : -0.004
    "tmYLMajorLengthF" : -0.008
    "tmYLMode" : "MANUAL"
    "tmYLTickSpacingF" : ytick(i)
    "tmYLMinorPerMajor" : ymtck(i)
    "tmYLLabelFont" : 21
    "tmYLLabelFontHeightF" : 0.011
    "tmYLFormat" : "@*+^sg"

    "tiMainFont" : 22
    "tiMainFontHeightF" : 0.015
    "tiMainString" : ""
    "tiXAxisFont" : 21
    "tiXAxisFontHeightF" : 0.014
    "tiXAxisString" : xstr(i)
    "tiYAxisFont" : 21
    "tiYAxisFontHeightF" : 0.014
    "tiYAxisString" : ystr(i)

    "vpYF" : yf(i)
    "vpXF" : xf(i)
    "vpHeightF" : he(i)
    "vpWidthF" : wi(i)
  end create
  getvalues plot0
    "xyCoordDataSpec" : spec
  end getvalues
  setvalues spec
    "xyMarkLineMode" : mark(i)
;    "xyDashPattern"    : ldsh(i)
;    "xyLineThicknessF" : lthc(i)
  end setvalues
  draw(plot0)
  delete(plot0)
  end do


  xf2 = (/0.1,0.1,0.55,0.55/)
  lcol = (/"black","blue","black","blue"/)
  data2 = new((/3,4/),graphic)
  do j=0, 3
  data2(0,j) = create "00" coordArraysClass defaultapp
    "caXArray" : lag
    "caYArray" : eof(j,:,0)
  end create
  data2(1,j) = create "00" coordArraysClass defaultapp
    "caXArray" : t2(:1116)
    "caYArray" : pc1(j,:1116)
  end create
  data2(2,j) = create "00" coordArraysClass defaultapp
    "caXArray" : t(:) 
    "caYArray" : rc1(j,:,0)
  end create
  end do


  do i=2, 4
  do j=0, 3
  plot0 = create "plot0" xyPlotClass wks
    "xyCoordData" : data2(i-2,j)

    "trYMaxF" : ymax(i)
    "trYMinF" : ymin(i)
    "trXMaxF" : xmax(i)
    "trXMinF" : xmin(i)

    "tmXBMode" : "manual"
    "tmXBTickSpacingF" : xtick(i)
    "tmXBMinorPerMajor" : xmtck(i)
    "tmXBLabelFont" : 21
    "tmXBLabelFontHeightF" : 0.011
    "tmXBMajorLengthF" :  -0.008
    "tmXBMinorLengthF" :  -0.004
    "tmXBFormat" : "@*+^sg"

    "tmYLMinorLengthF" : -0.004
    "tmYLMajorLengthF" : -0.008
    "tmYLMode" : "MANUAL"
    "tmYLTickSpacingF" : ytick(i)
    "tmYLMinorPerMajor" : ymtck(i)
    "tmYLLabelFont" : 21
    "tmYLLabelFontHeightF" : 0.011
    "tmYLFormat" : "@*+^sg"

    "tiMainFont" : 22
    "tiMainFontHeightF" : 0.015
    "tiMainString" : ""
    "tiXAxisFont" : 21
    "tiXAxisFontHeightF" : 0.014
    "tiXAxisString" : xstr(i)
    "tiYAxisFont" : 21
    "tiYAxisFontHeightF" : 0.014
    "tiYAxisString" : ystr(i)

    "vpYF" : yf(i)
    "vpXF" : xf2(j)
    "vpHeightF" : he(i)
    "vpWidthF" : wi(i)
  end create
  getvalues plot0
    "xyCoordDataSpec" : spec
  end getvalues
  setvalues spec
    "xyMarkLineMode" : mark(i)
    "xyLineColor" : lcol(j)
  end setvalues
  draw(plot0)
  delete(plot0)
  end do
  end do

;-----------------------------------------------

frame(wks)

end



