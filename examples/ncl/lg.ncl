begin


  wks = create "work" psWorkstationClass defaultapp
    "wkPSFileName" : "lg.ps"
    "wkColorMap"   : "WhBlGrYeRe"
  end create


  lg1 = create "lg1" legendClass wks
    "lgAutoManage"       : False
    "lgItemCount"        : 4
    "lgLabelStrings"     : (/"October","July","April","January"/)

    "lgMonoDashIndex"      : True
    "lgMonoLineThickness"  : True
    "lgLineThicknessF"     : 2.5
    "lgLineColors"         : (/95,1,48,25/)
    "lgLabelFontHeightF" : 0.013
    "lgLabelJust"        : "CenterLeft"
    "lgLabelOffsetF"     : 0.10
    "lgPerimOn"          : False

    "vpXF"        : 0.50
    "vpYF"        : 0.90
    "vpWidthF"    : 0.25
    "vpHeightF"   : 0.25
  end create
  draw(lg1)

  frame(wks)


end




