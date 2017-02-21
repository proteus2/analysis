;23456 
      begin
      
      MIS=1.e+32
 
      fname = "/export15/kjs/js_paper/ncl/exiall.nc"
      
      fn = addfile(fname,"R")
      ns = fn->NS
      val = fn->VAL

      ns2  = new(8,float)
      day2 = new(137,float)
      val2 = new((/8,137/),float)
      val3 = new((/8,137/),float)
 
      do j=0,7
        do i=0,136
          val2(j,i)=MIS
        end do
      end do
 
      do j=0,7
        ns2(j) = 1.0*j
      end do
      do i=0,136
        day2(i)=i*1.0 
      end do
      do j=1,6
        do i=0,67
          k=i*2
          val2(j,k+1)=val(j-1,i)
        end do
      end do 

      val3(0,:)=val2(0,:)
      val3(1,:)=val2(2,:)
      val3(2,:)=val2(3,:)
      val3(3,:)=val2(5,:)
      val3(4,:)=val2(4,:)
      val3(5,:)=val2(6,:)
      val3(6,:)=val2(1,:)
      val3(7,:)=val2(7,:)


     cmap = (/(/1.000,1.000,1.000/), (/0.000,0.000,0.000/), \
            (/0.750,1.000,0.000/), (/1.000,0.693,0.000/), \
            (/0.000,0.333,0.667/), (/0.000,0.667,0.333/), \
            (/0.864,0.000,0.136/), (/0.970,0.970,0.000/), \
            (/0.750,0.250,0.750/), (/0.000,0.500,0.500/) /)
     
;      cmap = (/(/1.000,1.000,1.000/), (/0.000,0.000,0.000/), \
;               (/0.075,0.075,0.075/), (/0.100,0.100,0.100/), \
;               (/0.320,0.320,0.320/), (/0.470,0.470,0.470/), \
;               (/0.620,0.620,0.620/), (/0.750,0.750,0.750/), \
;               (/0.875,0.875,0.875/), (/0.900,0.900,0.900/) /)

      title="SONDE-DATA USED FOR ANALYSIS"
      CN = (/1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0/)     
;      FV = (/2, 4, 6, 5, 7, 8/) 
;      FV = (/1, 8, 7, 5, 6, 4/) 
;      FV = (/2, 3, 4, 5, 6, 7/) 
;      FV = (/3, 4, 5, 6, 7, 8/) 
      FV = (/8,6,3, 5,7,4/)
      labelsx = (/"20", "21", "22", "23", "24", "25", "26", "27", "28", \
                 "29", "30", "31", "01", "02", "03", "04", "05"        /)
      labellocsx = (/1, 9, 17, 25, 33, 41, 49, 57, 65, 73, 81, 89, 97, 105, 113, 121, 129/)
;      labelsy = (/"B","G","K","O","P","S"/)
      labelsy = (/"G","K","P","O","S","B"/)
      labellocsy = (/1,2,3,4,5,6/)

      wks = create "WPWork" psWorkstationClass defaultapp
        "wkPSFileName"           :    "dataex_all_col.ps" 
        "wkColorMap"             :    cmap
      end create

      titleid = create "title" textItemClass wks
        "txFont"                 :  26
        "txPosXF"                :  0.5
        "txPosYF"                :  0.955
        "txFontHeightF"          :  0.022
        "txString"               :  title
      end create 

      dataid = create "dataid" scalarFieldClass defaultapp
        "sfDataArray"            :   val3
        "sfYArray"               :   ns2
        "sfXArray"               :   day2
        "sfYCStartV"             :   0.0
        "sfYCEndV"               :   7.0
        "sfMissingValueV"        :   MIS
      end create

      cnplot = create "cnplot" contourPlotClass wks
        "cnScalarFieldData"      :  dataid
        "vpXF"                   :  0.10
        "vpYF"                   :  0.88
        "vpWidthF"               :  0.80
        "vpHeightF"              :  0.40
        "cnLevelSelectionMode"   : "ExplicitLevels"
        "cnInfoLabelOn"          : False
        "cnLevels"               :  CN
        "cnFillColors"           :  FV
        "cnFillOn"               : True
        "cnLinesOn"              : False
        "cnLineLabelsOn"         : False
        "cnRasterModeOn"         : True
        "tmXBMode"               : "explicit"
        "tmXBValues"             : labellocsx
        "tmXBLabels"             : labelsx
        "tmXBFormat"             : "f"
        "tmXBLabelFontHeightF"   : 0.010
        "tmXBMinorPerMajor"      : 1
        "tmYLMode"               : "explicit"
        "tmYLMajorLengthF"       : 0.000
        "tmYLValues"             : labellocsy
        "tmYLLabels"             : labelsy
        "tmYLLabelFontHeightF"   : 0.014
        "tmYLMinorPerMajor"      : 1
        "tiXAxisString"          : "TIME (day)"
        "tiXAxisFontHeightF"     : 0.012
        "tiYAxisString"          : ""
        "tiYAxisFontHeightF"     : 0.012
        "tiMainString"           : ""
        "tiMainPosition"         : "Left"
        "tiMainFontHeightF"      : 0.012
      end create

      lblid = create "labelbar" labelBarClass wks
        "vpXF"                   : 0.25
        "vpYF"                   : 0.10
        "vpWidthF"               : 0.5
        "vpHeightF"              : 0.08
        "lbOrientation"          : "horizontal"
        "lbAutoManage"           : False
        "lbLabelFont"            : 22
        "lbLabelFontHeightF"     : 0.011
        "lbMonoFillPattern"      : True
        "lbLabelAlignment"       : "ExternalEdges"
        "lbPerimOn"              : False
        "lbLabelJust"            : "CenterCenter"
        "lbFillColors"           : FV
        "lbLabelStrings"         : ""
        "lbBoxCount"             : 7
        "lbTitlePosition"        : "Bottom"
        "lbTitleFont"            : 22
        "lbTitleFontHeightF"     : 0.018
      end create

;      draw(lblid)
;      draw(titleid)
      draw(cnplot)

      frame(wks)

      end 
