     begin

     lev1=34 ; 20
     lev2=40 ; 23
     lev3=46 ; 26
     lev4=52 ; 29
     lev5=58 ; 32
     lev6=64 ; 35
     lev7=70 ; 38
     lev8=76 ; 41
     lev9=82 ; 44

     load "nclscr/gsn_code.ncl"
     load "nclscr/gsn_csm.ncl"
     load "nclscr/contour.ncl"

     fname1="res/temp2/plotdata/mfx-up_kh_phi_plot-data_p1_f1.nc"
     file1=addfile(fname1,"R")
     kh=file1->kh
     z=file1->height
     ang=file1->direction
     mfx1=file1->mfx

     fname1="res/temp2/plotdata/mfy-up_kh_phi_plot-data_p1_f1.nc"
     file1=addfile(fname1,"R")
     mfy1=file1->mfy

     nang=dimsizes(ang)
     nx=dimsizes(kh)

     mf1_1=new((/nx,nang/),float)
     mf1_2=new((/nx,nang/),float)
     mf1_3=new((/nx,nang/),float)
     mf2_1=new((/nx,nang/),float)
     mf2_2=new((/nx,nang/),float)
     mf2_3=new((/nx,nang/),float)
     mf3_1=new((/nx,nang/),float)
     mf3_2=new((/nx,nang/),float)
     mf3_3=new((/nx,nang/),float)

     do j=0, nang-1
       do i=0, nx-1
         mf1_1(i,j)=sqrt(mfx1(lev1,j,i)*mfx1(lev1,j,i)+mfy1(lev1,j,i)*mfy1(lev1,j,i))
         mf1_2(i,j)=sqrt(mfx1(lev2,j,i)*mfx1(lev2,j,i)+mfy1(lev2,j,i)*mfy1(lev2,j,i))
         mf1_3(i,j)=sqrt(mfx1(lev3,j,i)*mfx1(lev3,j,i)+mfy1(lev3,j,i)*mfy1(lev3,j,i))
         mf2_1(i,j)=sqrt(mfx1(lev4,j,i)*mfx1(lev4,j,i)+mfy1(lev4,j,i)*mfy1(lev4,j,i))
         mf2_2(i,j)=sqrt(mfx1(lev5,j,i)*mfx1(lev5,j,i)+mfy1(lev5,j,i)*mfy1(lev5,j,i))
         mf2_3(i,j)=sqrt(mfx1(lev6,j,i)*mfx1(lev6,j,i)+mfy1(lev6,j,i)*mfy1(lev6,j,i))
         mf3_1(i,j)=sqrt(mfx1(lev7,j,i)*mfx1(lev7,j,i)+mfy1(lev7,j,i)*mfy1(lev7,j,i))
         mf3_2(i,j)=sqrt(mfx1(lev8,j,i)*mfx1(lev8,j,i)+mfy1(lev8,j,i)*mfy1(lev8,j,i))
         mf3_3(i,j)=sqrt(mfx1(lev9,j,i)*mfx1(lev9,j,i)+mfy1(lev9,j,i)*mfy1(lev9,j,i))
       end do
     end do

     wks = create "WPWork" psWorkstationClass defaultapp
       "wkPSFileName"          : "mf_kh_dir-up.eps"
       "wkColorMap"            : "div_h1_009"
       "wkPSFormat"            : "eps"
     end create

     plot=new(1,graphic)

       res=True

       res@vpXF                    =  0.05
       res@vpYF                    =  1.1
       res@vpWidthF                =  0.21

       res@mpMinLatF               =  70.    ; 0-60 m/s
;       res@mpMinLatF               =  77.5   ; 0-50 m/s
       res@mpCenterLonF            =  270. 
       res@mpFillOn                =  False
       res@mpOutlineOn             =  False
       res@mpGridLineDashPattern   =  1
       res@mpGridLonSpacingF       =  30 
       res@mpGridLatSpacingF       =  5. 
       res@mpGridAndLimbOn         =  True
       res@mpGridAndLimbDrawOrder  =  "PostDraw" 
       res@mpGridLineColor         =  1 
       res@mpGridLineThicknessF    =  1.0
       res@mpGridMaxLatF           =  90.

       res@gsnPolarLabelFontHeightF = 0.012
       res@gsnPolarLabelFont        = "helvetica"
       res@gsnPolarLabelDistance    =  1.1

       res@lbLabelBarOn            =  False

       res@cnInfoLabelOn           =  False
       res@cnHighLabelsOn          =  False
       res@cnLowLabelsOn           =  False
       res@cnLinesOn               =  False
       res@cnLineLabelsOn          =  False
       res@cnFillOn                =  True
       res@cnLevelSelectionMode    =  "ExplicitLevels"
       res@cnLevels                =  (/0.05,0.1,0.2,0.4,0.6,0.8,1.0/)*10
       res@cnFillColors            =  (/2,3,4,6,7,8,9,10/)

       res@cnFillDrawOrder         =  "draw"

       res@tiMainString            =  "z="+z(lev1)/1000.+" km (P1)"
       res@tiMainFontHeightF       =  0.016
       res@tiMainFont              =  "helvetica-bold"
       res@tiMainOffsetYF          =  0.01
       res@tiMainOffsetXF          =  0.00
       res@gsnFrame                =  False       

       plot(0)=contour_angle(wks,mf1_1,res)


       res@vpYF                    =  0.8
       res@tiMainString            =  "z="+z(lev2)/1000.+" km (P1)"
       plot(0)=contour_angle(wks,mf1_2,res)

       res@vpYF                    =  0.5
       res@tiMainString            =  "z="+z(lev3)/1000.+" km (P1)"
       plot(0)=contour_angle(wks,mf1_3,res)



       res@vpXF                    =  0.35
       res@vpYF                    =  1.1
       res@tiMainString            =  "z="+z(lev4)/1000.+" km (P2)"
       plot(0)=contour_angle(wks,mf2_1,res)

       res@vpYF                    =  0.8
       res@tiMainString            =  "z="+z(lev5)/1000.+" km (P2)"
       plot(0)=contour_angle(wks,mf2_2,res)

       res@vpYF                    =  0.5
       res@tiMainString            =  "z="+z(lev6)/1000.+" km (P2)"
       plot(0)=contour_angle(wks,mf2_3,res)


       res@vpXF                    =  0.65
       res@vpYF                    =  1.1
       res@tiMainString            =  "z="+z(lev7)/1000.+" km (P3)"
       plot(0)=contour_angle(wks,mf3_1,res)

       res@vpYF                    =  0.8
       res@tiMainString            =  "z="+z(lev8)/1000.+" km (P3)"
       plot(0)=contour_angle(wks,mf3_2,res)
       
       res@vpYF                    =  0.5
       res@tiMainString            =  "z="+z(lev9)/1000.+" km (P3)"
       plot(0)=contour_angle(wks,mf3_3,res)


     labelid1 = create "LabelBar" labelBarClass wks
       "vpXF"                    :  0.25
       "vpYF"                    :  0.05
       "vpWidthF"                :  0.4
       "vpHeightF"               :  0.04
       "lbBoxCount"              :  dimsizes(res@cnFillColors)
       "lbLabelStrings"          :  res@cnLevels * 1.e6
       "lbFillColors"            :  res@cnFillColors
       "lbLabelAlignment"        :  "InteriorEdges"
       "lbMonoFillPattern"       :  True
       "lbOrientation"           :  "horizontal"
       "lbBoxLinesOn"            :  True
       "lbLabelsOn"              :  True
       "lbPerimOn"               :  False
       "lbAutoManage"            :  False
       "lbLabelFontHeightF"      :  0.011
       "lbLabelFont"             :  "helvetica"
       "lbTitleString"           :  "MF [x10:S:-6:N: N m:S:-2:N:/(m s:S:-1:N:)/deg]"
       "lbTitleAngleF"           :  0
       "lbTitlePosition"         :  "bottom"
       "lbTitleDirection"        :  "Across"
       "lbTitleFont"             :  "helvetica"
       "lbTitleFontHeightF"      :  0.011
       "lbTitleOffsetF"          :  0.1
     end  create
     draw(labelid1)

     frame(wks)
     end 
