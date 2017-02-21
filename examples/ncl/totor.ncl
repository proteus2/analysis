load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin
;-------------------------------------------  
  stid = 47185 
  dz =50
  tconst =10  
  wconst =20 
;-------------------------------------------  
  fnam = "~miok/KEOP05/"+sprinti("%0.2i",stid)+"/per"+sprinti("%0.2i",stid)
  data = asciiread(fnam,-1,"string")
  date = asciiread(fnam,-1,"integer")
  nn = dimsizes(date)
  time = new(nn,integer)
  hr = new(nn,integer)
  day = new(nn,integer)
  time = date - 60000
do i=0,nn-1
  if(time(i) .lt. 10000) then 
    hr(i) = time(i) - floattoint(time(i)/100.)*100
    day(i) = floattoint(time(i)/100.)*24
    time(i) = day(i)+hr(i)
  else
    time(i) = time(i) - 10000
    hr(i) = time(i) - floattoint(time(i)/100.)*100
    day(i) = (floattoint(time(i)/100.)+30)*24
    time(i) = day(i) +hr(i)
  end if 
end do
time = (time - 360)/6
  fdat = "/usr/users/miok/KEOP05/"+stid+"/ort"
  gdat = "/usr/users/miok/KEOP05/"+stid+"/orw"
;--------------------------------------------

  cmap = (/   (/1.000,1.000,1.000/), \
              (/0.000,0.000,0.000/), \
              (/1.000,0.000,0.000/), \
              (/0.000,0.353,0.765/), \
              (/0.000,0.706,0.294/), \
              (/0.294,1.000,0.000/), \
              (/0.588,1.000,0.000/), \
              (/0.765,1.000,0.000/), \
              (/0.882,1.000,0.000/), \
              (/1.000,0.941,0.000/), \
              (/0.750,0.750,0.750/)  /)

  wks = create "1sttry" psWorkstationClass defaultapp
    "wkPSFileName" : "or_"+stid+".ps"
    "wkColorMap"   : cmap
    end create

do j=0,nn-1
;------------------------------------------

   f  = addfile(fdat+data(j)+".nc","r")
   g  = addfile(gdat+data(j)+".nc","r")
   T  = f->T -273.
   U  = g->U
   V  = g->V
   tz = f->Z/1.e3
   wz = g->Z/1.e3
   nt = dimsizes(tz)
   nw = dimsizes(wz)
;------------------------------------------

   var1 = new((/1,nt/),"float")
   var1(0,:) = T(:)+time(j)*tconst
   var2 = new((/2,nw/),"float")
   var2(0,:) = U(:)+time(j)*wconst
   var2(1,:) = V(:)+time(j)*wconst

   xydata = new(3,graphic)
   xydata(0)= create "xydata" coordArraysClass defaultapp
      "caXArray"  : var1(0,:)
      "caYArray"  : tz
     end create
  do i=1,2
   xydata(i)= create "xydata" coordArraysClass defaultapp
      "caXArray"  : var2(i-1,:)
      "caYArray"  : wz
     end create
  end do
;-----------------------------------------
  xstr = (/"T (:F35:J:F22:C)","U (m s:S:-1:N:)","V (m s:S:-1:N:)"/)
  ystr = (/"height (km)","height (km)","height (km)"/)

  xmin = (/-80, -40,-40/)
  xmax = (/40+tconst*(time(nn-1)-1),70+wconst*(time(nn-1)-1),60+wconst*(time(nn-1)-1)/)
  xtick = (/4*tconst,4*wconst,4*wconst/)
  vy = (/0.95,0.62,0.29/)
  ticksp = (/50,100,100/)
  minortick = (/4,9,9/)
;------------------------------------------
lbx = (/"15","16","17","18","19","20","21","22","23","24",\
        "25","26","27","28","29",  "30","1", "2","3", "4", \
        "5", "6", "7", "8", "9","10","11","12","13","14","15","16","17","18"/)
lbxf = new((/3,dimsizes(lbx)/),float)
do i=0,dimsizes(lbx)-1
  lbxf(0,i) = i*4*tconst -50 
  lbxf(1,i) = i*4*wconst
  lbxf(2,i) = i*4*wconst
end do

;------------------------------------------
  do i=0,2

   plot =  create "lplot" xyPlotClass wks

   "xyCoordData" : xydata(i)

   "trYMaxF" : 34 
   "trYMinF" : 0 
   "trXMaxF" : xmax(i)
   "trXMinF" : xmin(i)

   "tmXBMode" : "explicit"
   "tmXBLabelFont" : 21
;   "tmXBTickSpacingF" :1
;   "tmXBMinorPerMajor" :1

   "tmXBLabelFontHeightF" : 0.012
   "tmXBMajorLengthF" : -0.006
   "tmXBMinorLengthF" : -0.003
;   "tmXBLabelsOn"      :True
   "tmXBLabels"     :lbx 
   "tmXBValues"      :lbxf(i,:)

   "tmXUseBottom"  : False
   "tmXTLabelsOn"   : True
   "tmXTMode"      : "Manual"
   "tmXTTickSpacingF" : ticksp(i)
   "tmXTMinorPerMajor": minortick(i) 
   "tmXTLabelFont"    : 21
   "tmXTLabelFontHeightF" : 0.0065
   "tmXTFormat"  : "@*+^sg"

   "tmYLMinorLengthF" : -0.004
   "tmYLMajorLengthF" : -0.008
   "tmYLMode" : "MANUAL"
   "tmYLTickSpacingF" : 10 
   "tmYLMinorPerMajor" : 4
   "tmYLLabelFont" : 21
   "tmYLLabelFontHeightF" : 0.011
   "tmYLFormat" : "@*+^sg"

   "tiMainFont" : 22
   "tiMainFontHeightF" : 0.015
   "tiMainString" :xstr(i) 
   "tiXAxisFont" : 22
   "tiXAxisFontHeightF" : 0.014
   "tiXAxisString" :"Day"  
   "tiYAxisFont" : 22
   "tiYAxisFontHeightF" : 0.014
   "tiYAxisString" : ystr(i)
 
   "vpYF" : vy(i)
   "vpXF" : 0.08 
   "vpHeightF" : 0.2
   "vpWidthF" : 0.9

    end create

    if (j .eq. 0) then
      xbox = (/xmin(i),xmin(i),xmax(i),xmax(i)/)
      ybox = (/30,17,17,30/)
      resp = True
      resp@gsFillColor = 10
      gsn_polygon(wks,plot,xbox,ybox,resp)
    end if

    draw(plot)

;------------------------------------------------
   end do  ; i
;------------------------------------------------

  delete(f)
  delete(g)
  delete(T)
  delete(U)
  delete(V)
  delete(tz)
  delete(wz)
  delete(var1)
  delete(var2)
  delete(xydata)
  delete(plot)
;------------------------------------------------
end do ;j

;do j=0,2
;   tid = create "Text" textItemClass wks
;     "txString" : "June"
;     "txPosXF"  : 0.10 + 0.02
;     "txPosYF"  : vy(j)-0.23-0.02
;    "txFont"   : 21
;     "txFontHeightF" : 0.014
;   end create
;   draw(tid)

;   setvalues tid
;     "txString" : "July"
;     "txPosXF"  : 0.10+0.85 - 0.02
;   end setvalues
;   draw(tid)
;   delete(tid)
;end do
   frame(wks)
  end


