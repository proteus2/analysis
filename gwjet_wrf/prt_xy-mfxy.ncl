load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

;axisfontratio = 0.9

begin

;--------------------------------------- parameters
 fdir = "$DATD/gwjet_wrf"
 var1 = "W"
 var2 = "U"
 var3 = "V"
 lev = 300
 nave = 27 ;12*12
; n1s = ispan(nave/2-1,1151,nave)
; n1s = ispan(nave-1,1151,nave)
 n1s = ispan(719,1151-1,nave)
 var1name = "prt_"+var1+"_"+lev
 var2name = "prt_"+var2+"_"+lev
 var3name = "prt_"+var3+"_"+lev
;	Time = UNLIMITED ; // (1152 currently)
 
; wks = cr_wks("fig/prt16_xy-"+var1name+"-"+var2name,600,0)
 wks = cr_wks("fig/prt16_xy-prt_mfxy_"+lev,600,0)
; gsn_define_colormap(wks,"StepSeq25")
 gsn_define_colormap(wks,"cmp_b2r")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 fname1 = fdir+"/x1/prt_d300km/"+var1name+"__x1__08_00_05__12_00_00.nc"
 fname2 = fdir+"/x1/prt_d300km/"+var2name+"__x1__08_00_05__12_00_00.nc"
 fname3 = fdir+"/x1/prt_d300km/"+var3name+"__x1__08_00_05__12_00_00.nc"
 f1 = addfile(fname1,"r")
 f2 = addfile(fname2,"r")
 f3 = addfile(fname3,"r")
 nx = 432
 ny = 648+2
 dx = 1./10.8
 dy = 1./10.8
 ys = 20. - 0.5*dy  ; yn = 80. + 0.5*dy
 v1 = f1->$var1name$(n1s,:,:)
 v2 = f2->$var2name$(n1s,:,:)
 v3 = f3->$var3name$(n1s,:,:)

 vv = v1*v2
; vv = v1*v3 
; vv = sqrt( (v1*v2)^2 + (v1*v3)^2 )

;--------------------------------------- process 0
 nv = dimsizes(v1)
 nt = nv(0)
 ndim = dimsizes(nv)

 x = ispan(0,nx-1,1)*dx
 y = ys + ispan(0,ny-1,1)*dy

; xc = ispan(0,nx,1)*dx
; nvc = nv
; nvc(ndim-1) = nvc(ndim-1) + 1
; vc = new(nvc,float)
; vc(:,:,:nx-1) = v0(:,:,:nx-1)
; vc(:,:,nx   ) = v0(:,:,0    )

 ev = (lev/1.e3 - 0.252)*3.141592/2.
 sinlat = sin(y*3.141592/180.)
 coslat = cos(y*3.141592/180.)
 tz = 288.*(lev/1.e3)^(287.*5.e-3/9.80616)

 ty = tz + 0.75*(lev/1.e3)*3.141592*35./287.*sin(ev)*sqrt(cos(ev))* \
      ( (-2*sinlat^6*(coslat^2 + 1./3.) + 10./63.)*2.*35.*cos(ev)^1.5 + \
        (1.6*coslat^3*(sinlat^2 + 2./3.) - 3.141592/4.)*6371.229e3*7.29212e-5 )

 rho = conform( vv, lev*1.e2/287./ty, 1 )

 vv = vv*rho

;--------------------------------------- cn level
 cnlev1 = fspan(-10,10,21)*2.
 cnlev2 = fspan(-10,10,21)*3.
; cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)
 cnfil = ispan(2,65,3)

;--------------------------------------- contour
 wi = 32
 he = 20
 xf = 12 + (/0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3/)*(wi+3)
 yf = 10 + (/0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3/)*(he+6)

 mstr  = new(nt,string)
 xstr  = new(nt,string)
 ystr  = new(nt,string)
 lbopt = new(nt,string)

 mstr(:) = "Day "+(/"6  12::00","6  14::15","6  16::30","6  18::45", \
                    "6  21::00","6  23::15","7  01::30","7  03::45", \
                    "7  06::00","7  08::15","7  10::30","7  12::45", \
                    "7  15::00","7  17::15","7  19::30","7  21::45"/)
 xstr(:) = ""
; xstr(nt/2:) = "lon [deg]"
 xstr(12:) = "lon [deg]"
 ystr(:) = ""
 ystr(ispan(0,nt-1,nt/2)) = "lat [deg]"
 lbopt(:) = ""
 lbopt(nt-1) = "[mPa]"
 subl = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)"/)

do il=0, nt-1

 contour1 = cnshadeplot(wks,x,y,vv(il,:,:)*1.e3,cnlev1,cnfil,lbopt(il))
; cnaddline(contour1,wks,x,y,v0(il,:,:),cnlev2,1,"D",(/0,0,14/),"T",0.25,"C",(/1,-1,1/))
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.08
  "pmLabelBarParallelPosF"   : -1.20;+0.05 ;+ 0.55
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 1.5
  "pmLabelBarHeightF"        : 0.12
;  "lbLabelStride"      : 2
  "lbLabelFontHeightF" : 0.03
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.06
  "lbTitleFontHeightF" : 0.04
  "tmXBLabelsOn" : True
  "tmYLLabelsOn" : True
 end setvalues
 if (mod(il,4) .ne. 0) then
  setvalues contour1
   "tmYLLabelsOn" : False
  end setvalues
 end if
 if (il .lt. 12) then
  setvalues contour1
   "tmXBLabelsOn" : False
  end setvalues
 end if

 set_size(contour1,xf(il),yf(il),wi,he)
 axis_range(contour1,0,40,20,80,False,False)
 axis_str(contour1,"",xstr(il),ystr(il))

 cnfinalize(contour1,wks)

 draw(contour1)

 delete(contour1)

 dr_txt(wks,mstr(il),xf(il)+wi/2.,yf(il)-2.5,22,5.5,0.)

end do

 frame(wks)

end



