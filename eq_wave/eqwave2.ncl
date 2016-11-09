load "$UTIL/ncl/header.ncl"

xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 xmax = 20
 ymax = 1.0
 lzmin = 2.4
 lzmax = 50.

 beta = 2.*(7.292116e-5)/6371229.
 nbv = sqrt(6.)*1.e-2

 wks = cr_wks("eqwave2",600,0)
 colind = NhlNewColor(wks,0.4,0.4,0.4)

 setvalues wks
   "wkBackgroundColor" : (/-1.,-1.,-1./)
 end setvalues

;--------------------------------------- read data

;--------------------------------------- process 0
 x = fspan(-20,20,41*10)
 y = fspan(0,ymax,toint(ymax*90.)*10)

 nx = dimsizes(x)
 ny = dimsizes(y)

 k = x/6371229.
 o = y*2.*3.141592/86400.
 okbi = new((/ny,nx/),float)
 okbi = conform_dims((/ny,nx/),o,0)*conform_dims((/ny,nx/),k,1)/beta

 lzigw1 = new((/ny,nx/),float,1.e30)
 lzigw1 = 1.e30
 lzigw2 = lzigw1
 lzrw1  = lzigw1
 lzrw2  = lzigw1
 lzrgw  = lzigw1
 lzkw   = lzigw1

 do j=1, ny-1
   lzigw1(j,:) = 2.*3.141592/1.e3/ \
      (nbv*beta/o(j)^2*(1.5+sqrt(1.5^2+okbi(j,:)*(1.+okbi(j,:)))))
   lzigw2(j,:) = 2.*3.141592/1.e3/ \
      (nbv*beta/o(j)^2*(2.5+sqrt(2.5^2+okbi(j,:)*(1.+okbi(j,:)))))
 end do

 do j=1, ny-1
 do i=0, nx/2-1
   if (okbi(j,i) .gt. -1) then
     lzrw1(j,i) = 2.*3.141592/1.e3/ \
        (nbv*beta/o(j)^2*(1.5-sqrt(1.5^2+okbi(j,i)*(1.+okbi(j,i)))))
     lzrw2(j,i) = 2.*3.141592/1.e3/ \
        (nbv*beta/o(j)^2*(2.5-sqrt(2.5^2+okbi(j,i)*(1.+okbi(j,i)))))
   end if
 end do
 end do

 do j=1, ny-1
 do i=0, nx-1
   if (okbi(j,i) .gt. -1) then
     lzrgw(j,i) = 2.*3.141592/1.e3/(nbv*beta/o(j)^2*(1.+okbi(j,i)))
   end if
 end do
 end do

 do j=1, ny-1
 do i=nx/2+1, nx-1
   lzkw(j,i) = 2.*3.141592/1.e3/ (nbv*k(i)/o(j))
 end do
 end do

 lz1 = (/lzigw1,lzrw1,lzkw/)
 lz2 = (/lzigw2,lzrw2,lzrgw/)

 dome = 1./90.

;--------------------------------------- cn level
 cnlev = (/lzmin,lzmax/)
 lcol1 = (/1,2,3/)
 lcol2 = (/1,2,4/)

;--------------------------------------- contour
 wi = 32
 he = 32
 xf = 10 + (/0,0/)
 yf = 10 + (/1,2/)*40

do il=0, 2

 contour1 = cnlineplot(wks,x,y,lz1(il,:,:),cnlev,1)
 cnline_att(contour1,"D",0,"T",1.,"C",lcol1(il))
 setvalues contour1
  "cnLineLabelFontHeightF" : 0.022
  "tmXBOn" : False
  "tmXTOn" : False
  "tmYLOn" : False
  "tmYROn" : False
 end setvalues
 cnline_att(contour1,"D",0,"T",1.,"C",lcol1(il))
 set_size(contour1,xf(0),yf(0),wi,he)
 axis_range(contour1,-xmax-0.5,xmax+0.5,0,ymax+0.5*dome,False,False)
; axis_str(contour1,"Symmetric component","","intrinsic frequency [cyc day:S:-1:N:]")
 cnfinalize(contour1,wks)
 draw(contour1)

 contour2 = cnlineplot(wks,x,y,lz2(il,:,:),cnlev,1)
 cnline_att(contour2,"D",0,"T",1.,"C",lcol2(il))
 setvalues contour2
  "cnLineLabelFontHeightF" : 0.022
  "tmXBOn" : False
  "tmXTOn" : False
  "tmYLOn" : False
  "tmYROn" : False
 end setvalues
 set_size(contour2,xf(1),yf(1),wi,he)
 axis_range(contour2,-xmax-0.5,xmax+0.5,0,ymax+0.5*dome,False,False)
; axis_str(contour2,"Anti-symmetric component","zonal wavenumber","intrinsic frequency [cyc day:S:-1:N:]")
 cnfinalize(contour2,wks)
 draw(contour2)

end do

 frame(wks)

end

