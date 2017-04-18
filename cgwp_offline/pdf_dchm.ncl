load "$UTIL/ncl/header.ncl"
xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"
axisfontratio = 0.9

begin

;--------------------------------------- parameters
 expname = "L60CGW"
 expcode = "uanuj"
 y1 = 1979
 y2 = 2006
 lata = 10
 vars = (/"zcta","zcba","rho_ct"/)
 varo = (/"z:B:ct:N: ,  z:B:cb:N: [km]","",":F33:r:F21: [kg m:S:-3:N:]"/)
 ymin0 = 1.e-4
 ymin = (/0 ,0 ,-999,-999/)
 ymax = (/10,10, 999, 999/)

; ncol = (y2-y1+1)*12*(30*8)*(192*(toint(lata*2/1.25)+1))
 ncol = (y2-y1+1)*12*(30*8-1)*(192*(toint(lata*2/1.25)+1))  ; fail for 00 UTC 1

 wks = cr_wks("fig/pdf_dchm",600,0)
 gsn_define_colormap(wks,"default")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 fdir = "$DATD/"+expname+"/dchm_pdf"
 fname = expcode+".dchm_pdf."+y1+"-"+y2+".01-12.nc"
 f0 = addfile(fdir+"/"+fname,"r")

 nv = dimsizes(vars)

 x0 = f0->x_pdf
; dx = log10(x0(4)) - log10(x0(3))
 dx = round((log10(x0(4)) - log10(x0(3)))*1.e6,1)*1.e-6
 npop0 = f0->N_pop / tofloat(ncol) * 100.
 ix = ind( npop0 .ne. 0. )

 x = x0(ix)
 npop = npop0(ix)

 nx = dimsizes(x)

 vv = new((/nv,nx/),float)
 do iv=0, nv-1
   vv(iv,:) = f0->$vars(iv)$(0,ix)
 end do
 vv(0:1,:) = vv(0:1,:)*1.e-3  ; [km]

; do iv=0, nv-1
;   vv(iv,:) = where( npop .lt. ymin0, vv@_FillValue, vv(iv,:) )
; end do

 perc_cgw = sum(npop)
print("columns for CGWP : "+perc_cgw+" %")

;--------------------------------------- cn level
 lcol = (/1,1,2,3,6,4/)
 ldsh = (/0,1,0,0,0,0/)
 ltck = 0.8

;--------------------------------------- plot
 wi = 50
 he = 40
 xf = 10 + (/0,0,0,0/)*(wi+5)
 yf = 10 + (/0,0,0,0/)*(he+5)

 xmin = 1.e-3  ; x(1)

 xpg = new(nx*2+2,float)
 ypg = new(nx*2+2,float)
 do i=0, nx-1
   ii = i*2
   xpg(ii) = 10.^(log10(x(i)) - 0.5*dx)
   xpg(ii+1) = 10.^(log10(x(i)) + 0.5*dx)
   ypg(ii:ii+1) = npop(i)
 end do
 xpg(nx*2) = xpg(nx*2-1)
 ypg(nx*2) = 1.e-10  ; 0.
 xpg(nx*2+1) = xpg(0)
 ypg(nx*2+1) = ypg(nx*2)

 respg = True
 respg@gsFillColor = colind

 mstr = ""
 xstr = "Max. convective heating  [K hr:S:-1:N:]"
 ystr = "PDF [%]"
 
 txtpdfsum = ":F34:e:F21: PDF = "+(round(perc_cgw*10.,1)*0.1)+"%"


 plot0 = xylineplot(wks,x,npop,"D",0,"T",ltck,"C",-1)
 dum = gsn_add_polygon(wks,plot0,xpg,ypg,respg)
 setvalues plot0
  "tmXMajorGrid" : True
  "tmXMajorGridThicknessF" : 0.4
 end setvalues

 set_size(plot0,xf(0),yf(0),wi,he)
 axis_range(plot0,xmin,999,ymin0,999,True,True)
 axis_str(plot0,mstr,xstr,ystr)

 dr_txt(wks,txtpdfsum,xf(0)+wi-7,yf(0)-2.5,21,4,0.) 

 draw(plot0)


 xr = xf(0)+wi+3

do il=0, nv-1

 plot1 = xylineplot(wks,x,vv(il,:),"D",ldsh(il),"T",ltck,"C",lcol(il))

 set_size(plot1,xf(il),yf(il),wi,he)
 axis_range(plot1,xmin,999,ymin(il),ymax(il),True,False)
 axis_str(plot1,"","","")
 setvalues plot1
  "tmXBOn" : False
  "tmXTOn" : False
  "tmYLOn" : False
  "tmYROn" : False
 end setvalues

 draw(plot1)

 if (vars(il) .eq. "zcba") then
   continue
 end if

 plotx = xylineplot(wks,x,vv(il,:),"D",0,"T",0.001,"C",-1)
 set_size(plotx,xr,yf(il),1,he)
 axis_range(plotx,xmin,999,ymin(il),ymax(il),True,False)
 axis_str(plotx,"","",varo(il))
 setvalues plotx
  "tmXBOn" : False
  "tmXTOn" : False
  "tmYLOn" : False
  "tmYROn" : True
  "tmXBBorderOn" : False
  "tmXTBorderOn" : False
  "tmYLBorderOn" : False
;  "tmYLMode" : "Explicit"
;  "tmYLValues" : zs
;  "tmYLLabels" : ps
  "tmYUseLeft"   : False
  "tmBorderLineColor"  : lcol(il)
  "tmYRMajorLineColor" : lcol(il)
  "tmYRMinorLineColor" : lcol(il)
  "tmYRLabelsOn"     : True
  "tmYRLabelFontColor" : lcol(il)
  "tiYAxisSide"      : "Right"
  "tiYAxisAngleF"    : 270
  "tiYAxisFontColor" : lcol(il)
 end setvalues

 draw(plotx)

 xr = xr + 10
end do

 frame(wks)

end

