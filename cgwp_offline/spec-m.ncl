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
; fmiddle = "dchm"
 fmiddle = "dchm-midlev"
; fmiddle = "dchm-nonmidlev"
; dir0 = 0
 dir0 = 45
 y1 = 1979
 y2 = 2006

 wks = cr_wks("fig/spec-m_"+fmiddle+"-"+dir0,600,0)
; gsn_define_colormap(wks,"default")
 gsn_define_colormap(wks,"so4_21")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 fdir = "$DATD/"+expname+"/dchm_pdf"
 fname = expcode+"."+fmiddle+"_pdf."+y1+"-"+y2+".01-12.nc"
 f0 = addfile(fdir+"/"+fname,"r")
 xpdf0 = f0->x_pdf({0.1:10})

 fdir = "."
 fname = "res-"+dir0+"-m.nc"
 f1 = addfile(fdir+"/"+fname,"r")

 spct0 = f1->mflx_ct_spec
 c = f1->c_ph
 dir = f1->dir

 dc = c(1) - c(0)

 ipdf = 5
 xpdf = xpdf0(::ipdf)
 if (dir0 .eq. 0) then
   spec = spct0(::ipdf,0,:) - spct0(::ipdf,2,::-1)

   sptot = dim_sum(abs(spec))*dc

 else if (dir0 .eq. 45) then
   cosdir0 = cos(dir0/180.*3.141592)
;cosdir0 = 1.

   c = c/cosdir0
   dc2 = c(1) - c(0)

   spec = ( (spct0(::ipdf,0,:) - spct0(::ipdf,2,::-1)) +  \
            (spct0(::ipdf,3,:) - spct0(::ipdf,1,::-1)) )*(dc/dc2*cosdir0)

   sptot = dim_sum(abs(spec))*dc2

 end if
 end if
print(sptot)

 specn = spec/conform_dims(dimsizes(spec),sptot,0)

;--------------------------------------- cn level
 lcol = (/16,14,12,10,8/)
 ldsh = 0
 ltck = 0.8

;--------------------------------------- plot
 wi = 60
 he = 30
 xf = 10 + (/0,0,0,0,0/)*(wi+5)
 yf = 10 + (/0,0,0,0,0/)*(he+5)

 mstr = ""
 xstr = "phase speed [m s:S:-1:N:]"
 ystr = "normalized Mom. flux [ / m s:S:-1:N:]"
 
 norder = 10^(toint(log10(xpdf)+10)-10)
 xpdfr = round(xpdf/norder,1)*norder
 norder = 10^(toint(log10(sptot)+10)-10)
 sptotr = round(sptot/norder*10.,1)*norder/10.
 lgtxt = "Q:B:max:N: = "+xpdfr+" K/hr  ( M = "+(sptotr*1.e3)+" mPa )"

 plot1 = xylineplot(wks,c,specn,"D",ldsh,"T",ltck,"C",lcol)
 setvalues plot1
  "tmXMajorGrid" : True
  "tmXMajorGridThicknessF" : 0.4
 end setvalues

 set_size(plot1,xf(0),yf(0),wi,he)
 axis_range(plot1,-20,20,-0.5,0.5,False,False)
if (dir0 .eq. 45) then
 axis_range(plot1,-20,20,-0.3,0.3,False,False)
end if
 axis_str(plot1,mstr,xstr,ystr)

 lg = lgplot(wks,lgtxt,"D",ldsh,"T",ltck,"C",lcol,2)
 set_size_lg(lg,xf(0)+wi*0.5,yf(0)+he*0.6,wi*0.20,12.5)
 setvalues lg
  "lgLabelFontHeightF" : 0.005
 end setvalues
 draw(lg)

; dr_txt(wks,txtpdfsum,xf(0)+wi-7,yf(0)-2.5,21,4,0.) 

 draw(plot1)

 frame(wks)

end

