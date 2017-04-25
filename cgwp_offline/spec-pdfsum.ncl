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
; fmiddle = "dchm-midlev"
 fmiddle = "dchm-nonmidlev"
; dir0 = 0
 dir0 = 45
 y1 = 1979
 y2 = 2006
 lata = 10
 z0 = 16.25e3  ; 100 hPa

; ncol = (y2-y1+1)*12*(30*8)*(192*(toint(lata*2/1.25)+1))
 ncol = (y2-y1+1)*12*(30*8-1)*(192*(toint(lata*2/1.25)+1))  ; fail for 00 UTC 1

 wks = cr_wks("fig/spec-pdfsum_"+fmiddle+"-"+dir0,600,0)
 gsn_define_colormap(wks,"default")
; gsn_define_colormap(wks,"so4_21")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 fdir = "$DATD/"+expname+"/dchm_pdf"
 fname = expcode+"."+fmiddle+"_pdf."+y1+"-"+y2+".01-12.nc"
 f0 = addfile(fdir+"/"+fname,"r")
 xpdf0 = f0->x_pdf({0.1:10})

 npop0 = f0->N_pop({0.1:10}) / tofloat(ncol) * 100.

 fdir = "."
 fname = "res-"+dir0+".nc"
 f1 = addfile(fdir+"/"+fname,"r")

 spct0 = f1->mflx_ct_spec
 sp0 = f1->mflx_spec(:,{z0},:,:)
 c = f1->c_ph
 dir = f1->dir

 dc = c(1) - c(0)

 nc = dimsizes(c)
 nx = dimsizes(xpdf0)

 xpdf = xpdf0
 if (dir0 .eq. 0) then
   spec_pos = dim_sum_n(spct0(:,0,:)*conform_dims((/nx,nc/),npop0/100.,0),0)
   spec_neg = dim_sum_n(-spct0(:,2,::-1)*conform_dims((/nx,nc/),npop0/100.,0),0)

   spec_pos2 = dim_sum_n(sp0(:,0,:)*conform_dims((/nx,nc/),npop0/100.,0),0)
   spec_neg2 = dim_sum_n(-sp0(:,2,::-1)*conform_dims((/nx,nc/),npop0/100.,0),0)

   sptot = (sum(abs(spec_pos))+sum(abs(spec_neg)))*dc

   sptot2 = (sum(abs(spec_pos2))+sum(abs(spec_neg2)))*dc

 else if (dir0 .eq. 45) then
   cosdir0 = cos(dir0/180.*3.141592)
;cosdir0 = 1.

   c = c/cosdir0
   dc2 = c(1) - c(0)

   spec_pos = dim_sum_n((spct0(:,0,:) + spct0(:,3,:))*(dc/dc2*cosdir0)*  \
                        conform_dims((/nx,nc/),npop0/100.,0),0)
   spec_neg = dim_sum_n(-(spct0(:,1,::-1) + spct0(:,2,::-1))*(dc/dc2*cosdir0)*  \
                        conform_dims((/nx,nc/),npop0/100.,0),0)

   spec_pos2 = dim_sum_n((sp0(:,0,:) + sp0(:,3,:))*(dc/dc2*cosdir0)*  \
                        conform_dims((/nx,nc/),npop0/100.,0),0)
   spec_neg2 = dim_sum_n(-(sp0(:,1,::-1) + sp0(:,2,::-1))*(dc/dc2*cosdir0)*  \
                        conform_dims((/nx,nc/),npop0/100.,0),0)

   sptot = (sum(abs(spec_pos))+sum(abs(spec_neg)))*dc2

   sptot2 = (sum(abs(spec_pos2))+sum(abs(spec_neg2)))*dc2

 end if
 end if
print(sptot)
print(sptot2)
print(sum(abs(spec_pos2))*dc2+"  "+sum(abs(spec_neg2))*dc2)

;--------------------------------------- cn level
 lcol = (/1,1,2,2/)
 ldsh = (/0,1,0,1/)
 ltck = 0.8

;--------------------------------------- plot
 wi = 60
 he = 30
 xf = 10 + (/0,0,0,0,0/)*(wi+5)
 yf = 10 + (/0,0,0,0,0/)*(he+5)

 mstr = ""
 xstr = "phase speed [m s:S:-1:N:]"
 ystr = "Mom. flux [mPa / m s:S:-1:N:]"
 
 norder = 10^(toint(log10(sptot2)+10)-10)
 sptotr = round(sptot2/norder*10.,1)*norder/10.
 txtstr = "|M:B:East:N:| + |M:B:West:N:| = "+(sptotr*1.e3)+" mPa"

 plot1 = xylineplot(wks,c,(/spec_pos,spec_neg,spec_pos2,spec_neg2/)*1.e3,"D",ldsh,"T",ltck,"C",lcol)
 xyaddline(plot1,c,0,"D",0,"T",0.5,"C",colind)
 setvalues plot1
  "tmXMajorGrid" : True
  "tmXMajorGridThicknessF" : 0.4
 end setvalues

 set_size(plot1,xf(0),yf(0),wi,he)
 axis_range(plot1,-25,25,-999,999,False,False)
 axis_str(plot1,mstr,xstr,ystr)

 dr_txt(wks,txtstr,xf(0)+wi-12,yf(0)-2.5,21,4,0.) 

 draw(plot1)

 frame(wks)

end

