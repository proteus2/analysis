load "$UTIL/ncl/header.ncl"

xspace = 83   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

;axisfontratio = 1.0

begin

;--------------------------------------- parameters
 expname = "era-int_f"
 dirname = "ERA-I-nr"
 uname   = "u"
; expname = "jra55c_f"
; dirname = "JRA55C-nr"
; uname   = "U_GRD_GDS4_HYBL"

 lat0 = 25.
 lat0avg = 15
 p0    = 100  ; 54, 66, 80, 96
if (expname .eq. "era-int_f") then
 p0u   = 23   ; 22, 23, 24, 25
else
 p0u   = 38   ; 39, 38, 37, 36
end if
 ymax0 = 0.5

 x_wmrgw = 1 ;3

 wks = cr_wks("fig/epf-kw-p"+numeric2int(p0,1)+"-"+expname(0),600,0)
 gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam0 = "$DATD/clim/"+dirname(0)+"/"+expname(0)+".u_ypt.nc"
   fnam1 = "$DATD/tem/"+dirname(0)+"/epf/"+expname(0)+".epfsa_yz_recon0.nc"
   fnam2 = "$DATD/tem/"+dirname(0)+"/epf/"+expname(0)+".epf_yz_recon0.nc"
print(fnam1+"")

   f0 = addfile(fnam0,"r")
   f1 = addfile(fnam1,"r")
   f2 = addfile(fnam2,"r")
   lat1 = f1->lat({-lat0:lat0})

   ny = dimsizes(lat1)
   nt = 1 + 407 + 36
if (expname .eq. "jra55c_f") then
  nt = 1 + 407 - 12
end if
   nw = 8

   u = f0->$uname(0)$(:nt-1,{p0u},{-lat0:lat0})

   rw_i = new((/nt-1,nw,ny/),float)
   rw_i(:,:,:) = f1->f_z(:nt-2,:nw-1,{p0},{-lat0:lat0})
   rw_i(:,nw-2:,:) = rw_i(:,nw-2:,:) + \
                     f2->f_z(:nt-2,8:9,{p0},{-lat0:lat0})
   rw_i(:,nw-1,:) = rw_i(:,nw-1,:) + \
                    ( f2->f_z(:nt-2,2,{p0},{-lat0:lat0}) + \
                      f2->f_z(:nt-2,3,{p0},{-lat0:lat0}) )

   rw_i(:,3,:) = rw_i(:,3,:) + f1->f_z(:nt-2,13,{p0},{-lat0:lat0})
   rw_i(:,6:7,:) = rw_i(:,6:7,:) + f1->f_z(:nt-2,16:17,{p0},{-lat0:lat0})

   rw_i = rw_i / 6371.e3 * 1.e3

   rw0 = new((/nt,nw,ny/),float,1.e20)
   rw0!2 = "lat"
   rw0&lat = lat1
   rw0(1:,:,:) = rw_i(:,:,:)
;   rw0(0,:,:) = rw0(1,:,:)
   rw0(0,:,:) = 1.e20

   rw = (/ rw0(:,2:3,:), rw0(:,4:5,:), rw0(:,6:7,:)/)

   fz = new((/3,2,ny,nt/),float,1.e20)
   fz!2 = "lat"
   fz&lat = lat1
   do n=0, nt-1
     fz(:,:,:,n) = (/ rw(:,n,:,:) /)
   end do

   t = fspan(0,13,14)

   fzss = new((/3,2,ny,14/),float,1.e20)
   fzss!2 = "lat"
   fzss&lat = lat1
   do n=1, 12
     fzss(:,:,:,n) = dim_avg(fz(:,:,:,ispan(n-1,nt-1,12)))
   end do
   fzss(:,:,:,0 ) = fzss(:,:,:,12)
   fzss(:,:,:,13) = fzss(:,:,:,1 )

;--------------------------------------- process 0
 um = tofloat( dim_avg(u(:,{-lat0avg:lat0avg})) )
 fzm = dim_avg_n(fz(:,:,{-lat0avg:lat0avg},:),2)
 fzssm = dim_avg_n(fzss(:,:,{-lat0avg:lat0avg},:),2)

print(min(um)+"  "+max(um))

; ph0 = fspan(-11,3,8)  ; -12 ~ 4  ; U80
 ph0 = fspan(-15.5,5.5,8)  ; -17 ~ 7  ; U70

; ph = new(nt,integer)
; ph = round((um+11.)/2.,3)  ; -12 ~ 4  ; U80
 ph = round((um+15.5)/3.,3)  ; -17 ~ 7  ; U70

 ph = where(ph .lt. 0, 0, ph)
 ph = where(ph .gt. 7, 7, ph)

;--------------------------------------- cn level
 imrk = (/1,1,1,1,1/)*2

;--------------------------------------- contour
 wi = 60
 he = 50
 xf = 15
 yf = 10

 ystr = "-F:B:z:N: [mPa]"
 xstr = "month"

nc = 12
ne = 1
npt = nt/nc
nph = 8

 txtstr = ispan(1,nc,1)
 lcol0 = round(fspan(2,195,nph),3)
 lcol = lcol0
 lcol(1:) = lcol0(0:dimsizes(lcol0)-2)
 lcol(0) = lcol0(dimsizes(lcol0)-1)

 lcol2 = new((/nc*ne,npt/),integer)
 dumm = new(nt,graphic)

print(max(fzm(0,0,:))+"  "+min(fzm(0,0,:)))

do ic=0, nc-1

 do ie=0, ne-1

   if (ie+ic .eq. 0) then
     plot = xylineplot(wks,(-ne/2)*0.15,-fzm(0,0,ispan(ic,nt-1,nc)),"D",0,"T",1,"C",-1)
   else
     xyaddline(plot,ic+(ie-ne/2)*0.15,-fzm(0,0,ispan(ic,nt-1,nc)),"D",0,"T",1,"C",-1)
   end if
   lcol2(ic*ne+ie,:) = lcol(ph(ispan(ic,nt-1,nc)))

   polyres = True
   polyres@gsMarkerIndex      = imrk
   polyres@gsMarkerSizeF      = 4.
   polyres@gsMarkerThicknessF = 0.5
   do n=0, npt-1
     polyres@gsMarkerColor = lcol(ph(ic+n*nc))
     dumm(ic+n*nc) = gsn_add_polymarker(wks,plot,ic,-fzm(0,0,ic+n*nc),polyres)
   end do

 end do

end do

;  xyaddmarker(plot,imrk,"S",0.8,"T",5.0,"C",lcol2)

  do ic=0, nc-2
    xyaddline(plot,ic+0.5,(/-999,999/),"D",0,"T",0.5,"C",colind)
  end do
;  xyaddline(plot,ispan(-1,nc,1),0,"D",0,"T",0.5,"C",1)

  xyaddline(plot,ispan(-1,12,1),-fzssm(0,0,0:13),"D",0,"T",0.5,"C",1)

 set_size(plot,xf,yf,wi,he)
 axis_range(plot,-0.5,nc-0.5,0,1.,False,False)
 axis_str(plot,"",xstr,ystr)
 setvalues plot
  "tmXBMode"   : "Explicit"
  "tmXBValues" : ispan(0,nc-1,1)
  "tmXBLabels" : txtstr
;  "tmXBLabelAngleF" : 30
  "tmXBMinorPerMajor" : 0
  "tmXTOn" : False
;  "tmXBMajorLengthF" : -0.002
 end setvalues

 draw(plot)

 lg = lgplot(wks,ph0+"","D",0,"T",0.5,"C",lcol,2)
 set_size_lg(lg,xf,yf+he+2.5,wi*0.2,16)
 setvalues lg
  "lgLabelFontHeightF" : 0.0045
 end setvalues
 draw(lg)

; dr_txt2(wks,txtstr(ic),xf+1,yf(ic)+he-2,22,6,0.,"CenterLeft",1)

 frame(wks)

end



