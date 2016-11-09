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

 lat0 = 25.
 lat0avg = 10
 p0    = 96.  ; 54, 66, 80, 96
 p0u0  = 23   ; 22, 23, 24, 25
 p0u = (/23,24,25/)
 ymax0 = 0.5

 x_wmrgw = 1 ;3

 wks = cr_wks("fig/epf-kw_t-ush-p"+numeric2int(p0,1)+"-"+expname(0),600,0)
; gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
 colind = NhlNewColor(wks,0.7,0.7,0.7)

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
   nt = 1 + 407
   nw = 8

   u0 = f0->u(:nt-1,{p0u0},{-lat0:lat0})
   u = f0->u(:nt-1,{p0u},{-lat0:lat0})
   u = u - f0->u(:nt-1,{p0u+1},{-lat0:lat0})

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

   t = 1979 + ( ispan(0,nt-1,1)+0.5 )/12.

   fzss = new((/3,2,ny,14/),float,1.e20)
   fzss!2 = "lat"
   fzss&lat = lat1
   do n=1, 12
     fzss(:,:,:,n) = dim_avg(fz(:,:,:,ispan(n-1,nt-1,12)))
   end do
   fzss(:,:,:,0 ) = fzss(:,:,:,12)
   fzss(:,:,:,13) = fzss(:,:,:,1 )

;--------------------------------------- process 0
 u0m = tofloat( dim_avg(u0(:,{-lat0avg:lat0avg})) )
 um = tofloat( dim_avg_n(u(:,:,{-lat0avg:lat0avg}),2) )
 fzm = dim_avg_n(fz(:,:,{-lat0avg:lat0avg},:),2)
 fzssm = dim_avg_n(fzss(:,:,{-lat0avg:lat0avg},:),2)

 fzdsm = fzm
 fzssmp = fzm
 do n=1, 12
   nn = ispan(n-1,nt-1,12)
   fzdsm(:,:,nn) = fzm(:,:,nn) -  \
       conform_dims((/3,2,dimsizes(nn)/),fzssm(:,:,n),(/0,1/))
   fzssmp(:,:,nn) = conform_dims((/3,2,dimsizes(nn)/),fzssm(:,:,n),(/0,1/))
   delete(nn)
 end do

; ph0 = fspan(-11,3,8)  ; -12 ~ 4  ; U80
 ph0 = fspan(-15.5,5.5,8)  ; -17 ~ 7  ; U70

; ph = new(nt,integer)
; ph = round((u0m+11.)/2.,3)  ; -12 ~ 4  ; U80
 ph = round((u0m+15.5)/3.,3)  ; -17 ~ 7  ; U70

 ph = where(ph .lt. 0, 0, ph)
 ph = where(ph .gt. 7, 7, ph)

 ncor = 12
 ccorr = new(ncor*2+1,float)
 do k=0, dimsizes(p0u)-1
   ccorr(ncor:0) = esccr(-fzdsm(0,0,1:),um(1:,k),ncor)
   ccorr(ncor:) = esccr(um(1:,k),-fzdsm(0,0,1:),ncor)
   i_maxcc = ind(max(abs(ccorr)) .eq. abs(ccorr))
print( ccorr(i_maxcc)+"  "+(i_maxcc-12))
 end do

;--------------------------------------- cn level
 imrk = (/1,1,1,1,1/)*2
 icol = (/2,3,4,5/)

;--------------------------------------- contour
 wi = 160
 he = 40
 xf = 15
 yf = 10 + ispan(0,3,1)*(he+10)

 ystr = "[mPa]"
 xstr = "year"

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

do ie=0, ne-1

 plot = xylineplot(wks,t,-fzm(0,0,:),"D",0,"T",1,"C",1)
 xyaddline(plot,t,avg(-fzm(0,0,:)),"D",0,"T",1,"C",colind)

 set_size(plot,xf,yf(0),wi,he)
 axis_range(plot,1979,2013,-999,999,False,False)
 axis_str(plot,"",xstr,ystr)
 setvalues plot
  "tmXBMaxTicks" : 8
 end setvalues

 draw(plot)

 plot = xylineplot(wks,t,-fzdsm(0,0,:),"D",0,"T",1,"C",1)
 xyaddline(plot,t,-fzssmp(0,0,:)-avg(-fzssmp(0,0,:)),"D",0,"T",0.8,"C",2)

 set_size(plot,xf,yf(1),wi,he)
 axis_range(plot,1979,2013,-999,999,False,False)
 axis_str(plot,"",xstr,ystr)
 setvalues plot
  "tmXBMaxTicks" : 8
 end setvalues

 draw(plot)

 plot = xylineplot(wks,t,-fzdsm(0,0,:),"D",0,"T",1,"C",1)
 do il=0, dimsizes(p0u)-1
   xyaddline(plot,t,(um(:,il)-avg(um(:,il)))*0.05,"D",0,"T",0.8,"C",icol(il))
 end do

 set_size(plot,xf,yf(2),wi,he)
 axis_range(plot,1979,2013,-999,999,False,False)
 axis_str(plot,"",xstr,ystr)
 setvalues plot
  "tmXBMaxTicks" : 8
 end setvalues

 draw(plot)

end do

 frame(wks)

end

