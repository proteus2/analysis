load "$UTIL/ncl/header.ncl"

xspace = 83   ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 datname    = (/"ERA-I"    ,"MERRA2"   ,"JRA-55"  ,"JRA-55C"  /)
 expname    = (/"era-int_f","merra2_f" ,"jra55_f" ,"jra55c_f" /)
 dirname    = (/"ERA-I-nr" ,"MERRA2-nr","JRA55-nr","JRA55C-nr"/)
 varname_u  = (/"u"        ,"U"        ,"U_GRD_GDS4_HYBL","U_GRD_GDS4_HYBL"/)
 fvarname_u = (/"u"        ,"u0"       ,"ugrd"    ,"ugrd"     /)
; iexps = (/0,1,2/)
; icol = (/1,14,2,21/)
 iexps = (/2,3/)
 icol = (/2,1/)
; iexps = 0

 lat0  = 15.
 p0    = 100.  ; 54, 66, 80, 96 - ERA
 p0u0  = 23    ; 22, 23, 24, 25 - ERA

 y0 = (/1979,1980,1979,1979/)
 y9 = 2016
 y1 = 1980
 y2 = 2016

 if ( dimsizes(iexps) .eq. 1 ) then
   f_exp = expname(iexps)
 else
   f_exp = ""+dimsizes(iexps)
 end if

 wks = cr_wks("fig/epfs_t2-p"+numeric2int(p0,1)+"-"+f_exp,600,0)
; gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
 colind = NhlNewColor(wks,0.7,0.7,0.7)

;--------------------------------------- read data
 ne = dimsizes(iexps)
 nt = (y2 - y1 + 1)*12
 t = y1 + ( ispan(0,nt-1,1)+0.5 )/12.

 fnam0 = "$DATD/clim/"+dirname+"/"+expname+"."+fvarname_u+"_ypt."+y0+"-"+y9+".nc"
 fnam1 = "$DATD/tem/"+dirname+"/epf/"+expname+".epfsa_yz_recon2.nc"
 fnam2 = "$DATD/tem/"+dirname+"/epf/"+expname+".epf_yz_recon0.nc"
print(fnam1+"")

 u0 = new((/ne,nt/),float,1.e20)

 nw = 8
 rw = new((/ne,nw,nt/),float,1.e20)

 do ie=0, ne-1

   iexp = iexps(ie)

if (iexp .ne. 3) then
   f0 = addfile(fnam0(iexp),"r")
end if
   f1 = addfile(fnam1(iexp),"r")
;y   f2 = addfile(fnam2(iexp),"r")
   print( f1->p({p0})+"" )

   it1u = (y1 - y0(iexp))*12
if (iexp .ne. 3) then
   it9u = min( (/ dimsizes( f0->$varname_u(iexp)$(:,0,0) )-1, it1u+nt-1 /) )
end if
   ntu = it9u - it1u + 1

   i1mis = 1  ; The 1st is Feb.
   it1r = it1u - i1mis
   it9r = min( (/ dimsizes( f1->f_z(:,0,0,0) )-1, it1r+nt-1 /) )
   ntr = it9r - it1r + 1
   n1 = 0
   n2 = ntr-1
   if (it1r .lt. 0) then
     n1 = -it1r
     ntr = ntr + it1r
     it1r = 0
   end if
 
;   u0(ie,0:ntu-1) = dim_avg(tofloat(  \
;           f0->$varname_u(iexp)$(it1u:it9u,{p0u0},{-lat0:lat0}) ))
;u0 = conform_dims(dimsizes(u0),u0(ie,:),1)

   rw(ie,:,n1:n2) = transpose(dim_avg( f1->f_z(it1r:it9r,:nw-1,{p0},{-lat0:lat0}) ))
;y   rw(ie,nw-2:,n1:n2) = rw(ie,nw-2:,n1:n2) +  \
;y           transpose(dim_avg( f2->f_z(it1r:it9r,8:9,{p0},{-lat0:lat0}) ))
;y   rw(ie,nw-1,n1:n2) = rw(ie,nw-1,n1:n2) +  \
;y           dim_avg( f2->f_z(it1r:it9r,2,{p0},{-lat0:lat0}) +  \
;y                    f2->f_z(it1r:it9r,3,{p0},{-lat0:lat0}) )

   rw(ie,3,n1:n2) = rw(ie,3,n1:n2) +  \
           dim_avg( f1->f_z(it1r:it9r,13,{p0},{-lat0:lat0}) )
   rw(ie,6:7,n1:n2) = rw(ie,6:7,n1:n2) +  \
           transpose( dim_avg( f1->f_z(it1r:it9r,16:17,{p0},{-lat0:lat0}) ))

   rw(ie,:,:) = rw(ie,:,:) / 6371.e3 * 1.e3

 end do

 fz = rw

 fzss = new((/ne,nw,14/),float,1.e20)
 do n=1, 12
   fzss(:,:,n) = dim_avg(fz(:,:,ispan(n-1,nt-1,12)))
 end do
 fzss(:,:,0 ) = fzss(:,:,12)
 fzss(:,:,13) = fzss(:,:,1 )

;--------------------------------------- process 0
 fzds = fz
 fzssp = fz
 do n=1, 12
   nn = ispan(n-1,nt-1,12)
   fzds(:,:,nn) = fz(:,:,nn) -  \
       conform_dims((/ne,nw,dimsizes(nn)/),fzss(:,:,n),(/0,1/))
   fzssp(:,:,nn) = conform_dims((/ne,nw,dimsizes(nn)/),fzss(:,:,n),(/0,1/))
   delete(nn)
 end do

; ph0 = fspan(-11,3,8)  ; -12 ~ 4  ; U80
 ph0 = fspan(-15.5,5.5,8)  ; -17 ~ 7  ; U70

; ph = new(nt,integer)
; ph = round((u0+11.)/2.,3)  ; -12 ~ 4  ; U80
;; ph = round((u0+15.5)/3.,3)  ; -17 ~ 7  ; U70

;; ph = where(ph .lt. 0, 0, ph)
;; ph = where(ph .gt. 7, 7, ph)

;--------------------------------------- cn level

;--------------------------------------- contour
 nl = 5

 wi = 150
 he = 25
 xf = 15
 yf = 10 + ispan(0,nl-1,1)*(he+6)

 ig = (/1,3,3,2,1/)
 id = (/0,0,1,1,1/)
 iw = ig*2+id
 fsgn = where(id .eq. 0, -1, 1)

 ystr = "F:B:z:N: [mPa]"
 xstr = where(yf .eq. max(yf), "year", "")

do il=0, nl-1

 plot = xylineplot(wks,t,fz(:,iw(il),:),"D",0,"T",1,"C",icol)
 xyaddline(plot,t,conform_dims(dimsizes(fz(:,0,:)),dim_avg(fz(:,iw(il),:)),0),"D",0,"T",0.5,"C",icol)

 set_size(plot,xf,yf(il),wi,he)
 axis_range(plot,y1,2017,-999*fsgn(il),999*fsgn(il),False,False)
 axis_str(plot,"",xstr(il),ystr)
 setvalues plot
  "tmXBMaxTicks" : 8
 end setvalues

 draw(plot)

end do

 lg = lgplot(wks,datname(iexps),"D",0,"T",0.5,"C",icol,2)
 set_size_lg(lg,xf,max(yf)+he+3.0,wi*0.2,15)
; setvalues lg
;  "lgLabelFontHeightF" : 0.0045
; end setvalues
 draw(lg)


 frame(wks)

end

