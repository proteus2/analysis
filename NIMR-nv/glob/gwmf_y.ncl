load "$UTIL/ncl/header.ncl"

xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.0

begin

;--------------------------------------- parameters
; jobname = (/"ui","uj"/)
; expname = (/"L60CTL","L60CGW-t"/)
 jobname = (/"tk","tc"/)
 expname = (/"AOL60CTL","AOL60CGW"/)
 ppid_bg = "g"
 ppid_cg = "d"
 varname_bg = "Fp_"+(/"e","w"/)+"_gwussp"
 varname_cg = "F_"+(/"e","w"/)+"_gwconv"

 y0 = 1951
 y1 = 1951
 y2 = 2009
 mms = "ann"

 z = (/16,48/)

 wks = cr_wks("fig/gwmf_y-"+mms+"_"+y1+"-"+y2+"-"+jobname(1),600,0)
 gsn_define_colormap(wks,"Default")
; colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 nv = 4
 nyr = y2 - y1 + 1
 nz = dimsizes(z)
 ndr = 2

 mm3s = (/"djf","mam","jja","son"/)
 mm3i = (/ (/1,2,12/), (/3,4,5/), (/6,7,8/), (/9,10,11/) /)
 do is=0, dimsizes(mm3s)-1
   if ( mms .eq. mm3s(is) ) then
     mm = mm3i(is,:)
   end if
 end do
 if (mms .eq. "ann") then
   mm = ispan(1,12,1)
 end if
 nm = dimsizes(mm)

do ie=0, 1

 fname_bg = "../dat/"+expname(ie)+"/mmxm.uan"+jobname(ie)+"a.p"+ppid_bg+".nc"
 fname_cg = "../dat/"+expname(ie)+"/mmxm.uan"+jobname(ie)+"a.p"+ppid_cg+".nc"
 f_bg = addfile(fname_bg,"r")
 tag_cg = False
 if ( isfilepresent(fname_cg) ) then
   tag_cg = True
   f_cg = addfile(fname_cg,"r")
 end if
 n0 = (y1 - y0)*12

 ; get y
 tmpbe = f_bg->$varname_bg(0)$(n0:n0+nyr*12-1,:,:)
 tmpbw = f_bg->$varname_bg(1)$(n0:n0+nyr*12-1,:,:)
 if ( tag_cg ) then
   tmpce = f_cg->$varname_cg(0)$(n0:n0+nyr*12-1,:,:)
   tmpcw = f_cg->$varname_cg(1)$(n0:n0+nyr*12-1,:,:)
 end if
 lat = tmpbe&$tmpbe!2$
 ny = dimsizes(lat)

 if (ie .eq. 0) then
   vara = new((/nv,ndr,nm,nyr,nz,ny/),float)
 end if
 iv = 0
 if ( tag_cg ) then
   iv = 1
 end if
 do im=0, nm-1
   vara(iv,0,im,:,:,:) = tmpbe(ispan(mm(im)-1,nyr*12-1,12),{z*1.e3},:)
   vara(iv,1,im,:,:,:) = tmpbw(ispan(mm(im)-1,nyr*12-1,12),{z*1.e3},:)
 end do
 if ( tag_cg ) then
   do im=0, nm-1
     vara(2,0,im,:,:,:) = (/ tmpce(ispan(mm(im)-1,nyr*12-1,12),{z*1.e3},:) /)
     vara(2,1,im,:,:,:) = (/ tmpcw(ispan(mm(im)-1,nyr*12-1,12),{z*1.e3},:) /)* \
                          (-1.)
   end do
   vara(3,:,:,:,:,:) = vara(1,:,:,:,:,:) + vara(2,:,:,:,:,:)
 end if

 delete(tmpbe)
 delete(tmpbw)

end do  ; ie

;--------------------------------------- process 0
 var = dim_avg_n(dim_avg_n(vara,3),2)*1.e3

;--------------------------------------- cn level
 ldsh = (/14,0,0,0/)
 ltck = (/1.0,1.0,1.0,1.8/)
 lcol = (/1,1,2,1/)

;--------------------------------------- contour
 wi = 32
 he = 20
 xf = 10 + (/0,1,2/)*40
 yf = 10 + ispan(0,nz-1,1)*25

 mstr = new((/nz,ndr/),string,"")
 xstr = new(nz ,string,"")
 ystr = new(ndr,string,"")

 mstr(0,:) = (/"Eastward-MF","Westward-MF"/)
 xstr(nz-1) = "latitude (deg)"
 ystr(0) = "M (mPa)"
 lbopt = (/"(m s:S:-1:N:)","","(m s:S:-1:N:)"/)

do k=0, nz-1
do idr=0, ndr-1

 plot = xylineplot(wks,lat,var(:,idr,k,:),"D",ldsh,"T",ltck,"C",lcol)

 set_size(plot,xf(idr),yf(k),wi,he)
 axis_range(plot,-90,90,-999,999,False,False)
 axis_str(plot,mstr(k,idr),xstr(k),ystr(idr))

 draw(plot)

end do
end do

 frame(wks)

end

