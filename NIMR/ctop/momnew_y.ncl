load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = (/"tk","tl","tm","ui"/)
 ppid    = (/"e" ,"e" ,"e" ,"g" /)
 varname = (/"field420","field422"/)
 varnamec= (/"unspecified_1","unspecified_2"/)

 y1 = 1952
 y2 = 1957
 kz = (/28,35,50/)

 wks = cr_wks("fig/momann_y-"+expname(0),600,0)
; gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read base data
 ne = dimsizes(expname)
 nt0 = (y2-y1+1)*12

do ie=0, ne-1

 fname = new(nt0,string)
 fnamec= new(nt0,string)
 fdir  = "../dat/xm"
 n = 0
 do iy=y1, y2
 do im=1, 12
   umtail5 = systemfunc("ex_date_umstr.sh "+iy+" "+im)
   fname(n) = fdir+"/"+expname(ie)+"/mfb_xm.uan"+expname(ie)+"a.p"+ppid(ie)+umtail5+".nc"
   fnamec(n) = fdir+"/"+expname(ie)+"/mfc_xm.uan"+expname(ie)+"a.p"+ppid(ie)+umtail5+".nc"
   n = n + 1
 end do
 end do
 if (expname(ie) .eq. "ui") then
   fexist = isfilepresent(fname)
 else
   fexist = isfilepresent(fname) .and. isfilepresent(fnamec)
 end if

 n1 = -1
 n2 = nt0-1
 do n=0, nt0-1
   if ( fexist(n) ) then
     if (n1 .eq. -1) then
       n1 = n
     end if
   else
     if ( n1 .ne. -1 .and. n2 .eq. nt0-1 ) then
       n2 = n - 1
     end if
   end if
 end do
 if (n1 .eq. -1) then
   n2 = -1
 end if

 ; get y, z
 if (ie .eq. 0) then
   f0 = addfile(fname(n1),"r")
   tmp = f0->$varname(ie)$
   lat = tmp&$tmp!2$
   z   = tmp&$tmp!1$(kz) / 1.e3
   ny = dimsizes(lat)
   nz = dimsizes(z  )
   delete(tmp)
 end if

 if (ie .eq. 0) then
   mf = new((/ne*2-1,2,nz,ny/),float,1.e20)
 end if

 if (n1 .ne. -1) then
print(expname(ie)+"  "+((n2-n1+1)/12.))
   f1 = addfiles(fname(n1:n2),"r")
   ListSetType(f1,"cat")
   mfs = (/ dim_avg_n( f1[:]->$varname(0)$(:,kz,:), 0 ), \
            dim_avg_n( f1[:]->$varname(1)$(:,kz,:), 0 ) /)
   if (expname(ie) .ne. "ui") then
     f1c= addfiles(fnamec(n1:n2),"r")
     ListSetType(f1c,"cat")
     mfc = (/ dim_avg_n( f1c[:]->$varnamec(0)$(:,kz,:), 0 ), \
              dim_avg_n( f1c[:]->$varnamec(1)$(:,kz,:), 0 )*(-1.) /)
     mfs = mfs + mfc
     mf(ie*2:ie*2+1,:,:,:) = (/ mfc, mfs /)
   else
     mf(ie*2,:,:,:) = mfs
   end if
 end if

end do

;--------------------------------------- process 0
 mf = mf*1.e3

;--------------------------------------- plot

 xf = 20 + (/0,1,2/)*70
 yf = 10 + (/0,1/)*42
 wi = 60
 he = 30
 mstr = (/ sprintf("%2.0f",z)+" km", (/"","",""/) /)
 ystr = (/ (/"E","W"/)+"-MF (mPa)", (/"",""/), (/"",""/) /)
 llbl = (/"CGW","CGW+BGW","BGW (CTL)"/)
 ymax = (/2,2,0.8/)

; ldsh = (/0,0, 1/)
; lcol = (/2,1, 1/)
 ldsh = (/0,0,0,0,0,0, 1/)
 lcol = (/1,1,2,2,4,4, 1/)

 do k=0, dimsizes(kz)-1
 do idr=0, 1

   plot = xylineplot(wks,lat,mf(:,idr,k,:),"D",ldsh,"T",1.0,"C",lcol)

   set_size(plot,xf(k),yf(idr),wi,he)
   axis_range(plot,-90,90,0,ymax(k),False,False)
   axis_str(plot,mstr(idr,k),"lat (deg)",ystr(k,idr))

   draw(plot)

 end do
 end do

; lg = lgplot(wks,llbl,"D",ldsh,"T",1.0,"C",lcol,2)
 lg = lgplot(wks,expname,"D",ldsh((/0,2,4,6/)),"T",1.0,"C",lcol((/0,2,4,6/)),2)
 set_size_lg(lg,xf,yf(0),wi*0.25,he*0.4)
 draw(lg)


 frame(wks)

end

