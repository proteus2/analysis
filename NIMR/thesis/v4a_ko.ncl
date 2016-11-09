load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = "uj"
 varname = (/"u","rhov","rhodz_dt","theta"/)

 yyyymm1 = 197901
 yyyymm2 = 200611
 lat0 = 15.  ; 15. ; 5.
 z0 = 18.59  ; 17.41 ; 18.59

 kmax = 15;20
 omax = 0.75;1.0

 beta = 2.*(7.292116e-5)/6371229.
 nbv = 2.5e-2

 wks = cr_wks("fig/v4a_ko-z"+numeric2int(z0,1)+"-y"+numeric2int(lat0,1)+"-"+expname(0),600,0)
 gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 nfig = 2
 nv = dimsizes(varname)

if ( yyyymm1 .eq. 195301 .and. yyyymm2 .eq. 200611 ) then

 fdir  = "/prime0/kyh/dat/L60CGW-t/fcoef_ko/avg"
 fname = fdir+"/uan"+expname+".fc2_"+varname(0)+"_wave_ko.nc"
 if ( varname(0) .eq. "rhodz_dt" ) then
   fname = fdir+"/uan"+expname+".fc2_"+varname(0)+"-wc_wave_ko.nc"
 end if
 f1 = addfile(fname,"r")
 no = dimsizes( f1->ome_fr )
 nk = kmax + 1
 iz = minind( abs( f1->z - z0*1.e3 ) )
 print(f1->z(iz))

 kwn = fspan(-kmax,kmax,kmax*2+1)
 ome = fspan(0,no/2,no/2+1) / 90.
 dome = 1./90.

 vspec = new((/2,nv,no,nk/),float)
 do iv=0, nv-1
   fname = fdir+"/uan"+expname+".fc2_"+varname(iv)+"_wave_ko.nc"
   if ( varname(iv) .eq. "rhodz_dt" ) then
     fname = fdir+"/uan"+expname+".fc2_"+varname(iv)+"-wc_wave_ko.nc"
   end if
   f1 = addfile(fname,"r")

   tmps = f1->$("fc2_"+varname(iv)+"_s")$(iz,{-lat0:lat0},:,:kmax)
   tmpa = f1->$("fc2_"+varname(iv)+"_a")$(iz,{-lat0:lat0},:,:kmax)
   tmp = dim_avg_n( (/ tmps, tmpa /), 1 )*2.
   if ( varname(iv) .eq. "theta" ) then
     tmpsl = f1->$("fc2_"+varname(iv)+"_s")$(iz-1,{-lat0:lat0},:,:kmax)
     tmpal = f1->$("fc2_"+varname(iv)+"_a")$(iz-1,{-lat0:lat0},:,:kmax)
     tmpl = dim_avg_n( (/ tmpsl, tmpal /), 1 )*2.
     tmp = 0.5*(tmp + tmpl)
   end if
   vspec(:,iv,:,:) = tmp(:,:,:)
 end do

else

 y1 = yyyymm1/100
 y2 = yyyymm2/100

 fdir  = "/prime0/kyh/dat/L60CGW-t/fcoef_ko"
 fname = fdir+"/"+y1+"/uan"+expname+".fft_"+varname(0)+"_ko."+y1+".01.nc"
 if ( varname(0) .eq. "rhodz_dt" ) then
   fname = fdir+"/"+y1+"/uan"+expname+".fft_"+varname(0)+"-wc_ko."+y1+".01.nc"
 end if
 f1 = addfile(fname,"r")
 no = dimsizes( f1->ome_fr )
 nk = kmax + 1
 iz = minind( abs( f1->ht - z0*1.e3 ) )
 print(f1->ht(iz))

 kwn = fspan(-nk+1,nk-1,nk*2-1)
 ome = fspan(0,no/2,no/2+1) / 90.
 dome = 1./90.

 vspec = new((/2,nv,no,nk/),float)
 vspec = 0.
 n = 0
 do iy=y1, y2
   im1 = 1
   im2 = 12
   if ( iy .eq. y1 ) then
     im1 = yyyymm1 - y1*100
   end if
   if ( iy .eq. y2 ) then
     im2 = yyyymm2 - y2*100
   end if
   do im=im1, im2
     do iv=0, nv-1
       fname = fdir+"/"+iy+"/uan"+expname+".fft_"+varname(iv)+"_ko."+ \
               iy+"."+sprinti("%2.2i",im)+".nc"
       if ( varname(iv) .eq. "rhodz_dt" ) then
         fname = fdir+"/"+iy+"/uan"+expname+".fft_"+varname(iv)+"-wc_ko."+ \
               iy+"."+sprinti("%2.2i",im)+".nc"
       end if

       f1 = addfile(fname,"r")

       tmp1r = f1->$("fcr_"+varname(iv))$(iz,{-lat0:lat0},:,:kmax)
       tmp2r = tmp1r({lat0:-lat0},:,:)
       tmp1i = f1->$("fci_"+varname(iv))$(iz,{-lat0:lat0},:,:kmax)
       tmp2i = tmp1i({lat0:-lat0},:,:)

       tmp = dim_avg_n( \
             (/ (0.5*(tmp1r + tmp2r))^2 + (0.5*(tmp1i + tmp2i))^2, \
                (0.5*(tmp1r - tmp2r))^2 + (0.5*(tmp1i - tmp2i))^2 /), 1 )
       tmp(:,:,:) = tmp(:,:,:)*2.

       if ( varname(iv) .eq. "theta" ) then
         tmp1r = f1->$("fcr_"+varname(iv))$(iz-1,{-lat0:lat0},:,:kmax)
         tmp2r = tmp1r({lat0:-lat0},:,:)
         tmp1i = f1->$("fci_"+varname(iv))$(iz-1,{-lat0:lat0},:,:kmax)
         tmp2i = tmp1i({lat0:-lat0},:,:)

         tmpl = dim_avg_n( \
                (/ (0.5*(tmp1r + tmp2r))^2 + (0.5*(tmp1i + tmp2i))^2, \
                   (0.5*(tmp1r - tmp2r))^2 + (0.5*(tmp1i - tmp2i))^2 /), 1 )
         tmpl(:,:,:) = tmpl(:,:,:)*2.

         tmp = 0.5*(tmp + tmpl)
       end if

       vspec(:,iv,:,:) = vspec(:,iv,:,:) + tmp(:,:,:)
     end do
     n = n + 1
   end do
 end do
 nt = n
 vspec(:,:,:,:) = vspec(:,:,:,:)/nt

end if

 ; switch for v'
 do iv=0, nv-1
   if (varname(iv) .eq. "rhov") then
     vstmp = vspec(0,iv,:,:)
     vspec(0,iv,:,:) = vspec(1,iv,:,:)
     vspec(1,iv,:,:) = vstmp
     delete( vstmp )
   end if
 end do

 fname2 = "$DATD/L60CGW-t/cmxm.uan"+expname+"a.pj.nc"
 f2 = addfile(fname2,"r")

 rho0 = dim_avg(f2->rho({z0*1.e3},{-lat0:lat0}))

 k0 = minind(abs( f2->$( f2->theta!0 )$ - z0*1.e3 ))
 k1 = k0 + 1
 if ( f2->$( f2->theta!0 )$(k0) .gt. z0*1.e3 ) then
   k1 = k0 - 1
 end if
 th0 = dim_avg(f2->theta(k0,{-lat0:lat0}))
 th1 = dim_avg(f2->theta(k1,{-lat0:lat0}))
 thr0 = exp(0.5*(log(th0) + log(th1)))

 do iv=0, nv-1
   if ( varname(iv) .eq. "rhov" ) then
     vspec(:,iv,:,:) = vspec(:,iv,:,:) / rho0^2
   end if
   if ( varname(iv) .eq. "rhodz_dt" ) then
     vspec(:,iv,:,:) = vspec(:,iv,:,:) / rho0^2 * 1.e6
   end if
   if ( varname(iv) .eq. "theta" ) then
     vspec(:,iv,:,:) = vspec(:,iv,:,:) / thr0^2
   end if
 end do

;--------------------------------------- process 0
 vsp = new((/nfig,nv,no/2+1,nk*2-1/),float)
 vsp = 0.
 do i=0, nk-2
   vsp(:,:,:,i) = vspec(:,:,0:no/2,nk-1-i)
 end do
 vsp(:,:,0,nk:) = vspec(:,:,0,1:);*0.5
 do n=1, no/2
   vsp(:,:,n,nk:) = vspec(:,:,no-n,1:)
 end do
 vsp(:,:,no/2,:) = vsp(:,:,no/2,:)*0.5

 vsp = vsp/dome  ; PS density

 vsp(:,:,:,nk-1) = vsp@_FillValue  ; k = 0
print(max(vsp))

 nk10 = nk*2*10
 no10 = no/2*10

 nw = 3
 dispc = new((/2,nw,no10+1,nk10+1/),float)
 dispc = 1.e3

 kwn2 = fspan(-nk,nk,nk10+1)
 ome2 = fspan(0,no/2,no10+1) / 90.

 two_pi = 6.283185
 k = kwn2/6371229.
 o = ome2*two_pi/86400.
 okbi = conform_dims((/no10+1,nk10+1/),o,0)* \
        conform_dims((/no10+1,nk10+1/),k,1)/beta

 ; K, IG / MRG, IG
 do n=0, no10
   dispc(0,0,n,nk10/2+1:) = two_pi/1.e3*o(n)/(nbv*k(nk10/2+1:))
   dispc(0,1,n,:) = two_pi/1.e3*o(n)^2/  \
        (nbv*beta*(1.5+sqrt(1.5*1.5+okbi(n,:)*(1.+okbi(n,:)))))
   dispc(1,0,n,:) = two_pi/1.e3*o(n)^2/(nbv*beta*(1.+okbi(n,:)))
   dispc(1,1,n,:) = two_pi/1.e3*o(n)^2/  \
        (nbv*beta*(2.5+sqrt(2.5*2.5+okbi(n,:)*(1.+okbi(n,:)))))
 end do
 dispc(0,0,:,nk10/2) = 1.e3
 dispc(1,0,:,:) = where(okbi .gt. -1., dispc(1,0,:,:), 1.e30)
 ; R
 do n=1, no10
   dispc(0,2,n,:nk10/2-1) = two_pi/1.e3*o(n)^2/  \
        (nbv*beta*(1.5-sqrt(1.5*1.5+okbi(n,:nk10/2-1)*(1.+okbi(n,:nk10/2-1)))))
   dispc(1,2,n,:nk10/2-1) = two_pi/1.e3*o(n)^2/  \
        (nbv*beta*(2.5-sqrt(2.5*2.5+okbi(n,:nk10/2-1)*(1.+okbi(n,:nk10/2-1)))))
 end do
 dispc(:,2,0,:) = 0.
 dispc(0,2,:,:) = where(okbi .gt. -1., dispc(0,2,:,:), 1.e30)
 dispc(1,2,:,:) = where(okbi .gt. -1., dispc(1,2,:,:), 1.e30)

;--------------------------------------- cn level
 vars = (/":F10:u:F13:'",":F8:(r:F10:v:F8:):F13:'/:F8:r:F13B:0:N:",":F8:(r:F10:w:F8:):F13:'/:F8:r:F13B:0:N:",":F8:q:F13:'/:F8:q:F13::B:0:N:"/)
 facts = (/0,0,0,-5/)
 cnlevs = (/ (/0.02,0.05,0.1,0.2,0.5,1,2,5,10,20/), \
             (/0.02,0.05,0.1,0.2,0.5,1,2,5,10,20/), \
             (/0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10/), \
             (/0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10/) /)

 cnfil = (/5,6,7,10,12,14,15,16,17,19,21/)

 var   = new(nv,string )
 fact  = new(nv,integer)
 cnlev = new((/nv,dimsizes(cnlevs(0,:))/),float)
 lbtit = new(nv,string )
 do iv=0, nv-1
   if ( varname(iv) .eq. "u" ) then
     var  (iv) = vars (0)
     fact (iv) = facts(0)
	 cnlev(iv,:) = cnlevs(0,:)
     lbtit(iv) = "m:S:2:N: s:S:-2:N: / cyc day:S:-1:N: ]"
   end if
   if ( varname(iv) .eq. "rhov" ) then
     var  (iv) = vars (1)
     fact (iv) = facts(1)
     cnlev(iv,:) = cnlevs(1,:)
     lbtit(iv) = "m:S:2:N: s:S:-2:N: / cyc day:S:-1:N: ]"
   end if
   if ( varname(iv) .eq. "rhodz_dt" ) then
     var  (iv) = vars (2)
     fact (iv) = facts(2)
     cnlev(iv,:) = cnlevs(2,:)
     lbtit(iv) = "mm:S:2:N: s:S:-2:N: / cyc day:S:-1:N: ]"
   end if
   if ( varname(iv) .eq. "theta" ) then
     var  (iv) = vars (3)
     fact (iv) = facts(3)
     cnlev(iv,:) = cnlevs(3,:)
     lbtit(iv) = "/ cyc day:S:-1:N: ]"
   end if
   lbt_hd = "[ 10:S:"+fact(iv)+":N: "
   if (fact(iv) .eq. 0) then
     lbt_hd = "[ "
   end if
   lbtit(iv) = lbt_hd+lbtit(iv)
 end do

 lzs = (/ (/ (/2.,5.,15./),    \
             (/2.5,5.,5./),    \
             (/2.,10.,60./) /), \
          (/ (/2.,5.,30./),    \
             (/2.5,5.,5./),    \
             (/2.,10.,60./) /) /)
 ldsh = (/0,16,14/)

;--------------------------------------- ctr
 wi = 32
 he = 32
 xf = 12 + (/0,1,2,3/)*(wi+8)
 yf = 10 + (/0,1/)*(/he+10/)

 mstr  = new((/nfig,nv/),string)
 xstr  = new(nfig,string)
 ystr  = new(nv,string)

 mstr = ""
 mstr(0,:) = var(:)
 mstr(1,:) = var(:)
 xstr = ""
; xstr(nfig-1) = "zonal wavenumber"
; xstr = "zonal wavenumber"
 ystr = ""
 ystr(0) = "frequency [cyc day:S:-1:N:]"

do ip=1, 1 ;0, nfig-1
do iv=0, nv-1

 ctr1 = cnshadeplot(wks,kwn,ome,vsp(ip,iv,:,:)/(10.^fact(iv)),cnlev(iv,:),cnfil,"")
; cnaddline(ctr1,wks,kwn,ome,vsp(ip,iv,:,:),cnlev(iv,:),0,"D",0,"T",1,"C",1)
 do iw=0, nw-1
   cnaddline(ctr1,wks,kwn2,ome2,dispc(ip,iw,:,:),lzs(ip,iw,:),0,"D",ldsh(iw),"T",1.2,"C",1)
 end do
 setvalues ctr1
  "cnFillMode" : "RasterFill"
 end setvalues

 set_size(ctr1,xf(iv),yf(ip),wi,he)
; axis_range(ctr1,-kmax-0.5,kmax+0.5,0,omax+0.5*dome,False,False)
 axis_range(ctr1,-kmax-0.5,kmax+0.5,0,omax,False,False)
 axis_str(ctr1,mstr(ip,iv),xstr(ip),ystr(iv))

 cnfinalize(ctr1,wks)

 draw(ctr1)

end do
end do

 do iv=0, 0 ;nv-1
   label = create "lab" labelBarClass wks
     "lbOrientation" : "Horizontal"
     "lbBoxCount"   : dimsizes(cnfil)
     "lbFillColors" : cnfil
     "lbMonoFillPattern" : True
     "lbLabelsOn" : True
     "lbLabelStrings" : cnlev(iv,:)
     "lbLabelStride" : 1
;     "lbAutoManage" : False
     "lbLabelAlignment" : "InteriorEdges"
     "lbPerimOn" : False
     "lbTitleString" : lbtit(iv)
     "lbTitlePosition" : "Bottom"
   end create
   set_size_lb(label,xf(max((/0,nv-2/))),yf(nfig-1)+he+18+iv*8,wi*2+8,3)
;   draw(label)
 end do

 frame(wks)

end

