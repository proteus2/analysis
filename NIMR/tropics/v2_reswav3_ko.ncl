load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = "uj"
 varname = (/"u","rhodz_dt"/)

 yyyymm1 = 197902
 yyyymm2 = 200611
 lat0 = 5.
 z0 = 18.59 ;26.9 ;17.41 ;18.59

 kmax = 20
 omax = 1.0
 ocutp = 0.75
 ocutn = 0.75;0.4

 wks = cr_wks("fig/v2_reswav3_ko-z"+numeric2int(z0,1)+"-"+expname(0),600,0)
 gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 nfig = 2
 nv = dimsizes(varname)

 y1 = yyyymm1/100
 y2 = yyyymm2/100

 nt = (y2-y1+1)*12 - (yyyymm1-y1*100-1) + (yyyymm2-y2*100-12)

 fdir  = "/prime0/kyh/dat/L60CGW-t/fcoef_ko"
 fname = fdir+"/"+y1+"/uan"+expname+".fft_"+varname(0)+"_ko."+y1+".01.nc"
 f1 = addfile(fname,"r")
 no = dimsizes( f1->ome_fr )
 nk = kmax + 2
 iz = minind( abs( f1->ht - z0*1.e3 ) )
 print(f1->ht(iz))
 ny = dimsizes( f1->latitude({-lat0:lat0}) )

 kwn = fspan(-nk+1,nk-1,nk*2-1)
 ome = fspan(0,no/2,no/2+1) / 90.
 dome = 1./90.
 kr = kwn(nk-1:)/(2.*3.141592*6371229.)

 iocutn = numeric2int(ocutn*90.,1)
 iocutp = no - numeric2int(ocutp*90.,1)

 ub = new((/nt,ny/),float)
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
     fnameu = "$DATD/L60CGW-t/xm/xm.uan"+expname+"a.pj_"+ \
              iy+sprinti("%2.2i",im)+".nc"
     fu = addfile(fnameu,"r")
     ub(n,:) = fu->u(0,{z0*1.e3},{-lat0:lat0})
     n = n + 1
   end do
 end do

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
     ub0 = avg(ub(n,:))
     sht = numeric2int(ub0*abs(kr(:))*86400./dome,1)
     isht = abs(sht)
     i1 = nk
     do i=1, nk-1
       if (isht(i) .gt. 0) then
         i1 = i
         break
       end if
     end do
     do iv=0, nv-1
       fname = fdir+"/"+iy+"/uan"+expname+".fft_"+varname(iv)+"_ko."+ \
               iy+"."+sprinti("%2.2i",im)+".nc"

       f1 = addfile(fname,"r")

       tmp1r = f1->$("fcr_"+varname(iv))$(iz,{-lat0:lat0},:,:kmax+1)
       tmp1i = f1->$("fci_"+varname(iv))$(iz,{-lat0:lat0},:,:kmax+1)

       tmp2r = tmp1r({lat0:-lat0},:,:)
       tmp2i = tmp1i({lat0:-lat0},:,:)

       tmpsr = 0.5*(tmp1r + tmp2r)
       tmpsi = 0.5*(tmp1i + tmp2i)
       tmpar = 0.5*(tmp1r - tmp2r)
       tmpai = 0.5*(tmp1i - tmp2i)

       tmp4 = (/ tmpsr^2 + tmpsi^2, tmpar^2 + tmpai^2 /)  ; (sym,y,o,k)

       if ( varname(iv) .eq. "rhodz_dt" .or. varname(iv) .eq. "theta" ) then
         tmp1r = f1->$("fcr_"+varname(iv))$(iz-1,{-lat0:lat0},:,:kmax+1)
         tmp1i = f1->$("fci_"+varname(iv))$(iz-1,{-lat0:lat0},:,:kmax+1)

         tmp2r = tmp1r({lat0:-lat0},:,:)
         tmp2i = tmp1i({lat0:-lat0},:,:)

         tmpsr = 0.5*(tmp1r + tmp2r)
         tmpsi = 0.5*(tmp1i + tmp2i)
         tmpar = 0.5*(tmp1r - tmp2r)
         tmpai = 0.5*(tmp1i - tmp2i)

         tmp4l = (/ tmpsr^2 + tmpsi^2, tmpar^2 + tmpai^2 /)  ; (sym,y,o,k)

         tmp4 = 0.5*(tmp4 + tmp4l)
       end if

       tmp4(:,:,iocutn:iocutp,:) = 0.

       tmp4b = tmp4(:,:,:,:)

       if (sht(nk-1) .gt. 0) then
         do i=i1, nk-1
           tmp4(:,:,isht(i):no-1,i) = tmp4b(:,:,0:no-1-isht(i),i)
           tmp4(:,:,isht(i):no-1,i) = tmp4b(:,:,0:no-1-isht(i),i)
           tmp4(:,:,0:isht(i)-1,i) = tmp4b(:,:,no-isht(i):no-1,i)
           tmp4(:,:,0:isht(i)-1,i) = tmp4b(:,:,no-isht(i):no-1,i)
         end do
       end if
       if (sht(nk-1) .lt. 0) then
         do i=i1, nk-1
           tmp4(:,:,0:no-1-isht(i),i) = tmp4b(:,:,isht(i):no-1,i)
           tmp4(:,:,0:no-1-isht(i),i) = tmp4b(:,:,isht(i):no-1,i)
           tmp4(:,:,no-isht(i):no-1,i) = tmp4b(:,:,0:isht(i)-1,i)
           tmp4(:,:,no-isht(i):no-1,i) = tmp4b(:,:,0:isht(i)-1,i)
         end do
       end if

       tmp = dim_avg_n( tmp4, 1 )  ; lat. avg.
       tmp(:,:,:) = tmp(:,:,:)*2.

       vspec(:,iv,:,:) = vspec(:,iv,:,:) + tmp(:,:,:)
     end do
     n = n + 1
   end do
 end do
 nt = n
 vspec(:,:,:,:) = vspec(:,:,:,:)/nt

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
   if ( varname(iv) .eq. "rhov" .or. varname(iv) .eq. "rhodz_dt" ) then
     vspec(:,iv,:,:) = vspec(:,iv,:,:) / rho0^2
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

;--------------------------------------- cn level
 vars = ":F14:"+(/"u'",":F8:(r:F10:w:F8:):F14:'/:F8:r:F13B:0:N:"/)
 facts = (/-1,-7/)
 cnlevs = (/ (/0.1,0.2,0.5,1,2,5,10,20,50,100,200/), \
             (/0.1,0.2,0.5,1,2,5,10,20,50,100,200/) /)

 cnfil = (/5,6,7,8,10,12,14,15,16,17,19,21/)

 var   = new(nv,string )
 fact  = new(nv,integer)
 cnlev = new((/nv,dimsizes(cnlevs(0,:))/),float)
 lbtit = new(nv,string )
 do iv=0, nv-1
   if ( varname(iv) .eq. "u" ) then
     var  (iv) = vars (0)
     fact (iv) = facts(0)
	 cnlev(iv,:) = cnlevs(0,:)
     lbtit(iv) = "PSD:B: "+var(iv)+":F22N:  [ 10:S:"+fact(iv)+":N: m:S:2:N: s:S:-2:N: / cyc day:S:-1:N: ]"
   end if
   if ( varname(iv) .eq. "rhodz_dt" ) then
     var  (iv) = vars (1)
     fact (iv) = facts(1)
     cnlev(iv,:) = cnlevs(1,:)
     lbtit(iv) = "PSD:B: "+var(iv)+":F22N:  [ 10:S:"+fact(iv)+":N: m:S:2:N: s:S:-2:N: / cyc day:S:-1:N: ]"
   end if
 end do

;--------------------------------------- ctr
 wi = 32
 he = 32
 xf = 12 + (/0,1,2,3/)*(wi+8)
 yf = 10 + (/0,1/)*(/he+8/)

 mstr  = new((/nfig,nv/),string)
 xstr  = new(nfig,string)
 ystr  = new(nv,string)

 mstr = ""
 mstr(0,:) = vars(:)
 xstr = ""
 xstr(nfig-1) = "zonal wavenumber"
 ystr = ""
 ystr(0) = "intrinsic freq. [cyc day:S:-1:N:]"

do ip=0, nfig-1
do iv=0, nv-1

 ctr1 = cnshadeplot(wks,kwn,ome,vsp(ip,iv,:,:)/(10.^fact(iv)),cnlev(iv,:),cnfil,"")
; cnaddline(ctr1,wks,kwn,ome,vsp(ip,iv,:,:),cnlev(iv,:),0,"D",0,"T",1,"C",1)
 setvalues ctr1
  "cnFillMode" : "RasterFill"
 end setvalues

 set_size(ctr1,xf(iv),yf(ip),wi,he)
 axis_range(ctr1,-kmax-0.5,kmax+0.5,0,omax+0.5*dome,False,False)
 axis_str(ctr1,mstr(ip,iv),xstr(ip),ystr(iv))

 cnfinalize(ctr1,wks)

 draw(ctr1)

end do
end do

 do iv=0, nv-1
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
   set_size_lb(label,xf(max((/0,nv-2/))),yf(nfig-1)+he+16+iv*8,wi*2+8,3)
   draw(label)
 end do

 frame(wks)

end

