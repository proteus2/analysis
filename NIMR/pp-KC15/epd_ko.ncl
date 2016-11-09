load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = "uj"
 varname = (/"f_z","epd_z","Fz/D"/)

 yyyymm1 = 195301
 yyyymm2 = 200611
 lat1 = 5 ; -20.  ; 15. ; 5.
 z0 = 32 ; 36  ; 17.41 ; 18.59

 kmax = 35;20
 omax = 2.4

 beta = 2.*(7.292116e-5)/6371229.
 nbv = 2.5e-2

 sl = toint(abs(lat1))+"N"
 if (lat1 .lt. 0.) then
   sl = toint(abs(lat1))+"S"
 end if
 if (lat1 .eq. 0.) then
   sl = "EQ"
 end if

 wks = cr_wks("fig/epd_ko-z"+numeric2int(z0,1)+"-"+sl+"-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 nfig = 1
 nv = dimsizes(varname)

if ( False .and. yyyymm1 .eq. 195301 .and. yyyymm2 .eq. 200611 ) then

 fdir  = "/prime0/kyh/dat/L60CGW-t/epf-wc/avg"
 fname = fdir+"/uan"+expname+".epf_koyz.nc"
 f1 = addfile(fname,"r")
 no = dimsizes( f1->ome_fr )
 nk = kmax + 1

 kwn = fspan(-kmax,kmax,kmax*2+1)
 ome = fspan(0,no/2,no/2+1) / 90.
 dome = 1./90.

 vspec = new((/nv,no,nk/),float)
 do iv=0, nv-1
   if ( varname(iv) .eq. "f_vT" .or. varname(iv) .eq. "Fz/D" ) then
     break
   end if
   fname = fdir+"/uan"+expname+".epf_koyz.nc"
   if ( varname(iv) .ne. "f_z" .and. varname(iv) .ne. "f_y" .and. \
        varname(iv) .ne. "epd" .and. varname(iv) .ne. "epd_z" ) then
     fname = fdir+"/uan"+expname+".epf2_koyz.nc"
   end if
   f1 = addfile(fname,"r")

   tmp1 = f1->$varname(iv)$({z0*1.e3},{lat1-5:lat1+5},:,:kmax)
   vspec(iv,:,:) = dim_avg_n( tmp1, 0 )*2.
 end do

else

 y1 = yyyymm1/100
 y2 = yyyymm2/100

 fdir  = "/prime0/kyh/dat/L60CGW-t/epf-wc"
 fname = fdir+"/"+y1+"/uan"+expname+".epf_koyz."+y1+".01.nc"
 f1 = addfile(fname,"r")
 no = dimsizes( f1->ome_fr )
 nk = kmax + 1

 kwn = fspan(-nk+1,nk-1,nk*2-1)
 ome = fspan(0,no/2,no/2+1) / 90.
 dome = 1./90.

 vspec = new((/nv,no,nk/),float)
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
if ( im .ge. 12 .or. im .le. 2 ) then
     do iv=0, nv-1
       if ( varname(iv) .eq. "f_vT" .or. varname(iv) .eq. "Fz/D" ) then
         break
       end if
       fname = fdir+"/"+iy+"/uan"+expname+".epf_koyz."+ \
               iy+"."+sprinti("%2.2i",im)+".nc"
       if ( varname(iv) .ne. "f_z" .and. varname(iv) .ne. "f_y" .and. \
            varname(iv) .ne. "epd" .and. varname(iv) .ne. "epd_z" ) then
         fname = fdir+"/"+iy+"/uan"+expname+".epf2_koyz."+ \
                 iy+"."+sprinti("%2.2i",im)+".nc"
       end if
       f1 = addfile(fname,"r")

       tmp1 = f1->$varname(iv)$({z0*1.e3},{lat1-5:lat1+5},:,:kmax)

       tmp = dim_avg_n( tmp1, 0 )*2.

       vspec(iv,:,:) = vspec(iv,:,:) + tmp(:,:)
     end do
     n = n + 1
end if
   end do
 end do
 nt = n
 vspec(:,:,:) = vspec(:,:,:)/nt

end if

 if ( varname(nv-1) .eq. "f_vT" ) then
   do iv=0, nv-2
     if ( varname(iv) .eq. "f_z" ) then
       iv_z = iv
     end if
     if ( varname(iv) .eq. "f_uw" ) then
       iv_uw = iv
     end if
   end do
   vspec(nv-1,:,:) = vspec(iv_z,:,:) - vspec(iv_uw,:,:)
 end if
 if ( varname(nv-1) .eq. "Fz/D" ) then
   do iv=0, nv-2
     if ( varname(iv) .eq. "f_z" ) then
       iv_z = iv
     end if
     if ( varname(iv) .eq. "epd" .or. varname(iv) .eq. "epd_z" ) then
       iv_d = iv
     end if
   end do
   vspec(nv-1,:,:) = vspec(iv_z,:,:)/where( vspec(iv_d,:,:) .ne. 0., \
                     vspec(iv_d,:,:), -1.e-20 )*(-1.)
 end if

 fname2 = "$DATD/L60CGW-t/cmxm.uan"+expname+"a.pj.nc"
 f2 = addfile(fname2,"r")

 rho0 = dim_avg(f2->rho({z0*1.e3},{lat1-5:lat1+5}))

;--------------------------------------- process 0
 do iv=0, nv-1
   if ( varname(iv) .eq. "epd" .or. varname(iv) .eq. "epd_z" .or. \
        varname(iv) .eq. "epd_y" ) then
     vspec(iv,:,:) = vspec(iv,:,:)*30.
   else
     if ( varname(iv) .eq. "Fz/D" ) then
       vspec(iv,:,:) = vspec(iv,:,:)/(6371.e3*rho0*cos(lat1*0.01745329)) \
                       * 86400. * 1.e-3  ; [km]
     else
       vspec(iv,:,:) = vspec(iv,:,:)/6371.e3 * 1.e3
     end if
   end if
 end do

 vsp = new((/nv,no/2+1,nk*2-1/),float)
 vsp = 0.
 do i=0, nk-2
   vsp(:,:,i) = vspec(:,0:no/2,nk-1-i)
 end do
 vsp(:,0,nk:) = vspec(:,0,1:);*0.5
 do n=1, no/2
   vsp(:,n,nk:) = vspec(:,no-n,1:)
 end do
 vsp(:,no/2,:) = vsp(:,no/2,:)*0.5

 vsp = vsp/dome  ; PS density

 vsp(:,:,nk-1) = vsp@_FillValue  ; k = 0
print(max(vsp))

;--------------------------------------- cn level
 vars = ":F10:"+(/"F:B:z:N:","F:B:z,M:N:","F:B:z,H:N:",":F21:EPD",":F21:EPDz",":F21:-F/D"/)
 facts = (/-3,-3,-3,-2,-2,3/)
 cnlevs = (/ (/-20,-10,-5,-2,-1,-0.5,-0.2,-0.1, \
               0.1,0.2,0.5,1,2,5,10,20/), \
             (/-20,-10,-5,-2,-1,-0.5,-0.2,-0.1, \
               0.1,0.2,0.5,1,2,5,10,20/), \
             (/-20,-10,-5,-2,-1,-0.5,-0.2,-0.1, \
               0.1,0.2,0.5,1,2,5,10,20/), \
             (/-20,-10,-5,-2,-1,-0.5,-0.2,-0.1, \
               0.1,0.2,0.5,1,2,5,10,20/), \
             (/-20,-10,-5,-2,-1,-0.5,-0.2,-0.1, \
               0.1,0.2,0.5,1,2,5,10,20/), \
             (/-20,-10,-5,-2,-1,-0.5,-0.2,-0.1, \
               0.1,0.2,0.5,1,2,5,10,20/) /)
 cnfil = (/17,19,20,21,23,24,25,26,0,11,10,9,8,6,5,4,2/)
 cnfil = cnfil(dimsizes(cnfil)-1:0)

 var   = new(nv,string )
 fact  = new(nv,integer)
 cnlev = new((/nv,dimsizes(cnlevs(0,:))/),float)
 lbtit = new(nv,string )
 do iv=0, nv-1
   if ( varname(iv) .eq. "f_z" ) then
     var  (iv) = vars (0)
     fact (iv) = facts(0)
	 cnlev(iv,:) = cnlevs(0,:)
     lbtit(iv) = " mPa / cyc day:S:-1:N: ]"
   end if
   if ( varname(iv) .eq. "f_uw" ) then
     var  (iv) = vars (1)
     fact (iv) = facts(1)
     cnlev(iv,:) = cnlevs(1,:)
     lbtit(iv) = " mPa / cyc day:S:-1:N: ]"
   end if
   if ( varname(iv) .eq. "f_vT" ) then
     var  (iv) = vars (2)
     fact (iv) = facts(2)
     cnlev(iv,:) = cnlevs(2,:)
     lbtit(iv) = " mPa / cyc day:S:-1:N: ]"
   end if
   if ( varname(iv) .eq. "epd" ) then
     var  (iv) = vars (3)
     fact (iv) = facts(3)
     cnlev(iv,:) = cnlevs(3,:)
     lbtit(iv) = " m s:S:-1:N: month:S:-1:N: / cyc day:S:-1:N: ]"
   end if
   if ( varname(iv) .eq. "epd_z" ) then
     var  (iv) = vars (4)
     fact (iv) = facts(4)
     cnlev(iv,:) = cnlevs(4,:)
     lbtit(iv) = " m s:S:-1:N: month:S:-1:N: / cyc day:S:-1:N: ]"
   end if
   if ( varname(iv) .eq. "Fz/D" ) then
     var  (iv) = vars (5)
     fact (iv) = facts(5)
     cnlev(iv,:) = cnlevs(5,:)
     lbtit(iv) = " km ]"
   end if
   if (fact(iv) .eq. 0) then
     lbtit(iv) = "["+lbtit(iv)
   else
     lbtit(iv) = "[ 10:S:"+fact(iv)+":N:"+lbtit(iv)
   end if
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
 xstr = ""
 xstr(nfig-1) = "zonal wavenumber"
 ystr = ""
 ystr(0) = "frequency [cyc day:S:-1:N:]"

do ip=0, nfig-1
do iv=0, nv-1

 ctr1 = cnshadeplot(wks,kwn,ome,vsp(iv,:,:)/(10.^fact(iv)),cnlev(iv,:),cnfil,"")
; cnaddline(ctr1,wks,kwn,ome,vsp(iv,:,:),cnlev(iv,:),0,"D",0,"T",1,"C",1)
 setvalues ctr1
  "cnFillMode" : "RasterFill"
 end setvalues

 set_size(ctr1,xf(iv),yf(ip),wi,he)
; axis_range(ctr1,-kmax-0.5,kmax+0.5,0.,omax+0.5*dome,False,False)
 axis_range(ctr1,-kmax-0.5,kmax+0.5,0.,omax,False,False)
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
   set_size_lb(label,xf(0)+wi*0.25,yf(nfig-1)+he+18+iv*8,wi*3,3)
   draw(label)
 end do

 frame(wks)

end

