load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
; expname = "merra"
; cexp = "MERRA"
 expname = "era-int"
 cexp = "ERA-I"
 varname = (/"f_z","f_uw","f_vT"/)

 yyyymm1 = 197902
 yyyymm2 = 201012
 lat0 = 15.  ; 15. ; 5.
 p0 = 50.  ; 17.41 ; 18.59

 kmax = 15
 omax = 0.75;1.0

 beta = 2.*(7.292116e-5)/6371229.
 nbv = 2.5e-2

 wks = cr_wks("fig/fzsa_ko-p"+numeric2int(p0,1)+"-y"+numeric2int(lat0,1)+"-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 nfig = 2
 nv = dimsizes(varname)

if ( yyyymm1 .eq. 197902 .and. yyyymm2 .eq. 200612 ) then

 fdir  = "$DATD/tem/"+cexp+"/epf/avg"
 fname = fdir+"/"+expname+".epfsa0_koyz.nc"
 f1 = addfile(fname,"r")
 no = dimsizes( f1->ome_fr )
 nk = kmax + 1

 kwn = fspan(-kmax,kmax,kmax*2+1)
 ome = fspan(0,no/2,no/2+1) / 90.
 dome = 1./90.

 vspec = new((/2,nv,no,nk/),float)
 do iv=0, nv-1
   if ( varname(iv) .eq. "f_vT" ) then
     break
   end if
   fname = fdir+"/"+expname+".epfsa0_koyz.nc"
   if ( varname(iv) .ne. "f_z" .and. varname(iv) .ne. "f_y" ) then
     fname = fdir+"/"+expname+".epf2sa0_koyz.nc"
   end if
   f1 = addfile(fname,"r")

   tmps = f1->$(varname(iv)+"_s")$({p0},{-lat0:lat0},:,:kmax)
   tmpa = f1->$(varname(iv)+"_a")$({p0},{-lat0:lat0},:,:kmax)
   vspec(:,iv,:,:) = dim_avg_n( (/ tmps, tmpa /), 1 )*2.
 end do

else

 y1 = yyyymm1/100
 y2 = yyyymm2/100

 fdir  = "$DATD/tem/"+cexp+"/epf"
 fname = fdir+"/"+y1+"/"+expname+".epfsa0_koyz."+y1+".02.nc"
 f1 = addfile(fname,"r")
 no = dimsizes( f1->ome_fr )
 nk = kmax + 1

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
       if ( varname(iv) .eq. "f_vT" ) then
         break
       end if
       fname = fdir+"/"+iy+"/"+expname+".epfsa0_koyz."+ \
               iy+"."+sprinti("%2.2i",im)+".nc"
       if ( varname(iv) .ne. "f_z" .and. varname(iv) .ne. "f_y" ) then
         fname = fdir+"/"+iy+"/"+expname+".epf2sa0_koyz."+ \
                 iy+"."+sprinti("%2.2i",im)+".nc"
       end if
       f1 = addfile(fname,"r")

       tmps = f1->$(varname(iv)+"_s")$({p0},{-lat0:lat0},:,:kmax)
       tmpa = f1->$(varname(iv)+"_a")$({p0},{-lat0:lat0},:,:kmax)

       tmp = dim_avg_n( (/ tmps, tmpa /), 1 )*2.

       vspec(:,iv,:,:) = vspec(:,iv,:,:) + tmp(:,:,:)
     end do
     n = n + 1
   end do
 end do
 nt = n
 vspec(:,:,:,:) = vspec(:,:,:,:)/nt

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
   vspec(:,nv-1,:,:) = vspec(:,iv_z,:,:) - vspec(:,iv_uw,:,:)
 end if

;--------------------------------------- process 0
 vspec = vspec/6371.e3 * 1.e3

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
 vars = ":F10:"+(/"F:B:z:N:","F:B:z,M:N:","F:B:z,H:N:"/)
 facts = (/-3,-3,-3/)
; cnlevs = (/ (/-100,-50,-20,-10,-5,-2,-1,-0.5, \
;               0.5,1,2,5,10,20,50,100/), \
;             (/-100,-50,-20,-10,-5,-2,-1,-0.5, \
;               0.5,1,2,5,10,20,50,100/), \
;             (/-100,-50,-20,-10,-5,-2,-1,-0.5, \
;               0.5,1,2,5,10,20,50,100/) /)
 cnlevs = (/ (/-50,-20,-10,-5,-2,-1,-0.5,-0.2, \
               0.2,0.5,1,2,5,10,20,50/), \
             (/-50,-20,-10,-5,-2,-1,-0.5,-0.2, \
               0.2,0.5,1,2,5,10,20,50/), \
             (/-50,-20,-10,-5,-2,-1,-0.5,-0.2, \
               0.2,0.5,1,2,5,10,20,50/) /)
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
     lbtit(iv) = "[ 10:S:"+fact(iv)+":N: mPa / cyc day:S:-1:N: ]"
   end if
   if ( varname(iv) .eq. "f_uw" ) then
     var  (iv) = vars (1)
     fact (iv) = facts(1)
     cnlev(iv,:) = cnlevs(1,:)
     lbtit(iv) = "[ 10:S:"+fact(iv)+":N: mPa / cyc day:S:-1:N: ]"
   end if
   if ( varname(iv) .eq. "f_vT" ) then
     var  (iv) = vars (2)
     fact (iv) = facts(2)
     cnlev(iv,:) = cnlevs(2,:)
     lbtit(iv) = "[ 10:S:"+fact(iv)+":N: mPa / cyc day:S:-1:N: ]"
   end if
   if (fact(iv) .eq. 0) then
     lbtit(iv) = "[ mPa / cyc day:S:-1:N: ]"
   end if
 end do

 lzs = (/ (/ (/2.,5.,15./),    \
             (/2.5,5.0,10./),    \
             (/2.,10.,60./) /), \
          (/ (/2.,5.,30./),    \
             (/2.5,5.0,10./),    \
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
   set_size_lb(label,xf(0)+wi*0.25,yf(nfig-1)+he+18+iv*8,wi*3,3)
   draw(label)
 end do

 frame(wks)

end

