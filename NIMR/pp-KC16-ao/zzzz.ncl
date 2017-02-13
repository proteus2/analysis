load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = "tc"
 varname = (/"f_z","f_uw","f_vT"/)

 yyyymm1 = 200501
 yyyymm2 = 204412
 ydif = 50
 lat0 = 15.  ; 15. ; 5.
 z0 = 17.41  ; 17.41 ; 18.59

 kmax = 15;20
 omax = 0.7

 beta = 2.*(7.292116e-5)/6371229.
 nbv = 2.5e-2

; wks = cr_wks("fig/fzsa_ko-2-z"+numeric2int(z0,1)+"-y"+numeric2int(lat0,1)+"-"+expname(0),600,0)
 wks = cr_wks("zzzz",600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 nv = dimsizes(varname)

 wcode = (/"03","05","04"/)
 ng = dimsizes(wcode)
 nfig = ng

do nn=0, 1

 y1 = yyyymm1/100 + nn*ydif
 y2 = yyyymm2/100 + nn*ydif

 fdir  = "$DATD/AOL60CGW/epf-wc_sepa/avg"
 fname = fdir+"/uan"+expname+".epf_koyz."+y1+"-"+y2+"_sepa0_"+wcode(0)+".nc"

if ( isfilepresent(fname) ) then

 f1 = addfile(fname,"r")
 no = dimsizes( f1->ome_fr )
 nk = kmax + 1

 kwn = fspan(-kmax,kmax,kmax*2+1)
 ome = fspan(0,no/2,no/2+1) / 90.
 dome = 1./90.

 if ( .not. isdefined("vspec") ) then
   vspec = new((/2,ng,nv,no,nk/),float)
 end if
 do ig=0, ng-1
 do iv=0, nv-1
   if ( varname(iv) .eq. "f_vT" ) then
     break
   end if
   fname = fdir+"/uan"+expname+".epf_koyz."+y1+"-"+y2+"_sepa0_"+wcode(ig)+".nc"
   f1 = addfile(fname,"r")

   vspec(nn,ig,iv,:,:) = dim_avg_n(  \
      f1->$varname(iv)$({z0*1.e3},{-lat0:lat0},:,:kmax), 0 )*2.

   if (wcode(ig) .eq. "04") then
     fname2 = fdir+"/uan"+expname+".epf_koyz."+y1+"-"+y2+"_sepa0_14.nc"
     fname3 = fdir+"/uan"+expname+".epf_koyz."+y1+"-"+y2+"_sepa0_07.nc"
     fname4 = fdir+"/uan"+expname+".epf_koyz."+y1+"-"+y2+"_sepa0_17.nc"
     f2 = addfile(fname2,"r")
     f3 = addfile(fname3,"r")
     f4 = addfile(fname4,"r")

     vspec(nn,ig,iv,:,:) = vspec(nn,ig,iv,:,:) + dim_avg_n(  \
        f2->$varname(iv)$({z0*1.e3},{-lat0:lat0},:,:kmax) +  \
        f3->$varname(iv)$({z0*1.e3},{-lat0:lat0},:,:kmax) +  \
        f4->$varname(iv)$({z0*1.e3},{-lat0:lat0},:,:kmax) , 0 )*2.
   end if

 end do
 end do

else

 fdir  = "$DATD/AOL60CGW/epf-wc_sepa"
 fname = fdir+"/"+y1+"/uan"+expname+".epf_koyz."+y1+".01_sepa0_"+wcode(0)+".nc"
 f1 = addfile(fname,"r")
 no = dimsizes( f1->ome_fr )
 nk = kmax + 1

 kwn = fspan(-nk+1,nk-1,nk*2-1)
 ome = fspan(0,no/2,no/2+1) / 90.
 dome = 1./90.

 if ( .not. isdefined("vspec") ) then
   vspec = new((/2,ng,nv,no,nk/),float)
 end if
 vspec(nn,:,:,:,:) = 0.
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
     do ig=0, ng-1
     do iv=0, nv-1
       if ( varname(iv) .eq. "f_vT" ) then
         break
       end if
       fname = fdir+"/"+iy+"/uan"+expname+".epf_koyz."+ \
               iy+"."+sprinti("%2.2i",im)+"_sepa0_"+wcode(ig)+".nc"
       f1 = addfile(fname,"r")

       tmp = dim_avg_n( f1->$varname(iv)$({z0*1.e3},{-lat0:lat0},:,:kmax),  \
                        0 )*2.

       vspec(nn,ig,iv,:,:) = vspec(nn,ig,iv,:,:) + tmp(:,:)

       if (wcode(ig) .eq. "04") then
         fname2 = fdir+"/"+iy+"/uan"+expname+".epf_koyz."+ \
                 iy+"."+sprinti("%2.2i",im)+"_sepa0_14.nc"
         fname3 = fdir+"/"+iy+"/uan"+expname+".epf_koyz."+ \
                 iy+"."+sprinti("%2.2i",im)+"_sepa0_07.nc"
         fname4 = fdir+"/"+iy+"/uan"+expname+".epf_koyz."+ \
                 iy+"."+sprinti("%2.2i",im)+"_sepa0_17.nc"
         f2 = addfile(fname2,"r")
         f3 = addfile(fname3,"r")
         f4 = addfile(fname4,"r")

         tmp = dim_avg_n( f2->$varname(iv)$({z0*1.e3},{-lat0:lat0},:,:kmax) +  \
                          f3->$varname(iv)$({z0*1.e3},{-lat0:lat0},:,:kmax) +  \
                          f4->$varname(iv)$({z0*1.e3},{-lat0:lat0},:,:kmax)    \
                          , 0 )*2.

         vspec(nn,ig,iv,:,:) = vspec(nn,ig,iv,:,:) + tmp(:,:)
       end if

     end do
     end do
     n = n + 1
   end do
 end do
 nt = n
 vspec(nn,:,:,:,:) = vspec(nn,:,:,:,:)/nt

end if

end do  ; nn

 if ( varname(nv-1) .eq. "f_vT" ) then
   do iv=0, nv-2
     if ( varname(iv) .eq. "f_z" ) then
       iv_z = iv
     end if
     if ( varname(iv) .eq. "f_uw" ) then
       iv_uw = iv
     end if
   end do
   vspec(:,:,nv-1,:,:) = vspec(:,:,iv_z,:,:) - vspec(:,:,iv_uw,:,:)
 end if

; fname2 = "$DATD/L60CGW-t/cmxm.uan"+expname+"a.pj.nc"
; f2 = addfile(fname2,"r")
;
; rho0 = dim_avg(f2->rho({z0*1.e3},{-lat0:lat0}))

;--------------------------------------- process 0
 vspec = vspec/6371.e3 * 1.e3

 vsp = new((/3,nfig,nv,no/2+1,nk*2-1/),float)
 vsp = 0.
 do i=0, nk-2
   vsp(:1,:,:,:,i) = vspec(:,:,:,0:no/2,nk-1-i)
 end do
 vsp(:1,:,:,0,nk:) = vspec(:,:,:,0,1:);*0.5
 do n=1, no/2
   vsp(:1,:,:,n,nk:) = vspec(:,:,:,no-n,1:)
 end do
 vsp(:1,:,:,no/2,:) = vsp(:1,:,:,no/2,:)*0.5

 vsp(2,:,:,:,:) = vsp(1,:,:,:,:) - vsp(0,:,:,:,:)

 vsp = vsp/dome  ; PS density

 vsp(:,:,:,:,nk-1) = vsp@_FillValue  ; k = 0
print(max(vsp))

 vspd = vsp(1,:,:,:,:) - vsp(0,:,:,:,:)
 vsp0 = vsp(0,:,:,:,:)
 do iter=0, 3
   tmpd = vspd
   tmp0 = vsp0
   do n=1, no/2-1
     vspd(:,:,n,:) = 0.5*tmpd(:,:,n,:) +  \
                     0.25*(tmpd(:,:,n-1,:) + tmpd(:,:,n+1,:))
     vsp0(:,:,n,:) = 0.5*tmp0(:,:,n,:) +  \
                     0.25*(tmp0(:,:,n-1,:) + tmp0(:,:,n+1,:))
   end do
   vspd(:,:,0,:) = 0.5*tmpd(:,:,0,:) + 0.25*(tmpd(:,:,1,nk*2-2:0) + tmpd(:,:,1,:))
   vsp0(:,:,0,:) = 0.5*tmp0(:,:,0,:) + 0.25*(tmp0(:,:,1,nk*2-2:0) + tmp0(:,:,1,:))
 end do
 vsp@_FillValue = 1.e20
 vsp0@_FillValue = 1.e20
 vsp = where( vsp .eq. 0., vsp@_FillValue, vsp )
 vsp0 = where( vsp(0,:,:,:,:) .eq. 0. , vsp0@_FillValue, vsp0 )
 vspd = where( abs(vsp0) .lt. 0.5e-3, vspd@_FillValue, vspd )
 vspd = vspd/vsp0 * 100.

 nk10 = nk*2*10
 no10 = no/2*10

 nw = 6 ;4
 dispc = new((/nw,no10+1,nk10+1/),float)
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
   dispc(0,n,nk10/2+1:) = two_pi/1.e3*o(n)/(nbv*k(nk10/2+1:))
   dispc(3,n,:) = two_pi/1.e3*o(n)^2/  \
        (nbv*beta*(1.5+sqrt(1.5*1.5+okbi(n,:)*(1.+okbi(n,:)))))
   dispc(1,n,:) = two_pi/1.e3*o(n)^2/(nbv*beta*(1.+okbi(n,:)))
   dispc(5,n,:) = two_pi/1.e3*o(n)^2/  \
        (nbv*beta*(2.5+sqrt(2.5*2.5+okbi(n,:)*(1.+okbi(n,:)))))
 end do
 dispc(0,:,nk10/2) = 1.e3
 dispc(1,:,:) = where(okbi .gt. -1., dispc(1,:,:), 1.e30)
 ; R
 do n=1, no10
   dispc(2,n,:nk10/2-1) = two_pi/1.e3*o(n)^2/  \
        (nbv*beta*(1.5-sqrt(1.5*1.5+okbi(n,:nk10/2-1)*(1.+okbi(n,:nk10/2-1)))))
   dispc(4,n,:nk10/2-1) = two_pi/1.e3*o(n)^2/  \
        (nbv*beta*(2.5-sqrt(2.5*2.5+okbi(n,:nk10/2-1)*(1.+okbi(n,:nk10/2-1)))))
 end do
 dispc(2,0,:) = 0.
 dispc(4,0,:) = 0.
 dispc(2,:,:) = where(okbi .gt. -1., dispc(2,:,:), 1.e30)
 dispc(4,:,:) = where(okbi .gt. -1., dispc(4,:,:), 1.e30)

;--------------------------------------- cn level
 vars = ":F10:"+(/"F:B:z:N:","F:B:z,M:N:","F:B:z,H:N:"/)
 facts = (/-3,-3,-3/)
 cnlevs = (/ (/-256,-128,-64,-32,-16,-8,-4,-2,-1, \
               1,2,4,8,16,32,64,128,256/), \
             (/-256,-128,-64,-32,-16,-8,-4,-2,-1, \
               1,2,4,8,16,32,64,128,256/), \
              (/-256,-128,-64,-32,-16,-8,-4,-2,-1, \
               1,2,4,8,16,32,64,128,256/) /)
 cnfil = (/17,18,19,20,21,23,24,25,26,0,11,10,9,8,6,5,4,3,2/)
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

 lzs = (/2.,5.,15./)
 ldsh = (/0,16,14/)

;--------------------------------------- ctr
 np = 3

 wi = 32
 he = 32
 xf = 12 + (/0,1,2/)*(wi+8) + (/0,0,4/)
 yf = 10 + (/0,1,2,3,4,5/)*(/he+10/)

 mstr  = new((/nfig,np/),string)
 xstr  = new(nfig,string)
 ystr  = new(np,string)

 mstr = ""
; mstr(0,:) = var(:)
 mstr(0,:) = (/ (yyyymm1/100)+"-"+(yyyymm2/100),  \
                (yyyymm1/100+ydif)+"-"+(yyyymm2/100+ydif),  \
                "Diff." /)
 xstr = ""
 xstr(nfig-1) = "zonal wavenumber"
 ystr = ""
 ystr(0) = "frequency [cyc day:S:-1:N:]"

do ip=0, 2
do ig=0, nfig-1
 iv = 0

 ctr1 = cnshadeplot(wks,kwn,ome,vsp(2,ig,ip,:,:)/(10.^fact(iv)),cnlev(iv,:),cnfil,"")
; cnaddline(ctr1,wks,kwn,ome,vsp(0,ig,iv,:,:),cnlev(iv,:),0,"D",0,"T",1,"C",1)
; cnaddline(ctr1,wks,kwn,ome,vspd(ig,iv,:,:),(/-100,-50,-30,-10,0,10,30,50,100/),1,"D",(/0,0,1/),"T",1.2,"C",(/1,0,1/))
 if (ig .lt. 2) then
   cnaddline(ctr1,wks,kwn2,ome2,dispc(ig,:,:),lzs,0,"D",ldsh(0),"T",1.2,"C",1)
 end if
 setvalues ctr1
  "cnFillMode" : "RasterFill"
 end setvalues

 set_size(ctr1,xf(ip),yf(ig),wi,he)
; axis_range(ctr1,-kmax-0.5,kmax+0.5,0,omax+0.5*dome,False,False)
 axis_range(ctr1,-kmax-0.5,kmax+0.5,0,omax,False,False)
 axis_str(ctr1,mstr(ig,ip),xstr(ig),ystr(ip))

 cnfinalize(ctr1,wks)

 draw(ctr1)

end do
end do

 do ip=0, 0 ;nv-1
   iv = 0
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
   set_size_lb(label,xf(0)+wi*0.25,yf(nfig-1)+he+18+ip*8,wi*3,3)
   draw(label)
 end do

 frame(wks)

end

