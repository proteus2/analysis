load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = "uj"
 varname = (/"dT_conv_nsc_1","dT_conv_nsc"/)

 mm1 = 1
 mm2 = 12
 lat0 = 20.
 z1 = 1.9
 z2 = 16.9

 kmax = 20
 omax = 1.0
 kmax0 = 40
 omax0 = 3.2

 wks = cr_wks("fig/ch_ko-z"+numeric2int(z1,1)+"-"+numeric2int(z2,1)+"-"+expname(0),600,0)
 gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 if ( mm1 .ne. 1 .and. mm2 .eq. 12 ) then
   print("Take care for Dec. (1953-2005, not 2006), STOP")
   status_exit(12)
 end if

 nfig = 2
 nv = dimsizes(varname)

 fdir  = "/prime0/kyh/dat/L60CGW/fcoef_ko/avg"
 fname = fdir+"/uan"+expname+".fc2_"+varname(0)+"_ko.nc"
 f1 = addfile(fname,"r")
 no = dimsizes( f1->ome_fr )
 nk = dimsizes( f1->k_wn   )
 nz = dimsizes( f1->ht({z1*1.e3:z2*1.e3}) )
print( f1->ht({z1*1.e3:z2*1.e3})/1.e3 )

 kwn = fspan(-nk+1,nk-1,nk*2-1)
 ome = fspan(0,no/2,no/2+1) / 90.
 dome = 1./90.

 vspec4 = new((/nv,nz,no,nk/),float)
 vspec4 = 0.
 if ( mm1 .eq. 1 .and. mm2 .eq. 12 ) then
   do iv=0, nv-1
     fname = fdir+"/uan"+expname+".fc2_"+varname(iv)+"_ko.nc"
     f1 = addfile(fname,"r")
     tmp = dim_avg_n( \
             f1->$("fc2_"+varname(iv))$({z1*1.e3:z2*1.e3},{-lat0:lat0},:,:), 1 )
     tmp(:,:,:nk-2) = tmp(:,:,:nk-2)*2.
     vspec4(iv,:,:,:) = tmp
   end do
 else
   n = 0
   do im=mm1, mm2
     do iv=0, nv-1
       fname = fdir+"/uan"+expname+".fc2_"+varname(iv)+"_ko."+ \
               sprinti("%2.2i",im)+".nc"

       f1 = addfile(fname,"r")

       tmp = dim_avg_n( \
               f1->$("fc2_"+varname(iv))$({z1*1.e3:z2*1.e3},{-lat0:lat0},:,:), 1 )

       tmp(:,:,:nk-2) = tmp(:,:,:nk-2)*2.

       vspec4(iv,:,:,:) = vspec4(iv,:,:,:) + tmp
     end do
     n = n + 1
   end do
   nt = n
   vspec4(:,:,:,:) = vspec4(:,:,:,:)/nt
 end if

 fname2 = "$DATD/L60CGW/cmxm.uan"+expname+"a.pj.nc"
 f2 = addfile(fname2,"r")

 th0 = dim_avg(f2->theta({z1*1.e3:z2*1.e3},{-lat0:lat0}))

 k0 = minind(abs( f2->$( f2->rho!0 )$ - z1*1.e3 ))
 k2 = minind(abs( f2->$( f2->rho!0 )$ - z2*1.e3 ))
 k1 = k0 + 1
 k3 = k2 + 1
 if ( f2->$( f2->rho!0 )$(k0) .gt. z1*1.e3 ) then
   k1 = k0 - 1
 end if
 if ( f2->$( f2->rho!0 )$(k2) .gt. z2*1.e3 ) then
   k3 = k2 - 1
 end if
 if ( k0-k1 .ne. k2-k3 ) then
   k2 = k3
   k3 = k2 - k0 + k1
 end if
 rho0 = dim_avg(f2->rho(k0:k2,{-lat0:lat0}))
 rho1 = dim_avg(f2->rho(k1:k3,{-lat0:lat0}))
 rho0t = exp(0.5*(log(rho0) + log(rho1)))
 p0 = dim_avg(f2->p(k0:k2,{-lat0:lat0}))
 p1 = dim_avg(f2->p(k1:k3,{-lat0:lat0}))
 p0t = exp(0.5*(log(p0) + log(p1)))
 t2 = (th0*(p0t/1.e5)^0.286)^2

 norm = rho0t / t2

 do k=0, nz-1
   vspec4(:,k,:,:) = vspec4(:,k,:,:) * norm(k)
 end do

 vspec = dim_avg_n(vspec4,1)

;--------------------------------------- process 0
 vsp = new((/nfig,nv,no/2+1,nk*2-1/),float)
 vsp = 0.
 do i=0, nk-2
   vsp(0,:,:,i) = vspec(:,0:no/2,nk-1-i)
 end do
 vsp(0,:,0,nk:) = vspec(:,0,1:);*0.5
 do n=1, no/2
   vsp(0,:,n,nk:) = vspec(:,no-n,1:)
 end do
 vsp(0,:,no/2,:) = vsp(0,:,no/2,:)*0.5

 vsp(1,:,:,:) = vsp(0,:,:,:)

 if (kmax0 .ne. 999) then
;   vsp(0,:,:,nk-1-kmax0) = dim_sum( vsp(0,:,:,:nk-1-kmax0) )
;   vsp(0,:,:,nk-1+kmax0) = dim_sum( vsp(0,:,:,nk-1+kmax0:) )
   vsp(0,:,:,:nk-1-kmax0-1) = 0.
   vsp(0,:,:,nk-1+kmax0+1:) = 0.
   kmax0 = kmax0 + 1
 end if
 if (omax0 .ne. 999) then
   nn = toint(omax0*90.)
;   vsp(0,:,nn,:) = dim_sum_n( vsp(0,:,nn:,:), 1 )
   vsp(0,:,nn+1:,:) = 0.
   omax0 = omax0 + dome
 end if

 vsp = vsp/dome  ; PS density

 vsp(:,:,:,nk-1) = vsp@_FillValue  ; k = 0
print(max(vsp))

;--------------------------------------- cn level
 vars = ":F10:"+(/"H'/T:B:0:N:","H'/T:B:0:N:"/)
 facts = (/-11,-11/)
 cnlevs = (/ (/0.1,0.2,0.5,1,2,5,10,20,50,100,200/), \
             (/0.1,0.2,0.5,1,2,5,10,20,50,100,200/) /)

 cnfil = (/5,6,7,8,10,12,14,15,16,17,19,21/)

 var   = new(nv,string )
 fact  = new(nv,integer)
 cnlev = new((/nv,dimsizes(cnlevs(0,:))/),float)
 lbtit = new(nv,string )
 do iv=0, nv-1
;   if ( varname(iv) .eq. "theta" ) then
     var  (iv) = vars (iv)
     fact (iv) = facts(iv)
     cnlev(iv,:) = cnlevs(iv,:)
     lbtit(iv) = "[ 10:S:"+fact(iv)+":N: s:S:-2:N: / cyc day:S:-1:N: ]"
;   end if
 end do

;--------------------------------------- ctr
 wi = 32
 he = 32
 xf = 12 + (/0,1,2,3/)*(wi+8)
 yf = 10 + (/0,1/)*(/he+15/)

 mstr  = new((/nfig,nv/),string)
 xstr  = new(nfig,string)
 ystr  = new(nv,string)

 mstr = ""
 mstr(0,:) = vars(:)
 xstr = "zonal wavenumber"
 ystr = ""
 ystr(0) = "frequency [cyc day:S:-1:N:]"

do ip=0, nfig-1
do iv=0, nv-1

 ctr1 = cnshadeplot(wks,kwn,ome,vsp(ip,iv,:,:)/(10.^fact(iv)),cnlev(iv,:),cnfil,"")
; cnaddline(ctr1,wks,kwn,ome,vsp(ip,iv,:,:),cnlev(iv,:),0,"D",0,"T",1,"C",1)
 setvalues ctr1
  "cnFillMode" : "RasterFill"
 end setvalues

 set_size(ctr1,xf(iv),yf(ip),wi,he)
 axis_range(ctr1,tofloat(-kmax0),tofloat(kmax0),0,omax0+0.5*dome,False,False)
 if (ip .eq. 1) then
   axis_range(ctr1,-kmax-0.5,kmax+0.5,0,omax+0.5*dome,False,False)
 end if
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
   set_size_lb(label,xf(max((/0,nv-2/))),yf(nfig-1)+he+20+iv*8,wi*2+8,3)
   draw(label)
 end do

 frame(wks)

end

