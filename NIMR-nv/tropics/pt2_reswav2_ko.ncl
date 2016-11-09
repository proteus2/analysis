load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = "uj"
 varname = "theta"

 yyyymm1 = 200202
 yyyymm2 = 200603
 lat0 = 14.
 z0 = 21.0 ;17.41 ;18.59

 kmax = 10

 wks = cr_wks("fig/pt2_reswav2_ko-z"+numeric2int(z0,1)+"-"+expname(0),600,0)
 gsn_define_colormap(wks,"rainbow+gray")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 nfig = 2

 fdir  = "/prime0/kyh/dat/L60CGW/fcoef_ko"
 fname = fdir+"/"+(yyyymm1/100)+"/uan"+expname+".fft_theta_ko."+(yyyymm1/100)+".01.nc"
 f1 = addfile(fname,"r")
 no = dimsizes( f1->ome_fr )
 nk = dimsizes( f1->k_wn(:kmax+1))

 vspec = new((/nfig,no,nk/),float)
 vspec = 0.
 n = 0
 do iy=yyyymm1/100, yyyymm2/100
   im1 = 1
   im2 = 12
   if ( iy .eq. yyyymm1/100 ) then
     im1 = yyyymm1 - iy*100
   end if
   if ( iy .eq. yyyymm2/100 ) then
     im2 = yyyymm2 - iy*100
   end if
   do im=im1, im2
print(iy+"  "+im)
     fname = fdir+"/"+iy+"/uan"+expname+".fft_theta_ko."+iy+"."+sprinti("%2.2i",im)+".nc"

     f1 = addfile(fname,"r")

       tmp1r = f1->$("fcr_"+varname)$({z0*1.e3},{-lat0:lat0},:,:kmax+1)
       tmp2r = tmp1r({lat0:-lat0},:,:)
       tmp1i = f1->$("fci_"+varname)$({z0*1.e3},{-lat0:lat0},:,:kmax+1)
       tmp2i = tmp1i({lat0:-lat0},:,:)

       tmp = dim_avg_n( \
             (/ (0.5*(tmp1r + tmp2r))^2 + (0.5*(tmp1i + tmp2i))^2, \
                (0.5*(tmp1r - tmp2r))^2 + (0.5*(tmp1i - tmp2i))^2 /), 1 )
       tmp(:,:,:) = tmp(:,:,:)*2.

       vspec(:,:,:) = vspec(:,:,:) + tmp(:,:,:)

     n = n + 1
   end do
 end do
 nt = n
 vspec(:,:,:) = vspec(:,:,:)/nt

 fname2 = "$DATD/L60CGW/cmxm.uan"+expname+"a.pj.nc"
 f2 = addfile(fname2,"r")

 k0 = minind(abs( f2->hybrid_ht - z0*1.e3 ))
 k1 = k0 + 1
 if ( f2->hybrid_ht(k0) .gt. z0*1.e3 ) then
   k1 = k0 - 1
 end if
k1 = k0  ; temp
 p0 = dim_avg(f2->p(k0,{-lat0:lat0}))
 p1 = dim_avg(f2->p(k1,{-lat0:lat0}))
 pi = exp(0.5*(log(p0) + log(p1)))

 vspec = vspec * ((pi/1.e5)^0.287)^2

;--------------------------------------- process 0
 kwn = fspan(0,nk-1,nk)
 ome = fspan(-no/2,no/2,no+1) / 90.
 dome = 1./90.
print(ome)

 vsp = new((/nfig,no+1,nk/),float)
 vsp = 0.
 do n=1, no/2
   vsp(:,no/2-1+n,:) = vspec(:,no-n,:)
 end do
 do n=0, no/2-1
   vsp(:,n,:) = vspec(:,no/2-n,:)
 end do

 vsp(:,:,:) = vsp(:,:,:)/dome

; vsp(:,:,nk-1) = vsp@_FillValue
print(max(vsp))
print(sum(vsp(0,no/2-no/16:no/2+no/16,:7))*dome)
print(sum(vsp(1,no/2-no/16:no/2+no/16,:7))*dome)

;--------------------------------------- cn level
 cnlev = new((/2,30/),float)
 cnlev(0,:) = fspan(0.03333333333333,1.,30)
 cnlev(1,:) = cnlev(0,:)*0.3
 cnfil = ispan(1,211,7)
 cnfil(16:) = cnfil(16:) + 20
 cnfil(25:) = cnfil(25:) + 5
;print(dimsizes(cnfil))

 cnlev2= (/1,2,4,8,16/)

 cnpc = (/5,10,20,30,40,50,60,70,80,90,100/)

;--------------------------------------- ctr
 wi = 50
 he = 70
 xf = 34 + (/0,1,1,1,2,2/)*(wi+20)
 yf = 10 + (/0,0,0,1,0,1/)*55

 mstr  = new(nfig,string)
 xstr  = new(nfig,string)
 ystr  = new(nfig,string)
 lbopt = new(nfig,string)

 mstr(:) = ""
 xstr(:) = "zonal wavenumber"
; xstr(nfig-1) = "year"
 ystr(:) = "frequency [cyc day:S:-1:N:]"

do ip=0, nfig-1

 ctr1 = cnshadeplot(wks,kwn,ome,vsp(ip,:,:),cnlev(ip,:),cnfil,"")
; cnaddline(ctr1,wks,kwn,ome,vsp(ip,:,:),cnlev2,0,"D",0,"T",1,"C",1)
; setvalues ctr1
;  "cnFillMode" : "RasterFill"
; end setvalues

 set_size(ctr1,xf(ip),yf(ip),wi,he)
 axis_range(ctr1,0,7,-1,1,False,False)
 axis_str(ctr1,mstr(ip),xstr(ip),ystr(ip))
 setvalues ctr1
  "tmYLMode" : "Manual"
  "tmYLTickSpacingF"  : 0.5
;  "tmXBTickStartF"    : 1942
  "tmYLMinorPerMajor" : 4
 end setvalues

 cnfinalize(ctr1,wks)

 draw(ctr1)

   label = create "lab" labelBarClass wks
;     "lbAutoManage" : False
     "lbOrientation" : "Vertical"
     "lbBoxCount"   : dimsizes(cnfil)
     "lbFillColors" : cnfil
     "lbMonoFillPattern" : True
     "lbLabelsOn" : False ;True
     "lbLabelStrings" : cnlev(ip,:)
     "lbLabelStride" : 3
;     "lbAutoManage" : False
     "lbLabelAlignment" : "InteriorEdges"
     "lbPerimOn" : False
;     "lbTitleString" : "[ K:S:2:N: / cyc day:S:-1:N: ]"
;     "lbTitlePosition" : "Bottom"
   end create
   set_size_lb(label,xf(ip)+wi+2,yf(ip),5,he)
   draw(label)

end do

 frame(wks)

end



