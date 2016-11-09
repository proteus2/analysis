load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = "uj"
 varname = "epd"

 y1 = 1953
 y2 = 2005
 lat0 = 5.
 z0 = 27.0 ;18.59

 wks = cr_wks("fig/epd_reswav3_ko-z"+numeric2int(z0,1)+"-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 nfig = 2

 fdir  = "/data4/kyh/dat/L60CGW-t/epf"
 fname = fdir+"/"+y1+"/uan"+expname+".epf_wave_koyz."+y1+".01.nc"
 f1 = addfile(fname,"r")
 no = dimsizes( f1->ome_fr )
 nk = dimsizes( f1->k_wn   )

 nt = (y2-y1+1)*12
 ub = new(nt,float)
 vspec = new((/nfig,nt,no,nk/),float)
 vspec = 0.
 n = 0
 do iy=y1, y2
 do im=1, 12
   fname = fdir+"/"+iy+"/uan"+expname+".epf_wave_koyz."+iy+"."+sprinti("%2.2i",im)+".nc"

   f1 = addfile(fname,"r")

   tmp = dim_avg_n( f1->$(varname+"_s")$({z0*1.e3},{-lat0:lat0},:,:), 0 )

   vspec(0,n,:,:) = tmp*2.

   tmp = dim_avg_n( f1->$(varname+"_a")$({z0*1.e3},{-lat0:lat0},:,:), 0 )

   vspec(1,n,:,:) = tmp*2.

   fname2 = "$DATD/L60CGW-t/xm/xm.uan"+expname+"a.pj_"+iy+sprinti("%2.2i",im)+".nc"
   f2 = addfile(fname2,"r")

   ub(n) = avg(f2->u(0,{z0*1.e3},{-lat0:lat0}))

   n = n + 1
 end do
 end do
print(ub)

;--------------------------------------- process 0
 kwn = fspan(-nk+1,nk-1,nk*2-1)
 ome = fspan(0,no/2,no/2+1) / 90.
 dome = 1./90.

 vsp = new((/nfig,nt,no/2+1,nk*2-1/),float)
 vsp = 0.
 do i=0, nk-2
   vsp(:,:,:,i) = vspec(:,:,0:no/2,nk-1-i)
 end do
 vsp(:,:,0,nk:) = vspec(:,:,0,1:);*0.5
 do n=1, no/2
   vsp(:,:,n,nk:) = vspec(:,:,no-n,1:)
 end do
; vsp(:,:,no/2,:) = vsp(:,:,no/2,:)*0.5

 vsp(:,:,:,:) = vsp(:,:,:,:)/dome

 vsp(:,:,:,nk-1) = vsp@_FillValue

 delete( vspec )

 vspi = new((/nfig,nt,no/2+1,nk*2-1/),float)
 vspi = 0.
 kr = kwn/(2.*3.141592*6371229.)
 do n=0, nt-1
 do i=0, nk-2
   sht = numeric2int(ub(n)*abs(kr(i))*86400./dome,1)
   isht = abs(sht)
   if (sht .eq. 0) then
     vspi(:,n,:,:) = vsp(:,n,:,:)
   end if
   if (sht .gt. 0) then
     vspi(:,n,isht:no/2,i) = vsp(:,n,0:no/2-isht,i)
     vspi(:,n,0:no/2-isht,nk*2-2-i) = vsp(:,n,isht:no/2,nk*2-2-i)
     vspi(:,n,0:isht-1,i) = vsp(:,n,isht-1:0,nk*2-2-i)
   end if
   if (sht .lt. 0) then
     vspi(:,n,0:no/2-isht,i) = vsp(:,n,isht:no/2,i)
     vspi(:,n,isht:no/2,nk*2-2-i) = vsp(:,n,0:no/2-isht,nk*2-2-i)
     vspi(:,n,0:isht-1,nk*2-2-i) = vsp(:,n,isht-1:0,i)
   end if
 end do
 end do

 vspa = dim_avg_n(vspi,1)
print(max(vspa))

;--------------------------------------- cn level
 cnlev = (/-200,-100,-50,-20,-10,-5,-2,-1,-0.5,-0.2,0,0.2,0.5,1,2,5,10,20,50,100,200/)
; cnlev = fspan(-0.5,0.5,21)
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)
 cnlev2= (/-16,-8,-4,-2,-1,1,2,4,8,16/)

 cnpc = (/5,10,20,30,40,50,60,70,80,90,100/)

;--------------------------------------- ctr
 wi = 40
 he = 40
 xf = 34 + (/0,1,0,1,0,1/)*(wi+10)
 yf = 10 + (/0,0,1,1,2,2/)*45+(/0,0,0,0,1,1/)*18

 mstr  = new(nfig,string)
 xstr  = new(nfig,string)
 ystr  = new(nfig,string)
 lbopt = new(nfig,string)

 mstr(:) = ""
 xstr(:) = "zonal wavenumber"
; xstr(nfig-1) = "year"
 ystr(:) = ""
 ystr(0) = "intrinsic frequency [cyc day:S:-1:N:]"

do ip=0, nfig-1

 ctr1 = cnshadeplot(wks,kwn,ome,vspa(ip,:,:)*1.e3,cnlev,cnfil,"")
; cnaddline(ctr1,wks,kwn,ome,vspa(ip,:,:),cnlev2,0,"D",0,"T",1,"C",1)
 setvalues ctr1
  "cnFillMode" : "RasterFill"
 end setvalues

 set_size(ctr1,xf(ip),yf(ip),wi,he)
 axis_range(ctr1,-20.5,20.5,0,1.0-0.5*dome,False,False)
 axis_str(ctr1,mstr(ip),xstr(ip),ystr(ip))

 cnfinalize(ctr1,wks)

 draw(ctr1)

end do

   label = create "lab" labelBarClass wks
     "lbOrientation" : "Horizontal"
     "lbBoxCount"   : dimsizes(cnfil)
     "lbFillColors" : cnfil
     "lbMonoFillPattern" : True
     "lbLabelsOn" : True
     "lbLabelStrings" : cnlev
     "lbLabelStride" : 1
;     "lbAutoManage" : False
     "lbLabelAlignment" : "InteriorEdges"
     "lbPerimOn" : False
     "lbTitleString" : "[ 10:S:-3:N: m s:S:-1:N: day:S:-1:N: / cyc day:S:-1:N: ]"
     "lbTitlePosition" : "Bottom"
   end create
   set_size_lb(label,xf(0),yf(0)+he+15,wi*2+10,3)
   draw(label)


 frame(wks)

end



