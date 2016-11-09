load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = "uj"
 varname = "f_y"

 y1 = 1953
 y2 = 2005
 lat0 = 5.
 z0 = 17.41 ;18.59

 kmax = 40
 omax = 3.2

 wks = cr_wks("fig/fy_reswav_ko-z"+numeric2int(z0,1)+"-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 nfig = 2

 fdir  = "/data4/kyh/dat/L60CGW-t/epf"
 fname = fdir+"/"+y1+"/uan"+expname+".epf_koyz."+y1+".01.nc"
 f1 = addfile(fname,"r")
 no = dimsizes( f1->ome_fr )
 nk = dimsizes( f1->k_wn   )

 vspec = new((/nfig,no,nk/),float)
 vspec = 0.
 n = 0
 do iy=y1, y2
 do im=1, 12
   fname = fdir+"/"+iy+"/uan"+expname+".epf_koyz."+iy+"."+sprinti("%2.2i",im)+".nc"

   f1 = addfile(fname,"r")
   tmp = dim_avg_n( f1->$varname$({z0*1.e3},{-lat0:lat0},:,:), 0 )*(-1.)
   tmp(:,:nk-2) = tmp(:,:nk-2)*2.
;   tmp = abs(tmp)

   vspec(0,:,:) = vspec(0,:,:) + tmp

   n = n + 1
 end do
 end do
 nt = n
 vspec(0,:,:) = vspec(0,:,:)/nt

;--------------------------------------- process 0
 kwn = fspan(-nk+1,nk-1,nk*2-1)
 ome = fspan(0,no/2,no/2+1) / 90.
 dome = 1./90.

 vsp = new((/nfig,no/2+1,nk*2-1/),float)
 vsp = 0.
 do i=0, nk-2
   vsp(0,:,i) = vspec(0,0:no/2,nk-1-i)
 end do
 vsp(0,0,nk:) = vspec(0,0,1:);*0.5
 do n=1, no/2
   vsp(0,n,nk:) = vspec(0,no-n,1:)
 end do
 vsp(0,no/2,:) = vsp(0,no/2,:)*0.5

 vsp(1,:,:) = vsp(0,:,:)

 if (kmax .ne. 999) then
   vsp(0,:,nk-1-kmax) = dim_sum( vsp(0,:,:nk-1-kmax) )
   vsp(0,:,nk-1+kmax) = dim_sum( vsp(0,:,nk-1+kmax:) )
   vsp(0,:,:nk-1-kmax-1) = 0.
   vsp(0,:,nk-1+kmax+1:) = 0.
   kmax = kmax + 1
 end if
 if (omax .ne. 999) then
   nn = toint(omax*90.)
   vsp(0,nn,:) = dim_sum_n( vsp(0,nn:,:), 0 )
   vsp(0,nn+1:,:) = 0.
   omax = omax + 1./90.
 end if

 vsp = vsp/dome

 vsp(:,:,nk-1) = vsp@_FillValue
print(max(vsp))

;--------------------------------------- cn level
 cnlev = (/-100,-50,-20,-10,-5,-2,-1,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,1,2,5,10,20,50,100/)
; cnlev = fspan(-0.5,0.5,21)
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)
 cnlev2= (/-16,-8,-4,-2,-1,1,2,4,8,16/)

;--------------------------------------- ctr
 wi = 40
 he = 40
 xf = 34 + (/0,1/)*(wi+10)
 yf = 10 + (/0,0/)*28

 mstr  = new(nfig,string)
 xstr  = new(nfig,string)
 ystr  = new(nfig,string)
 lbopt = new(nfig,string)

 mstr(:) = ""
 xstr(:) = "zonal wavenumber"
; xstr(nfig-1) = "year"
 ystr(:) = ""
 ystr(0) = "frequency [cyc day:S:-1:N:]"

do ip=0, nfig-1

 ctr1 = cnshadeplot(wks,kwn,ome,vsp(ip,:,:)/1.e4,cnlev,cnfil,"")
; cnaddline(ctr1,wks,kwn,ome,vsp(ip,:,:),cnlev2,0,"D",0,"T",1,"C",1)
 setvalues ctr1
  "cnFillMode" : "RasterFill"
 end setvalues

 set_size(ctr1,xf(ip),yf(ip),wi,he)
 axis_range(ctr1,tofloat(-kmax),tofloat(kmax),0,omax+0.5*dome,False,False)
 if (ip .eq. 1) then
   axis_range(ctr1,-20.5,20.5,0,1.0+0.5*dome,False,False)
 end if
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
     "lbTitleString" : "[ 10:S:4:N: Pa m / cyc day:S:-1:N: ]"
     "lbTitlePosition" : "Bottom"
   end create
   set_size_lb(label,xf(0),yf(0)+he+20,wi*2+10,3)
   draw(label)

 frame(wks)

end



