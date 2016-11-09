load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = "uj"
 varname = (/"f_z"/);,"f_y"/)

 yyyymm1 = 195301
 yyyymm2 = 200611
 lat0 = 5.
 z0 = 17.41 ;18.59

 kmax = 20
 omax = 1.0
 kmax0 = 40
 omax0 = 3.2

 wks = cr_wks("fig/fz_reswav_ko-z"+numeric2int(z0,1)+"-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 nfig = 2
 nv = dimsizes(varname)

 y1 = yyyymm1/100
 y2 = yyyymm2/100

 fdir  = "/prime0/kyh/dat/L60CGW/epf"
 fname = fdir+"/"+y1+"/uan"+expname+".epf_koyz."+y1+".01.nc"
 f1 = addfile(fname,"r")
 no = dimsizes( f1->ome_fr )
 nk = dimsizes( f1->k_wn   )

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
     do iv=0, nv-1
       fname = fdir+"/"+iy+"/uan"+expname+".epf_koyz."+ \
               iy+"."+sprinti("%2.2i",im)+".nc"

       f1 = addfile(fname,"r")

       tmp = dim_avg_n( \
             ( f1->$(varname(iv))$({z0*1.e3},{-lat0:lat0},:,:) ), 0 )*(-1.)
       tmp(:,:nk-2) = tmp(:,:nk-2)*2.

       vspec(iv,:,:) = vspec(iv,:,:) + tmp
     end do
     n = n + 1
   end do
 end do
 nt = n
 vspec(:,:,:) = vspec(:,:,:)/nt

 fname2 = "$DATD/L60CGW/cmxm.uan"+expname+"a.pj.nc"
 f2 = addfile(fname2,"r")

 rho0 = dim_avg(f2->rho({z0*1.e3},{-lat0:lat0}))

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
 vars = ":F10:"+(/"-F:B:z:N:","-F:B:y:N:"/)
 facts = (/1,4/)

 cnlevs = (/ (/-200,-100,-50,-20,-10,-5,-2,-1,-0.5,-0.2, \
               0,0.2,0.5,1,2,5,10,20,50,100,200/), \
             (/-100,-50,-20,-10,-5,-2,-1,-0.5,-0.2,-0.1, \
               0,0.1,0.2,0.5,1,2,5,10,20,50,100/) /)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

 var   = new(nv,string )
 fact  = new(nv,integer)
 cnlev = new((/nv,dimsizes(cnlevs(0,:))/),float)
 lbtit = new(nv,string )
 do iv=0, nv-1
   if ( varname(iv) .eq. "f_z" ) then
     var  (iv) = vars (0)
     fact (iv) = facts(0)
     cnlev(iv,:) = cnlevs(0,:)
     lbtit(iv) = "[ 10:S:"+fact(iv)+":N: Pa m / cyc day:S:-1:N: ]"
   end if
   if ( varname(iv) .eq. "f_y" ) then
     var  (iv) = vars (1)
     fact (iv) = facts(1)
     cnlev(iv,:) = cnlevs(1,:)
     lbtit(iv) = "[ 10:S:"+fact(iv)+":N: Pa m / cyc day:S:-1:N: ]"
   end if
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
   set_size_lb(label,xf(max((/0,nv-2/))),yf(nfig-1)+he+20+iv*8,wi*2+8,3)
   draw(label)
 end do

 frame(wks)

end

