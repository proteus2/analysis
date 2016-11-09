load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 load "./param0.ncl"
 load "./param1.ncl"

 raname  = param0_ranames (param0_i_exp)
 expname = param0_expnames(param0_i_exp)
 dirname = param0_dirnames(param0_i_exp)
 var    = param0_vars(param0_i_var(0))
 vname  = param0_vnames (param0_i_exp,param0_i_var(0))
 fvname = param0_fvnames(param0_i_exp,param0_i_var(0))
 vstr   = param0_vstrs(param0_i_var(0))

 yyyymm1 = param0_yyyymm1
 yyyymm2 = param0_yyyymm2
 lat_b = param0_lat_b
 p     = param0_p

 nbv     = param0_nbv
 h_scale = param0_h_scale

 kmax = param1_kmax_draw
 omax = param1_omax_draw

 lzs = param1_lz_draw

 l_sym_asym = False  ; True
 l_normalize = False

 wks = cr_wks("fig/var-"+var+"_ko-p"+numeric2int(p,1)+"-y"+numeric2int(lat_b,1),600,0)
; gsn_define_colormap(wks,"prcp_3")
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.6,0.6,0.6)

;--------------------------------------- read data
 nfig = 1
 ne = dimsizes(vname)

 y1 = yyyymm1/100
 y2 = yyyymm2/100
 f_ko = "_ko"
 if ( l_sym_asym ) then
   f_ko = "_wave_ko"
   nfig = 2
 end if

 l_precal = False
 fdir   = "$DATD/fspec/"+dirname+"/fcoef_ko"
 fname0 = fdir+"/"+y1+"/"+expname+".fft_"+fvname+f_ko+"."+y1+".02.nc"
; if ( yyyymm1 .eq. 197902 .and. yyyymm2 .eq. 200612 ) then
 if ( yyyymm1 .eq. 198101 .and. yyyymm2 .eq. 201012 ) then
   l_precal = True
   fdir   = "$DATD/fspec/"+dirname+"/fcoef_ko/avg_81-10"
   fname0 = fdir+"/"+expname+".fc2_"+fvname+f_ko+".nc"
 end if

 no0 = new(ne,integer)
 do ie=0, ne-1
   f0 = addfile(fname0(ie),"r")
   no0(ie) = dimsizes( f0->ome_fr )
 end do

 no = min( no0 )
 nk = kmax + 1

 kwn = fspan(-kmax,kmax,kmax*2+1)
 ome = fspan(0,no/2,no/2+1) / 90.
 dome = 1./90.


do ie=0, ne-1

 f0 = addfile(fname0(ie),"r")

if ( l_precal ) then
 vname0f = "fc2_"+vname(ie)
else
 vname0f = "fcr_"+vname(ie)
end if
if ( l_sym_asym ) then
 vname0f = vname0f+"_s"
end if
 stmp = f0->$vname0f$!0
 print( "p asked: "+p+"  /  read: "+f0->$stmp$({p}) )

 vspec = new((/nfig,no0(ie),nk/),float)

if ( l_precal ) then

 fname = fdir(ie)+"/"+expname(ie)+".fc2_"+fvname(ie)+f_ko+".nc"
 f1 = addfile(fname,"r")
if ( .not. l_sym_asym ) then
 vspec = dim_avg_n(  \
     f1->$("fc2_"+vname(ie))$({p},{-lat_b:lat_b},:,:kmax), 0 )*2.
else 
 tmps = f1->$("fc2_"+vname(ie)+"_s")$({p},{-lat_b:lat_b},:,:kmax)
 tmpa = f1->$("fc2_"+vname(ie)+"_a")$({p},{-lat_b:lat_b},:,:kmax)
 vspec = dim_avg_n( (/ tmps, tmpa /), 1 )*2.
 delete( [/tmps,tmpa/] )
end if
; end do

else

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
     fname = fdir(ie)+"/"+iy+"/"+expname(ie)+".fft_"+fvname(ie)+f_ko+"."+ \
             iy+"."+sprinti("%2.2i",im)+".nc"
     f1 = addfile(fname,"r")

     tmp1r = f1->$("fcr_"+vname(ie))$({p},{-lat_b:lat_b},:,:kmax)
     tmp1i = f1->$("fci_"+vname(ie))$({p},{-lat_b:lat_b},:,:kmax)
if ( .not. l_sym_asym ) then
     vspec(0,:,:) = vspec(0,:,:) + dim_avg_n( tmp1r^2 + tmp1i^2, 0 )*2.
else
     tmp2r = tmp1r({lat_b:-lat_b},:,:)
     tmp2i = tmp1i({lat_b:-lat_b},:,:)
     tmp = dim_avg_n( \
           (/ (0.5*(tmp1r + tmp2r))^2 + (0.5*(tmp1i + tmp2i))^2, \
              (0.5*(tmp1r - tmp2r))^2 + (0.5*(tmp1i - tmp2i))^2 /), 1 )*2.
     vspec = vspec + tmp
     delete( [/tmp2r,tmp2i,tmp/] )
end if
     delete( [/tmp1r,tmp1i/] )
     n = n + 1
   end do
 end do
 nt = n
 vspec = vspec/tofloat(nt)

end if

if ( l_sym_asym ) then
 ; switch for v'
 if ( var .eq. "V" ) then
   vstmp = vspec(0,:,:)
   vspec(0,:,:) = vspec(1,:,:)
   vspec(1,:,:) = vstmp
   delete( vstmp )
 end if
end if

 if ( l_normalize .and. var .eq. "T" ) then
   fname2 = "/data18/kyh/dat/clim/"+dirname+"/"+expname+".t_yp.1979-2006.nc"
   f2 = addfile(fname2(ie),"r")
   t0 = dim_avg(f2->$vname(ie)$({p},{-lat_b:lat_b}))
   vspec = vspec / t0^2
   vstr = vstr+"/:F10:T:B:0:N:"
 end if
 if ( var .eq. "W" ) then
   vspec = vspec*(h_scale/(p*1.e2))^2 * 1.e6
   vstr = vstr+":F8:(:F10:H:B:s:N::F13:/:F10:p:F8:)"
 end if
 
;--------------------------------------- process 0
 if ( .not. isdefined("vsp") ) then
   vsp = new((/nfig,ne,no/2+1,nk*2-1/),float)
   vsp = 0.
 end if
 vsp(:,ie,:,0:nk-2) = vspec(:,0:no/2,nk-1:1)
 vsp(:,ie,0,nk:) = vspec(:,0,1:);*0.5
 vsp(:,ie,1:no/2,nk:) = vspec(:,no0(ie)-1:no0(ie)-no/2,1:)
 vsp(:,ie,no/2,:) = vsp(:,ie,no/2,:)*0.5

 ; k = 0
 vsp(:,ie,:,nk-1) = 0.125*(vspec(:,0:no/2,0)+vspec(:,no0(ie)-1:no0(ie)-no/2-1,0))
 vsp(:,:,0,nk-1) = vsp@_FillValue

 delete(vspec)

end do  ; ie


 vsp = vsp/dome  ; PS density

; vsp(:,:,:,nk-1) = vsp@_FillValue  ; k = 0

print(max(vsp))

 nw = 3

 nk10 = nk*2*10
 no10 = no/2*10

 kwn2 = fspan(-nk,nk,nk10+1)
 ome2 = fspan(0,no/2,no10+1) / 90.

 dispc = disp_eq(kwn2,ome2,False,nbv)

;--------------------------------------- cn level
; cnfil = (/5,6,7,10,12,14,15,16,17,19,21/)
 cnfil = (/26,24,21,19,17,16,14,12,2,5,11/)

 fact = 0
 if ( var .eq. "T" ) then
   cnlev = (/0.02,0.05,0.1,0.2,0.5,1,2,5,10,20/)*0.1
   lbtit = "K:S:2:N:"
   if ( l_normalize ) then
     fact = -6
     cnlev = (/0.1,0.2,0.5,1,2,5,10,20,50,100/)
     lbtit = ""
   end if
 end if
 if ( var .eq. "U" ) then
   cnlev = (/0.02,0.05,0.1,0.2,0.5,1,2,5,10,20/)
   lbtit = "m:S:2:N: s:S:-2:N:"
 end if
 if ( var .eq. "V" ) then
   cnlev = (/0.02,0.05,0.1,0.2,0.5,1,2,5,10,20/)
   lbtit = "m:S:2:N: s:S:-2:N:"
 end if
 if ( var .eq. "W" ) then
   cnlev = (/0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10/)
   lbtit = "mm:S:2:N: s:S:-2:N:"
 end if
 lbt_hd = ""
 if (fact .ne. 0) then
   lbt_hd = "10:S:"+fact+":N: "
 end if
 lbtit = "[ "+lbt_hd+lbtit+" / cyc day:S:-1:N: ]"

 ldsh = (/0,16,14/)

;--------------------------------------- ctr
 wi = 32
 he = 32
 xf = 12 + (/0,1,2,0,1,2/)*(wi+8)
 yf = 10 + (/0,0,0,1,1,1/)*(he+20)

 mstr = new((/nfig,ne/),string)
 xstr = new(nfig,string)
 ystr = new(ne,string)

 mstr = ""
 mstr(0,:) = raname
 xstr = ""
 xstr(nfig-1) = "zonal wavenumber"
 ystr = ""
 ystr = where(xf(:ne-1) .eq. xf(0), "frequency [cyc day:S:-1:N:]", ystr)

 l0 = 1.0
 resp = True
 resp@gsFillColor = 1

do ip=0, nfig-1
do ie=0, ne-1

 ctr1 = cnshadeplot(wks,kwn,ome,vsp(ip,ie,:,:)/(10.^fact),cnlev,cnfil,"")
; cnaddline(ctr1,wks,kwn,ome,vsp(ip,ie,:,:)/(10.^fact),cnlev,0,"D",0,"T",1,"C",1)
if ( l_sym_asym ) then
 do iw=0, nw-1
   cnaddline(ctr1,wks,kwn2,ome2,dispc(ip,iw,:,:),lzs(ip,iw,:),0,"D",ldsh(iw),"T",1.2,"C",1)
 end do
end if
; setvalues ctr1
;  "cnFillMode" : "RasterFill"
; end setvalues
 dum = gsn_add_polygon(wks,ctr1,(/-l0,l0,l0,-l0/),(/0,0,omax,omax/),resp)

 set_size(ctr1,xf(ie),yf(ie),wi,he)
 axis_range(ctr1,-kmax,kmax,0,omax,False,False)
 axis_str(ctr1,mstr(ip,ie),xstr(ip),ystr(ie))

 cnfinalize(ctr1,wks)

 draw(ctr1)

end do
end do

 label = create "lab" labelBarClass wks
   "lbOrientation" : "Horizontal"
   "lbBoxCount"   : dimsizes(cnfil)
   "lbFillColors" : cnfil
   "lbMonoFillPattern" : True
   "lbAutoManage"   : False
   "lbLabelsOn"     : True
   "lbLabelStrings" : cnlev
   "lbLabelStride"  : 1
   "lbLabelAlignment" : "InteriorEdges"
   "lbPerimOn"        : False
   "lbTitleString"    : lbtit
   "lbTitlePosition"  : "Right"
   "lbTitleDirection" : "Across"
 end create
 set_size_lb(label,xf(0),yf(ne-1)+he+15,wi*5,4)
 draw(label)

 frame(wks)

end

