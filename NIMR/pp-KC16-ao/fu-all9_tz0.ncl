load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.85

begin

;--------------------------------------- parameters
 cexp    = (/"CC"/)
 expname = (/"tc"/)
 expdir  = (/"AOL60CGW"/)
 ppid    = (/"j" /)
 varname = "u"

 y0 = 1951
 y1 = 1970 ;2070 ;1970
 y2 = 1989 ;2089 ;1989
 lat0 = 10.

 wks = cr_wks("fig/fu-all9_tz0-"+y1+"-"+y2+"-"+expname(0),600,0)
 gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "$DATD/"+expdir(0)+"/epf-wc/uan"+expname(0)+".epfsa_yz_recon0.nc"
   fnam2 = "$DATD/"+expdir(0)+"/epf-wc/uan"+expname(0)+".epf_yz_recon0.nc"

   f1 = addfile(fnam1,"r")
   f2 = addfile(fnam2,"r")
   z  = f1->z / 1.e3

   nz = dimsizes(z)
   nw = 8

   nt = (y2 - y1 + 1)*12 ;647 + 1
   n1 = (y1 - y0)*12 - 1  ; Jan. of y0 is absent.
   n2 = n1 + nt - 1

;   rw_i = new((/nt-1,nw,nz/),float)
   rw_i = new((/nt,nw,nz/),float)
   rw_i(:,:,:) = dim_avg( f1->epd(n1:n2,:nw-1,:,{-lat0:lat0}) )
   rw_i(:,nw-2:,:) = rw_i(:,nw-2:,:) + \
                     dim_avg( f2->epd(n1:n2,8:9,:,{-lat0:lat0}) )

; migrating tide
;   rw_i(:,nw-1,:) = rw_i(:,nw-1,:) + \
;                    ( dim_avg( f2->epd(n1:n2,2,:,{-lat0:lat0}) ) + \
;                      dim_avg( f2->epd(n1:n2,3,:,{-lat0:lat0}) ) )

   rw_i(:,3,:) = rw_i(:,3,:) + dim_avg( f1->epd(n1:n2,13,:,{-lat0:lat0}) )
   rw_i(:,6:7,:) = rw_i(:,6:7,:) + dim_avg( f1->epd(n1:n2,16:17,:,{-lat0:lat0}) )

;   rw0 = new((/nt,nw,nz/),float)
;   rw0(:nt-2,:,:) = rw_i(:nt-2,:,:)
;   rw0(nt-1,:,:) = rw0(nt-2,:,:)
   rw0 = rw_i

   n1 = (y1 - y0)*12
   n2 = n1 + nt - 1

   z1 = f1->z(0)
   fnamg = "$DATD/"+expdir(0)+"/xm.uan"+expname(0)+"a.pe.nc"
   fg = addfile(fnamg,"r")
   zg = fg->hybrid_ht_1 / 1.e3

   gw = new((/nt,nz/),float)
   gw(:,:) = dim_avg( fg->du_gw      (n1:n2,{z1:},{-lat0:lat0}) )*72. + \
             dim_avg( fg->dudt_gwconv(n1:n2,{z1:},{-lat0:lat0}) )*86400.


   fnam1 = "$DATD/"+expdir(0)+"/xm.uan"+expname(0)+"a.pj.nc"

   f1 = addfile(fnam1,"r")
   zu = f1->hybrid_ht / 1.e3

   nzu = dimsizes(zu)

   u0 = transpose( dim_avg( f1->u(n1:n2,:,{-lat0:lat0}) ) )

   fnam2 = "$DATD/"+expdir(0)+"/tem/uan"+expname(0)+".tem-wc_yzt.nc"
   f2 = addfile(fnam2,"r")
   za = f2->z / 1.e3

   ut = transpose( dim_avg_Wrap( f2->u_tend(n1:n2,:,{-lat0:lat0}) ) )

   advz = transpose( dim_avg_Wrap( f2->uadv_z(n1:n2,:,{-lat0:lat0}) ) )

   uf = transpose( dim_avg_Wrap( f2->u_force(n1:n2,:,{-lat0:lat0}) ) )

   ue = transpose( dim_avg_Wrap( f2->epd(n1:n2,:,{-lat0:lat0}) ) )

   advy = transpose( dim_avg_Wrap( f2->uadv_y(n1:n2,:,{-lat0:lat0}) ) )
   advy = advy + transpose( dim_avg_Wrap( f2->cor(n1:n2,:,{-lat0:lat0}) ) )

   ureq = ut

   ureq = (/ ut - advz - advy /)

   wf = (/ rw0(:,2,:), rw0(:,3,:), rw0(:,4,:)+rw0(:,5,:), \
           rw0(:,6,:)+rw0(:,7,:), gw, gw /)

   fu = new((/6,nz,nt/),float)
   do n=0, nt-1
     fu(:,:,n) = wf(:,n,:)
   end do
print(dimsizes(fu))
print(dimsizes(ureq))

   fu!1 = "z"
   fu&z = z*1.e3

   fu(5,{12.e3:},:) = ureq({12.e3:},:) - ue({12.e3:},:) - fu(5,{12.e3:},:)

   t = y1 + (fspan(1,nt,nt)-0.5)/12.

;   tt = y0 + (fspan(1,nt-1,nt-1))/12.

;--------------------------------------- process 0
 fu = fu*30.
 ut = ut*30.
 advz = advz*30.
 ureq = ureq*30.

 u3 = (/ut, advz, ureq/)

;   uss = new((/2,nz,12/),float)
;   do n=0, 11
;     uss(:,:,n) = dim_avg(u(:,:,ispan(n,nt-1,12)))
;   end do
;
;   udes = new(dimsizes(u),float)
;   do n=0, nt-1
;     udes(:,:,n) = u(:,:,n) - uss(:,:,mod(n,12))
;   end do

;--------------------------------------- cn level
 cscale = (/4,4,8,8,1,1/)

 cnlev = fspan(-10,10,21)*2.
; cnfil = (/2,3,4,5,6,7,8,9,10,11,0,0,26,25,24,23,22,21,20,19,18,17/)
 cnfil = (/17,18,19,20,21,22,23,24,25,26,0,0,11,10,9,8,7,6,5,4,3,2/)

 cnlevu = fspan(-10,10,3)

;--------------------------------------- contour
 nl = 9

 wi = 120
 he = 15
 xf = 15 + (/0,1/)*(wi+5.)
 yf = 10 + ispan(0,nl-1,1)*(he+5)

 mstr  = new(nl,string)
 xstr  = new(nl,string)
 lbopt = new(nl,string)

 mstr(:) = ""
 mstr(0) = cexp
 xstr(:) = ""
 xstr(nl-1) = "year"
 ystr = (/"height (km)",""/)
 lbopt(:) = ""
 lbopt(nl-1) = "(m s:S:-1:N: month:S:-1:N:)"
 txtstr = (/":F18:s:F22:U/:F18:s:F22:t","ADVz","REQ","Kelvin","Rossby","MRG","IG","CGW+BGW", \
            "D:B:u:H-10V15:_:N:"/)
 subl = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)"/)
 do il=3, nl-1
   if (cscale(il-3) .ne. 1) then
     txtstr(il) = txtstr(il)+" (:F34:4:F22: "+cscale(il-3)+")"
   end if
 end do

do ie=0, 0
do il=0, nl-1

 if (il .lt. 3) then
   contour1 = cnshadeplot(wks,t,za,u3(il,:,:),cnlev,cnfil,"")
 else
   contour1 = cnshadeplot(wks,t,z,fu(il-3,:,:)*cscale(il-3),cnlev,cnfil,"")
 end if
 setvalues contour1
  "pmLabelBarOrthogonalPosF" : 0.12 + 0.2
  "pmLabelBarParallelPosF"   : 0.50 ;+ 0.55
  "pmLabelBarKeepAspect"     : False
  "pmLabelBarWidthF"         : 0.55
  "pmLabelBarHeightF"        : 0.20
  "lbLabelStride"      : 2
  "lbLabelFontHeightF" : 0.012
  "lbTitlePosition"    : "Right"
  "lbTitleDirection"   : "Across"
  "lbTitleOffsetF"     : 0.05 + 0.05
  "lbTitleFontHeightF" : 0.015
 end setvalues

 set_size(contour1,xf(ie),yf(il),wi,he)
 axis_range(contour1,y1,y2+1,19,33,False,False)
 axis_str(contour1,mstr(il),xstr(il),ystr(ie))
 setvalues contour1
  "tmXBMode" : "Manual"
  "tmXBTickSpacingF"  : 1
  "tmXBMinorPerMajor" : 0
  "tiMainFontHeightF" : sz_mstr*axisfontratio*1.25
 end setvalues

 contouru = cnlineplot(wks,t,zu,u0,cnlevu,0)
 cnline_att(contouru,"D",(/0,0,1/),"T",0.6,"C",1)

 set_size(contouru,xf(ie),yf(il),wi,he)
 axis_range(contouru,y1,y2+1,19,33,False,False)
 axis_str(contouru,"","","")
 setvalues contouru
  "tmXBOn" : False
  "tmXTOn" : False
 end setvalues

 cnfinalize(contour1,wks)

 draw(contour1)

 cnfinalize(contouru,wks)

 draw(contouru)

 dr_txt2(wks,txtstr(il),xf(ie)+1,yf(il)+he-2,22,6,0.,"CenterLeft",1)
 if (ie .eq. 0) then
   dr_txt(wks,subl(il),xf(0)-10,yf(il),25,8,0.)
 end if

 delete(contour1)
 delete(contouru)

end do
end do


 frame(wks)

end

