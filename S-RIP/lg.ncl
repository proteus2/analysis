load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 0.8

begin

;--------------------------------------- parameters
 load "./param0.ncl"
 load "./ftn_get_color_S-RIP.ncl"
 load "./ftn_read_spec.ncl"

 ra = param0_ra_p
; ra = array_append_record(param0_ra_p,param0_ra_m,0)

 ranames = get_ra_names(ra)
 
 wks = cr_wks("fig/lg",600,0)
 icol_s = get_col_srip(wks)

;--------------------------------------- read data
 ltck = (/1,1,1,1,1,1,1,1/)*0.8
 ldsh = (/0,0,0,0,0,14,14,14/)
 lcol = (/ icol_s@ei , icol_s@m , icol_s@m2 , icol_s@c , icol_s@j5 ,  \
           icol_s@ei , icol_s@m2 , icol_s@j5 /)

;--------------------------------------- ctr
 lg = lgplot(wks,ranames(:,1),"D",ldsh,"T",ltck,"C",lcol,2)
 set_size_lg(lg,10,10,14,2.8*5)
 draw(lg)

 frame(wks)

end

