load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 fvarname = (/"u","v","omega","t"/)

 lat0 = 15.
 p0 = 70

 wks = cr_wks("fig/v4_ko-shading-p"+numeric2int(p0,1)+"-y"+numeric2int(lat0,1),600,0)
 gsn_define_colormap(wks,"prcp_3")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
 nv = dimsizes(fvarname)

;--------------------------------------- cn level
 vars = (/":F10:u:F13:'",":F10:v:F13:'",":F10:-:F8:w:F13:':F8:(:F10:H:B:s:N::F13:/:F10:p:F8:)",":F10:T:F13:'/:F10:T:B:0:N:"/)
 facts = (/0,0,0,-5/)
 cnlevs = (/ (/0.02,0.05,0.1,0.2,0.5,1,2,5,10,20/), \
             (/0.02,0.05,0.1,0.2,0.5,1,2,5,10,20/), \
             (/0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10/), \
             (/0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10/) /)

 cnfil = (/5,6,7,10,12,14,15,16,17,19,21/)

 var   = new(nv,string )
 fact  = new(nv,integer)
 cnlev = new((/nv,dimsizes(cnlevs(0,:))/),float)
 lbtit = new(nv,string )
 do iv=0, nv-1
   if ( fvarname(iv) .eq. "u" ) then
     var  (iv) = vars (0)
     fact (iv) = facts(0)
	 cnlev(iv,:) = cnlevs(0,:)
     lbtit(iv) = "m:S:2:N: s:S:-2:N: / cyc day:S:-1:N: ]"
   end if
   if ( fvarname(iv) .eq. "v" ) then
     var  (iv) = vars (1)
     fact (iv) = facts(1)
     cnlev(iv,:) = cnlevs(1,:)
     lbtit(iv) = "m:S:2:N: s:S:-2:N: / cyc day:S:-1:N: ]"
   end if
   if ( fvarname(iv) .eq. "omega" .or. fvarname(iv) .eq. "ome" ) then
     var  (iv) = vars (2)
     fact (iv) = facts(2)
     cnlev(iv,:) = cnlevs(2,:)
     lbtit(iv) = "mm:S:2:N: s:S:-2:N: / cyc day:S:-1:N: ]"
   end if
   if ( fvarname(iv) .eq. "t" ) then
     var  (iv) = vars (3)
     fact (iv) = facts(3)
     cnlev(iv,:) = cnlevs(3,:)
     lbtit(iv) = "/ cyc day:S:-1:N: ]"
   end if
   lbt_hd = "[ 10:S:"+fact(iv)+":N: "
   if (fact(iv) .eq. 0) then
     lbt_hd = "[ "
   end if
   lbtit(iv) = lbt_hd+lbtit(iv)
 end do

;--------------------------------------- ctr
 wi = 32
 he = 32
 xf = 12 + (/0,1,2,3/)*(wi+8)
 yf = 10

 lbtit(2) = lbtit(2)+" ,  "+lbtit(3)
 do iv=0, nv-1, 2
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
;   set_size_lb(label,xf(0),yf+he+18+iv*8,wi*2+8,3)
   set_size_lb(label,xf(iv),yf,wi*2+8,3)
   draw(label)
 end do

 frame(wks)

end

