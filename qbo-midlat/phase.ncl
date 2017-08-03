load "$UTIL/ncl/header.ncl"
xspace =  83  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"
axisfontratio = 1.2
begin

;--------------------------------------- parameters
 load "./param-qbodata-obs.ncl"
 load "./param-1-ra.ncl"
 load "./ftn_qbo_phase.ncl"

 season = par1_season
 y1 = par1_y1
 y2 = par1_y2

 p1 = par0_p_q(0)
 p2 = par0_p_q(1:)
 if (par0_case_q .ne. "obs") then
   lat0 = par0_latbnd_q
 end if

 y0 = par0_y0_q

;-------------------------------------------------
 nma = par1_nma

;--------------------------------------- read data
 nt = (y2-y1+1)*nma

 it0 = (y1-y0)*12
 it9 = (y2+1-y0)*12 - 1

 its = new(nt,integer)
 do n=0, nma-1
   its(n::nma) = ispan(it0+par1_imon0(par1_iseason,n),it9,12)
 end do

 f1 = addfile(par0_fname_q,"r")
 tmpi = dimsizes( f1->$par0_varname_q$ )
 if (dimsizes(tmpi) .eq. 2) then
   if (tmpi(0) .lt. tmpi(1)) then
     u = f1->$par0_varname_q$(:,its)
   else
     u = transpose( f1->$par0_varname_q$(its,:) )
   end if
 else
   u = transpose( dim_avg( f1->$par0_varname_q$(its,:,{-lat0:lat0}) ) )
 end if
 p = u&$u!0$
 z = 7.*log(1.e3/p)

;--------------------------------------- process 0
 it_phs = qbo_phase_ind(season,par0_fname_q,par0_varname_q,  \
                        p1,p2,10,par0_y0_q,y1,y2)
;                        p1,p2,10,par0_y0_q,par0_y1_q,par0_y2_q)
 nph = it_phs@n_phase
 nt_phs = it_phs@nt_phase
 phs = it_phs@phases
 ua1 = it_phs@ua1
 ua2 = it_phs@ua2
print(nt_phs)

;--------------------------------------- wks
 wks = cr_wks("fig/phase_"+season+"_p"+p1+"-"+p2(0)+"_"+y1+"-"+y2,600,0)
 colind = NhlNewColor(wks,0.5,0.5,0.5)

;--------------------------------------- cn level
 mcol = (/colind,2,3,4,6,7,8,13,5/)

;--------------------------------------- plot
 mstr = (/"",""/)
 xstr = (/"",""/)
 ystr = (/"",""/)

 wi = 30
 he = 30
 xf = 15 + ispan(0,1,1)*(wi+15)
 yf = 10 + ispan(0,1,1)*(he+15)

 xmnmx = (/ toint(min(ua1))-2, toint(max(ua1))+2 /) + (/-4,3/)
 ymnmx = (/ toint(min(ua2))-2, toint(max(ua2))+2 /) + (/-5,5/)


 xstr(0) = "U"+p1+"a [m s:S:-1:N:]"
 ystr(0) = "U"+p2(0)+"a [m s:S:-1:N:]"
 if (dimsizes(p2) .eq. 2) then
   ystr(0) = "U"+p2(0)+" - U"+p2(1)+"a [m s:S:-1:N:]"
 end if

 plot = xylineplot(wks,ua1,ua2,"D",0,"T",0.5,"C",-1)
 xyaddmarker(plot,16,"S",0.2,"T",1.0,"C",1)

 dumm = new(nt,graphic)
 plres = True
 plres@gsMarkerIndex = 16
 plres@gsMarkerSizeF = 0.001
 do n=0, nt-1
   plres@gsMarkerColor = mcol(phs(n))
   dumm(n) = gsn_add_polymarker(wks,plot,ua1(n),ua2(n),plres)
 end do

 set_size(plot,xf(0),yf(0),wi,he)
; axis_range(plot,-999,999,-999,999,False,False)
 axis_range(plot,xmnmx(0),xmnmx(1),ymnmx(0),ymnmx(1),False,False)
 axis_str(plot,mstr(0),xstr(0),ystr(0))

; dr_txt(wks,p(ip)+" hPa",xf(0)+8,yf(ip)+3,22,8,0.)

 draw(plot)


 xstr(1) = "U [m s:S:-1:N:]"
 ystr(1) = "p [hPa]"

 plot2 = xylineplot(wks,0,z,"D",0,"T",0.5,"C",colind)
 do il=0, nph-1
   u_z = dim_avg(u(:,it_phs(il,:nt_phs(il)-1)))
   xyaddline(plot2,u_z,z,"D",0,"T",1.0,"C",mcol(il+1))
 end do

 set_size(plot2,xf(1),yf(0),wi,he)
; axis_range(plot2,-999,999,-999,999,False,False)
 axis_range(plot2,-40,30,min(z),max(z),False,False)
 axis_str(plot2,mstr(1),xstr(1),ystr(1))
 setvalues plot2
  "tmYLMode"   : "Explicit"
  "tmYLValues" : z
  "tmYLLabels" : p
 end setvalues

 draw(plot2)


 frame(wks)

end

