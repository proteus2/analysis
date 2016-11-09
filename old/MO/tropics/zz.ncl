load "$UTIL/ncl/header.ncl"

xspace = 170  ; 83, 170  -  QJRMS
yspace = 240  ; 240
load "$UTIL/ncl/figsize.ncl"
load "$UTIL/ncl/attribute.ncl"

axisfontratio = 1.2

begin

;--------------------------------------- parameters
 expname = (/"fd","y"/)

 nmon = 4+144
 lat0 = 5.
 lat1 = 20.
 kz  = (/59,55/)
;  100,50,30,10, 3, 1,0.3
;   47,53,55,59,62,64,66

 wks = cr_wks("zz_t-"+expname(0),600,0)
; gsn_define_colormap(wks,"StepSeq25")
 colind = NhlNewColor(wks,0.4,0.4,0.4)

;--------------------------------------- read data
   fnam1 = "../dat/u_"+expname(0)+".nc"

   f1 = addfile(fnam1,"r")
   u1 = dim_avg( f1->u(4:nmon-1,kz,{-lat0:lat0},0) )
   z1 = f1->hybrid_ht(kz)/1.e3

   fnam2 = "../dat/u_"+expname(1)+".nc"

   f2 = addfile(fnam2,"r")
   u2 = dim_avg( f2->u(4:nmon-1,kz,{-lat0:lat0},0) )

   fnam1 = "../dat/cg_"+expname(0)+".nc"
   f1 = addfile(fnam1,"r")
   g1 = dim_avg( f1->unspecified(4:nmon-1,kz,{-lat1:lat1},0) )*86400.*30.
   fnam1 = "../dat/ng_"+expname(0)+".nc"
   f1 = addfile(fnam1,"r")
   g1 = g1 + dim_avg( f1->field424(4:nmon-1,kz,{-lat1:lat1},0) )*86400.*30.

   fnam2 = "../dat/ng_"+expname(1)+".nc"
   f2 = addfile(fnam2,"r")
   g2 = dim_avg( f2->field424(4:nmon-1,kz,{-lat1:lat1},0) )*86400.*30.


   nz = dimsizes(z1)
   nt = nmon-4

   u = new((/4,nz,nt/),float)
   do n=0, nt-1
     u(:,:,n) = (/u2(n,:),u1(n,:),g2(n,:),g1(n,:)/)
   end do

   t = fspan(1,nt,nt)

   uss = new((/4,nz,12/),float)
   do n=0, 11
     uss(:,:,n) = dim_avg(u(:,:,ispan(n,nt-1,12)))
   end do
   umean = dim_avg(uss)

   udes = new(dimsizes(u),float)
   do n=0, nt-1
     udes(:,:,n) = u(:,:,n) - uss(:,:,mod(n,12))
   end do

   ifilt = nt/24-1

   fftu = ezfftf(udes)
   fftu(:,:,:,ifilt+1:) = 0.
   ulow = ezfftb(fftu,0.)

   udes = udes - ulow
iev = 6
print(esccr(udes(3,:,:)-runave(udes(3,:,:),iev,0),udes(1,:,:)-runave(udes(1,:,:),iev,0),3))

print(esccr(udes(2,:,:)-runave(udes(2,:,:),iev,0),udes(0,:,:)-runave(udes(0,:,:),iev,0),3))

;--------------------------------------- process 0

;--------------------------------------- plot
 xmax = nmon-4
 ymin = (/-45,-45/)
 ymax = (/45,45/)
 wi = 72
 he = 30
 xf = 12 + (/0,1,0/)*78
 yf = 10 + (/0,1,2,3,4/)*36

 mstr = new((/nz,2/),string)
 mstr(:,:) = ""
 mstr(0,:) = (/"CTL","SC"/)
 xstr = new(nz,string)
 xstr(:) = ""
 xstr(nz-1) = "month"
 ystr = (/"(m s:S:-1:N:), (m s:S:-1:N: mon:S:-1:N:)",""/)


 do k=0, nz-1
 do iv=0, 1

     plot = xylineplot(wks,t,0.,"D",0,"T",0.5,"C",1)
     xyaddline(plot,t,udes(iv,k,:),"D",0,"T",0.5,"C",1)
     xyaddline(plot,t,udes(iv+2,k,:),"D",0,"T",0.5,"C",2)

   set_size(plot,xf(iv),yf(k),wi,he)
   axis_range(plot,1,xmax,ymin(k),ymax(k),False,False)
   axis_str(plot,mstr(k,iv),xstr(k),ystr(iv))
   setvalues plot
    "tmXBMode" : "Manual"
    "tmXBTickSpacingF"  : 12
    "tmXBTickStartF"    : 0
    "tmXBMinorPerMajor" : 11
   end setvalues

   dr_txt(wks,floattoint(z1(k))+" km",xf(0)+8,yf(k)+3,22,8,0.)

   draw(plot)

  end do
  end do


 frame(wks)

end



