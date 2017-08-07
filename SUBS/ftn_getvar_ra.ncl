
;  getvar_ra( raname:string, varname:string, year, month, hour,  \
;            p_rng:numeric, lat_rng:numeric, lon_rng:numeric,  \
;            opt_havg:integer, opt_tavg:integer )
;
;    raname :  ERA-I, ...
;    varname :  T, U, V, W, Z, ...
;    year, month, hour :  array or scalar, integer or string
;    p_rng, XXX_rng :  p, lat, and lon array or scalar to read. If dimsizes = 2,
;                      the two values are used as upper/lower bounds.
;    opt_havg :  (/opt_lat,opt_lon/), options for horizontal averaging.
;                opt_lon = 1 or 0 (average or not). opt_lat can be 1 or 2,
;                (simple or cos[lat]-weighted averaging), only when opt_lon = 1.
;    opt_tavg :  option for temporal averaging.
;                0 : non-averaging
;                1 : daily averaging ( e.g., X(00)/4 + X(12)/2 + X(24)/4 )
;                2 : daily averaging ( e.g., [X(00) + X(12)]/2 )

;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

load "$ANAL/SUBS/info_ra.ncl"
; get the reference name of variables (par0_ra_var) and each RA's
;     directory (par0_XXXX_dir), variable names (par0_XXXX_var), and
;     file name info for variables (par0_XXXX_varfile)
; e.g., par0_ra_var = (/"T","U","V","W","Z"/)
;       par0_erai_var = (/"T","U","V","W","Z"/)+"_GDS0_ISBL"

;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

undef("read_dim_xyp")

function read_dim_xyp(ff:file, pname:string, latname:string, lonname:string,  \
                      p_rng:numeric, lat_rng:numeric, lon_rng:numeric)

local vout, ff, pname, latname, lonname, p_rng, lat_rng, lon_rng

begin

 vout = new(3,integer)

 if (dimsizes(p_rng) .eq. 2) then
   vout@p = ff->$pname$({p_rng(0):p_rng(1)})
 else
   vout@p = ff->$pname$({p_rng})
 end if
 if (dimsizes(lat_rng) .eq. 2) then
   vout@lat = ff->$latname$({lat_rng(0):lat_rng(1)})
 else
   vout@lat = ff->$latname$({lat_rng})
 end if
 if (dimsizes(lon_rng) .eq. 2) then
   vout@lon = ff->$lonname$({lon_rng(0):lon_rng(1)})
 else
   vout@lon = ff->$lonname$({lon_rng})
 end if

 vout = (/ dimsizes(vout@p),dimsizes(vout@lat),dimsizes(vout@lon) /)

 return(vout)

end

;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

undef("read_xypt_havg")

function read_xypt_havg(ff:list, vname:string, p_rng:numeric,  \
                        lat_rng:numeric, lon_rng:numeric, opt_havg:integer)

local vout, ff, vname, p_rng, lat_rng, lon_rng, opt_havg,  \
      nlatrng, nlonrng, pname, p, latname, coslat

begin

 nlatrng = dimsizes(lat_rng)
 nlonrng = dimsizes(lon_rng)
 pname = ff[0]->$vname$!1
 if (dimsizes(p_rng) .eq. 2) then
   p = ff[0]->$pname$({p_rng(0):p_rng(1)})
 else
   p = ff[0]->$pname$({p_rng})
 end if
 if (opt_havg(0) .eq. 2) then
   latname = ff[0]->$vname$!2
   coslat = ff[0]->$latname$({lat_rng(0):lat_rng(1)})
   coslat = cos(coslat*3.141592/180.)
 end if

 if ( opt_havg(1) .eq. 1 .and. opt_havg(0) .eq. 1 ) then
   vout = dim_avg_Wrap( dim_avg_Wrap(  \
       ff[:]->$vname$(:,{p},{lat_rng(0):lat_rng(1)},{lon_rng(0):lon_rng(1)}) ) )
 else if ( opt_havg(1) .eq. 1 .and. opt_havg(0) .eq. 2 ) then
   vout = dim_avg_wgt_Wrap( dim_avg_Wrap(  \
       ff[:]->$vname$(:,{p},{lat_rng(0):lat_rng(1)},{lon_rng(0):lon_rng(1)})  \
       ), coslat, 0 )
 else if ( opt_havg(1) .eq. 1 .and. opt_havg(0) .eq. 0 ) then
   if (nlatrng .eq. 2) then
     vout = dim_avg_Wrap( ff[:]->$vname$(:,{p},{lat_rng(0):lat_rng(1)},  \
                                          {lon_rng(0):lon_rng(1)}) )
   else
     vout = dim_avg_Wrap( ff[:]->$vname$(:,{p},{lat_rng},                \
                                          {lon_rng(0):lon_rng(1)}) )
   end if
 else if ( opt_havg(1) .eq. 0 .and. opt_havg(0) .eq. 0 ) then
   if ( nlatrng .eq. 2 .and. nlonrng .eq. 2 ) then
     vout = ff[:]->$vname$(:,{p},{lat_rng(0):lat_rng(1)},  \
                            {lon_rng(0):lon_rng(1)})
   else if ( nlatrng .eq. 2 .and. nlonrng .ne. 2 ) then
     vout = ff[:]->$vname$(:,{p},{lat_rng(0):lat_rng(1)},{lon_rng})
   else if ( nlatrng .ne. 2 .and. nlonrng .eq. 2 ) then
     vout = ff[:]->$vname$(:,{p},{lat_rng},{lon_rng(0):lon_rng(1)})
   else
     vout = ff[:]->$vname$(:,{p},{lat_rng},{lon_rng})
   end if
   end if
   end if
 end if
 end if
 end if
 end if

 return(vout)

end

;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

undef("read_xypt_havg_j")

function read_xypt_havg_j(ff:list, vname:string, p_rng:numeric,  \
                          lat_rng:numeric, lon_rng:numeric, opt_havg:integer,  \
                          it:integer)

local vout, ff, vname, p_rng, lat_rng, lon_rng, opt_havg, it,  \
      nlatrng, nlonrng, pname, p, latname, coslat

begin

 nlatrng = dimsizes(lat_rng)
 nlonrng = dimsizes(lon_rng)
 pname = ff[0]->$vname$!1
 if (dimsizes(p_rng) .eq. 2) then
   p = ff[0]->$pname$({p_rng(0):p_rng(1)})
 else
   p = ff[0]->$pname$({p_rng})
 end if
 if (opt_havg(0) .eq. 2) then
   latname = ff[0]->$vname$!2
   coslat = ff[0]->$latname$({lat_rng(0):lat_rng(1)})
   coslat = cos(coslat*3.141592/180.)
 end if
 
 if ( opt_havg(1) .eq. 1 .and. opt_havg(0) .eq. 1 ) then
   vout = dim_avg_Wrap( dim_avg_Wrap(  \
       ff[:]->$vname$(:,it,{p},{lat_rng(0):lat_rng(1)},  \
                       {lon_rng(0):lon_rng(1)}) ) )
 else if ( opt_havg(1) .eq. 1 .and. opt_havg(0) .eq. 2 ) then
   vout = dim_avg_wgt_Wrap( dim_avg_Wrap(  \
       ff[:]->$vname$(:,it,{p},{lat_rng(0):lat_rng(1)},  \
                       {lon_rng(0):lon_rng(1)}) ), coslat, 0 )
 else if ( opt_havg(1) .eq. 1 .and. opt_havg(0) .eq. 0 ) then
   if (nlatrng .eq. 2) then
     vout = dim_avg_Wrap( ff[:]->$vname$(:,it,{p},{lat_rng(0):lat_rng(1)},  \
                                          {lon_rng(0):lon_rng(1)}) )
   else
     vout = dim_avg_Wrap( ff[:]->$vname$(:,it,{p},{lat_rng},                \
                                          {lon_rng(0):lon_rng(1)}) )
   end if
 else if ( opt_havg(1) .eq. 0 .and. opt_havg(0) .eq. 0 ) then
   if ( nlatrng .eq. 2 .and. nlonrng .eq. 2 ) then
     vout = ff[:]->$vname$(:,it,{p},{lat_rng(0):lat_rng(1)},  \
                            {lon_rng(0):lon_rng(1)})
   else if ( nlatrng .eq. 2 .and. nlonrng .ne. 2 ) then
     vout = ff[:]->$vname$(:,it,{p},{lat_rng(0):lat_rng(1)},{lon_rng})
   else if ( nlatrng .ne. 2 .and. nlonrng .eq. 2 ) then
     vout = ff[:]->$vname$(:,it,{p},{lat_rng},{lon_rng(0):lon_rng(1)})
   else
     vout = ff[:]->$vname$(:,it,{p},{lat_rng},{lon_rng})
   end if
   end if
   end if
 end if
 end if
 end if
 end if

 return(vout)

end

;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

undef("getvar_erai")

function getvar_erai(vname:string, vfname:string, year_i:integer,  \
                     mon_i:integer, hour_i:integer, p_rng:numeric,  \
                     lat_rng:numeric, lon_rng:numeric, opt_havg:integer,  \
                     opt_tavg:integer)

local vout, vname, vfname, year_i, mon_i, hour_i, p_rng, lat_rng, lon_rng,  \
      opt_havg, opt_tavg,  \
      nmon, nhr, mm, hh, wgt, fname, iyr, yyyy, imn, mm2, ff, yyyymmdd00,  \
      nt,dim,nz,ny,nx,var00,var,varhh,dd00,itl,var0100n,var24,ndim0,  \
      tmp1,tmp2, n

begin

 nmon = dimsizes(year_i)*dimsizes(mon_i)
 nhr = dimsizes(hour_i)

 if (nmon .gt. 1024) then
   print("ERA-I: Too many files to read at once: "+nmon+" > 1024")
   exit
 end if

 mm = where(mon_i  .lt. 10, "0"+mon_i , ""+mon_i )
 hh = where(hour_i .lt. 10, "0"+hour_i, ""+hour_i)

 wgt = 1./tofloat(nhr)

 if (opt_tavg .eq. 1) then
   fname = new((/nhr+1,nmon/),string)
 else
   fname = new((/nhr,nmon/),string)
 end if

 n = 0
 do iyr=0, dimsizes(year_i)-1
   yyyy = year_i(iyr)+""
   do imn=0, dimsizes(mon_i)-1
     fname(:nhr-1,n) = par0_erai_dir+"/"+yyyy+"/"+mm(imn)+"/"+  \
                       vfname+hh+".pl."+yyyy+mm(imn)+".nc"
     n = n + 1
   end do
 end do

 if (opt_tavg .eq. 1) then
   mm2 = where(mon_i+1 .lt. 10, "0"+(mon_i+1), ""+(mon_i+1))
   if (mm2 .eq. "13") then
     mm2 = "01"
   end if
   n = 0
   do iyr=0, dimsizes(year_i)-1
     do imn=0, dimsizes(mon_i)-1
       yyyy = year_i(iyr)+""
       if (mm2 .eq. "01") then
         yyyy = (year_i(iyr)+1)+""
       end if
       fname(nhr,n) = par0_erai_dir+"/"+yyyy+"/"+mm2(imn)+"/"+  \
                      vfname+"00.pl."+yyyy+mm2(imn)+".nc"
       n = n + 1
     end do
   end do
 end if

 ff = addfiles(fname(0,:),"r")
 ListSetType(ff, "cat")
 yyyymmdd00 = ff[:]->initial_time0_encoded
 nt = dimsizes( yyyymmdd00 )  ; number of days

 dim = read_dim_xyp(ff[0],"lv_ISBL1","g0_lat_2","g0_lon_3",  \
                    p_rng,lat_rng,lon_rng)
 nz = dim(0)
 ny = dim(1)
 nx = dim(2)
 if (opt_havg(1) .gt. 0) then
   nx = 1
   if (opt_havg(0) .gt. 0) then
     ny = 1
   end if
 end if

 ; 1st HH
 var00 = read_xypt_havg(ff,vname,p_rng,lat_rng,lon_rng,opt_havg)
 if (opt_tavg .eq. 0) then
   var = new((/nt,nhr,nz,ny,nx/),float,1.e20)
   var(:,0,:,:,:) = reshape( var00, (/nt,nz,ny,nx/) )
 else
   if (opt_tavg .eq. 1) then
     var = 0.5*wgt*var00
   else
     var = wgt*var00
   end if
   copy_VarMeta(var00,var)
   var@yyyymmdd = floor(yyyymmdd00 / 100)
 end if
 
 if (nhr .gt. 1) then
   ; other HH
   do n=1, nhr-1
     ff = addfiles(fname(n,:),"r")
     ListSetType(ff, "cat")
     varhh = read_xypt_havg(ff,vname,p_rng,lat_rng,lon_rng,opt_havg)
     if (opt_tavg .eq. 0) then
       var(:,n,:,:,:) = reshape( varhh, (/nt,nz,ny,nx/) )
     else
       var = var + wgt*varhh
     end if
   end do
   delete(varhh)
 end if

 if (opt_tavg .eq. 1) then
   ; 00 UTC on the next day
   dd00 = yyyymmdd00 - floor(yyyymmdd00 / 10000)*10000
   ff = addfiles(fname(nhr,:),"r")
   ListSetType(ff, "join")
   itl = new(nmon,integer)
   itl(:nmon-2) = ind( dd00(1:) .eq. 100 )  ; the last days of each month
   itl(nmon-1) = nt-1
   var0100n = read_xypt_havg_j(ff,vname,p_rng,lat_rng,lon_rng,opt_havg,0)
   var24 = var00
   if (dimsizes(dimsizes(var00)) .eq. 4) then
     var24(:nt-2,:,:,:) = (/ var00(1:,:,:,:) /)
     var24(itl,:,:,:) = (/ var0100n /)
   else if (dimsizes(dimsizes(var00)) .eq. 3) then
     var24(:nt-2,:,:) = (/ var00(1:,:,:) /)
     var24(itl,:,:) = (/ var0100n /)
   else if (dimsizes(dimsizes(var00)) .eq. 2) then
     var24(:nt-2,:) = (/ var00(1:,:) /)
     var24(itl,:) = (/ var0100n /)
   else
     var24(:nt-2) = (/ var00(1:) /)
     var24(itl) = (/ var0100n /)
   end if
   end if
   end if
   var = var + 0.5*wgt*var24
   delete([/var24,var0100n/])
 end if

 delete(ff)

 ; ERA-I only: convert GP to GPH
 if (vname .eq. "Z_GDS0_ISBL") then
   var = var/9.81
 end if

 if (opt_tavg .eq. 0) then
   ndim0 = (/nt*nhr,nz,ny,nx/)
   vout = reshape( var, ndim0( ind(ndim0 .ne. 1) ) )
   copy_VarAtts(var00,vout)
   do n=1, dimsizes(dimsizes(var00))-1
     vout!n = var00!n
     vout&$vout!n$ = var00&$var00!n$
   end do
   tmp1 = conform_dims((/dimsizes(var00&$var00!0$),nhr/),var00&$var00!0$,0)
   tmp2 = conform_dims((/dimsizes(yyyymmdd00),nhr/),yyyymmdd00,0)
   do n=1, nhr-1
     tmp1(:,n) = tmp1(:,0) + (hour_i(n) - hour_i(0))
     tmp2(:,n) = tmp2(:,0) + (hour_i(n) - hour_i(0))
   end do
   vout!0 = var00!0
   vout&$vout!0$ = ndtooned(tmp1)
   vout@yyyymmddhh = ndtooned(tmp2)
   delete([/tmp1,tmp2/])
 else
   vout = var
 end if
 delete(var00)
printVarSummary(vout)

 return(vout)

end

;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

undef("getvar_ra")

function getvar_ra(raname:string, varname:string, year, month, hour,  \
                   p_rng:numeric, lat_rng:numeric, lon_rng:numeric,  \
                   opt_havg:integer, opt_tavg:integer)

local vout, raname, varname, year, month, hour, p_rng, lat_rng, lon_rng,  \
      opt_havg, opt_tavg,  \
      iv_reg, vname, vfname, year_i, mon_i, hour_i

begin

 iv_reg = ind(str_upper(varname) .eq. par0_ra_var)
 vname = par0_erai_var(iv_reg)
 vfname = par0_erai_varfile(iv_reg)

 year_i = toint(year)
 mon_i = toint(month)
 hour_i = toint(hour)

 if ( opt_tavg .eq. 1 .and. hour_i(0) .ne. 0 ) then
   print("the 1st hour should be 00 UTC for opt_tavg = 1.")
   exit
 end if

 if (raname .eq. "ERA-I") then
   vout = getvar_erai(vname,vfname,year_i,mon_i,hour_i,p_rng,lat_rng,lon_rng,  \
                      opt_havg,opt_tavg)
 end if

 return(vout)

end

;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


