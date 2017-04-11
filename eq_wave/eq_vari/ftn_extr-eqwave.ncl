undef("read_eqwave_rc0")
function read_eqwave_rc0(wav_name:string, var_name:string, f_name:string,  \
                     lat:numeric, p:numeric)
local var, wav_name, var_name, f_name, lat, p,  \
      iw_k, iw_me, iw_mw, iw_r, iw_ie, iw_iw,  \
      iw_ie_ss, iw_iw_ss, iw_t_dm, iw_t_sm,  \
      f1, f2, s_axis, lats, ps, nw, nv, ny, nz, nt, wav_name_1, wav_name_9,  \
      itmp, iw, nstr, var5d, dim5, iv, tmp, tmp2, iw_m, iw_i, iw_i_ss
begin

;=========  in file 1 : reconstr_epfsa  ==========
; Kelvin -----------------------------------------
 iw_k = 2
; Mixed Rossby-gravity ---------------------------
 iw_me = 4
 iw_mw = 5
; Rossby -----------------------------------------
 iw_r = 3 + (/0,10/)
; iw_r = 19
; Inertio-gravity --------------------------------
 iw_ie = 6 + (/0,10/)
 iw_iw = 7 + (/0,10/)
; iw_ie = 20
; iw_iw = 21
;=========  in file 2 : reconstr_epf  ============
; smaller-scale waves ( Lx < 2000 km or T < 4/3 days )
 iw_ie_ss = 8
 iw_iw_ss = 9
; diurnal and semi-diurnal migrating tides
 iw_t_dm = 2
 iw_t_sm = 3
;=================================================

 f1 = addfile(f_name(0),"r")
 if (dimsizes(f_name) .gt. 1) then
   f2 = addfile(f_name(1),"r")
 end if

 s_axis = f1->$var_name(0)$!3
 lats = f1->$s_axis$({-lat:lat})
 s_axis = f1->$var_name(0)$!2
 ps = f1->$s_axis$({p})

 nw = dimsizes(wav_name)
 nv = dimsizes(var_name)
 ny = dimsizes(lats)
 nz = dimsizes(ps)
 nt = dimsizes( f1->$var_name(0)$(:,0,0,0) )

 wav_name_1 = str_upper( str_get_cols( wav_name, 0, 0 ) )
 wav_name_9 = new(nw,string)
 itmp = dimsizes( tochar("abc") ) - 3   ; = 1 for old NCL versions
 do iw=0, nw-1
   nstr = dimsizes( tochar(wav_name(iw)) ) - itmp
   wav_name_9(iw) = str_upper( str_get_cols( wav_name(iw), nstr-1, nstr-1 ) )
 end do

 var5d = new((/nv,nw,nt,nz,ny/),float)
 dim5 = dimsizes(var5d)

 do iw=0, nw-1

   if ( wav_name_1(iw) .eq. "K" ) then
     do iv=0, nv-1
       tmp = f1->$var_name(iv)$(:,iw_k,{p},{-lat:lat})
       var5d(iv,iw,:,:,:) = (/ tmp /)
     end do
   end if
   if ( wav_name_1(iw) .eq. "M" ) then
     if ( wav_name_9(iw) .eq. "E" ) then
       iw_m = iw_me
     else if ( wav_name_9(iw) .eq. "W" ) then
       iw_m = iw_mw
     end if
     end if
     if ( .not. isdefined("iw_m") ) then
       print("Error in the argument (0) in the function read_eqwave_rc0.")
       exit
     end if
     do iv=0, nv-1
       tmp = f1->$var_name(iv)$(:,iw_m,{p},{-lat:lat})
       var5d(iv,iw,:,:,:) = (/ tmp /)
     end do
     delete(iw_m)
   end if
   if ( wav_name_1(iw) .eq. "R" ) then
     do iv=0, nv-1
       tmp  = f1->$var_name(iv)$(:,iw_r(0),{p},{-lat:lat})
       tmp2 = f1->$var_name(iv)$(:,iw_r(1),{p},{-lat:lat})
       tmp = tmp + tmp2
       var5d(iv,iw,:,:,:) = (/ tmp /)
     end do
   end if
   if ( wav_name_1(iw) .eq. "I" ) then
     if ( wav_name_9(iw) .eq. "E" ) then
       iw_i    = iw_ie
       iw_i_ss = iw_ie_ss
     else if ( wav_name_9(iw) .eq. "W" ) then
       iw_i    = iw_iw
       iw_i_ss = iw_iw_ss
     end if
     end if
     if ( .not. isdefined("iw_i") ) then
       print("Error in the argument (0) in the function read_eqwave_rc0.")
       exit
     end if
     do iv=0, nv-1
       tmp  = f1->$var_name(iv)$(:,iw_i(0),{p},{-lat:lat})
       tmp2 = f1->$var_name(iv)$(:,iw_i(1),{p},{-lat:lat})
       tmp = tmp + tmp2
       ; smaller-scale waves
       tmp2 = f2->$var_name(iv)$(:,iw_i_ss,{p},{-lat:lat})
       tmp = tmp + tmp2
       if ( wav_name_9(iw) .eq. "W" ) then
         ; migrating tides
         tmp2 = f2->$var_name(iv)$(:,iw_t_dm,{p},{-lat:lat}) +  \
                f2->$var_name(iv)$(:,iw_t_sm,{p},{-lat:lat})
         tmp = tmp + tmp2
       end if
       var5d(iv,iw,:,:,:) = (/ tmp /)
     end do
     delete(iw_i)
   end if
   if ( wav_name_1(iw) .eq. "T" ) then
     do iv=0, nv-1
       tmp = f2->$var_name(iv)$(:,iw_t_dm,{p},{-lat:lat})
       if ( wav_name_9(iw) .ne. "1" ) then
         tmp = tmp + f2->$var_name(iv)$(:,iw_t_sm,{p},{-lat:lat})
       end if
       var5d(iv,iw,:,:,:) = (/ tmp /)
     end do
   end if

   if ( .not. isdefined("tmp") ) then
     print("Error in the argument (0) in the function read_eqwave_rc0.")
     exit
   end if
 
   delete(tmp)

 end do

 if ( any( dim5 .eq. 1 ) ) then
   var = onedtond( ndtooned(var5d), dim5( ind( dim5 .ne. 1 ) ) )
 else
   var = var5d
 end if
 var!(dimsizes(dimsizes(var))-1) = "lat"
 var&lat = lats
 if (nz .ne. 1) then
   var!(dimsizes(dimsizes(var))-2) = "p"
   var&p = ps
 else
   var@p = ps
 end if
 if (nv .ne. 1) then
   var!0 = "var"
   var&var = var_name
 end if
 if (nw .ne. 1) then
   if (nv .eq. 1) then
     var!0 = "wave"
   else
     var!1 = "wave"
   end if
   var&wave = wav_name
 end if

return(var)
end

