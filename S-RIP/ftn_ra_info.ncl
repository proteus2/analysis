
load "./param_ra.ncl"

;-------------------------------------------------------------------------------

undef("get_ra_ind")
function get_ra_ind(ra:string)
local ira, ra, nr, ra_123, ra_9, ra_89, ra_789, itmp, nstr, ir,  \
      nr0, ra0_123, ra0_9, ra0_89, ra0_789, ii, ii2
begin

 nr = dimsizes(ra)
 ra_123 = str_upper( str_get_cols( ra, 0, 2 ) )
 ra_9   = new(nr,string)
 ra_89  = new(nr,string)
 ra_789 = new(nr,string)
 itmp = dimsizes( tochar("abc") ) - 3   ; = 1 for old NCL versions
 do ir=0, nr-1
   nstr = dimsizes( tochar(ra(ir)) ) - itmp
   ra_9(ir) = str_upper( str_get_cols( ra(ir), nstr-1, nstr-1 ) )
   ra_89(ir) = str_upper( str_get_cols( ra(ir), nstr-2, nstr-1 ) )
   ra_789(ir) = str_upper( str_get_cols( ra(ir), nstr-3, nstr-1 ) )
 end do

 nr0 = dimsizes(param_ra_ranames)
 ra0_123 = str_upper( str_get_cols( param_ra_ranames, 0, 2 ) )
 ra0_9   = new(nr0,string)
 ra0_89  = new(nr0,string)
 ra0_789 = new(nr0,string)
 do ir=0, nr0-1
   nstr = dimsizes( tochar(param_ra_ranames(ir)) ) - itmp
   ra0_9(ir) = str_upper( str_get_cols( param_ra_ranames(ir),  \
                   nstr-1, nstr-1 ) )
   ra0_89(ir) = str_upper( str_get_cols( param_ra_ranames(ir),  \
                    nstr-2, nstr-1 ) )
   ra0_789(ir) = str_upper( str_get_cols( param_ra_ranames(ir),  \
                     nstr-3 , nstr-1 ) )
 end do

 ira = new(nr,integer,999)

 do ir=0, nr-1
   ii = ind( ra0_123 .eq. ra_123(ir) .and. ra0_9 .eq. ra_9(ir) )
   if ( dimsizes(ii) .eq. 1 ) then
     ira(ir) = ii
   else
     ii2 = ind( ra0_123 .eq. ra_123(ir) .and. ra0_89 .eq. ra_89(ir) )
     if ( dimsizes(ii2) .eq. 1 ) then
       ira(ir) = ii2
     else
       ira(ir) = ind( ra0_123 .eq. ra_123(ir) .and. ra0_789 .eq. ra_789(ir) )
     end if
     delete(ii2)
   end if
   delete(ii)
 end do

return(ira)
end

;-------------------------------------------------------------------------------

undef("get_ra_names_att")
procedure get_ra_names_att(var)
local var
begin

 var@n_info = 4
 var@dim1_0 = "index based on the parameter file"
 var@dim1_1 = "dataset name"
 var@dim1_2 = "file header"
 var@dim1_3 = "directory name"
 var@file_h   = 2
 var@dirname  = 3

return(var)
end

;-------------------------------------------------------------------------------

undef("get_ra_names_2d")
function get_ra_names_2d(ra:string)
local sra2d, ra, ira, n_info
begin

 ira = get_ra_ind(ra)

 n_info = 4
 sra2d = new((/dimsizes(ra),n_info/),string)

 sra2d(:,0) = ira
 sra2d(:,1) = param_ra_ranames (ira)
 sra2d(:,2) = param_ra_expnames(ira)
 sra2d(:,3) = param_ra_dirnames(ira)

 get_ra_names_att(sra2d)  ; should be matching

 if (sra2d@n_info .ne. n_info) then
   print("Error in the function get_ra_names_3d.")
   exit
 end if

return(sra2d)
end

;-------------------------------------------------------------------------------

undef("get_ra_names")
function get_ra_names(ra:string)
local sra, ra, ira, n_info, sra2d
begin

 sra2d = get_ra_names_2d(ra)

 if (dimsizes(ra) .eq. 1) then
   sra = sra2d(0,:)
   get_ra_names_att(sra)
 else
   sra = sra2d
 end if

return(sra)
end

;-------------------------------------------------------------------------------

undef("get_ra_varind")
function get_ra_varind(var:string)
local ivar, var, nv, var_ren, iv
begin

 nv = dimsizes(var)

 var_ren = where(str_upper(str_get_cols(var,0,2)) .eq. "OME", "w", var)

 ivar = new(nv,integer)
 do iv=0, nv-1
   ivar(iv) = ind( str_upper(param_ra_vars) .eq. str_upper(var_ren(iv)) )
 end do

return(ivar)
end

;-------------------------------------------------------------------------------

undef("get_ra_varnames_att")
procedure get_ra_varnames_att(var)
local var
begin

 var@n_info = 4
 var@dim1_0 = "index based on the parameter file"
 var@dim1_1 = "variable name"
 var@dim1_2 = "variable name as saved in files"
 var@dim1_3 = "variable name part of file names"
 var@varname  = 2
 var@fvarname = 3

return(var)
end

;-------------------------------------------------------------------------------

undef("get_ra_varnames_3d")
function get_ra_varnames_3d(ra:string, var:string)
local vra, ra, var, nr, nv, ira, ivar, n_info, vra3d, ir
begin

 nr = dimsizes(ra)
 nv = dimsizes(var)

 ira = get_ra_ind(ra)
 ivar = get_ra_varind(var)

 n_info = 4
 vra3d = new((/nr,n_info,nv/),string)

 do ir=0, nr-1
   vra3d(ir,0,:) = ivar
   vra3d(ir,1,:) = param_ra_vars   (ivar)
   vra3d(ir,2,:) = param_ra_vnames (ira(ir),ivar)
   vra3d(ir,3,:) = param_ra_fvnames(ira(ir),ivar)
 end do
 
 get_ra_names_att(vra3d)  ; should be matching

 if (vra3d@n_info .ne. n_info) then
   print("Error in the function get_ra_varnames_3d.")
   exit
 end if

return(vra3d)
end

;-------------------------------------------------------------------------------

undef("get_ra_varnames")
function get_ra_varnames(ra:string, var:string)
local vra, ra, var, vra3d, dim3d, dim
begin

 vra3d = get_ra_varnames_3d(ra,var)

 dim3d = dimsizes(vra3d)
 dim = dim3d( ind( dim3d .ne. 1 ) )
 vra = onedtond( ndtooned(vra3d), dim )
 get_ra_varnames_att(vra)

return(vra)
end

