!234
   module netcdfio

!--------------------------------------------------------------
!
!  Authour In-Sun Song
!
!  Purpose
!
!  To be easier in reading and writing output file in NetCDF format.
!
!  July 5, 2001
!  This module is first written
!
!  July 9, 2001
!  Included title and missing value attribute
!  Missing value is 1.e+32
!
!  August 28, 2001
!  Included the documentation for the use of module outvar
!
!  Usage
!
!  opennc    : Open NetCDF file
!
!  opennc (fname,ncid)
!
!     fname  : File name
!     ncid   : NetCDF file id number
!
!     opennc automatically open the specfied file on NF_NOWRITE 
!     (=readonly) mode
!
!
!  get1d, get2d, get3d, get4d
!
!            : Extract variable from a NetCDF file whose id is
!              be already known
!
!  get1d (ncid,varname,nx,var)
!
!     ncid   : NetCDF file id number ( known )
!     varname: Variable name defined in the NetCDF file
!     nx     : The number of dimension    
!     var    : Variable that is already defined
!
!  get2d (ncid,varname,nx,ny,var)
!
!     ncid   : NetCDF file id number ( known )
!     varname: Variable name defined in the NetCDF file
!     nx     : The number of the first dimension    
!     ny     : The number of the second dimension    
!     var    : Variable that is already defined
!
!  get3d (ncid,varname,nx,ny,nz,var)
!
!     ncid   : NetCDF file id number ( known )
!     varname: Variable name defined in the NetCDF file
!     nx     : The number of the first dimension    
!     ny     : The number of the second dimension    
!     nz     : The number of the third dimension    
!     var    : Variable that is already defined
!
!  get4d (ncid,varname,nx,ny,nz,nt,var)
!
!     ncid   : NetCDF file id number ( known )
!     varname: Variable name defined in the NetCDF file
!     nx     : The number of the first dimension ( usually logitude )    
!     ny     : The number of the second dimension ( usually latitude )   
!     nz     : The number of the third dimension ( usually height or pressre )
!     nt     : The number of the fourth dimension ( usually time) 
!     var    : Variable that is already defined
!
!     get?d reads variable from NetCDF file and store the variable
!     to the predefined variable in Fortran program
!
!     The dimension of variable in NetCDF file should be the same
!     as that of predefined variable in Fortran program
!
!  dget1d, dget2d, dget3d, and dget4d are double-precision version of get?d
!  iget1d, iget2d, iget3d, and iget4d are integer version of get?d
!
!  closenc   : Close NetCDF file
!
!  closenc (ncid)
!
!     ncid   : NetCDF file id number ( known )
!     
!  
!  out1d     : Dump 1-dimensional data in NetCDF format
!
!  out1d (wname,name,var,d1name,dim1,axis1,title)
!
!     wname  : File name
!     name   : The name of variables to be dumped
!     var    : Varialble to be dumped          
!     d1name : The name of dimension
!     dim1   : The number of element of dimension array
!     axis1  : Dimension array
!     title  : File title 
!  
!  out2d     : Dump 2-dimensional data in NetCDF format
!
!  out2d (wname,name,var,d1name,dim1,axis1,d2name,dim2,axis2,title)
!
!     wname  : File name
!     name   : The name of variables to be dumped
!     var    : Varialble to be dumped          
!     d1name : The name of the first dimension
!     dim1   : The number of element of the first dimension array
!     axis1  : The first dimension array
!     d2name : The name of the second dimension
!     dim2   : The number of element of the second dimension array
!     axis2  : The second dimension array
!     title  : File title 
!  
!  out3d     : Dump 3-dimensional data in NetCDF format
!
!  out3d (wname,name,var,d1name,dim1,axis1,d2name,dim2,axis2,
!         d3name,dim3,axis3,title)
!
!     wname  : File name
!     name   : The name of variables to be dumped
!     var    : Varialble to be dumped          
!     d1name : The name of the first dimension
!     dim1   : The number of element of the first dimension array
!     axis1  : The first dimension array
!     d2name : The name of the second dimension
!     dim2   : The number of element of the second dimension array
!     axis2  : The second dimension array
!     d3name : The name of the third dimension
!     dim3   : The number of element of the third dimension array
!     axis3  : The third dimension array
!     title  : File title 
!  
!  out4d     : Dump 4-dimensional data in NetCDF format
!
!  out4d (wname,name,var,d1name,dim1,axis1,d2name,dim2,axis2,
!         d3name,dim3,axis3,d4name,dim4,axis4,title)
!
!     wname  : File name
!     name   : The name of variables to be dumped
!     var    : Varialble to be dumped          
!     d1name : The name of the first dimension
!     dim1   : The number of element of the first dimension array
!     axis1  : The first dimension array
!     d2name : The name of the second dimension
!     dim2   : The number of element of the second dimension array
!     axis2  : The second dimension array
!     d3name : The name of the third dimension
!     dim3   : The number of element of the third dimension array
!     axis3  : The third dimension array
!     d4name : The name of the fourth dimension
!     dim4   : The number of element of the fourth dimension array
!     axis4  : The fourth dimension array
!     title  : File title 
!  
!---------------------------------------------------------------

   implicit none

   include 'netcdf.inc'

   contains

!---------------------------------------------------------------

   subroutine opennc(fname,ncid)

   implicit none

   character(len=*), intent(in) :: fname
   integer, intent(out) :: ncid

   integer :: st,l

   write(*,*) trim(fname)
 
   st = nf_open(trim(fname),nf_nowrite,ncid)

   write(*,*) 'NetCDF file id = ',ncid

   return
   end subroutine opennc

!---------------------------------------------------------------

   subroutine closenc(ncid)     

   implicit none
 
   integer, intent(in) :: ncid
   integer :: st

   st = nf_close(ncid)
   call errhdr(st)

   return
   end subroutine closenc

!---------------------------------------------------------------

   subroutine getdimlen(ncid,dimname,nx)

   implicit none

   integer, intent(in) :: ncid
   character(len=*), intent(in) :: dimname
   integer, intent(out) :: nx

   integer :: dimid, st

   st = nf_inq_dimid(ncid,trim(dimname),dimid)
   call errhdr(st)
   st = nf_inq_dimlen(ncid,dimid,nx)
   call errhdr(st)
 
   return
   end subroutine getdimlen

!---------------------------------------------------------------

   subroutine get1d(ncid,varname,nx,var)

   implicit none

   integer, intent(in) :: ncid,nx
   character(len=*), intent(in) :: varname
   real, dimension(nx), intent(out) :: var

   integer :: varid,st

   st = nf_inq_varid(ncid,trim(varname),varid)
   call errhdr(st)
   st = nf_get_var_real(ncid,varid,var)
   call errhdr(st)

   write(*,*) trim(varname)
   write(*,*) 'Max : ',maxval(var)
   write(*,*) 'Min : ',minval(var)

   return
   end subroutine get1d 

!---------------------------------------------------------------

   subroutine get2d(ncid,varname,nx,ny,var)

   implicit none

   integer, intent(in) :: ncid,nx,ny
   character(len=*), intent(in) :: varname
   real, dimension(nx,ny), intent(out) :: var

   integer :: varid,st

   st = nf_inq_varid(ncid,trim(varname),varid)
   call errhdr(st)
   st = nf_get_var_real(ncid,varid,var)
   call errhdr(st)

   write(*,*) trim(varname)
   write(*,*) 'Max : ',maxval(var)
   write(*,*) 'Min : ',minval(var)

   return
   end subroutine get2d

!---------------------------------------------------------------

   subroutine get3d(ncid,varname,nx,ny,nz,var)

   implicit none

   integer, intent(in) :: ncid,nx,ny,nz
   character(len=*), intent(in) :: varname
   real, dimension(nx,ny,nz), intent(out) :: var

   integer :: varid,st

   st = nf_inq_varid(ncid,trim(varname),varid)
   call errhdr(st)
   st = nf_get_var_real(ncid,varid,var)
   call errhdr(st)

   write(*,*) trim(varname)
   write(*,*) 'Max : ',maxval(var)
   write(*,*) 'Min : ',minval(var)

   return
   end subroutine get3d

!---------------------------------------------------------------

   subroutine get4d(ncid,varname,nx,ny,nz,nt,var)

   implicit none

   integer, intent(in) :: ncid,nx,ny,nz,nt
   character(len=*), intent(in) :: varname
   real, dimension(nx,ny,nz,nt), intent(out) :: var

   integer :: varid,st

   st = nf_inq_varid(ncid,trim(varname),varid)
   call errhdr(st)
   st = nf_get_var_real(ncid,varid,var)
   call errhdr(st)

   write(*,*) trim(varname)
   write(*,*) 'Max : ',maxval(var)
   write(*,*) 'Min : ',minval(var)

   return
   end subroutine get4d

!---------------------------------------------------------------

   subroutine dget1d(ncid,varname,nx,var)

   implicit none   

   integer, intent(in) :: ncid,nx 
   character(len=*), intent(in) :: varname
   double precision, dimension(nx), intent(out) :: var

   integer :: varid,st

   st = nf_inq_varid(ncid,trim(varname),varid)
   call errhdr(st)
   st = nf_get_var_double(ncid,varid,var)
   call errhdr(st)

   write(*,*) trim(varname)
   write(*,*) 'Max : ',maxval(var)
   write(*,*) 'Min : ',minval(var)

   return
   end subroutine dget1d

!---------------------------------------------------------------

   subroutine dget2d(ncid,varname,nx,ny,var)

   implicit none

   integer, intent(in) :: ncid,nx,ny
   character(len=*), intent(in) :: varname
   double precision, dimension(nx,ny), intent(out) :: var

   integer :: varid, st

   st = nf_inq_varid(ncid,trim(varname),varid)
   call errhdr(st)
   st = nf_get_var_double(ncid,varid,var)
   call errhdr(st)

   write(*,*) trim(varname)
   write(*,*) 'Max : ',maxval(var)
   write(*,*) 'Min : ',minval(var)

   return
   end subroutine dget2d

!---------------------------------------------------------------

   subroutine dget3d(ncid,varname,nx,ny,nz,var)

   implicit none

   integer, intent(in) :: ncid,nx,ny,nz
   character(len=*), intent(in) :: varname
   double precision, dimension(nx,ny,nz), intent(out) :: var

   integer :: varid, st

   st = nf_inq_varid(ncid,trim(varname),varid)
   call errhdr(st)
   st = nf_get_var_double(ncid,varid,var)
   call errhdr(st)

   write(*,*) trim(varname)
   write(*,*) 'Max : ',maxval(var)
   write(*,*) 'Min : ',minval(var)

   return
   end subroutine dget3d

!---------------------------------------------------------------

   subroutine dget4d(ncid,varname,nx,ny,nz,nt,var)

   implicit none

   integer, intent(in) :: ncid,nx,ny,nz,nt
   character(len=*), intent(in) :: varname
   double precision, dimension(nx,ny,nz,nt), intent(out) :: var

   integer :: varid, st

   st = nf_inq_varid(ncid,trim(varname),varid)
   call errhdr(st)
   st = nf_get_var_double(ncid,varid,var)
   call errhdr(st)

   write(*,*) trim(varname)
   write(*,*) 'Max : ',maxval(var)
   write(*,*) 'Min : ',minval(var)

   return
   end subroutine dget4d

!---------------------------------------------------------------

   subroutine iget1d(ncid,varname,nx,var)

   implicit none 

   integer, intent(in) :: ncid,nx
   character(len=*), intent(in) :: varname
   integer, dimension(nx), intent(out) :: var 

   integer :: varid, st

   st = nf_inq_varid(ncid,trim(varname),varid)
   call errhdr(st)
   st = nf_get_var_int(ncid,varid,var)
   call errhdr(st)

   write(*,*) trim(varname)
   write(*,*) 'Max : ',maxval(var)
   write(*,*) 'Min : ',minval(var)

   return
   end subroutine iget1d

!---------------------------------------------------------------

   subroutine iget2d(ncid,varname,nx,ny,var)

   implicit none

   integer, intent(in) :: ncid,nx,ny
   character(len=*), intent(in) :: varname
   integer, dimension(nx,ny), intent(out) :: var

   integer :: varid, st

   st = nf_inq_varid(ncid,trim(varname),varid)
   call errhdr(st)
   st = nf_get_var_int(ncid,varid,var)
   call errhdr(st)

   write(*,*) trim(varname)
   write(*,*) 'Max : ',maxval(var)
   write(*,*) 'Min : ',minval(var)

   return
   end subroutine iget2d

!---------------------------------------------------------------

   subroutine iget3d(ncid,varname,nx,ny,nz,var)

   implicit none

   integer, intent(in) :: ncid,nx,ny,nz
   character(len=*), intent(in) :: varname
   integer, dimension(nx,ny,nz), intent(out) :: var

   integer :: varid, st

   st = nf_inq_varid(ncid,trim(varname),varid)
   call errhdr(st)
   st = nf_get_var_int(ncid,varid,var)
   call errhdr(st)

   write(*,*) trim(varname)
   write(*,*) 'Max : ',maxval(var)
   write(*,*) 'Min : ',minval(var)

   return
   end subroutine iget3d

!---------------------------------------------------------------

   subroutine iget4d(ncid,varname,nx,ny,nz,nt,var)

   implicit none

   integer, intent(in) :: ncid,nx,ny,nz,nt
   character(len=*), intent(in) :: varname
   integer, dimension(nx,ny,nz,nt), intent(out) :: var

   integer :: varid, st

   st = nf_inq_varid(ncid,trim(varname),varid)
   call errhdr(st)
   st = nf_get_var_int(ncid,varid,var)
   call errhdr(st)

   write(*,*) trim(varname)
   write(*,*) 'Max : ',maxval(var)
   write(*,*) 'Min : ',minval(var)

   return
   end subroutine iget4d

!---------------------------------------------------------------

   subroutine out1d(wname,name,var,d1name,dim1,axis1,title)

   implicit none

   character(len=*), intent(in) :: wname,name
   character(len=*), intent(in) :: d1name 
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim1
   real, dimension(dim1), intent(in) :: axis1
   real, dimension(dim1), intent(in) :: var

   integer :: i,j,k,it,st
   integer :: len1, len2, len3, len4, len5, len6, len7
   integer :: ncid
   integer :: d1id, d2id, d3id, d4id
   integer :: a1id, a2id, a3id, a4id
   integer :: varid
   integer, dimension(1) :: vardimid

   len1 = len_trim(wname)
   len2 = len_trim(name)
   len3 = len_trim(d1name)
   len4 = len_trim(title)

   st = nf_create(wname(1:len1),nf_clobber,ncid)

   st = nf_def_dim(ncid,d1name(1:len3),dim1,d1id)

   st = nf_def_var(ncid,d1name(1:len3),nf_real,1,d1id,a1id)

   vardimid(1) = d1id
   st = nf_def_var(ncid,name(1:len2),nf_real,1,vardimid,varid) 

   st = nf_put_att_real(ncid,varid,'_FillValue',nf_real,1,1.e+32)
   st = nf_put_att_text(ncid,nf_global,'title',len4,title)

   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a1id,axis1)
   st = nf_put_var_real(ncid,varid,var)

   st = nf_close(ncid)

   end subroutine out1d

!---------------------------------------------------------------

   subroutine out1d2(wname,name1,name2,  &
                     var1,var2,d1name,dim1,axis1,title)

   implicit none

   character(len=*), intent(in) :: wname,name1,name2
   character(len=*), intent(in) :: d1name
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim1
   real, dimension(dim1), intent(in) :: axis1
   real, dimension(dim1), intent(in) :: var1,var2

   integer :: i,j,k,it,st
   integer :: len1, len2, len3, len4, len5, len6, len7, len8
   integer :: ncid
   integer :: d1id
   integer :: a1id
   integer :: var1id,var2id,var3id
   integer, dimension(1) :: vardimid

   len1 = len_trim(wname)
   len2 = len_trim(name1)
   len3 = len_trim(name2)
   len6 = len_trim(d1name)
   len7 = len_trim(title)

   st = nf_create(wname(1:len1),nf_clobber,ncid)

   st = nf_def_dim(ncid,d1name(1:len6),dim1,d1id)

   st = nf_def_var(ncid,d1name(1:len6),nf_real,1,d1id,a1id)

   vardimid(1) = d1id

   st = nf_def_var(ncid,name1(1:len2),nf_real,1,vardimid,var1id)
   st = nf_def_var(ncid,name2(1:len3),nf_real,1,vardimid,var2id)

   st = nf_put_att_real(ncid,var1id,'_FillValue',nf_real,1,1.e+32)
   st = nf_put_att_real(ncid,var2id,'_FillValue',nf_real,1,1.e+32)

   st = nf_put_att_text(ncid,nf_global,'title',len7,title)

   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a1id,axis1)

   st = nf_put_var_real(ncid,var1id,var1)
   st = nf_put_var_real(ncid,var2id,var2)

   st = nf_close(ncid)

   end subroutine out1d2

!---------------------------------------------------------------

   subroutine out1d3(wname,name1,name2,name3,  &
                     var1,var2,var3,d1name,dim1,axis1,title)

   implicit none

   character(len=*), intent(in) :: wname,name1,name2,name3
   character(len=*), intent(in) :: d1name
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim1
   real, dimension(dim1), intent(in) :: axis1
   real, dimension(dim1), intent(in) :: var1,var2,var3

   integer :: i,j,k,it,st
   integer :: len1, len2, len3, len4, len5, len6, len7, len8
   integer :: ncid
   integer :: d1id
   integer :: a1id
   integer :: var1id,var2id,var3id
   integer, dimension(1) :: vardimid

   len1 = len_trim(wname)
   len2 = len_trim(name1)
   len3 = len_trim(name2)
   len4 = len_trim(name3)
   len6 = len_trim(d1name)
   len7 = len_trim(title)

   st = nf_create(wname(1:len1),nf_clobber,ncid)

   st = nf_def_dim(ncid,d1name(1:len6),dim1,d1id)

   st = nf_def_var(ncid,d1name(1:len6),nf_real,1,d1id,a1id)

   vardimid(1) = d1id

   st = nf_def_var(ncid,name1(1:len2),nf_real,1,vardimid,var1id)
   st = nf_def_var(ncid,name2(1:len3),nf_real,1,vardimid,var2id)
   st = nf_def_var(ncid,name3(1:len4),nf_real,1,vardimid,var3id)

   st = nf_put_att_real(ncid,var1id,'_FillValue',nf_real,1,1.e+32)
   st = nf_put_att_real(ncid,var2id,'_FillValue',nf_real,1,1.e+32)
   st = nf_put_att_real(ncid,var3id,'_FillValue',nf_real,1,1.e+32)

   st = nf_put_att_text(ncid,nf_global,'title',len7,title)

   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a1id,axis1)

   st = nf_put_var_real(ncid,var1id,var1)
   st = nf_put_var_real(ncid,var2id,var2)
   st = nf_put_var_real(ncid,var3id,var3)

   st = nf_close(ncid)

   end subroutine out1d3

!---------------------------------------------------------------

   subroutine out1d4(wname,name1,name2,name3,name4,  &
                     var1,var2,var3,var4,d1name,dim1,axis1,title)

   implicit none

   character(len=*), intent(in) :: wname,name1,name2,name3,name4
   character(len=*), intent(in) :: d1name
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim1
   real, dimension(dim1), intent(in) :: axis1
   real, dimension(dim1), intent(in) :: var1,var2,var3,var4

   integer :: i,j,k,it,st
   integer :: len1, len2, len3, len4, len5, len6, len7, len8
   integer :: ncid
   integer :: d1id
   integer :: a1id
   integer :: var1id,var2id,var3id,var4id
   integer, dimension(1) :: vardimid

   len1 = len_trim(wname)
   len2 = len_trim(name1)
   len3 = len_trim(name2)
   len4 = len_trim(name3)
   len5 = len_trim(name4)
   len6 = len_trim(d1name)
   len7 = len_trim(title)

   st = nf_create(wname(1:len1),nf_clobber,ncid)

   st = nf_def_dim(ncid,d1name(1:len6),dim1,d1id)

   st = nf_def_var(ncid,d1name(1:len6),nf_real,1,d1id,a1id)

   vardimid(1) = d1id

   st = nf_def_var(ncid,name1(1:len2),nf_real,1,vardimid,var1id)
   st = nf_def_var(ncid,name2(1:len3),nf_real,1,vardimid,var2id)
   st = nf_def_var(ncid,name3(1:len4),nf_real,1,vardimid,var3id)
   st = nf_def_var(ncid,name4(1:len5),nf_real,1,vardimid,var4id)

   st = nf_put_att_real(ncid,var1id,'_FillValue',nf_real,1,1.e+32)
   st = nf_put_att_real(ncid,var2id,'_FillValue',nf_real,1,1.e+32)
   st = nf_put_att_real(ncid,var3id,'_FillValue',nf_real,1,1.e+32)
   st = nf_put_att_real(ncid,var4id,'_FillValue',nf_real,1,1.e+32)

   st = nf_put_att_text(ncid,nf_global,'title',len7,title)

   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a1id,axis1)

   st = nf_put_var_real(ncid,var1id,var1)
   st = nf_put_var_real(ncid,var2id,var2)
   st = nf_put_var_real(ncid,var3id,var3)
   st = nf_put_var_real(ncid,var4id,var4)

   st = nf_close(ncid)

   end subroutine out1d4

!---------------------------------------------------------------

   subroutine out1d8(wname,name1,name2,name3,name4,name5,name6,name7,name8, &
                     var1,var2,var3,var4,var5,var6,var7,var8,d1name,dim1,axis1,title)

   implicit none

   character(len=*), intent(in) :: wname,name1,name2,name3,name4
   character(len=*), intent(in) :: name5,name6,name7,name8
   character(len=*), intent(in) :: d1name
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim1
   real, dimension(dim1), intent(in) :: axis1
   real, dimension(dim1), intent(in) :: var1,var2,var3,var4
   real, dimension(dim1), intent(in) :: var5,var6,var7,var8

   integer :: i,j,k,it,st
   integer :: len1, len2, len3, len4, len5, len6
   integer :: len7, len8, len9, len10,len11,len12
   integer :: ncid
   integer :: d1id
   integer :: a1id
   integer :: var1id,var2id,var3id,var4id
   integer :: var5id,var6id,var7id,var8id
   integer, dimension(1) :: vardimid

   len1 = len_trim(wname)
   len2 = len_trim(name1)
   len3 = len_trim(name2)
   len4 = len_trim(name3)
   len5 = len_trim(name4)
   len6 = len_trim(name5)
   len7 = len_trim(name6)
   len8 = len_trim(name7)
   len9 = len_trim(name8)
   len10= len_trim(d1name)
   len11= len_trim(title)

   st = nf_create(wname(1:len1),nf_clobber,ncid)

   st = nf_def_dim(ncid,d1name(1:len10),dim1,d1id)

   st = nf_def_var(ncid,d1name(1:len10),nf_real,1,d1id,a1id)

   vardimid(1) = d1id

   st = nf_def_var(ncid,name1(1:len2),nf_real,1,vardimid,var1id)
   st = nf_def_var(ncid,name2(1:len3),nf_real,1,vardimid,var2id)
   st = nf_def_var(ncid,name3(1:len4),nf_real,1,vardimid,var3id)
   st = nf_def_var(ncid,name4(1:len5),nf_real,1,vardimid,var4id)
   st = nf_def_var(ncid,name5(1:len6),nf_real,1,vardimid,var5id)
   st = nf_def_var(ncid,name6(1:len7),nf_real,1,vardimid,var6id)
   st = nf_def_var(ncid,name7(1:len8),nf_real,1,vardimid,var7id)
   st = nf_def_var(ncid,name8(1:len9),nf_real,1,vardimid,var8id)

   st = nf_put_att_real(ncid,var1id,'_FillValue',nf_real,1,1.e+32)
   st = nf_put_att_real(ncid,var2id,'_FillValue',nf_real,1,1.e+32)
   st = nf_put_att_real(ncid,var3id,'_FillValue',nf_real,1,1.e+32)
   st = nf_put_att_real(ncid,var4id,'_FillValue',nf_real,1,1.e+32)
   st = nf_put_att_real(ncid,var5id,'_FillValue',nf_real,1,1.e+32)
   st = nf_put_att_real(ncid,var6id,'_FillValue',nf_real,1,1.e+32)
   st = nf_put_att_real(ncid,var7id,'_FillValue',nf_real,1,1.e+32)
   st = nf_put_att_real(ncid,var8id,'_FillValue',nf_real,1,1.e+32)

   st = nf_put_att_text(ncid,nf_global,'title',len11,title)

   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a1id,axis1)

   st = nf_put_var_real(ncid,var1id,var1)
   st = nf_put_var_real(ncid,var2id,var2)
   st = nf_put_var_real(ncid,var3id,var3)
   st = nf_put_var_real(ncid,var4id,var4)
   st = nf_put_var_real(ncid,var5id,var5)
   st = nf_put_var_real(ncid,var6id,var6)
   st = nf_put_var_real(ncid,var7id,var7)
   st = nf_put_var_real(ncid,var8id,var8)

   st = nf_close(ncid)

   end subroutine out1d8

!---------------------------------------------------------------

   subroutine out2d(wname,name,var,d1name,dim1,axis1,d2name,dim2,axis2,  &
                    title)

   implicit none

   character(len=*), intent(in) :: wname,name
   character(len=*), intent(in) :: d1name, d2name 
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim1, dim2
   real, dimension(dim1), intent(in) :: axis1
   real, dimension(dim2), intent(in) :: axis2
   real, dimension(dim1,dim2), intent(in) :: var

   integer :: i,j,k,it,st
   integer :: len1, len2, len3, len4, len5, len6, len7
   integer :: ncid
   integer :: d1id, d2id, d3id, d4id
   integer :: a1id, a2id, a3id, a4id
   integer :: varid
   integer, dimension(2) :: vardimid

   len1 = len_trim(wname)
   len2 = len_trim(name)
   len3 = len_trim(d1name)
   len4 = len_trim(d2name)
   len5 = len_trim(title)

   st = nf_create(wname(1:len1),nf_clobber,ncid)

   st = nf_def_dim(ncid,d1name(1:len3),dim1,d1id)
   st = nf_def_dim(ncid,d2name(1:len4),dim2,d2id)

   st = nf_def_var(ncid,d1name(1:len3),nf_real,1,d1id,a1id)
   st = nf_def_var(ncid,d2name(1:len4),nf_real,1,d2id,a2id)

   vardimid(1) = d1id
   vardimid(2) = d2id
   st = nf_def_var(ncid,name(1:len2),nf_real,2,vardimid,varid) 

   st = nf_put_att_real(ncid,varid,'_FillValue',nf_real,1,1.e+32)
   st = nf_put_att_text(ncid,nf_global,'title',len5,title)
   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a1id,axis1)
   st = nf_put_var_real(ncid,a2id,axis2)
   st = nf_put_var_real(ncid,varid,var)

   st = nf_close(ncid)
   end subroutine out2d

!---------------------------------------------------------------

   subroutine out2d2(wname ,name1,name2,var1,var2,  &
                     d1name,dim1,axis1,d2name,dim2,axis2,title)

   implicit none

   character(len=*), intent(in) :: wname,name1,name2
   character(len=*), intent(in) :: d1name, d2name
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim1, dim2
   real, dimension(dim1), intent(in) :: axis1
   real, dimension(dim2), intent(in) :: axis2
   real, dimension(dim1,dim2), intent(in) :: var1,var2

   integer :: i,j,k,it,st
   integer :: len1, len2, len3, len4, len5, len6, len7, len8, len9
   integer :: ncid
   integer :: d1id, d2id
   integer :: a1id, a2id
   integer :: var1id,var2id
   integer, dimension(2) :: vardimid

   len1 = len_trim(wname)
   len2 = len_trim(name1)
   len3 = len_trim(name2)
   len6 = len_trim(d1name)
   len7 = len_trim(d2name)
   len8 = len_trim(title)

   st = nf_create(wname(1:len1),nf_clobber,ncid)

   st = nf_def_dim(ncid,d1name(1:len6),dim1,d1id)
   st = nf_def_dim(ncid,d2name(1:len7),dim2,d2id)

   st = nf_def_var(ncid,d1name(1:len6),nf_real,1,d1id,a1id)
   st = nf_def_var(ncid,d2name(1:len7),nf_real,1,d2id,a2id)

   vardimid(1) = d1id
   vardimid(2) = d2id

   st = nf_def_var(ncid,name1(1:len2),nf_real,2,vardimid,var1id)
   st = nf_def_var(ncid,name2(1:len3),nf_real,2,vardimid,var2id)

   st = nf_put_att_real(ncid,var1id,'_FillValue',nf_real,1,1.e+32)
   st = nf_put_att_real(ncid,var2id,'_FillValue',nf_real,1,1.e+32)

   st = nf_put_att_text(ncid,nf_global,'title',len8,title)
   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a1id,axis1)
   st = nf_put_var_real(ncid,a2id,axis2)
   st = nf_put_var_real(ncid,var1id,var1)
   st = nf_put_var_real(ncid,var2id,var2)

   st = nf_close(ncid)
   end subroutine out2d2

!---------------------------------------------------------------

   subroutine out2d3(wname ,name1,name2,name3 ,var1,var2,var3,  &
                     d1name,dim1 ,axis1,d2name,dim2,axis2,title)  

   implicit none

   character(len=*), intent(in) :: wname,name1,name2,name3
   character(len=*), intent(in) :: d1name, d2name
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim1, dim2
   real, dimension(dim1), intent(in) :: axis1
   real, dimension(dim2), intent(in) :: axis2
   real, dimension(dim1,dim2), intent(in) :: var1,var2,var3

   integer :: i,j,k,it,st
   integer :: len1, len2, len3, len4, len5, len6, len7, len8, len9
   integer :: ncid
   integer :: d1id, d2id
   integer :: a1id, a2id
   integer :: var1id,var2id,var3id
   integer, dimension(2) :: vardimid

   len1 = len_trim(wname)
   len2 = len_trim(name1)
   len3 = len_trim(name2)
   len4 = len_trim(name3)
   len6 = len_trim(d1name)
   len7 = len_trim(d2name)
   len8 = len_trim(title)

   st = nf_create(wname(1:len1),nf_clobber,ncid)

   st = nf_def_dim(ncid,d1name(1:len6),dim1,d1id)
   st = nf_def_dim(ncid,d2name(1:len7),dim2,d2id)

   st = nf_def_var(ncid,d1name(1:len6),nf_real,1,d1id,a1id)
   st = nf_def_var(ncid,d2name(1:len7),nf_real,1,d2id,a2id)

   vardimid(1) = d1id
   vardimid(2) = d2id

   st = nf_def_var(ncid,name1(1:len2),nf_real,2,vardimid,var1id)
   st = nf_def_var(ncid,name2(1:len3),nf_real,2,vardimid,var2id)
   st = nf_def_var(ncid,name3(1:len4),nf_real,2,vardimid,var3id)

   st = nf_put_att_real(ncid,var1id,'_FillValue',nf_real,1,1.e+32)
   st = nf_put_att_real(ncid,var2id,'_FillValue',nf_real,1,1.e+32)
   st = nf_put_att_real(ncid,var3id,'_FillValue',nf_real,1,1.e+32)

   st = nf_put_att_text(ncid,nf_global,'title',len8,title)
   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a1id,axis1)
   st = nf_put_var_real(ncid,a2id,axis2)
   st = nf_put_var_real(ncid,var1id,var1)
   st = nf_put_var_real(ncid,var2id,var2)
   st = nf_put_var_real(ncid,var3id,var3)

   st = nf_close(ncid)
   end subroutine out2d3

!---------------------------------------------------------------

   subroutine out2d4(wname ,name1,name2,name3,name4,var1,var2,var3,var4,  &
                     d1name,dim1 ,axis1,d2name,dim2,axis2,title)

   implicit none

   character(len=*), intent(in) :: wname,name1,name2,name3,name4
   character(len=*), intent(in) :: d1name, d2name
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim1, dim2
   real, dimension(dim1), intent(in) :: axis1
   real, dimension(dim2), intent(in) :: axis2
   real, dimension(dim1,dim2), intent(in) :: var1,var2,var3,var4

   integer :: i,j,k,it,st
   integer :: len1, len2, len3, len4, len5, len6, len7, len8, len9
   integer :: ncid
   integer :: d1id, d2id
   integer :: a1id, a2id
   integer :: var1id,var2id,var3id,var4id
   integer, dimension(2) :: vardimid

   len1 = len_trim(wname)
   len2 = len_trim(name1)
   len3 = len_trim(name2)
   len4 = len_trim(name3)
   len5 = len_trim(name4)
   len6 = len_trim(d1name)
   len7 = len_trim(d2name)
   len8 = len_trim(title)

   st = nf_create(wname(1:len1),nf_clobber,ncid)

   st = nf_def_dim(ncid,d1name(1:len6),dim1,d1id)
   st = nf_def_dim(ncid,d2name(1:len7),dim2,d2id)

   st = nf_def_var(ncid,d1name(1:len6),nf_real,1,d1id,a1id)
   st = nf_def_var(ncid,d2name(1:len7),nf_real,1,d2id,a2id)

   vardimid(1) = d1id
   vardimid(2) = d2id

   st = nf_def_var(ncid,name1(1:len2),nf_real,2,vardimid,var1id)
   st = nf_def_var(ncid,name2(1:len3),nf_real,2,vardimid,var2id)
   st = nf_def_var(ncid,name3(1:len4),nf_real,2,vardimid,var3id)
   st = nf_def_var(ncid,name4(1:len5),nf_real,2,vardimid,var4id)

   st = nf_put_att_real(ncid,var1id,'_FillValue',nf_real,1,1.e+32)
   st = nf_put_att_real(ncid,var2id,'_FillValue',nf_real,1,1.e+32)
   st = nf_put_att_real(ncid,var3id,'_FillValue',nf_real,1,1.e+32)
   st = nf_put_att_real(ncid,var4id,'_FillValue',nf_real,1,1.e+32)

   st = nf_put_att_text(ncid,nf_global,'title',len8,title)
   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a1id,axis1)
   st = nf_put_var_real(ncid,a2id,axis2)
   st = nf_put_var_real(ncid,var1id,var1)
   st = nf_put_var_real(ncid,var2id,var2)
   st = nf_put_var_real(ncid,var3id,var3)
   st = nf_put_var_real(ncid,var4id,var4)

   st = nf_close(ncid)
   end subroutine out2d4

!---------------------------------------------------------------

   subroutine out3d(wname,name,var,d1name,dim1,axis1,d2name,dim2,axis2, &
                                   d3name,dim3,axis3,title)

   implicit none

   character(len=*), intent(in) :: wname,name
   character(len=*), intent(in) :: d1name, d2name, d3name 
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim1, dim2, dim3
   real, dimension(dim1), intent(in) :: axis1
   real, dimension(dim2), intent(in) :: axis2
   real, dimension(dim3), intent(in) :: axis3
   real, dimension(dim1,dim2,dim3), intent(in) :: var

   integer :: i,j,k,it,st
   integer :: len1, len2, len3, len4, len5, len6, len7
   integer :: ncid
   integer :: d1id, d2id, d3id, d4id
   integer :: a1id, a2id, a3id, a4id
   integer :: varid
   integer, dimension(3) :: vardimid

   len1 = len_trim(wname)
   len2 = len_trim(name)
   len3 = len_trim(d1name)
   len4 = len_trim(d2name)
   len5 = len_trim(d3name)
   len6 = len_trim(title)

   st = nf_create(wname(1:len1),nf_clobber,ncid)

   st = nf_def_dim(ncid,d1name(1:len3),dim1,d1id)
   st = nf_def_dim(ncid,d2name(1:len4),dim2,d2id)
   st = nf_def_dim(ncid,d3name(1:len5),dim3,d3id)

   st = nf_def_var(ncid,d1name(1:len3),nf_real,1,d1id,a1id)
   st = nf_def_var(ncid,d2name(1:len4),nf_real,1,d2id,a2id)
   st = nf_def_var(ncid,d3name(1:len5),nf_real,1,d3id,a3id)

   vardimid(1) = d1id
   vardimid(2) = d2id
   vardimid(3) = d3id
   st = nf_def_var(ncid,name(1:len2),nf_real,3,vardimid,varid) 
   
   st = nf_put_att_real(ncid,varid,'_FillValue',nf_real,1,1.e+32)
   st = nf_put_att_text(ncid,nf_global,'title',len6,title)

   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a1id,axis1)
   st = nf_put_var_real(ncid,a2id,axis2)
   st = nf_put_var_real(ncid,a3id,axis3)
   st = nf_put_var_real(ncid,varid,var)

   st = nf_close(ncid)
   end subroutine out3d

   subroutine out4d(wname,name,var,d1name,dim1,axis1,d2name,dim2,axis2, &
                                   d3name,dim3,axis3,d4name,dim4,axis4,  &
                                   title)

   implicit none

   character(len=*), intent(in) :: wname,name
   character(len=*), intent(in) :: d1name, d2name, d3name, d4name 
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim1, dim2, dim3, dim4
   real, dimension(dim1), intent(in) :: axis1
   real, dimension(dim2), intent(in) :: axis2
   real, dimension(dim3), intent(in) :: axis3
   real, dimension(dim4), intent(in) :: axis4
   real, dimension(dim1,dim2,dim3,dim4), intent(in) :: var

   integer :: i,j,k,it,st
   integer :: len1, len2, len3, len4, len5, len6, len7
   integer :: ncid
   integer :: d1id, d2id, d3id, d4id
   integer :: a1id, a2id, a3id, a4id
   integer :: varid
   integer, dimension(4) :: vardimid

   len1 = len_trim(wname)
   len2 = len_trim(name)
   len3 = len_trim(d1name)
   len4 = len_trim(d2name)
   len5 = len_trim(d3name)
   len6 = len_trim(d4name)
   len7 = len_trim(title)

   st = nf_create(wname(1:len1),nf_clobber,ncid)

   st = nf_def_dim(ncid,d1name(1:len3),dim1,d1id)
   st = nf_def_dim(ncid,d2name(1:len4),dim2,d2id)
   st = nf_def_dim(ncid,d3name(1:len5),dim3,d3id)
   st = nf_def_dim(ncid,d4name(1:len6),dim4,d4id)

   st = nf_def_var(ncid,d1name(1:len3),nf_real,1,d1id,a1id)
   st = nf_def_var(ncid,d2name(1:len4),nf_real,1,d2id,a2id)
   st = nf_def_var(ncid,d3name(1:len5),nf_real,1,d3id,a3id)
   st = nf_def_var(ncid,d4name(1:len6),nf_real,1,d4id,a4id)

   vardimid(1) = d1id
   vardimid(2) = d2id
   vardimid(3) = d3id
   vardimid(4) = d4id
   st = nf_def_var(ncid,name(1:len2),nf_real,4,vardimid,varid) 

   st = nf_put_att_real(ncid,varid,'_FillValue',nf_real,1,1.e+32)
   st = nf_put_att_text(ncid,nf_global,'title',len7,title)

   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a1id,axis1)
   st = nf_put_var_real(ncid,a2id,axis2)
   st = nf_put_var_real(ncid,a3id,axis3)
   st = nf_put_var_real(ncid,a4id,axis4)
   st = nf_put_var_real(ncid,varid,var)

   st = nf_close(ncid)
   end subroutine out4d

   subroutine errhdr(ier,loc)

   implicit none

   include 'netcdf.inc'

   integer, intent(in) :: ier
   integer, intent(in), optional :: loc

   if (ier /= NF_NOERR) then
     print *,loc,':',NF_STRERROR(ier)
     stop
   end if

   return
   end subroutine errhdr

   end module netcdfio
