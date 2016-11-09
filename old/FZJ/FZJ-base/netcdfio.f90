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

   subroutine dilen(ncid,dimname,nx)

   implicit none
 
   integer, intent(in)           ::  ncid
   character(len=*), intent(in)  ::  dimname
   integer                       ::  dimid,st
   integer, intent(out)          ::  nx
  
   st = nf_inq_dimid(ncid,trim(dimname),dimid)
   call errhdr(st)
   st = nf_inq_dimlen(ncid,dimid,nx)
   call errhdr(st)
   
   write(*,*) trim(dimname)
   write(*,*) 'dimsize : ', nx

   return
   end subroutine dilen

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

   subroutine geta1d_yh(ncid,varname,startx,nx,var)

   implicit none

   integer, intent(in) :: ncid,nx
   integer, intent(in) :: startx
   character(len=*), intent(in) :: varname
   real, dimension(nx), intent(out) :: var

   integer :: varid,st

   st = nf_inq_varid(ncid,trim(varname),varid)
   call errhdr(st)
   st = nf_get_vara_real(ncid,varid,startx,nx,var)
   call errhdr(st)

   write(*,*) trim(varname)
   write(*,*) 'Max : ',maxval(var)
   write(*,*) 'Min : ',minval(var)

   return
   end subroutine geta1d_yh

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

   subroutine geta3d(ncid,varname,nx,ny,nz,start,count,var)
    
   implicit none

   integer, intent(in) :: ncid,nx,ny,nz
   character(len=*), intent(in) :: varname
   integer,dimension(3) :: count,start
   real, dimension(nx,ny,nz), intent(out) :: var
   

   integer :: varid,st

   st = nf_inq_varid(ncid,trim(varname),varid)
   call errhdr(st)
   st = nf_get_vara_real(ncid,varid,start,count,var)
   call errhdr(st)
   
   write(*,*) trim(varname)
   write(*,*) 'Max : ',maxval(var)
   write(*,*) 'Min : ',minval(var)

   return
   end subroutine geta3d

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

   subroutine geta4d(ncid,varname,nx,ny,nz,nt,start,count,var)
    
   implicit none

   integer, intent(in) :: ncid,nx,ny,nz,nt
   character(len=*), intent(in) :: varname
   integer,dimension(4) :: count,start
   real, dimension(nx,ny,nz,nt), intent(out) :: var
   

   integer :: varid,st

   st = nf_inq_varid(ncid,trim(varname),varid)
   call errhdr(st)
   st = nf_get_vara_real(ncid,varid,start,count,var)
   call errhdr(st)
   
   write(*,*) trim(varname)
   write(*,*) 'Max : ',maxval(var)
   write(*,*) 'Min : ',minval(var)

   return
   end subroutine geta4d

   subroutine geta4d_yh(ncid,varname,startx,nx,starty,ny,startz,nz, &
                     startt,nt,var)
   
   implicit none

   integer, intent(in) :: ncid,nx,ny,nz,nt
   integer, intent(in) :: startx, starty, startz, startt
   character(len=*), intent(in) :: varname
   real, dimension(nx,ny,nz,nt), intent(out) :: var


   integer :: varid,st

   st = nf_inq_varid(ncid,trim(varname),varid)
   call errhdr(st)
   st = nf_get_vara_real(ncid,varid,(/startx,starty,startz,startt/),&
                         (/nx,ny,nz,nt/),var)
   call errhdr(st)

   write(*,*) trim(varname)
   write(*,*) 'Max : ',maxval(var)
   write(*,*) 'Min : ',minval(var)

   return
   end subroutine geta4d_yh

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

   subroutine out1d(wname,nv,name,var,d1name,dim1,axis1,title)

   implicit none

   integer, intent(in) :: nv
   character(len=*), intent(in) :: wname,name(nv)
   character(len=*), intent(in) :: d1name
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim1
   real, dimension(dim1), intent(in) :: axis1
   real, dimension(dim1,nv), intent(in) :: var

   integer :: i,j,k,it,st
   integer :: len1, len2(nv), len3, len4, len5, len6, len7
   integer :: ncid
   integer :: d1id, d2id, d3id, d4id
   integer :: a1id, a2id, a3id, a4id
   integer :: varid(nv)
   integer, dimension(1) :: vardimid

   len1 = len_trim(wname)
   do i=1, nv
     len2(i) = len_trim(name(i))
   enddo
   len3 = len_trim(d1name)
   len4 = len_trim(title)

   st = nf_create(wname(1:len1),nf_share,ncid)

   st = nf_def_dim(ncid,d1name(1:len3),dim1,d1id)

   st = nf_def_var(ncid,d1name(1:len3),nf_real,1,d1id,a1id)

   vardimid(1) = d1id

   do i=1, nv
     st = nf_def_var(ncid,trim(name(i)),nf_real,1,vardimid,varid(i))
!     st = nf_def_var(ncid,name(1:len2(i),i),nf_real,1,vardimid,varid(i))
     st = nf_put_att_real(ncid,varid(i),'_FillValue',nf_real,1,1.e+32)
     st = nf_put_att_text(ncid,nf_global,'title',len4,title)
   enddo
   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a1id,axis1)
   do i=1, nv
     st = nf_put_var_real(ncid,varid(i),var(:,i))
   enddo

   st = nf_close(ncid)

   end subroutine out1d

   subroutine out1d_yh(wname,nv4,name4,var4, &
                       d4name1,dim41,axis41, title)

   implicit none

   integer, intent(in) :: nv4
   character(len=*), intent(in) :: wname, name4(nv4)
   character(len=*), intent(in) :: d4name1
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim41
   real, intent(in)    :: axis41(dim41)
   real, dimension(dim41,nv4), intent(in) :: var4

   integer :: i,j,k,it,st
   integer :: ncid
   integer :: d4id1
   integer :: a4id1
   integer :: varid4(nv4)
   integer :: vardimid4(1)


   st = nf_create(trim(wname),nf_share,ncid)

   st = nf_def_dim(ncid,trim(d4name1),dim41,d4id1)

   st = nf_def_var(ncid,trim(d4name1),nf_real,1,d4id1,a4id1)

   vardimid4 = (/d4id1/)

   do i=1, nv4
     st = nf_def_var(ncid,trim(name4(i)),nf_real,1,vardimid4,varid4(i))
     st = nf_put_att_real(ncid,varid4(i),'_FillValue',nf_real,1,1.e+32)
   enddo
   st = nf_put_att_text(ncid,nf_global,'title',len_trim(title),title)
   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a4id1,axis41)

   do i=1, nv4
     st = nf_put_var_real(ncid,varid4(i),var4(:,i))
   enddo

   st = nf_close(ncid)
   end subroutine out1d_yh


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

   st = nf_create(wname(1:len1),nf_share,ncid)

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

   subroutine out2d_yh(wname,nv4,name4,var4, &
                       d4name1,dim41,axis41, &
                       d4name2,dim42,axis42, title)

   implicit none

   integer, intent(in) :: nv4
   character(len=*), intent(in) :: wname, name4(nv4)
   character(len=*), intent(in) :: d4name1, d4name2
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim41, dim42
   real, intent(in)    :: axis41(dim41), axis42(dim42)
   real, dimension(dim41,dim42,nv4), intent(in) :: var4

   integer :: i,j,k,it,st
   integer :: ncid
   integer :: d4id1, d4id2
   integer :: a4id1, a4id2
   integer :: varid4(nv4)
   integer :: vardimid4(2)


   st = nf_create(trim(wname),nf_share,ncid)

   st = nf_def_dim(ncid,trim(d4name1),dim41,d4id1)
   st = nf_def_dim(ncid,trim(d4name2),dim42,d4id2)

   st = nf_def_var(ncid,trim(d4name1),nf_real,1,d4id1,a4id1)
   st = nf_def_var(ncid,trim(d4name2),nf_real,1,d4id2,a4id2)

   vardimid4 = (/d4id1, d4id2/)

   do i=1, nv4
     st = nf_def_var(ncid,trim(name4(i)),nf_real,2,vardimid4,varid4(i))
     st = nf_put_att_real(ncid,varid4(i),'_FillValue',nf_real,1,1.e+32)
   enddo
   st = nf_put_att_text(ncid,nf_global,'title',len_trim(title),title)
   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a4id1,axis41)
   st = nf_put_var_real(ncid,a4id2,axis42)

   do i=1, nv4
     st = nf_put_var_real(ncid,varid4(i),var4(:,:,i))
   enddo

   st = nf_close(ncid)
   end subroutine out2d_yh

   subroutine out2d2(wname,name1,name2,var1,var2,d1name,dim1,axis1,d2name,dim2,axis2,  &
                    title)

   implicit none

   character(len=*), intent(in) :: wname,name1,name2
   character(len=*), intent(in) :: d1name, d2name
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim1, dim2
   real, dimension(dim1), intent(in) :: axis1
   real, dimension(dim2), intent(in) :: axis2
   real, dimension(dim1,dim2), intent(in) :: var1, var2

   integer :: i,j,k,it,st
   integer :: len1, len2, len3, len4, len5, len6, len7
   integer :: ncid
   integer :: d1id, d2id, d3id, d4id
   integer :: a1id, a2id, a3id, a4id
   integer :: var1id, var2id
   integer, dimension(2) :: vardimid

   len1 = len_trim(wname)
   len2 = len_trim(name1)
   len3 = len_trim(name2)
   len4 = len_trim(d1name)
   len5 = len_trim(d2name)
   len6 = len_trim(title)

   st = nf_create(wname(1:len1),nf_share,ncid)

   st = nf_def_dim(ncid,d1name(1:len4),dim1,d1id)
   st = nf_def_dim(ncid,d2name(1:len5),dim2,d2id)

   st = nf_def_var(ncid,d1name(1:len4),nf_real,1,d1id,a1id)
   st = nf_def_var(ncid,d2name(1:len5),nf_real,1,d2id,a2id)

   vardimid(1) = d1id
   vardimid(2) = d2id
   st = nf_def_var(ncid,name1(1:len2),nf_real,2,vardimid,var1id)
   st = nf_def_var(ncid,name2(1:len3),nf_real,2,vardimid,var2id)

   st = nf_put_att_real(ncid,var1id,'_FillValue',nf_real,1,1.e+32)
   st = nf_put_att_real(ncid,var2id,'_FillValue',nf_real,1,1.e+32)
   st = nf_put_att_text(ncid,nf_global,'title',len6,title)
   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a1id,axis1)
   st = nf_put_var_real(ncid,a2id,axis2)
   st = nf_put_var_real(ncid,var1id,var1)
   st = nf_put_var_real(ncid,var2id,var2)

   st = nf_close(ncid)
   end subroutine out2d2

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

   st = nf_create(wname(1:len1),nf_share,ncid)

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

   subroutine out3d_yh(wname,nv4,name4,var4, &
                       d4name1,dim41,axis41, &
                       d4name2,dim42,axis42, &
                       d4name3,dim43,axis43, title)

   implicit none

   integer, intent(in) :: nv4
   character(len=*), intent(in) :: wname, name4(nv4)
   character(len=*), intent(in) :: d4name1, d4name2, d4name3
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim41, dim42, dim43
   real, intent(in)    :: axis41(dim41), axis42(dim42), axis43(dim43)
   real, dimension(dim41,dim42,dim43,nv4), intent(in) :: var4

   integer :: i,j,k,it,st
   integer :: ncid
   integer :: d4id1, d4id2, d4id3
   integer :: a4id1, a4id2, a4id3
   integer :: varid4(nv4)
   integer :: vardimid4(3)


   st = nf_create(trim(wname),nf_share,ncid)

   st = nf_def_dim(ncid,trim(d4name1),dim41,d4id1)
   st = nf_def_dim(ncid,trim(d4name2),dim42,d4id2)
   st = nf_def_dim(ncid,trim(d4name3),dim43,d4id3)

   st = nf_def_var(ncid,trim(d4name1),nf_real,1,d4id1,a4id1)
   st = nf_def_var(ncid,trim(d4name2),nf_real,1,d4id2,a4id2)
   st = nf_def_var(ncid,trim(d4name3),nf_real,1,d4id3,a4id3)

   vardimid4 = (/d4id1, d4id2, d4id3/)

   do i=1, nv4
     st = nf_def_var(ncid,trim(name4(i)),nf_real,3,vardimid4,varid4(i))
     st = nf_put_att_real(ncid,varid4(i),'_FillValue',nf_real,1,1.e+32)
   enddo
   st = nf_put_att_text(ncid,nf_global,'title',len_trim(title),title)
   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a4id1,axis41)
   st = nf_put_var_real(ncid,a4id2,axis42)
   st = nf_put_var_real(ncid,a4id3,axis43)

   do i=1, nv4
     st = nf_put_var_real(ncid,varid4(i),var4(:,:,:,i))
   enddo

   st = nf_close(ncid)
   end subroutine out3d_yh

   subroutine out3d2(wname,name1,name2,var1,var2,&
                    d1name,dim1,axis1,d2name,dim2,axis2,d3name,dim3,axis3,title)

   implicit none

   character(len=*), intent(in) :: wname,name1,name2
   character(len=*), intent(in) :: d1name, d2name, d3name
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim1, dim2, dim3
   real, dimension(dim1), intent(in) :: axis1
   real, dimension(dim2), intent(in) :: axis2
   real, dimension(dim3), intent(in) :: axis3
   real, dimension(dim1,dim2,dim3), intent(in) :: var1,var2

   integer :: i,j,k,it,st
   integer :: len1, len2, len3, len4, len5, len6, len7
   integer :: ncid
   integer :: d1id, d2id, d3id, d4id
   integer :: a1id, a2id, a3id, a4id
   integer :: var1id, var2id
   integer, dimension(3) :: vardimid

   len1 = len_trim(wname)
   len2 = len_trim(name1)
   len3 = len_trim(name2)
   len4 = len_trim(d1name)
   len5 = len_trim(d2name)
   len6 = len_trim(d3name)
   len7 = len_trim(title)

   st = nf_create(wname(1:len1),nf_share,ncid)

   st = nf_def_dim(ncid,d1name(1:len4),dim1,d1id)
   st = nf_def_dim(ncid,d2name(1:len5),dim2,d2id)
   st = nf_def_dim(ncid,d3name(1:len6),dim3,d3id)

   st = nf_def_var(ncid,d1name(1:len4),nf_real,1,d1id,a1id)
   st = nf_def_var(ncid,d2name(1:len5),nf_real,1,d2id,a2id)
   st = nf_def_var(ncid,d3name(1:len6),nf_real,1,d3id,a3id)

   vardimid(1) = d1id
   vardimid(2) = d2id
   vardimid(3) = d3id
   st = nf_def_var(ncid,name1(1:len2),nf_real,3,vardimid,var1id)
   st = nf_def_var(ncid,name2(1:len3),nf_real,3,vardimid,var2id)
  
   st = nf_put_att_real(ncid,var1id,'_FillValue',nf_real,1,1.e+32)
   st = nf_put_att_real(ncid,var2id,'_FillValue',nf_real,1,1.e+32)
   st = nf_put_att_text(ncid,nf_global,'title',len7,title)

   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a1id,axis1)
   st = nf_put_var_real(ncid,a2id,axis2)
   st = nf_put_var_real(ncid,a3id,axis3)
   st = nf_put_var_real(ncid,var1id,var1)
   st = nf_put_var_real(ncid,var2id,var2)

   st = nf_close(ncid)
   end subroutine out3d2

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

   st = nf_create(wname(1:len1),nf_share,ncid)

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

   subroutine out4d_yh(wname,nv4,name4,var4, &
                    d4name1,dim41,axis41, &
                    d4name2,dim42,axis42, &
                    d4name3,dim43,axis43, &
                    d4name4,dim44,axis44, title)

   implicit none

   integer, intent(in) :: nv4
   character(len=*), intent(in) :: wname, name4(nv4)
   character(len=*), intent(in) :: d4name1, d4name2, d4name3, d4name4
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim41, dim42, dim43, dim44
   real, intent(in)    :: axis41(dim41), axis42(dim42), axis43(dim43), axis44(dim44)
   real, dimension(dim41,dim42,dim43,dim44,nv4), intent(in) :: var4

   integer :: i,j,k,it,st
   integer :: ncid
   integer :: d4id1, d4id2, d4id3, d4id4
   integer :: a4id1, a4id2, a4id3, a4id4
   integer :: varid4(nv4)
   integer :: vardimid4(4)


   st = nf_create(trim(wname),nf_share,ncid)

   st = nf_def_dim(ncid,trim(d4name1),dim41,d4id1)
   st = nf_def_dim(ncid,trim(d4name2),dim42,d4id2)
   st = nf_def_dim(ncid,trim(d4name3),dim43,d4id3)
   st = nf_def_dim(ncid,trim(d4name4),dim44,d4id4)

   st = nf_def_var(ncid,trim(d4name1),nf_real,1,d4id1,a4id1)
   st = nf_def_var(ncid,trim(d4name2),nf_real,1,d4id2,a4id2)
   st = nf_def_var(ncid,trim(d4name3),nf_real,1,d4id3,a4id3)
   st = nf_def_var(ncid,trim(d4name4),nf_real,1,d4id4,a4id4)

   vardimid4 = (/d4id1, d4id2, d4id3, d4id4/)

   do i=1, nv4
     st = nf_def_var(ncid,trim(name4(i)),nf_real,4,vardimid4,varid4(i))
     st = nf_put_att_real(ncid,varid4(i),'_FillValue',nf_real,1,1.e+32)
   enddo
   st = nf_put_att_text(ncid,nf_global,'title',len_trim(title),title)
   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a4id1,axis41)
   st = nf_put_var_real(ncid,a4id2,axis42)
   st = nf_put_var_real(ncid,a4id3,axis43)
   st = nf_put_var_real(ncid,a4id4,axis44)

   do i=1, nv4
     st = nf_put_var_real(ncid,varid4(i),var4(:,:,:,:,i))
   enddo

   st = nf_close(ncid)
   end subroutine out4d_yh

   subroutine errhdr(ier,loc)

   implicit none

   include 'netcdf.inc'

   integer, intent(in) :: ier
   integer, intent(in), optional :: loc

   if (ier /= NF_NOERR) then
     print *,loc,':',NF_STRERROR(ier),' in reading daily mean data file'
     stop
   end if

   return
   end subroutine errhdr

   end module netcdfio
