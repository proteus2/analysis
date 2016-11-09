!234
   module netio

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
!---------------------------------------------------------------

   implicit none

   include 'netcdf.inc'

   type ::  vset
     real, dimension(:,:,:,:), allocatable ::  var_out
     real, dimension(:),       allocatable ::  axis1, axis2, axis3, axis4
     integer                               ::  nd(4)
     character(len=32)                     ::  vname, axis(4)
   end type vset

   type ::  vset0
     real, dimension(:),       allocatable ::  axis1, axis2, axis3, axis4
     integer                               ::  nd(4)
     character(len=32)                     ::  vname, axis(4)
   end type vset0


   contains

!---------------------------------------------------------------

   subroutine opennc(fname,ncid)

   implicit none

   character(len=*), intent(in) :: fname
   integer, intent(out) :: ncid

   integer :: st,l

   st = nf_open(trim(fname),nf_nowrite,ncid)

   write(*,'(i3,2x,a)') ncid,trim(fname)

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

   subroutine gettitle(ncid,title)

   implicit none

   integer, intent(in) :: ncid
   character(len=*), intent(out) :: title

   integer :: st

   st = nf_get_att_text(ncid,nf_global,'title',title)
   call errhdr(st)

   return
   end subroutine gettitle

!---------------------------------------------------------------

   subroutine dilen(ncid,dimname,nx,l_err)

   implicit none
 
   integer, intent(in)           ::  ncid
   character(len=*), intent(in)  ::  dimname
   integer                       ::  dimid,st
   integer, intent(out)          ::  nx
   logical, intent(out)          ::  l_err
 
   l_err = .FALSE.
   st = nf_inq_dimid(ncid,trim(dimname),dimid)
   call errhdr_nonstop(st)
   st = nf_inq_dimlen(ncid,dimid,nx)
   call errhdr_nonstop(st)

   if (st /= NF_NOERR) then
     l_err = .TRUE.
   else
     write(*,*) trim(dimname), ' size : ', nx
   end if

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

!   write(*,*) trim(varname)
!   write(*,*) '   Max : ',maxval(var)
!   write(*,*) '   Min : ',minval(var)

   return
   end subroutine get1d 

!---------------------------------------------------------------

   subroutine geta1d(ncid,varname,startx,nx,var)

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
   write(*,*) '   Max : ',maxval(var)
   write(*,*) '   Min : ',minval(var)

   return
   end subroutine geta1d

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
   write(*,*) '   Max : ',maxval(var)
   write(*,*) '   Min : ',minval(var)

   return
   end subroutine get2d

!---------------------------------------------------------------

   subroutine geta2d(ncid,varname,startx,nx,starty,ny,var)

   implicit none

   integer, intent(in) :: ncid,nx,ny
   integer, intent(in) :: startx, starty
   character(len=*), intent(in) :: varname
   real, dimension(nx,ny), intent(out) :: var

   integer :: varid,st

   st = nf_inq_varid(ncid,trim(varname),varid)
   call errhdr(st)
   st = nf_get_vara_real(ncid,varid,(/startx,starty/),(/nx,ny/),var)
   call errhdr(st)

   write(*,*) trim(varname)
   write(*,*) '   Max : ',maxval(var)
   write(*,*) '   Min : ',minval(var)

   return
   end subroutine geta2d

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
   write(*,*) '   Max : ',maxval(var)
   write(*,*) '   Min : ',minval(var)

   return
   end subroutine get3d

!---------------------------------------------------------------

   subroutine geta3d(ncid,varname,startx,nx,starty,ny,startz,nz,var)
    
   implicit none

   integer, intent(in) :: ncid,nx,ny,nz
   integer, intent(in) :: startx, starty, startz
   character(len=*), intent(in) :: varname
   real, dimension(nx,ny,nz), intent(out) :: var
   

   integer :: varid,st

   st = nf_inq_varid(ncid,trim(varname),varid)
   call errhdr(st)
   st = nf_get_vara_real(ncid,varid,(/startx,starty,startz/),(/nx,ny,nz/),var)
   call errhdr(st)
   
   write(*,*) trim(varname)
   write(*,*) '   Max : ',maxval(var)
   write(*,*) '   Min : ',minval(var)

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
   write(*,*) '   Max : ',maxval(var)
   write(*,*) '   Min : ',minval(var)

   return
   end subroutine get4d

!---------------------------------------------------------------

   subroutine geta4d(ncid,varname,startx,nx,starty,ny,startz,nz, &
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

!   write(*,*) trim(varname)
!   write(*,*) '   Max : ',maxval(var)
!   write(*,*) '   Min : ',minval(var)

   return
   end subroutine geta4d

!----------------------------------------------------------------

   subroutine get4ds(ncid,varname,nx,ny,nz,nt,var)

   implicit none

   integer, intent(in) :: ncid,nx,ny,nz,nt
   character(len=*), intent(in) :: varname
   integer*2, dimension(nx,ny,nz,nt)  :: var_temp
   real, dimension(nx,ny,nz,nt), intent(out) :: var

   integer :: varid,st 
   
   st = nf_inq_varid(ncid,trim(varname),varid)
   call errhdr(st)
   st = nf_get_var_int2(ncid,varid,var_temp)
   call errhdr(st)
   var = float(var_temp)

   write(*,*) trim(varname)
   write(*,*) '   Max : ',maxval(var)
   write(*,*) '   Min : ',minval(var)

   return
   end subroutine get4ds

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
   write(*,*) '   Max : ',maxval(var)
   write(*,*) '   Min : ',minval(var)

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
   write(*,*) '   Max : ',maxval(var)
   write(*,*) '   Min : ',minval(var)

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
   write(*,*) '   Max : ',maxval(var)
   write(*,*) '   Min : ',minval(var)

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
   write(*,*) '   Max : ',maxval(var)
   write(*,*) '   Min : ',minval(var)

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
   write(*,*) '   Max : ',maxval(var)
   write(*,*) '   Min : ',minval(var)

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
   write(*,*) '   Max : ',maxval(var)
   write(*,*) '   Min : ',minval(var)

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
   write(*,*) '   Max : ',maxval(var)
   write(*,*) '   Min : ',minval(var)

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
   write(*,*) '   Max : ',maxval(var)
   write(*,*) '   Min : ',minval(var)

   return
   end subroutine iget4d

!---------------------------------------------------------------

   subroutine igeta4d(ncid,varname,startx,nx,starty,ny,startz,nz, &
                     startt,nt,var)

   implicit none

   integer, intent(in) :: ncid,nx,ny,nz,nt
   integer, intent(in) :: startx, starty, startz, startt
   character(len=*), intent(in) :: varname
   integer, dimension(nx,ny,nz,nt), intent(out) :: var


   integer :: varid,st

   st = nf_inq_varid(ncid,trim(varname),varid)
   call errhdr(st)
   st = nf_get_vara_int(ncid,varid,(/startx,starty,startz,startt/),&
                        (/nx,ny,nz,nt/),var)
   call errhdr(st)

!   write(*,*) trim(varname)
!   write(*,*) '   Max : ',maxval(var)
!   write(*,*) '   Min : ',minval(var)

   return
   end subroutine igeta4d

!---------------------------------------------------------------

   subroutine sgeta4d(ncid,varname,startx,nx,starty,ny,startz,nz, &
                     startt,nt,var)

   implicit none
   
   integer, intent(in) :: ncid,nx,ny,nz,nt
   integer, intent(in) :: startx, starty, startz, startt
   character(len=*), intent(in) :: varname
   integer*2, dimension(nx,ny,nz,nt), intent(out) :: var

   integer :: varid,st
   
   st = nf_inq_varid(ncid,trim(varname),varid)
   call errhdr(st)               
   st = nf_get_vara_int2(ncid,varid,(/startx,starty,startz,startt/),&
                         (/nx,ny,nz,nt/),var)
   call errhdr(st)
   
!   write(*,*) trim(varname)
!   write(*,*) '   Max : ',maxval(var)
!   write(*,*) '   Min : ',minval(var)
   
   return
   end subroutine sgeta4d

!----------------------------------------------------------------

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

   subroutine out1d2(wname,nv1,name1,var1,d1name,dim1,axis1, &
                     nv2,name2,var2,d2name,dim2,axis2,title)

   implicit none

   integer, intent(in) :: nv1, nv2
   character(len=*), intent(in) :: wname,name1(nv1),name2(nv2)
   character(len=*), intent(in) :: d1name, d2name
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim1, dim2
   real, dimension(dim1), intent(in) :: axis1
   real, dimension(dim2), intent(in) :: axis2
   real, dimension(dim1,nv1), intent(in) :: var1
   real, dimension(dim2,nv2), intent(in) :: var2

   integer :: i,j,k,it,st
   integer :: len1,len21(nv1),len22(nv2),len31,len32,len4,len5,len6,len7
   integer :: ncid
   integer :: d1id, d2id, d3id, d4id
   integer :: a1id, a2id, a3id, a4id
   integer :: varid1(nv1), varid2(nv2)
   integer, dimension(1) :: vardimid1, vardimid2

   len1 = len_trim(wname)
   len21(:) = len_trim(name1(:))
   len22(:) = len_trim(name2(:))
   len31 = len_trim(d1name)
   len32 = len_trim(d2name)
   len4 = len_trim(title)

   st = nf_create(wname(1:len1),nf_share,ncid)

   st = nf_def_dim(ncid,d1name(1:len31),dim1,d1id)
   st = nf_def_var(ncid,d1name(1:len31),nf_real,1,d1id,a1id)
   vardimid1(1) = d1id

   st = nf_def_dim(ncid,d2name(1:len32),dim2,d2id)
   st = nf_def_var(ncid,d2name(1:len32),nf_real,1,d2id,a2id)
   vardimid2(1) = d2id

   do i=1, nv1
     st = nf_def_var(ncid,trim(name1(i)),nf_real,1,vardimid1,varid1(i))
     st = nf_put_att_real(ncid,varid1(i),'_FillValue',nf_real,1,1.e+32)
   enddo
   do i=1, nv2
     st = nf_def_var(ncid,trim(name2(i)),nf_real,1,vardimid2,varid2(i))
     st = nf_put_att_real(ncid,varid2(i),'_FillValue',nf_real,1,1.e+32)
   enddo
   st = nf_put_att_text(ncid,nf_global,'title',len4,title)
   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a1id,axis1)
   st = nf_put_var_real(ncid,a2id,axis2)
   do i=1, nv1
     st = nf_put_var_real(ncid,varid1(i),var1(:,i))
   enddo
   do i=1, nv2
     st = nf_put_var_real(ncid,varid2(i),var2(:,i))
   enddo

   st = nf_close(ncid)

   end subroutine out1d2

   subroutine out2d(wname,nv,name,var,d1name,dim1,axis1,d2name,dim2,axis2, &
                    title)

   implicit none

   integer, intent(in) :: nv
   character(len=*), intent(in) :: wname,name(nv)
   character(len=*), intent(in) :: d1name, d2name 
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim1, dim2
   real, dimension(dim1), intent(in) :: axis1
   real, dimension(dim2), intent(in) :: axis2
   real, dimension(dim1,dim2,nv), intent(in) :: var

   integer :: i,j,k,it,st
   integer :: len1, len2(nv), len3, len4, len5, len6, len7
   integer :: ncid
   integer :: d1id, d2id, d3id, d4id
   integer :: a1id, a2id, a3id, a4id
   integer :: varid(nv)
   integer, dimension(2) :: vardimid

   len1 = len_trim(wname)
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

   do i=1, nv
     st = nf_def_var(ncid,trim(name(i)),nf_real,2,vardimid,varid(i))
     st = nf_put_att_real(ncid,varid(i),'_FillValue',nf_real,1,1.e+32)
     st = nf_put_att_text(ncid,nf_global,'title',len5,title)
   enddo
   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a1id,axis1)
   st = nf_put_var_real(ncid,a2id,axis2)
   do i=1, nv
     st = nf_put_var_real(ncid,varid(i),var(:,:,i))
   enddo

   st = nf_close(ncid)
   end subroutine out2d

   subroutine out2d2(wname,nv1,name1,var1,d1name,dim1,axis1,d2name,dim2,axis2, &
                    nv2,name2,var2,d3name,dim3,axis3,d4name,dim4,axis4,title)

   implicit none

   integer, intent(in) :: nv1, nv2
   character(len=*), intent(in) :: wname, name1(nv1), name2(nv2)
   character(len=*), intent(in) :: d1name, d2name, d3name, d4name
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim1, dim2, dim3, dim4
   real, dimension(dim1), intent(in) :: axis1
   real, dimension(dim2), intent(in) :: axis2
   real, dimension(dim3), intent(in) :: axis3
   real, dimension(dim4), intent(in) :: axis4
   real, dimension(dim1,dim2,nv1), intent(in) :: var1
   real, dimension(dim3,dim4,nv2), intent(in) :: var2

   integer :: i,j,k,it,st
   integer :: len1,len21(nv1),len22(nv2),len3,len4,len32,len42,len5,len6,len7
   integer :: ncid
   integer :: d1id, d2id, d3id, d4id
   integer :: a1id, a2id, a3id, a4id
   integer :: varid1(nv1), varid2(nv2)
   integer, dimension(2) :: vardimid1, vardimid2

   len1 = len_trim(wname)
   len21(:) = len_trim(name1(:))
   len22(:) = len_trim(name2(:))
   len3 = len_trim(d1name)
   len4 = len_trim(d2name)
   len32= len_trim(d3name)
   len42= len_trim(d4name)
   len5 = len_trim(title)

   st = nf_create(wname(1:len1),nf_share,ncid)

   st = nf_def_dim(ncid,d1name(1:len3),dim1,d1id)
   st = nf_def_dim(ncid,d2name(1:len4),dim2,d2id)
   st = nf_def_dim(ncid,d3name(1:len32),dim3,d3id)
   st = nf_def_dim(ncid,d4name(1:len42),dim4,d4id)

   st = nf_def_var(ncid,d1name(1:len3),nf_real,1,d1id,a1id)
   st = nf_def_var(ncid,d2name(1:len4),nf_real,1,d2id,a2id)
   st = nf_def_var(ncid,d3name(1:len32),nf_real,1,d3id,a3id)
   st = nf_def_var(ncid,d4name(1:len42),nf_real,1,d4id,a4id)

   vardimid1(1) = d1id
   vardimid1(2) = d2id
   vardimid2(1) = d3id
   vardimid2(2) = d4id

   do i=1, nv1
     st = nf_def_var(ncid,trim(name1(i)),nf_real,2,vardimid1,varid1(i))
     st = nf_put_att_real(ncid,varid1(i),'_FillValue',nf_real,1,1.e+32)
   enddo
   do i=1, nv2
     st = nf_def_var(ncid,trim(name2(i)),nf_real,2,vardimid2,varid2(i))
     st = nf_put_att_real(ncid,varid2(i),'_FillValue',nf_real,1,1.e+32)
   enddo
   st = nf_put_att_text(ncid,nf_global,'title',len5,title)
   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a1id,axis1)
   st = nf_put_var_real(ncid,a2id,axis2)
   st = nf_put_var_real(ncid,a3id,axis3)
   st = nf_put_var_real(ncid,a4id,axis4)
   do i=1, nv1
     st = nf_put_var_real(ncid,varid1(i),var1(:,:,i))
   enddo
   do i=1, nv2
     st = nf_put_var_real(ncid,varid2(i),var2(:,:,i))
   enddo

   st = nf_close(ncid)
   end subroutine out2d2

   subroutine out2d1d(wname,nv2,name2,var2, &
                      d2name1,dim21,axis21, &
                      d2name2,dim22,axis22, &
                      nv1,name1,var1,       &
                      d1name1,dim11,axis11, title)
   
   implicit none

   integer, intent(in) :: nv2, nv1
   character(len=*), intent(in) :: wname, name2(nv2), name1(nv1)
   character(len=*), intent(in) :: d2name1, d2name2
   character(len=*), intent(in) :: d1name1
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim21, dim22, dim11
   real, intent(in)    :: axis21(dim21), axis22(dim22)
   real, intent(in)    :: axis11(dim11)
   real, dimension(dim21,dim22,nv2), intent(in) :: var2
   real, dimension(dim11,nv1),       intent(in) :: var1
   
   integer :: i,j,k,it,st
   integer :: ncid
   integer :: d2id1, d2id2, d1id1
   integer :: a2id1, a2id2, a1id1
   integer :: varid2(nv2), varid1(nv1)
   integer :: vardimid2(2), vardimid1(1)


   st = nf_create(trim(wname),nf_share,ncid)

   st = nf_def_dim(ncid,trim(d2name1),dim21,d2id1)
   st = nf_def_dim(ncid,trim(d2name2),dim22,d2id2)
   st = nf_def_dim(ncid,trim(d1name1),dim11,d1id1)

   st = nf_def_var(ncid,trim(d2name1),nf_real,1,d2id1,a2id1)
   st = nf_def_var(ncid,trim(d2name2),nf_real,1,d2id2,a2id2)
   st = nf_def_var(ncid,trim(d1name1),nf_real,1,d1id1,a1id1)


   vardimid2 = (/d2id1, d2id2/)
   vardimid1 = (/d1id1/)

   do i=1, nv2
     st = nf_def_var(ncid,trim(name2(i)),nf_real,2,vardimid2,varid2(i))
     st = nf_put_att_real(ncid,varid2(i),'_FillValue',nf_real,1,1.e+32)
   enddo
   do i=1, nv1
     st = nf_def_var(ncid,trim(name1(i)),nf_real,1,vardimid1,varid1(i))
     st = nf_put_att_real(ncid,varid1(i),'_FillValue',nf_real,1,1.e+32)
   enddo
   st = nf_put_att_text(ncid,nf_global,'title',len_trim(title),title)
   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a2id1,axis21)
   st = nf_put_var_real(ncid,a2id2,axis22)
   st = nf_put_var_real(ncid,a1id1,axis11)

   do i=1, nv2
     st = nf_put_var_real(ncid,varid2(i),var2(:,:,i))
   enddo
   do i=1, nv1
     st = nf_put_var_real(ncid,varid1(i),var1(:,i))
   enddo

   st = nf_close(ncid)
   end subroutine out2d1d

   subroutine out3d(wname,nv3,name3,var3, &
                    d3name1,dim31,axis31, &
                    d3name2,dim32,axis32, &
                    d3name3,dim33,axis33, title)

   implicit none

   integer, intent(in) :: nv3
   character(len=*), intent(in) :: wname, name3(nv3)
   character(len=*), intent(in) :: d3name1, d3name2, d3name3
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim31, dim32, dim33
   real, intent(in)    :: axis31(dim31), axis32(dim32), axis33(dim33)
   real, dimension(dim31,dim32,dim33,nv3), intent(in) :: var3

   integer :: i,j,k,it,st
   integer :: ncid
   integer :: d3id1, d3id2, d3id3
   integer :: a3id1, a3id2, a3id3
   integer :: varid3(nv3)
   integer :: vardimid3(3)


   st = nf_create(trim(wname),nf_share,ncid)

   st = nf_def_dim(ncid,trim(d3name1),dim31,d3id1)
   st = nf_def_dim(ncid,trim(d3name2),dim32,d3id2)
   st = nf_def_dim(ncid,trim(d3name3),dim33,d3id3)

   st = nf_def_var(ncid,trim(d3name1),nf_real,1,d3id1,a3id1)
   st = nf_def_var(ncid,trim(d3name2),nf_real,1,d3id2,a3id2)
   st = nf_def_var(ncid,trim(d3name3),nf_real,1,d3id3,a3id3)


   vardimid3 = (/d3id1, d3id2, d3id3/)

   do i=1, nv3
     st = nf_def_var(ncid,trim(name3(i)),nf_real,3,vardimid3,varid3(i))
     st = nf_put_att_real(ncid,varid3(i),'_FillValue',nf_real,1,1.e+32)
   enddo
   st = nf_put_att_text(ncid,nf_global,'title',len_trim(title),title)
   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a3id1,axis31)
   st = nf_put_var_real(ncid,a3id2,axis32)
   st = nf_put_var_real(ncid,a3id3,axis33)

   do i=1, nv3
     st = nf_put_var_real(ncid,varid3(i),var3(:,:,:,i))
   enddo

   st = nf_close(ncid)
   end subroutine out3d

   subroutine out3d1d(wname,nv3,name3,var3, &
                      d3name1,dim31,axis31, &
                      d3name2,dim32,axis32, &
                      d3name3,dim33,axis33, &
                      nv1,name1,var1,       &
                      d1name1,dim11,axis11, title)

   implicit none

   integer, intent(in) :: nv3, nv1
   character(len=*), intent(in) :: wname, name3(nv3), name1(nv1)
   character(len=*), intent(in) :: d3name1, d3name2, d3name3
   character(len=*), intent(in) :: d1name1
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim31, dim32, dim33, dim11
   real, intent(in)    :: axis31(dim31), axis32(dim32), axis33(dim33)
   real, intent(in)    :: axis11(dim11)
   real, dimension(dim31,dim32,dim33,nv3), intent(in) :: var3
   real, dimension(dim11,nv1),             intent(in) :: var1

   integer :: i,j,k,it,st
   integer :: ncid
   integer :: d3id1, d3id2, d3id3, d1id1
   integer :: a3id1, a3id2, a3id3, a1id1
   integer :: varid3(nv3), varid1(nv1)
   integer :: vardimid3(3), vardimid1(1)


   st = nf_create(trim(wname),nf_share,ncid)

   st = nf_def_dim(ncid,trim(d3name1),dim31,d3id1)
   st = nf_def_dim(ncid,trim(d3name2),dim32,d3id2)
   st = nf_def_dim(ncid,trim(d3name3),dim33,d3id3)
   st = nf_def_dim(ncid,trim(d1name1),dim11,d1id1)

   st = nf_def_var(ncid,trim(d3name1),nf_real,1,d3id1,a3id1)
   st = nf_def_var(ncid,trim(d3name2),nf_real,1,d3id2,a3id2)
   st = nf_def_var(ncid,trim(d3name3),nf_real,1,d3id3,a3id3)
   st = nf_def_var(ncid,trim(d1name1),nf_real,1,d1id1,a1id1)


   vardimid3 = (/d3id1, d3id2, d3id3/)
   vardimid1 = (/d1id1/)

   do i=1, nv3
     st = nf_def_var(ncid,trim(name3(i)),nf_real,3,vardimid3,varid3(i))
     st = nf_put_att_real(ncid,varid3(i),'_FillValue',nf_real,1,1.e+32)
   enddo
   do i=1, nv1
     st = nf_def_var(ncid,trim(name1(i)),nf_real,1,vardimid1,varid1(i))
     st = nf_put_att_real(ncid,varid1(i),'_FillValue',nf_real,1,1.e+32)
   enddo
   st = nf_put_att_text(ncid,nf_global,'title',len_trim(title),title)
   st = nf_enddef(ncid)
 
   st = nf_put_var_real(ncid,a3id1,axis31)
   st = nf_put_var_real(ncid,a3id2,axis32)
   st = nf_put_var_real(ncid,a3id3,axis33)
   st = nf_put_var_real(ncid,a1id1,axis11)
 
   do i=1, nv3
     st = nf_put_var_real(ncid,varid3(i),var3(:,:,:,i))
   enddo
   do i=1, nv1
     st = nf_put_var_real(ncid,varid1(i),var1(:,i))
   enddo
 
   st = nf_close(ncid)
   end subroutine out3d1d

   subroutine out3d2d(wname,nv3,name3,var3, &
                      d3name1,dim31,axis31, &
                      d3name2,dim32,axis32, &
                      d3name3,dim33,axis33, &
                      nv2,name2,var2,       &
                      d2name1,dim21,axis21, &
                      d2name2,dim22,axis22, title)

   implicit none

   integer, intent(in) :: nv3, nv2
   character(len=*), intent(in) :: wname, name3(nv3), name2(nv2)
   character(len=*), intent(in) :: d3name1, d3name2, d3name3
   character(len=*), intent(in) :: d2name1, d2name2
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim31, dim32, dim33, dim21, dim22
   real, intent(in)    :: axis31(dim31), axis32(dim32), axis33(dim33)
   real, intent(in)    :: axis21(dim21), axis22(dim22)
   real, dimension(dim31,dim32,dim33,nv3), intent(in) :: var3
   real, dimension(dim21,dim22,nv2),       intent(in) :: var2

   integer :: i,j,k,it,st
   integer :: ncid
   integer :: d3id1, d3id2, d3id3, d2id1, d2id2
   integer :: a3id1, a3id2, a3id3, a2id1, a2id2
   integer :: varid3(nv3), varid2(nv2)
   integer :: vardimid3(3), vardimid2(2)


   st = nf_create(trim(wname),nf_share,ncid)

   st = nf_def_dim(ncid,trim(d3name1),dim31,d3id1)
   st = nf_def_dim(ncid,trim(d3name2),dim32,d3id2)
   st = nf_def_dim(ncid,trim(d3name3),dim33,d3id3)
   st = nf_def_dim(ncid,trim(d2name1),dim21,d2id1)
   st = nf_def_dim(ncid,trim(d2name2),dim22,d2id2)

   st = nf_def_var(ncid,trim(d3name1),nf_real,1,d3id1,a3id1)
   st = nf_def_var(ncid,trim(d3name2),nf_real,1,d3id2,a3id2)
   st = nf_def_var(ncid,trim(d3name3),nf_real,1,d3id3,a3id3)
   st = nf_def_var(ncid,trim(d2name1),nf_real,1,d2id1,a2id1)
   st = nf_def_var(ncid,trim(d2name2),nf_real,1,d2id2,a2id2)


   vardimid3 = (/d3id1, d3id2, d3id3/)
   vardimid2 = (/d2id1, d2id2/)

   do i=1, nv3
     st = nf_def_var(ncid,trim(name3(i)),nf_real,3,vardimid3,varid3(i))
     st = nf_put_att_real(ncid,varid3(i),'_FillValue',nf_real,1,1.e+32)
   enddo
   do i=1, nv2
     st = nf_def_var(ncid,trim(name2(i)),nf_real,2,vardimid2,varid2(i))
     st = nf_put_att_real(ncid,varid2(i),'_FillValue',nf_real,1,1.e+32)
   enddo
   st = nf_put_att_text(ncid,nf_global,'title',len_trim(title),title)
   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a3id1,axis31)
   st = nf_put_var_real(ncid,a3id2,axis32)
   st = nf_put_var_real(ncid,a3id3,axis33)
   st = nf_put_var_real(ncid,a2id1,axis21)
   st = nf_put_var_real(ncid,a2id2,axis22)

   do i=1, nv3
     st = nf_put_var_real(ncid,varid3(i),var3(:,:,:,i))
   enddo
   do i=1, nv2
     st = nf_put_var_real(ncid,varid2(i),var2(:,:,i))
   enddo

   st = nf_close(ncid)
   end subroutine out3d2d

   subroutine out3d2(wname,nv3,name3,var3, &
                     d3name1,dim31,axis31, &
                     d3name2,dim32,axis32, &
                     d3name3,dim33,axis33, &
                     nv3b,name3b,var3b,       &
                     d3bname1,dim3b1,axis3b1, &
                     d3bname2,dim3b2,axis3b2, &
                     d3bname3,dim3b3,axis3b3, title)

   implicit none

   integer, intent(in) :: nv3, nv3b
   character(len=*), intent(in) :: wname, name3(nv3), name3b(nv3b)
   character(len=*), intent(in) :: d3name1, d3name2, d3name3
   character(len=*), intent(in) :: d3bname1, d3bname2, d3bname3
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim31, dim32, dim33, dim3b1, dim3b2, dim3b3
   real, intent(in)    :: axis31(dim31), axis32(dim32), axis33(dim33)
   real, intent(in)    :: axis3b1(dim3b1), axis3b2(dim3b2), axis3b3(dim3b3)
   real, dimension(dim31,dim32,dim33,nv3),     intent(in) :: var3
   real, dimension(dim3b1,dim3b2,dim3b3,nv3b), intent(in) :: var3b

   integer :: i,j,k,it,st
   integer :: ncid
   integer :: d3id1, d3id2, d3id3, d3bid1, d3bid2, d3bid3
   integer :: a3id1, a3id2, a3id3, a3bid1, a3bid2, a3bid3
   integer :: varid3(nv3), varid3b(nv3b)
   integer :: vardimid3(3), vardimid3b(3)


   st = nf_create(trim(wname),nf_share,ncid)

   st = nf_def_dim(ncid,trim(d3name1),dim31,d3id1)
   st = nf_def_dim(ncid,trim(d3name2),dim32,d3id2)
   st = nf_def_dim(ncid,trim(d3name3),dim33,d3id3)
   st = nf_def_dim(ncid,trim(d3bname1),dim3b1,d3bid1)
   st = nf_def_dim(ncid,trim(d3bname2),dim3b2,d3bid2)
   st = nf_def_dim(ncid,trim(d3bname3),dim3b3,d3bid3)

   st = nf_def_var(ncid,trim(d3name1),nf_real,1,d3id1,a3id1)
   st = nf_def_var(ncid,trim(d3name2),nf_real,1,d3id2,a3id2)
   st = nf_def_var(ncid,trim(d3name3),nf_real,1,d3id3,a3id3)
   st = nf_def_var(ncid,trim(d3bname1),nf_real,1,d3bid1,a3bid1)
   st = nf_def_var(ncid,trim(d3bname2),nf_real,1,d3bid2,a3bid2)
   st = nf_def_var(ncid,trim(d3bname3),nf_real,1,d3bid3,a3bid3)


   vardimid3 = (/d3id1, d3id2, d3id3/)
   vardimid3b = (/d3bid1, d3bid2, d3bid3/)

   do i=1, nv3
     st = nf_def_var(ncid,trim(name3(i)),nf_real,3,vardimid3,varid3(i))
     st = nf_put_att_real(ncid,varid3(i),'_FillValue',nf_real,1,1.e+32)
   enddo
   do i=1, nv3b
     st = nf_def_var(ncid,trim(name3b(i)),nf_real,3,vardimid3b,varid3b(i))
     st = nf_put_att_real(ncid,varid3b(i),'_FillValue',nf_real,1,1.e+32)
   enddo
   st = nf_put_att_text(ncid,nf_global,'title',len_trim(title),title)
   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a3id1,axis31)
   st = nf_put_var_real(ncid,a3id2,axis32)
   st = nf_put_var_real(ncid,a3id3,axis33)
   st = nf_put_var_real(ncid,a3bid1,axis3b1)
   st = nf_put_var_real(ncid,a3bid2,axis3b2)
   st = nf_put_var_real(ncid,a3bid3,axis3b3)

   do i=1, nv3
     st = nf_put_var_real(ncid,varid3(i),var3(:,:,:,i))
   enddo
   do i=1, nv3b
     st = nf_put_var_real(ncid,varid3b(i),var3b(:,:,:,i))
   enddo

   st = nf_close(ncid)
   end subroutine out3d2

   subroutine out3d2d1d(wname,nv3,name3,var3, &
                        d3name1,dim31,axis31, &
                        d3name2,dim32,axis32, &
                        d3name3,dim33,axis33, &
                        nv2,name2,var2,       &
                        d2name1,dim21,axis21, &
                        d2name2,dim22,axis22, &
                        nv1,name1,var1,       &
                        d1name1,dim11,axis11, title)

   implicit none

   integer, intent(in) :: nv3, nv2, nv1
   character(len=*), intent(in) :: wname, name3(nv3), name2(nv2), name1(nv1)
   character(len=*), intent(in) :: d3name1, d3name2, d3name3
   character(len=*), intent(in) :: d2name1, d2name2
   character(len=*), intent(in) :: d1name1
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim31, dim32, dim33, dim21, dim22, dim11
   real, intent(in)    :: axis31(dim31), axis32(dim32), axis33(dim33)
   real, intent(in)    :: axis21(dim21), axis22(dim22)
   real, intent(in)    :: axis11(dim11)
   real, dimension(dim31,dim32,dim33,nv3), intent(in) :: var3
   real, dimension(dim21,dim22,nv2),       intent(in) :: var2
   real, dimension(dim11,nv1),             intent(in) :: var1

   integer :: i,j,k,it,st
   integer :: ncid
   integer :: d3id1, d3id2, d3id3, d2id1, d2id2, d1id1
   integer :: a3id1, a3id2, a3id3, a2id1, a2id2, a1id1
   integer :: varid3(nv3), varid2(nv2), varid1(nv1)
   integer :: vardimid3(3), vardimid2(2), vardimid1(1)


   st = nf_create(trim(wname),nf_share,ncid)

   st = nf_def_dim(ncid,trim(d3name1),dim31,d3id1)
   st = nf_def_dim(ncid,trim(d3name2),dim32,d3id2)
   st = nf_def_dim(ncid,trim(d3name3),dim33,d3id3)
   st = nf_def_dim(ncid,trim(d2name1),dim21,d2id1)
   st = nf_def_dim(ncid,trim(d2name2),dim22,d2id2)
   st = nf_def_dim(ncid,trim(d1name1),dim11,d1id1)

   st = nf_def_var(ncid,trim(d3name1),nf_real,1,d3id1,a3id1)
   st = nf_def_var(ncid,trim(d3name2),nf_real,1,d3id2,a3id2)
   st = nf_def_var(ncid,trim(d3name3),nf_real,1,d3id3,a3id3)
   st = nf_def_var(ncid,trim(d2name1),nf_real,1,d2id1,a2id1)
   st = nf_def_var(ncid,trim(d2name2),nf_real,1,d2id2,a2id2)
   st = nf_def_var(ncid,trim(d1name1),nf_real,1,d1id1,a1id1)


   vardimid3 = (/d3id1, d3id2, d3id3/)
   vardimid2 = (/d2id1, d2id2/)
   vardimid1 = (/d1id1/)

   do i=1, nv3
     st = nf_def_var(ncid,trim(name3(i)),nf_real,3,vardimid3,varid3(i))
     st = nf_put_att_real(ncid,varid3(i),'_FillValue',nf_real,1,1.e+32)
   enddo
   do i=1, nv2
     st = nf_def_var(ncid,trim(name2(i)),nf_real,2,vardimid2,varid2(i))
     st = nf_put_att_real(ncid,varid2(i),'_FillValue',nf_real,1,1.e+32)
   enddo
   do i=1, nv1
     st = nf_def_var(ncid,trim(name1(i)),nf_real,1,vardimid1,varid1(i))
     st = nf_put_att_real(ncid,varid1(i),'_FillValue',nf_real,1,1.e+32)
   enddo 

   st = nf_put_att_text(ncid,nf_global,'title',len_trim(title),title)
   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a3id1,axis31)
   st = nf_put_var_real(ncid,a3id2,axis32)
   st = nf_put_var_real(ncid,a3id3,axis33)
   st = nf_put_var_real(ncid,a2id1,axis21)
   st = nf_put_var_real(ncid,a2id2,axis22)
   st = nf_put_var_real(ncid,a1id1,axis11)

   do i=1, nv3
     st = nf_put_var_real(ncid,varid3(i),var3(:,:,:,i))
   enddo
   do i=1, nv2
     st = nf_put_var_real(ncid,varid2(i),var2(:,:,i))
   enddo
   do i=1, nv1
     st = nf_put_var_real(ncid,varid1(i),var1(:,i))
   enddo

   st = nf_close(ncid)
   end subroutine out3d2d1d

   subroutine out3d2d2d(wname,nv3,name3,var3,    &
                        d3name1,dim31,axis31,    &
                        d3name2,dim32,axis32,    &
                        d3name3,dim33,axis33,    &
                        nv2,name2,var2,          &
                        d2name1,dim21,axis21,    &
                        d2name2,dim22,axis22,    &
                        nv2b,name2b,var2b,       &
                        d2name1b,dim21b,axis21b, &
                        d2name2b,dim22b,axis22b,title)

   implicit none

   integer, intent(in) :: nv3, nv2, nv2b
   character(len=*), intent(in) :: wname, name3(nv3), name2(nv2), name2b(nv2b)
   character(len=*), intent(in) :: d3name1, d3name2, d3name3
   character(len=*), intent(in) :: d2name1, d2name2
   character(len=*), intent(in) :: d2name1b, d2name2b
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim31, dim32, dim33, dim21, dim22, dim21b, dim22b
   real, intent(in)    :: axis31(dim31), axis32(dim32), axis33(dim33)
   real, intent(in)    :: axis21(dim21), axis22(dim22)
   real, intent(in)    :: axis21b(dim21b), axis22b(dim22b)
   real, dimension(dim31,dim32,dim33,nv3), intent(in) :: var3
   real, dimension(dim21,dim22,nv2),       intent(in) :: var2
   real, dimension(dim21b,dim22b,nv2b),    intent(in) :: var2b

   integer :: i,j,k,it,st
   integer :: ncid
   integer :: d3id1, d3id2, d3id3, d2id1, d2id2, d2id1b, d2id2b
   integer :: a3id1, a3id2, a3id3, a2id1, a2id2, a2id1b, a2id2b
   integer :: varid3(nv3), varid2(nv2), varid2b(nv2b)
   integer :: vardimid3(3), vardimid2(2), vardimid2b(2)


   st = nf_create(trim(wname),nf_share,ncid)

   st = nf_def_dim(ncid,trim(d3name1),dim31,d3id1)
   st = nf_def_dim(ncid,trim(d3name2),dim32,d3id2)
   st = nf_def_dim(ncid,trim(d3name3),dim33,d3id3)
   st = nf_def_dim(ncid,trim(d2name1),dim21,d2id1)
   st = nf_def_dim(ncid,trim(d2name2),dim22,d2id2)
   st = nf_def_dim(ncid,trim(d2name1b),dim21b,d2id1b)
   st = nf_def_dim(ncid,trim(d2name2b),dim22b,d2id2b)

   st = nf_def_var(ncid,trim(d3name1),nf_real,1,d3id1,a3id1)
   st = nf_def_var(ncid,trim(d3name2),nf_real,1,d3id2,a3id2)
   st = nf_def_var(ncid,trim(d3name3),nf_real,1,d3id3,a3id3)
   st = nf_def_var(ncid,trim(d2name1),nf_real,1,d2id1,a2id1)
   st = nf_def_var(ncid,trim(d2name2),nf_real,1,d2id2,a2id2)
   st = nf_def_var(ncid,trim(d2name1b),nf_real,1,d2id1b,a2id1b)
   st = nf_def_var(ncid,trim(d2name2b),nf_real,1,d2id2b,a2id2b)


   vardimid3 = (/d3id1, d3id2, d3id3/)
   vardimid2 = (/d2id1, d2id2/)
   vardimid2b = (/d2id1b, d2id2b/)

   do i=1, nv3
     st = nf_def_var(ncid,trim(name3(i)),nf_real,3,vardimid3,varid3(i))
     st = nf_put_att_real(ncid,varid3(i),'_FillValue',nf_real,1,1.e+32)
   enddo
   do i=1, nv2
     st = nf_def_var(ncid,trim(name2(i)),nf_real,2,vardimid2,varid2(i))
     st = nf_put_att_real(ncid,varid2(i),'_FillValue',nf_real,1,1.e+32)
   enddo
   do i=1, nv2b 
     st = nf_def_var(ncid,trim(name2b(i)),nf_real,2,vardimid2b,varid2b(i))
     st = nf_put_att_real(ncid,varid2b(i),'_FillValue',nf_real,1,1.e+32)
   enddo
   st = nf_put_att_text(ncid,nf_global,'title',len_trim(title),title)
   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a3id1,axis31)
   st = nf_put_var_real(ncid,a3id2,axis32)
   st = nf_put_var_real(ncid,a3id3,axis33)
   st = nf_put_var_real(ncid,a2id1,axis21)
   st = nf_put_var_real(ncid,a2id2,axis22)
   st = nf_put_var_real(ncid,a2id1b,axis21b)
   st = nf_put_var_real(ncid,a2id2b,axis22b)

   do i=1, nv3
     st = nf_put_var_real(ncid,varid3(i),var3(:,:,:,i))
   enddo
   do i=1, nv2
     st = nf_put_var_real(ncid,varid2(i),var2(:,:,i))
   enddo
   do i=1, nv2b
     st = nf_put_var_real(ncid,varid2b(i),var2b(:,:,i))
   enddo

   st = nf_close(ncid)
   end subroutine out3d2d2d

   subroutine out4d(wname,nv4,name4,var4, &
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
   end subroutine out4d

   subroutine out4d1d(wname,nv4,name4,var4, &
                      d4name1,dim41,axis41, &
                      d4name2,dim42,axis42, &
                      d4name3,dim43,axis43, &
                      d4name4,dim44,axis44, &
                      nv1,name1,var1,       &
                      d1name1,dim11,axis11, title)

   implicit none

   integer, intent(in) :: nv4, nv1
   character(len=*), intent(in) :: wname, name4(nv4), name1(nv1)
   character(len=*), intent(in) :: d4name1, d4name2, d4name3, d4name4
   character(len=*), intent(in) :: d1name1
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim41, dim42, dim43, dim44
   integer, intent(in) :: dim11
   real, intent(in)    :: axis41(dim41), axis42(dim42)
   real, intent(in)    :: axis43(dim43), axis44(dim44)
   real, intent(in)    :: axis11(dim11)
   real, dimension(dim41,dim42,dim43,dim44,nv4), intent(in) :: var4
   real, dimension(dim11,nv1),                   intent(in) :: var1

   integer :: i,j,k,it,st
   integer :: ncid
   integer :: d4id1, d4id2, d4id3, d4id4
   integer :: d1id1
   integer :: a4id1, a4id2, a4id3, a4id4
   integer :: a1id1
   integer :: varid4(nv4), varid1(nv1)
   integer :: vardimid4(4), vardimid1(1)


   st = nf_create(trim(wname),nf_share,ncid)

   st = nf_def_dim(ncid,trim(d4name1),dim41,d4id1)
   st = nf_def_dim(ncid,trim(d4name2),dim42,d4id2)
   st = nf_def_dim(ncid,trim(d4name3),dim43,d4id3)
   st = nf_def_dim(ncid,trim(d4name4),dim44,d4id4)
   st = nf_def_dim(ncid,trim(d1name1),dim11,d1id1)

   st = nf_def_var(ncid,trim(d4name1),nf_real,1,d4id1,a4id1)
   st = nf_def_var(ncid,trim(d4name2),nf_real,1,d4id2,a4id2)
   st = nf_def_var(ncid,trim(d4name3),nf_real,1,d4id3,a4id3)
   st = nf_def_var(ncid,trim(d4name4),nf_real,1,d4id4,a4id4)
   st = nf_def_var(ncid,trim(d1name1),nf_real,1,d1id1,a1id1)


   vardimid4 = (/d4id1, d4id2, d4id3,d4id4/)
   vardimid1 = (/d1id1/)

   do i=1, nv4
     st = nf_def_var(ncid,trim(name4(i)),nf_real,4,vardimid4,varid4(i))
     st = nf_put_att_real(ncid,varid4(i),'_FillValue',nf_real,1,1.e+32)
   enddo
   do i=1, nv1
     st = nf_def_var(ncid,trim(name1(i)),nf_real,1,vardimid1,varid1(i))
     st = nf_put_att_real(ncid,varid1(i),'_FillValue',nf_real,1,1.e+32)
   enddo

   st = nf_put_att_text(ncid,nf_global,'title',len_trim(title),title)
   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a4id1,axis41)
   st = nf_put_var_real(ncid,a4id2,axis42)
   st = nf_put_var_real(ncid,a4id3,axis43)
   st = nf_put_var_real(ncid,a4id4,axis44)
   st = nf_put_var_real(ncid,a1id1,axis11)

   do i=1, nv4
     st = nf_put_var_real(ncid,varid4(i),var4(:,:,:,:,i))
   enddo
   do i=1, nv1
     st = nf_put_var_real(ncid,varid1(i),var1(:,i))
   enddo

   st = nf_close(ncid)
   end subroutine out4d1d

   subroutine out4d3d2d(wname,nv4,name4,var4, &
                        d4name1,dim41,axis41, &
                        d4name2,dim42,axis42, &
                        d4name3,dim43,axis43, &
                        d4name4,dim44,axis44, &
                        nv3,name3,var3,       &
                        d3name1,dim31,axis31, &
                        d3name2,dim32,axis32, &
                        d3name3,dim33,axis33, &
                        nv2,name2,var2,       &
                        d2name1,dim21,axis21, &
                        d2name2,dim22,axis22, title)

   implicit none

   integer, intent(in) :: nv4, nv3, nv2
   character(len=*), intent(in) :: wname, name4(nv4), name3(nv3)
   character(len=*), intent(in) :: name2(nv2)
   character(len=*), intent(in) :: d4name1, d4name2, d4name3, d4name4
   character(len=*), intent(in) :: d3name1, d3name2, d3name3
   character(len=*), intent(in) :: d2name1, d2name2
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim41, dim42, dim43, dim44
   integer, intent(in) :: dim31, dim32, dim33, dim21, dim22
   real, intent(in)    :: axis41(dim41), axis42(dim42)
   real, intent(in)    :: axis43(dim43), axis44(dim44)
   real, intent(in)    :: axis31(dim31), axis32(dim32), axis33(dim33)
   real, intent(in)    :: axis21(dim21), axis22(dim22)
   real, dimension(dim41,dim42,dim43,dim44,nv4), intent(in) :: var4
   real, dimension(dim31,dim32,dim33,nv3),       intent(in) :: var3
   real, dimension(dim21,dim22,nv2),             intent(in) :: var2

   integer :: i,j,k,it,st
   integer :: ncid
   integer :: d4id1, d4id2, d4id3, d4id4
   integer :: d3id1, d3id2, d3id3, d2id1, d2id2
   integer :: a4id1, a4id2, a4id3, a4id4
   integer :: a3id1, a3id2, a3id3, a2id1, a2id2
   integer :: varid4(nv4), varid3(nv3), varid2(nv2)
   integer :: vardimid4(4), vardimid3(3), vardimid2(2)


   st = nf_create(trim(wname),nf_share,ncid)

   st = nf_def_dim(ncid,trim(d4name1),dim41,d4id1)
   st = nf_def_dim(ncid,trim(d4name2),dim42,d4id2)
   st = nf_def_dim(ncid,trim(d4name3),dim43,d4id3)
   st = nf_def_dim(ncid,trim(d4name4),dim44,d4id4)
   st = nf_def_dim(ncid,trim(d3name1),dim31,d3id1)
   st = nf_def_dim(ncid,trim(d3name2),dim32,d3id2)
   st = nf_def_dim(ncid,trim(d3name3),dim33,d3id3)
   st = nf_def_dim(ncid,trim(d2name1),dim21,d2id1)
   st = nf_def_dim(ncid,trim(d2name2),dim22,d2id2)

   st = nf_def_var(ncid,trim(d4name1),nf_real,1,d4id1,a4id1)
   st = nf_def_var(ncid,trim(d4name2),nf_real,1,d4id2,a4id2)
   st = nf_def_var(ncid,trim(d4name3),nf_real,1,d4id3,a4id3)
   st = nf_def_var(ncid,trim(d4name4),nf_real,1,d4id4,a4id4)
   st = nf_def_var(ncid,trim(d3name1),nf_real,1,d3id1,a3id1)
   st = nf_def_var(ncid,trim(d3name2),nf_real,1,d3id2,a3id2)
   st = nf_def_var(ncid,trim(d3name3),nf_real,1,d3id3,a3id3)
   st = nf_def_var(ncid,trim(d2name1),nf_real,1,d2id1,a2id1)
   st = nf_def_var(ncid,trim(d2name2),nf_real,1,d2id2,a2id2)


   vardimid4 = (/d4id1, d4id2, d4id3,d4id4/)
   vardimid3 = (/d3id1, d3id2, d3id3/)
   vardimid2 = (/d2id1, d2id2/)

   do i=1, nv4
     st = nf_def_var(ncid,trim(name4(i)),nf_real,4,vardimid4,varid4(i))
     st = nf_put_att_real(ncid,varid4(i),'_FillValue',nf_real,1,1.e+32)
   enddo
   do i=1, nv3
     st = nf_def_var(ncid,trim(name3(i)),nf_real,3,vardimid3,varid3(i))
     st = nf_put_att_real(ncid,varid3(i),'_FillValue',nf_real,1,1.e+32)
   enddo
   do i=1, nv2
     st = nf_def_var(ncid,trim(name2(i)),nf_real,2,vardimid2,varid2(i))
     st = nf_put_att_real(ncid,varid2(i),'_FillValue',nf_real,1,1.e+32)
   enddo

   st = nf_put_att_text(ncid,nf_global,'title',len_trim(title),title)
   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a4id1,axis41)
   st = nf_put_var_real(ncid,a4id2,axis42)
   st = nf_put_var_real(ncid,a4id3,axis43)
   st = nf_put_var_real(ncid,a4id4,axis44)
   st = nf_put_var_real(ncid,a3id1,axis31)
   st = nf_put_var_real(ncid,a3id2,axis32)
   st = nf_put_var_real(ncid,a3id3,axis33)
   st = nf_put_var_real(ncid,a2id1,axis21)
   st = nf_put_var_real(ncid,a2id2,axis22)

   do i=1, nv4
     st = nf_put_var_real(ncid,varid4(i),var4(:,:,:,:,i))
   enddo
   do i=1, nv3
     st = nf_put_var_real(ncid,varid3(i),var3(:,:,:,i))
   enddo
   do i=1, nv2
     st = nf_put_var_real(ncid,varid2(i),var2(:,:,i))
   enddo

   st = nf_close(ncid)
   end subroutine out4d3d2d

   subroutine out4d3d2d1d(wname,nv4,name4,var4, &
                          d4name1,dim41,axis41, &
                          d4name2,dim42,axis42, &
                          d4name3,dim43,axis43, &
                          d4name4,dim44,axis44, &
                          nv3,name3,var3,       &
                          d3name1,dim31,axis31, &
                          d3name2,dim32,axis32, &
                          d3name3,dim33,axis33, &
                          nv2,name2,var2,       &
                          d2name1,dim21,axis21, &
                          d2name2,dim22,axis22, &
                          nv1,name1,var1,       &
                          d1name1,dim11,axis11, title)

   implicit none

   integer, intent(in) :: nv4, nv3, nv2, nv1
   character(len=*), intent(in) :: wname, name4(nv4), name3(nv3)
   character(len=*), intent(in) :: name2(nv2), name1(nv1)
   character(len=*), intent(in) :: d4name1, d4name2, d4name3, d4name4
   character(len=*), intent(in) :: d3name1, d3name2, d3name3
   character(len=*), intent(in) :: d2name1, d2name2
   character(len=*), intent(in) :: d1name1
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim41, dim42, dim43, dim44
   integer, intent(in) :: dim31, dim32, dim33, dim21, dim22, dim11
   real, intent(in)    :: axis41(dim41), axis42(dim42)
   real, intent(in)    :: axis43(dim43), axis44(dim44)
   real, intent(in)    :: axis31(dim31), axis32(dim32), axis33(dim33)
   real, intent(in)    :: axis21(dim21), axis22(dim22)
   real, intent(in)    :: axis11(dim11)
   real, dimension(dim41,dim42,dim43,dim44,nv4), intent(in) :: var4
   real, dimension(dim31,dim32,dim33,nv3),       intent(in) :: var3
   real, dimension(dim21,dim22,nv2),             intent(in) :: var2
   real, dimension(dim11,nv1),                   intent(in) :: var1

   integer :: i,j,k,it,st
   integer :: ncid
   integer :: d4id1, d4id2, d4id3, d4id4
   integer :: d3id1, d3id2, d3id3, d2id1, d2id2, d1id1
   integer :: a4id1, a4id2, a4id3, a4id4
   integer :: a3id1, a3id2, a3id3, a2id1, a2id2, a1id1
   integer :: varid4(nv4), varid3(nv3), varid2(nv2), varid1(nv1)
   integer :: vardimid4(4), vardimid3(3), vardimid2(2), vardimid1(1)


   st = nf_create(trim(wname),nf_share,ncid)

   st = nf_def_dim(ncid,trim(d4name1),dim41,d4id1) 
   st = nf_def_dim(ncid,trim(d4name2),dim42,d4id2) 
   st = nf_def_dim(ncid,trim(d4name3),dim43,d4id3) 
   st = nf_def_dim(ncid,trim(d4name4),dim44,d4id4)
   st = nf_def_dim(ncid,trim(d3name1),dim31,d3id1)
   st = nf_def_dim(ncid,trim(d3name2),dim32,d3id2)
   st = nf_def_dim(ncid,trim(d3name3),dim33,d3id3)
   st = nf_def_dim(ncid,trim(d2name1),dim21,d2id1)
   st = nf_def_dim(ncid,trim(d2name2),dim22,d2id2)
   st = nf_def_dim(ncid,trim(d1name1),dim11,d1id1)

   st = nf_def_var(ncid,trim(d4name1),nf_real,1,d4id1,a4id1)
   st = nf_def_var(ncid,trim(d4name2),nf_real,1,d4id2,a4id2)
   st = nf_def_var(ncid,trim(d4name3),nf_real,1,d4id3,a4id3)
   st = nf_def_var(ncid,trim(d4name4),nf_real,1,d4id4,a4id4)
   st = nf_def_var(ncid,trim(d3name1),nf_real,1,d3id1,a3id1)
   st = nf_def_var(ncid,trim(d3name2),nf_real,1,d3id2,a3id2)
   st = nf_def_var(ncid,trim(d3name3),nf_real,1,d3id3,a3id3)
   st = nf_def_var(ncid,trim(d2name1),nf_real,1,d2id1,a2id1)
   st = nf_def_var(ncid,trim(d2name2),nf_real,1,d2id2,a2id2)
   st = nf_def_var(ncid,trim(d1name1),nf_real,1,d1id1,a1id1)


   vardimid4 = (/d4id1, d4id2, d4id3,d4id4/)
   vardimid3 = (/d3id1, d3id2, d3id3/)
   vardimid2 = (/d2id1, d2id2/)
   vardimid1 = (/d1id1/)

   do i=1, nv4
     st = nf_def_var(ncid,trim(name4(i)),nf_real,4,vardimid4,varid4(i)) 
     st = nf_put_att_real(ncid,varid4(i),'_FillValue',nf_real,1,1.e+32) 
   enddo
   do i=1, nv3
     st = nf_def_var(ncid,trim(name3(i)),nf_real,3,vardimid3,varid3(i))
     st = nf_put_att_real(ncid,varid3(i),'_FillValue',nf_real,1,1.e+32)
   enddo
   do i=1, nv2
     st = nf_def_var(ncid,trim(name2(i)),nf_real,2,vardimid2,varid2(i))
     st = nf_put_att_real(ncid,varid2(i),'_FillValue',nf_real,1,1.e+32)
   enddo
   do i=1, nv1
     st = nf_def_var(ncid,trim(name1(i)),nf_real,1,vardimid1,varid1(i))
     st = nf_put_att_real(ncid,varid1(i),'_FillValue',nf_real,1,1.e+32)
   enddo

   st = nf_put_att_text(ncid,nf_global,'title',len_trim(title),title)
   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a4id1,axis41)
   st = nf_put_var_real(ncid,a4id2,axis42)
   st = nf_put_var_real(ncid,a4id3,axis43)
   st = nf_put_var_real(ncid,a4id4,axis44)
   st = nf_put_var_real(ncid,a3id1,axis31)
   st = nf_put_var_real(ncid,a3id2,axis32)
   st = nf_put_var_real(ncid,a3id3,axis33)
   st = nf_put_var_real(ncid,a2id1,axis21)
   st = nf_put_var_real(ncid,a2id2,axis22)
   st = nf_put_var_real(ncid,a1id1,axis11)

   do i=1, nv4
     st = nf_put_var_real(ncid,varid4(i),var4(:,:,:,:,i))
   enddo
   do i=1, nv3
     st = nf_put_var_real(ncid,varid3(i),var3(:,:,:,i))
   enddo
   do i=1, nv2
     st = nf_put_var_real(ncid,varid2(i),var2(:,:,i))
   enddo
   do i=1, nv1
     st = nf_put_var_real(ncid,varid1(i),var1(:,i))
   enddo

   st = nf_close(ncid)
   end subroutine out4d3d2d1d


   subroutine out5d(wname,nv5,name5,var5, &
                    d5name1,dim51,axis51, &
                    d5name2,dim52,axis52, &
                    d5name3,dim53,axis53, &
                    d5name4,dim54,axis54, &
                    d5name5,dim55,axis55, title)

   implicit none

   integer, intent(in) :: nv5
   character(len=*), intent(in) :: wname, name5(nv5)
   character(len=*), intent(in) :: d5name1, d5name2, d5name3, d5name4, d5name5
   character(len=*), intent(in) :: title
   integer, intent(in) :: dim51, dim52, dim53, dim54, dim55
   real, intent(in)    :: axis51(dim51), axis52(dim52), axis53(dim53), axis54(dim54), axis55(dim55)
   real, dimension(dim51,dim52,dim53,dim54,dim55,nv5), intent(in) :: var5

   integer :: i,j,k,it,st
   integer :: ncid
   integer :: d5id1, d5id2, d5id3, d5id4, d5id5
   integer :: a5id1, a5id2, a5id3, a5id4, a5id5
   integer :: varid5(nv5)
   integer :: vardimid5(5)


   st = nf_create(trim(wname),nf_share,ncid)

   st = nf_def_dim(ncid,trim(d5name1),dim51,d5id1)
   st = nf_def_dim(ncid,trim(d5name2),dim52,d5id2)
   st = nf_def_dim(ncid,trim(d5name3),dim53,d5id3)
   st = nf_def_dim(ncid,trim(d5name4),dim54,d5id4)
   st = nf_def_dim(ncid,trim(d5name5),dim55,d5id5)

   st = nf_def_var(ncid,trim(d5name1),nf_real,1,d5id1,a5id1)
   st = nf_def_var(ncid,trim(d5name2),nf_real,1,d5id2,a5id2)
   st = nf_def_var(ncid,trim(d5name3),nf_real,1,d5id3,a5id3)
   st = nf_def_var(ncid,trim(d5name4),nf_real,1,d5id4,a5id4)
   st = nf_def_var(ncid,trim(d5name5),nf_real,1,d5id5,a5id5)


   vardimid5 = (/d5id1, d5id2, d5id3, d5id4, d5id5/)

   do i=1, nv5
     st = nf_def_var(ncid,trim(name5(i)),nf_real,5,vardimid5,varid5(i))
     st = nf_put_att_real(ncid,varid5(i),'_FillValue',nf_real,1,1.e+32)
   enddo
   st = nf_put_att_text(ncid,nf_global,'title',len_trim(title),title)
   st = nf_enddef(ncid)

   st = nf_put_var_real(ncid,a5id1,axis51)
   st = nf_put_var_real(ncid,a5id2,axis52)
   st = nf_put_var_real(ncid,a5id3,axis53)
   st = nf_put_var_real(ncid,a5id4,axis54)
   st = nf_put_var_real(ncid,a5id5,axis55)

   do i=1, nv5
     st = nf_put_var_real(ncid,varid5(i),var5(:,:,:,:,:,i))
   enddo

   st = nf_close(ncid)
   end subroutine out5d


SUBROUTINE outnc(wname,nvar,set,title,nmode)

  implicit none

  integer,                     intent(in) ::  nvar
  type(vset), dimension(nvar), intent(in) ::  set
  character(len=*),            intent(in) ::  wname, title
  integer, optional,           intent(in) ::  nmode

  integer ::  iv, ia, ncid, cnt, aid, ii, st, ndim
  integer ::  n_rec
  integer, dimension(nvar*4) ::  dimid, axid
  character(len=64) ::  ctemp(nvar*4)
  integer, dimension(4,nvar) ::  tag
  integer, dimension(4)      ::  vardimid
  integer, dimension(nvar)   ::  varid

  n_rec = 0
  if ( present(nmode) ) then
    if ( nmode >= 1 .and. nmode <= 4 )  n_rec = nmode
  end if

  st = nf_create(trim(wname),nf_share,ncid)

  cnt = 0
  do iv=1, nvar
    ndim = 0
    do ia=1, 4
      cnt = cnt + 1
      ctemp(cnt) = set(iv)%axis(ia)
      if (trim(ctemp(cnt)) == '') then
        tag(ia,iv) = 0
      else
        ndim = ndim + 1
        tag(ia,iv) = cnt
        do ii=cnt-1, 1, -1
          if (ctemp(cnt) == ctemp(ii))  tag(ia,iv) = ii
        enddo
        if (tag(ia,iv) == cnt) then
          if (ia == n_rec) then
            st = nf_def_dim(ncid,trim(ctemp(cnt)),nf_unlimited,dimid(cnt))
          else
            st = nf_def_dim(ncid,trim(ctemp(cnt)),set(iv)%nd(ia),dimid(cnt))
          end if
          st = nf_def_var(ncid,trim(ctemp(cnt)),nf_real,1,dimid(cnt),axid(cnt))
        end if
        vardimid(ndim) = dimid(tag(ia,iv))
      end if
    enddo
    st = nf_def_var(ncid,trim(set(iv)%vname),nf_real,ndim,vardimid(1:ndim),varid(iv))
    st = nf_put_att_real(ncid,varid(iv),'_FillValue',nf_real,1,1.e+32)

  enddo
  st = nf_put_att_text(ncid,nf_global,'title',len_trim(title),title)

  st = nf_enddef(ncid)

  cnt = 0
  do iv=1, nvar
  do ia=1, 4
    cnt = cnt + 1
    if (tag(ia,iv) == cnt) then
      if (ia == 1)  st = nf_put_vara_real(ncid,axid(cnt),1,set(iv)%nd(ia),set(iv)%axis1)
      if (ia == 2)  st = nf_put_vara_real(ncid,axid(cnt),1,set(iv)%nd(ia),set(iv)%axis2)
      if (ia == 3)  st = nf_put_vara_real(ncid,axid(cnt),1,set(iv)%nd(ia),set(iv)%axis3)
      if (ia == 4)  st = nf_put_vara_real(ncid,axid(cnt),1,set(iv)%nd(ia),set(iv)%axis4)
    end if
  enddo
  enddo

  do iv=1, nvar
    st = nf_put_var_real(ncid,varid(iv),set(iv)%var_out)
  enddo

  st = nf_close(ncid)

END subroutine outnc

SUBROUTINE outnc_blank(wname,nvar,set0,title,nmode)

  implicit none

  integer,                         intent(in) ::  nvar
  type(vset0), dimension(nvar),    intent(in) ::  set0
  character(len=*),                intent(in) ::  wname, title
  integer, optional,               intent(in) ::  nmode

  integer ::  iv, ia, ncid, cnt, aid, ii, st, ndim
  integer ::  n_rec
  integer, dimension(nvar*4) ::  dimid, axid
  character(len=64) ::  ctemp(nvar*4)
  integer, dimension(4,nvar) ::  tag
  integer, dimension(4)      ::  vardimid
  integer, dimension(nvar)   ::  varid

  n_rec = 0
  if ( present(nmode) ) then
    if ( nmode >= 1 .and. nmode <= 4 )  n_rec = nmode
  end if

  st = nf_create(trim(wname),nf_share,ncid)

  cnt = 0
  do iv=1, nvar
    ndim = 0
    do ia=1, 4
      cnt = cnt + 1
      ctemp(cnt) = set0(iv)%axis(ia)
      if (trim(ctemp(cnt)) == '') then
        tag(ia,iv) = 0
      else
        ndim = ndim + 1
        tag(ia,iv) = cnt
        do ii=cnt-1, 1, -1
          if (ctemp(cnt) == ctemp(ii))  tag(ia,iv) = ii
        enddo
        if (tag(ia,iv) == cnt) then
          if (ia == n_rec) then
            st = nf_def_dim(ncid,trim(ctemp(cnt)),nf_unlimited,dimid(cnt))
          else
            st = nf_def_dim(ncid,trim(ctemp(cnt)),set0(iv)%nd(ia),dimid(cnt))
          end if
          st = nf_def_var(ncid,trim(ctemp(cnt)),nf_real,1,dimid(cnt),axid(cnt))
        end if
        vardimid(ndim) = dimid(tag(ia,iv))
      end if
    enddo
    st = nf_def_var(ncid,trim(set0(iv)%vname),nf_real,ndim,vardimid(1:ndim),varid(iv))
    st = nf_put_att_real(ncid,varid(iv),'_FillValue',nf_real,1,1.e+32)

  enddo
  st = nf_put_att_text(ncid,nf_global,'title',len_trim(title),title)

  st = nf_enddef(ncid)

  cnt = 0
  do iv=1, nvar
  do ia=1, 4
    cnt = cnt + 1
    if (tag(ia,iv) == cnt) then
      if (ia == 1)  st = nf_put_vara_real(ncid,axid(cnt),1,set0(iv)%nd(ia),set0(iv)%axis1)
      if (ia == 2)  st = nf_put_vara_real(ncid,axid(cnt),1,set0(iv)%nd(ia),set0(iv)%axis2)
      if (ia == 3)  st = nf_put_vara_real(ncid,axid(cnt),1,set0(iv)%nd(ia),set0(iv)%axis3)
      if (ia == 4)  st = nf_put_vara_real(ncid,axid(cnt),1,set0(iv)%nd(ia),set0(iv)%axis4)
    end if
  enddo
  enddo

  st = nf_close(ncid)

END subroutine outnc_blank

SUBROUTINE outnc_ss(wname,nvar,set,subsetd)

  implicit none

  integer,                         intent(in) ::  nvar
  type(vset), dimension(nvar),     intent(in) ::  set
  character(len=*),                intent(in) ::  wname
  integer, dimension(3), optional, intent(in) ::  subsetd

  integer ::  iv, ia, ncid, st
  integer, dimension(nvar)   ::  varid, ndim

  st = nf_open(trim(wname),nf_write,ncid)

  ndim(:) = 0
  do iv=1, nvar
    do ia=1, 4
      if (trim(set(iv)%axis(ia)) /= '') then
        ndim(iv) = ndim(iv) + 1
      end if
    enddo
  enddo
  if ( sum(abs(ndim(:)-4)) /= 0 ) then
    print*, 'outnc_ss: not yet coded for ndim /= 4.'
    STOP
  end if
 
  do iv=1, nvar
    st = nf_inq_varid(ncid,trim(set(iv)%vname),varid(iv))
  enddo

  if ( present(subsetd) ) then

    select case ( subsetd(1) )
    case ( 1 )
      do iv=1, nvar
        st = nf_put_vara_real(ncid,varid(iv),  &
             (/subsetd(2),1,1,1/),  &
             (/subsetd(3),set(iv)%nd(2),set(iv)%nd(3),set(iv)%nd(4)/),  &
             set(iv)%var_out)
      enddo
    case ( 2 )
      do iv=1, nvar
        st = nf_put_vara_real(ncid,varid(iv),  &
             (/1,subsetd(2),1,1/),  &
             (/set(iv)%nd(1),subsetd(3),set(iv)%nd(3),set(iv)%nd(4)/),  &
             set(iv)%var_out)
      enddo
    case ( 3 )
      do iv=1, nvar
        st = nf_put_vara_real(ncid,varid(iv),  &
             (/1,1,subsetd(2),1/),  &
             (/set(iv)%nd(1),set(iv)%nd(2),subsetd(3),set(iv)%nd(4)/),  &
             set(iv)%var_out)
      enddo
    case ( 4 )
      do iv=1, nvar
        st = nf_put_vara_real(ncid,varid(iv),  &
             (/1,1,1,subsetd(2)/),  &
             (/set(iv)%nd(1),set(iv)%nd(2),set(iv)%nd(3),subsetd(3)/),  &
             set(iv)%var_out)
      enddo
    case default
      print*, 'Error: subsetd(1) should be [1,4].'
      STOP
    end select

  else

    do iv=1, nvar
      st = nf_put_var_real(ncid,varid(iv),set(iv)%var_out)
    enddo

  end if

  st = nf_close(ncid)

END subroutine outnc_ss


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

   subroutine errhdr_nonstop(ier)

   implicit none

   include 'netcdf.inc'

   integer, intent(in) :: ier

   if (ier /= NF_NOERR) then
     print *,NF_STRERROR(ier)
   end if

   return
   end subroutine errhdr_nonstop


   subroutine varndim(ncid,varname,nd)

   implicit none

   integer,          intent(in)  :: ncid
   character(len=*), intent(in)  :: varname
   integer,          intent(out) :: nd

   integer :: varid, st

   st = nf_inq_varid(ncid,trim(varname),varid)
   call errhdr(st)
   st = nf_inq_varndims(ncid,varid,nd)
   call errhdr(st)

   return
   end subroutine varndim

   subroutine varlen(ncid,varname,nd,dlen)

   implicit none

   integer,                intent(in)  :: ncid, nd
   character(len=*),       intent(in)  :: varname
   integer, dimension(nd), intent(out) :: dlen

   integer :: varid, st, dimid(nd), i

   st = nf_inq_varid(ncid,trim(varname),varid)
   call errhdr(st)
   st = nf_inq_vardimid(ncid,varid,dimid)
   call errhdr(st)
   do i=1, nd
     st = nf_inq_dimlen(ncid,dimid(i),dlen(i))
     call errhdr(st)
   enddo

   return
   end subroutine varlen

   end module netio

