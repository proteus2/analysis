
  MODULE INTERPOLATE
  
      Implicit None
  
      contains 
          
      subroutine interpol(A1,A2,A3,A4,A5,A6,A7,A8,dx,dy,dz,delx,dely,delz,Af)
  
!------------------------------------------------------------------------------
!  Obtain Interpolation of Variable
!  
!  Author  
!  ILLIAN 
!  
!  Date 
!  JULY 16, 2004
!  
!  INPUT AND OUTPUT DESCRIPTION
!  
!  dx   :  Data Interval of X-Direction
!  dy   :  Data Interval of Y-Direction
!  dz   :  Data Interval of Z-Direction
!  delx :  Interpolation Interval of X-Direction from A1
!  dely :  Interpolation Interval of Y-Direction from A1
!  delz :  Interpolation Interval of Z-Direction from A1
!  Ai   :  To Interpolate Variable A of Index i
!  Af   :  Interpolation for X-Y-Z with delx,dely,delz      (double)
!
!------------------------------------------------------------------------------



      Real,   intent(in)                     ::  A1,A2,A3,A4,A5,A6,A7,A8
      Real*8, intent(in)                     ::  dx,dy,dz,delx,dely,delz
      Real*8, intent(out)                    ::  Af
      Real*8                                 ::  Aa1,Aa2,Aa3,Aa4,Ab1,Ab2

      Aa1 = ( dble(A2*delx + A1*(dx-delx)) )  / dx
      Aa2 = ( dble(A4*delx + A3*(dx-delx)) )  / dx
      Aa3 = ( dble(A6*delx + A5*(dx-delx)) )  / dx
      Aa4 = ( dble(A8*delx + A7*(dx-delx)) )  / dx
      
      Ab1 = ( Aa2*dely + Aa1*(dy-dely) ) / dy
      Ab2 = ( Aa4*dely + Aa3*(dy-dely) ) / dy

      Af = ( Ab2*delz + Ab1*(dz-delz) ) / dz

      end subroutine interpol

  END MODULE INTERPOLATE
