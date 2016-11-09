
      Subroutine T_vert_interp_to_p(                                    &
     &                           T, theta, row_length, rows,            &
     &                           model_levels, desired_p,               &
     &                           off_x, off_y,                          &
     &                            halo_i, halo_j,                       &
     &                           p_theta_levels,                        &
     &                           lapse_rate, R, g,                      &
     &                           boundary_layer_levels,                 &
     &                           exner_theta_levels,                    &
     &                           r_theta_levels,                        &
     &                           kappa, p_zero, T_out )

! Purpose:
!          Performs vertical interpolation of temperature to a
!          desired p surface assuming that T varies linearly with
!          geopotential height, which itself is assumed to vary
!          linearly with Exner pressure. Hence T varies linearly with
!          exner pressure. Where the desired surface is above/below
!          the top/bottom data point extrapolation is done
!          respectively: isothermal or
!          from the first level above a fixed height (200m)
!
      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
                        ! number of points on a row
     &, rows                                                            &
                        ! number of rows of data
     &, model_levels                                                    &
                        ! number of levels of data
     &, boundary_layer_levels                                           &
     &, off_x, off_y                                                    &
                        ! halo sizes
     &, halo_i, halo_j  ! large halo sizes

      Real                                                              &
     &  desired_p                                                       &
     &, p_zero                                                          &
     &, kappa                                                           &
     &, lapse_rate                                                      &
     &, R                                                               &
     &, g                                                               &
     &,earth_radius

      Real                                                              &
     &  T (row_length, rows, model_levels)                              &
     &, theta (1-off_x:row_length+off_x,                                &
     &         1-off_y:rows+off_y, model_levels)                        &
     &, r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                 1-halo_j:rows+halo_j,0:model_levels)             &
     &, exner_theta_levels (1-off_x:row_length+off_x,                   &
     &                      1-off_y:rows+off_y, model_levels)           &
     &, p_theta_levels(1-off_x:row_length+off_x,                        &
     &                 1-off_y:rows+off_y, model_levels)

! Arguments with Intent OUT. ie: Output variables.
      Real                                                              &
     &  T_out (row_length, rows)

! Local variables

      Integer                                                           &
     & i,j,k                                                            &
     & , level_extrap (row_length, rows)

      Integer                                                           &
     &  level_below(row_length, rows)

      Real                                                              &
     &  desired_exner                                                   &
     &, power                                                           &
     & , extrap_height  !height to determine level_extrap

! ----------------------------------------------------------------------
! Section 0. Initialize arrays
! ----------------------------------------------------------------------
      Do j = 1, rows
        Do i = 1, row_length
          level_below(i,j) = 0
          level_extrap(i,j) = 1
        End Do
      End Do

! ----------------------------------------------------------------------
! Section 1. Find level below which desired surface is.
!         Use nearest level below extrap_height if extrapolation needed
!         200m is chosed as realistic physical value for atmosphere
!         If first theta level is above this first level is used
! ----------------------------------------------------------------------
      extrap_height = 200 !m
      Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            If ( (r_theta_levels(i,j,k)-r_theta_levels(i,j,0) )  <=     &
     &            extrap_height ) Then
              level_extrap(i,j) = k
            End If
          End Do
        End Do
      End Do


! change desired pressure into exner equivalent.

      desired_exner = (desired_p/p_zero) ** kappa


      Do k = 1, model_levels - 1
        Do j = 1, rows
          Do i = 1, row_length
            If ( exner_theta_levels(i,j,k)  >=  desired_exner ) Then
              level_below(i,j) = k
            End If
          End Do
        End Do
      End Do

! if requested level is above top of model, set level_below to -1,
! which will be converted to use isothermal extrapolation

      Do j = 1, rows
        Do i = 1, row_length
          If ( desired_exner  <                                         &
     &         exner_theta_levels(i,j,model_levels) ) Then
            level_below(i,j) = -1
          End If
        End Do
      End Do

! ----------------------------------------------------------------------
! Section 2. Vertical interpolation.
! ----------------------------------------------------------------------

      power = Lapse_Rate * R / g
      Do j = 1, rows
        Do i = 1, row_length

          If (level_below(i,j)  ==  -1) Then
! isothermal extrapolation above top level
            T_out(i,j) = T(i,j,model_levels)

          Else If (level_below(i,j)  ==  0) Then
! extrapolate
            T_out(i,j) = T(i,j,level_extrap(i,j)) *                     &
     &                   (desired_p/                                    &
     &                 p_theta_levels(i,j,level_extrap(i,j)))           &
     &                  ** power

          End If

        End Do
      End Do

! end of routine

      Return
      END SUBROUTINE T_vert_interp_to_p

#endif
