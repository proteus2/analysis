MODULE switch_dump

  implicit none

  logical, public ::  l_drag_u_o, l_drag_v_o
  logical, public ::  l_mflx_u_o, l_mflx_v_o
  logical, public ::  l_mflx_u_ctop_o, l_mflx_v_ctop_o
  logical, public ::  l_spec_o, l_spec_ctop_o
  logical, public ::  l_mflx_u_ctop_0_o, l_mflx_v_ctop_0_o
  logical, public ::  l_spec_ctop_0_o
  logical, public ::  l_diag_znwcq_o

CONTAINS


SUBROUTINE switch_defaults

  l_drag_u_o = .True.  ;  l_drag_v_o = .True.
  l_mflx_u_o = .True.  ;  l_mflx_v_o = .True.
  l_mflx_u_ctop_o = .true.  ;  l_mflx_v_ctop_o = .true.
  l_spec_o = .true.  ;  l_spec_ctop_o = .true.
  l_mflx_u_ctop_0_o = .true.  ;  l_mflx_v_ctop_0_o = .true.
  l_spec_ctop_0_o = .true.
  l_diag_znwcq_o = .true.

END subroutine switch_defaults

SUBROUTINE get_nv_output(nv)

  integer, intent(out) ::  nv

  nv = 0
  if ( l_mflx_u_ctop_o )  nv = nv + 2
  if ( l_mflx_v_ctop_o )  nv = nv + 2
  if ( l_mflx_u_o      )  nv = nv + 2
  if ( l_mflx_v_o      )  nv = nv + 2
  if ( l_drag_u_o      )  nv = nv + 1
  if ( l_drag_v_o      )  nv = nv + 1
  if ( l_spec_ctop_o   )  nv = nv + 1
  if ( l_spec_o        )  nv = nv + 1

END subroutine get_nv_output

END module switch_dump

