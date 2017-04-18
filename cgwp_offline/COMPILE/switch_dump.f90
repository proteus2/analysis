MODULE switch_dump

  implicit none

  logical, public ::  l_mflx_regrid
  logical, public ::  l_drag_u_o, l_drag_v_o
  logical, public ::  l_mflx_u_o, l_mflx_v_o
  logical, public ::  l_mflx_u_ctop_o, l_mflx_v_ctop_o
  logical, public ::  l_spec_o, l_spec_ctop_o

CONTAINS


SUBROUTINE switch_defaults

  l_mflx_regrid = .True.
  l_drag_u_o = .True.  ;  l_drag_v_o = .True.
  l_mflx_u_o = .True.  ;  l_mflx_v_o = .True.
  l_mflx_u_ctop_o = .True.  ;  l_mflx_v_ctop_o = .True.
  l_spec_o = .True.  ;  l_spec_ctop_o = .True.

END subroutine switch_defaults

END module switch_dump

