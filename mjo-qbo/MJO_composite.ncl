;==============================================================================
; return indices that meet the following condition for MJO composite
;------------------------------------------------------------------------------
; define MJO event as in L'Heureux and Higgins 2008
; i.    the index > 1 for consecutive pentads -> at least 10 days
; ii.   the phase in numerical order (i.e., 6, 7, 8, 1, 2)
; iii.  events must continue at least 30 days
; iv.   events cannot remain in a phase for more than 20 days
;------------------------------------------------------------------------------

undef("LH_mjo_composite")
function LH_mjo_composite (index, phase)
begin             
;------------------------------------------------------------------------------
time	= index&time
ntim	= dimsizes (time)
ymd	= ut_calendar (time, -2)
;------------------------------------------------------------------------------
; i.    the index > 1 for consecutive pentads -> at least 10 days
;------------------------------------------------------------------------------
iCon1	= ymd
iCon1	= 0
do nt=0, ntim-10
  if (all(index(nt:nt+9).gt.1)) then
    iCon1(nt:nt+9) = 1
  end if
end do
;------------------------------------------------------------------------------
; ii.   the phase in numerical order (i.e., 6, 7, 8, 1, 2)
;------------------------------------------------------------------------------
iCon2	= ymd
iCon2	= 1
do nt=1, ntim-1
  if (ismissing(phase(nt-1)).or.ismissing(phase(nt))) then
    iCon2(nt) = 0
  else
  if (phase(nt-1).eq.8) then
    if (phase(nt).eq.8.or.phase(nt).eq.1) then
    else
      iCon2(nt) = 0
    end if
  else
    if (phase(nt).eq.phase(nt-1).or.phase(nt).eq.phase(nt-1)+1) then
    else
      iCon2(nt) = 0
    end if
  end if
  end if
end do
;------------------------------------------------------------------------------
; iii.  events must continue at least 30 days
;------------------------------------------------------------------------------
iCon3	= ymd
iCon3	= 0
do nt=0, ntim-30
  if (all(iCon2(nt:nt+29).eq.1)) then
    iCon3(nt:nt+29) = 1
  end if
end do
;------------------------------------------------------------------------------
; iv.   events cannot remain in a phase for more than 20 days
;------------------------------------------------------------------------------
iCon4	= ymd
iCon4	= 1
do nt=0, ntim-20
  if (all(phase(nt:nt+19).eq.phase(nt))) then
    iCon4(nt:nt+19) = 0
  end if
end do
;------------------------------------------------------------------------------
; find cases that satisfy all conditions above
;------------------------------------------------------------------------------
iCond	= ymd
iCond	= 0
do nt=0, ntim-30
  if (all(iCon3(nt:nt+29).eq.1.and.iCon4(nt:nt+29).eq.1)) then
    do ntt=0, 20
      if (all(iCon1(nt+ntt:nt+ntt+9).eq.1)) then
        iCond(nt:nt+29) = 1
      end if
    end do
  end if
end do
;------------------------------------------------------------------------------
return(iCond)
end 

;==============================================================================
