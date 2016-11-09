undef("get_col_srip")
function get_col_srip(wks)
local rgb_m, rgb_m2, rgb_ei, rgb_j5, rgb_c
begin

 rgb_m  = (/226,31 ,38 /) / 256.
 rgb_m2 = (/246,153,153/) / 256.
 rgb_ei = (/32 ,120,180/) / 256.
 rgb_j5 = (/107,62 ,152/) / 256.
 rgb_c  = (/52 ,160, 72/) / 256.

 icol = 5
 icol@m  = NhlNewColor(wks,rgb_m (0),rgb_m (1),rgb_m (2))
 icol@m2 = NhlNewColor(wks,rgb_m2(0),rgb_m2(1),rgb_m2(2))
 icol@ei = NhlNewColor(wks,rgb_ei(0),rgb_ei(1),rgb_ei(2))
 icol@j5 = NhlNewColor(wks,rgb_j5(0),rgb_j5(1),rgb_j5(2))
 icol@c  = NhlNewColor(wks,rgb_c (0),rgb_c (1),rgb_c (2))

return(icol)
end
