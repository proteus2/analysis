
param_ra_ranames  = (/"ERA-I"  ,"MERRA","MERRA2","CFSR","JRA-55","JRA-55C",  \
                      "ERA-I_m"  ,"MERRA2_m" ,"JRA-55_m","JRA-55C_m"/)
param_ra_expnames = (/"era-int","merra","merra2","cfsr","jra55" ,"jra55c" ,  \
                      "era-int_f","merra2_f" ,"jra55_f" ,"jra55c_f" /)
param_ra_dirnames = (/"ERA-I"  ,"MERRA","MERRA2","CFSR","JRA55" ,"JRA55C" ,  \
                      "ERA-I-nr" ,"MERRA2-nr","JRA55-nr","JRA55C-nr"/)

param_ra_plevs = new((/dimsizes(param_ra_expnames),30/),float)
param_ra_plevs@_FillValue = -999
param_ra_plevs = param_ra_plevs@_FillValue

param_ra_vnames  = new((/dimsizes(param_ra_expnames),4/),string)
param_ra_fvnames = new((/dimsizes(param_ra_expnames),4/),string)
param_ra_vars = (/"T","u","v","w"/)
; ERA-I
param_ra_vnames (0,:) = (/"T","U","V","W"  /)+"_GDS0_ISBL"
param_ra_fvnames(0,:) = (/"t","u","v","ome"/)
param_ra_plevs(0,:14) = (/200,175,150,125,100,70,50,30,20,10,7,5,3,2,1/)
; MERRA
param_ra_vnames (1,:) = (/"t","u","v","omega"/)
param_ra_fvnames(1,:) = (/"t","u","v","omega"/)
param_ra_plevs(1,:14) = (/200,150,100,70,50,40,30,20,10,7,5,4,3,2,1/)
; MERRA2
param_ra_vnames (2,:) = (/"T","U","V","OMEGA"/)
param_ra_fvnames(2,:) = (/"t","u","v","omega"/)
param_ra_plevs(2,:14) = (/200,150,100,70,50,40,30,20,10,7,5,4,3,2,1/)
; CFSR
param_ra_vnames (3,:) = (/"TMP","UGRD","VGRD","VVEL" /)
param_ra_fvnames(3,:) = (/"t"  ,"u"   ,"v"   ,"omega"/)
param_ra_plevs(3,:14) = (/200,175,150,125,100,70,50,30,20,10,7,5,3,2,1/)
; JRA55
param_ra_vnames (4,:) = (/"TMP","U_GRD","V_GRD","V_VEL"/)+"_GDS0_ISBL"
param_ra_fvnames(4,:) = (/"tmp","ugrd" ,"vgrd" ,"vvel" /)
param_ra_plevs(4,:14) = (/200,175,150,125,100,70,50,30,20,10,7,5,3,2,1/)
; JRA55C
param_ra_vnames (5,:) = param_ra_vnames (4,:)
param_ra_fvnames(5,:) = param_ra_fvnames(4,:)
param_ra_plevs  (5,:) = param_ra_plevs  (4,:)
; ERA-I-nr
param_ra_vnames (6,:) = (/"t","u","v","w"  /)
param_ra_fvnames(6,:) = (/"t","u","v","ome"/)
param_ra_plevs(6,:25) = (/176.6,153.6,132.5,113.1,95.74,80.17,66.42,  \
     54.43,44.17,35.65,28.77,23.22,18.74,15.13,12.21,9.855,7.954,6.420,5.182,  \
     4.177,3.348,2.654,2.065,1.564,1.140,0.7856/)
; MERRA2-nr
param_ra_vnames (7,:) = (/"T" ,"U" ,"V" ,"OMEGA" /)
param_ra_fvnames(7,:) = (/"t0","u0","v0","omega0"/)
param_ra_plevs(7,:28) = (/192.5,163.7,139.1,118.3,100.5,85.44,72.56,  \
     61.50,52.02,43.91,36.99,31.09,26.05,21.76,18.12,15.05,12.46,10.28,8.456,  \
     6.918,5.632,4.562,3.677,2.948,2.353,1.868,1.476,1.160,0.9072/)
; JRA55-nr
param_ra_vnames (8,:) = (/"TMP","U_GRD","V_GRD","V_VEL"/)+"_GDS4_HYBL"
param_ra_fvnames(8,:) = (/"tmp","ugrd" ,"vgrd" ,"vvel" /)
param_ra_plevs(8,:25) = (/176.9,153.9,132.9,113.9,96.89,81.64,67.64,  \
     55.15,44.67,36.08,29.15,23.53,18.99,15.32,12.35,9.971,8.050,6.493,5.240,  \
     4.222,3.383,2.684,2.089,1.584,1.160,0.8047/)
; JRA55C-nr
param_ra_vnames (9,:) = param_ra_vnames (8,:)
param_ra_fvnames(9,:) = param_ra_fvnames(8,:)
param_ra_plevs  (9,:) = param_ra_plevs  (8,:)

