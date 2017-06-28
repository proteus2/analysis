
;-- head of get_file_indx ----------------------------------------------
undef("get_file_indx")
function get_file_indx(exp:string, dir:string, fileid:string, season:string,  \
                       year:integer)
local files, exp, dir, fileid, season, year,  \
      stmp, d_target, ftmp0, ftmp1, ftmp2, f0, ytime, yt_center, indx
;-- body of get_file_indx ----------------------------------------------
begin

 stmp = dir+"/"+exp+".cam."+fileid+"."
 d_target = new(2,double)
 if (str_upper(season) .eq. "DJF") then
   ftmp0 = systemfunc("ls "+stmp+sprinti("%4.4i",year-1)+"-11*.nc")
   ftmp1 = systemfunc("ls "+stmp+sprinti("%4.4i",year-1)+"-12*.nc "+  \
                            stmp+sprinti("%4.4i",year)+"-0[12]*.nc")
   ftmp2 = systemfunc("ls "+stmp+sprinti("%4.4i",year)+"-03*.nc")
   d_target(0) = (year-1)*10000 + 12*100 + 1
   d_target(1) = year*10000 + 3*100 + 1
 else if (str_upper(season) .eq. "MAM") then
   ftmp0 = systemfunc("ls "+stmp+sprinti("%4.4i",year)+"-02*.nc")
   ftmp1 = systemfunc("ls "+stmp+sprinti("%4.4i",year)+"-0[345]*.nc")
   ftmp2 = systemfunc("ls "+stmp+sprinti("%4.4i",year)+"-06*.nc")
   d_target(0) = year*10000 + 3*100 + 1
   d_target(1) = year*10000 + 6*100 + 1
 else if (str_upper(season) .eq. "JJA") then
   ftmp0 = systemfunc("ls "+stmp+sprinti("%4.4i",year)+"-05*.nc")
   ftmp1 = systemfunc("ls "+stmp+sprinti("%4.4i",year)+"-0[678]*.nc")
   ftmp2 = systemfunc("ls "+stmp+sprinti("%4.4i",year)+"-09*.nc")
   d_target(0) = year*10000 + 6*100 + 1
   d_target(1) = year*10000 + 9*100 + 1
 else if (str_upper(season) .eq. "SON") then
   ftmp0 = systemfunc("ls "+stmp+sprinti("%4.4i",year)+"-08*.nc")
   ftmp1 = systemfunc("ls "+stmp+sprinti("%4.4i",year)+"-09*.nc "+  \
                            stmp+sprinti("%4.4i",year)+"-1[01]*.nc")
   ftmp2 = systemfunc("ls "+stmp+sprinti("%4.4i",year)+"-12*.nc")
   d_target(0) = year*10000 + 9*100 + 1
   d_target(1) = year*10000 + 12*100 + 1
 else if (str_lower(season) .eq. "all") then  ; from DEC to NOV
   ftmp0 = systemfunc("ls "+stmp+sprinti("%4.4i",year-1)+"-11*.nc")
   ftmp1 = systemfunc("ls "+stmp+sprinti("%4.4i",year-1)+"-12*.nc "+  \
                            stmp+sprinti("%4.4i",year)+"-0*.nc "+  \
                            stmp+sprinti("%4.4i",year)+"-1[01]*.nc")
   ftmp2 = systemfunc("ls "+stmp+sprinti("%4.4i",year)+"-12*.nc")
   d_target(0) = (year-1)*10000 + 12*100 + 1
   d_target(1) = year*10000 + 12*100 + 1
 else if (str_lower(season) .eq. "ann") then  ; from JAN to DEC
   ftmp0 = systemfunc("ls "+stmp+sprinti("%4.4i",year-1)+"-12*.nc")
   ftmp1 = systemfunc("ls "+stmp+sprinti("%4.4i",year)+"-*.nc")
   ftmp2 = systemfunc("ls "+stmp+sprinti("%4.4i",year+1)+"-01*.nc")
   d_target(0) = year*10000 + 1*100 + 1
   d_target(1) = (year+1)*10000 + 1*100 + 1
 end if
 end if
 end if
 end if
 end if
 end if
 files = array_append_record(array_append_record(  \
         ftmp0(dimsizes(ftmp0)-1),ftmp1,0),ftmp2(0),0)

 files@nfiles = dimsizes(files)

 f0 = addfiles(files,"r")
 ytime = yyyymmdd_to_yyyyfrac( todouble(f0[:]->date), 0 )
 yt_center = ytime - 0.5*(ytime(1) - ytime(0))
 indx = ind( yt_center .ge. yyyymmdd_to_yyyyfrac(d_target(0), 0) .and.  \
             yt_center .lt. yyyymmdd_to_yyyyfrac(d_target(1), 0) )

 files@t_indx_2 = (/min(indx),max(indx)/)
 files@t_2 = f0[:]->time(files@t_indx_2)
 files@date_2 = f0[:]->date(files@t_indx_2)
 files@dt = f0[:]->time(1) - f0[:]->time(0)
 
 return(files)

end

