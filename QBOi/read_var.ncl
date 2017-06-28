
;-- head of getvarfile0_q ----------------------------------------------
undef("getvarfile0_q")
function getvarfile0_q(model:string, vhead:string, yr, mon, hour)
local files, model, vhead, yr, mon, hour, vname, mm, hh, vfname, sysout,  \
      i0, nfile
;-- body of getvarfile0_q ----------------------------------------------
begin

  vname = str_lower(vhead)+"a"
  if (model .eq. "ERA-I-nr") then
    vname = str_upper(vhead)+"_GDS0_HYBL"
  else if (model .eq. "ERA-I") then
    vname = str_upper(vhead)+"_GDS0_ISBL"
  end if
  end if

  if ( isstring(mon) ) then
    mm = mon
  else
    mm = sprinti("%2.2i",mon)
  end if
  if ( isstring(hour) ) then
    hh = hour
  else
    hh = sprinti("%2.2i",hour)
  end if

  if ( any(str_lower(vhead) .eq. (/"u","v","w","t"/)) ) then

    if (str_get_cols(model,0,4) .eq. "ERA-I") then
      vfname = str_lower(vhead)
      if (vhead .eq. "w") then
        vfname = "ome"
      end if
      if (model .eq. "ERA-I-nr") then
        files = (/"$HOME/portal/ERA-I-nr/6hourly/0.75deg/"+yr+"/"+mm+"/"+  \
                  "era-int_f."+vfname+".anal."+hh+".ml."+yr+""+mm+".nc"/)
      else
        files = (/"$HOME/portal/ERA-I/6hourly/1.5deg/"+yr+"/"+mm+"/"+  \
                  "era-int."+vfname+".anal."+hh+".pl."+yr+""+mm+".nc"/)
      end if
      files@nfile = 1
    else
      sysout = systemfunc("env_qboi.sh 5 "+vname+" 6hr XYPT r1i1p1 "+yr+mm)

      i0 = ind( sysout .eq. "start_"+model ) + 2
      nfile = stringtoint(sysout(i0-1))
      files = sysout(i0:i0+nfile-1)
      files@nfile = nfile
    end if

    files@vname = vname

  else

    print("Invalid variable head in getfile_q")

  end if

  return(files)

end

;-- head of getvarfile_q -----------------------------------------------
undef("getvarfile_q")
function getvarfile_q(model:string, vhead:string, yr, mon, hour)
local files, model, vhead, yr, mon, hour, file0, file1
;-- body of getvarfile_q -----------------------------------------------
begin

  if ( any(str_lower(vhead) .eq. (/"u","v","w","t"/)) ) then

    files = getvarfile0_q(model,vhead,yr,mon,hour)

  else

    if (vhead .eq. "D") then
      file0 = getvarfile0_q(model,"u",yr,mon,hour)
      file1 = getvarfile0_q(model,"v",yr,mon,hour)
    else if (vhead .eq. "uw") then
      file0 = getvarfile0_q(model,"u",yr,mon,hour)
      file1 = getvarfile0_q(model,"w",yr,mon,hour)
    end if
    end if
    files = new((/2,file0@nfile/),string)
    files(0,:) = (/ file0 /)
    files(1,:) = (/ file1 /)
    files@nfile = file0@nfile
    files@vname = (/file0@vname,file1@vname/)

  end if

  return(files)

end

