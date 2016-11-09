
;-- head of axis_str ---------------------------------------------------

undef("axis_str")
procedure axis_str(obj:graphic, mstr:string, xstr:string, ystr:string)
local obj, mstr, xstr, ystr

;-- body of axis_str ---------------------------------------------------

begin

  setvalues obj
   "tiMainString"  : mstr
   "tiXAxisString" : xstr
   "tiYAxisString" : ystr
  end setvalues

end

;-- head of axis_range -------------------------------------------------

undef("axis_range")
procedure axis_range(obj:graphic, xl:numeric, xr:numeric, yb:numeric, \
                     yt:numeric, xlog:logical, ylog:logical)
local obj, xl, xr, yb, yt, xlog, ylog, \
      classname, autoreset_mm, autoreset_log

;-- body of axis_range -------------------------------------------------

begin

  obj@xleft   = xl
  obj@xright  = xr
  obj@ybottom = yb
  obj@ytop    = yt
  obj@xlog    = xlog
  obj@ylog    = ylog

  classname = NhlClassName(obj)

  autoreset_mm  = 0
  if (classname .eq. "contourPlotClass") then
    autoreset_mm  = 1
  end if
  autoreset_log = 0
  if ( classname .eq. "contourPlotClass" .or. \
       classname .eq. "logLinPlotClass" ) then
    autoreset_log = 1
  end if

  if ( autoreset_mm .eq. 0 .or. \
       ( min((/xl,xr/)) .ge. obj@xarrminmax(0) .and. \
         max((/xl,xr/)) .le. obj@xarrminmax(1) ) ) then

    if ( abs(xl) .ne. 999 .and. abs(xr) .ne. 999 ) then
      if (xl .lt. xr) then
        setvalues obj
         "trXMinF" : xl
         "trXMaxF" : xr
        end setvalues
      else
        setvalues obj
         "trXMinF"    : xr
         "trXMaxF"    : xl
         "trXReverse" : True
        end setvalues
      end if
    else
      if (abs(xl) .ne. 999) then
        if (xr .eq. 999) then
          setvalues obj
           "trXMinF" : xl
          end setvalues
        else
          setvalues obj
           "trXMaxF"    : xl
           "trXReverse" : True
          end setvalues
        end if
      end if
      if (abs(xr) .ne. 999) then  
        if (xl .eq. -999) then
          setvalues obj
           "trXMaxF" : xr
          end setvalues
        else
          setvalues obj
           "trXMinF"    : xr
           "trXReverse" : True
          end setvalues
        end if
      end if
      if ( xl .eq. 999 .and. xr .eq. -999 ) then
        setvalues obj
         "trXReverse" : True
        end setvalues
      end if
    end if

  end if

  if ( autoreset_mm .eq. 0 .or. \
       ( min((/yb,yt/)) .ge. obj@yarrminmax(0) .and. \
         max((/yb,yt/)) .le. obj@yarrminmax(1) ) ) then

    if ( abs(yb) .ne. 999 .and. abs(yt) .ne. 999 ) then
      if (yb .lt. yt) then
        setvalues obj
         "trYMinF" : yb
         "trYMaxF" : yt
        end setvalues
      else
        setvalues obj
         "trYMinF"    : yt
         "trYMaxF"    : yb
         "trYReverse" : True
        end setvalues
      end if
    else
      if (abs(yb) .ne. 999) then
        if (yt .eq. 999) then
          setvalues obj
           "trYMinF" : yb
          end setvalues
        else
          setvalues obj
           "trYMaxF"    : yb
           "trYReverse" : True
          end setvalues
        end if
      end if
      if (abs(yt) .ne. 999) then  
        if (yb .eq. -999) then
          setvalues obj
           "trYMaxF" : yt
          end setvalues
        else
          setvalues obj
           "trYMinF"    : yt
           "trYReverse" : True
          end setvalues
        end if
      end if
      if ( yb .eq. 999 .and. yt .eq. -999 ) then
        setvalues obj
         "trYReverse" : True
        end setvalues
      end if
    end if

  end if

  if (autoreset_log .eq. 0) then
    if ( xlog ) then
      setvalues obj
       "trXLog" : True
       "tmXBMinorPerMajor" : 8
       "tmXTMinorPerMajor" : 8
      end setvalues
    end if
    if ( ylog ) then
      setvalues obj
       "trYLog" : True
       "tmYLMinorPerMajor" : 8
       "tmYRMinorPerMajor" : 8
      end setvalues
    end if
  end if

end

;-- head of tick_exp ---------------------------------------------------

undef("tick_exp")
procedure tick_exp(obj:graphic, pos:string, vmaj:numeric, smaj:string, \
                   vmin:numeric)
local obj, pos, vmaj, smaj, vmin, \
      resmod, resvmaj, ressmaj, resvmin, tag_min

;-- body of tick_exp ---------------------------------------------------

begin

  resmod  = "tm"+pos+"Mode"
  resvmaj = "tm"+pos+"Values"
  ressmaj = "tm"+pos+"Labels"
  resvmin = "tm"+pos+"MinorValues"

  setvalues obj
   resmod  : "Explicit"
   resvmaj : vmaj
   ressmaj : smaj
  end setvalues

  tag_min = True
  if ( dimsizes(dimsizes(vmin)) .eq. 1 .and. dimsizes(vmin) .eq. 1 ) then
    if (abs(vmin) .eq. 999) then
      tag_min = False
    end if
  end if
  if ( tag_min ) then
    setvalues obj
     resvmin : vmin
    end setvalues
  end if

end

;-- head of axisvalue_form ---------------------------------------------

undef("axisvalue_form")
procedure axisvalue_form(res:logical)
local res

;-- body of axisvalue_form ---------------------------------------------

begin

  res@tmXBFormat = "@*+^sg"
  res@tmXTFormat = "@*+^sg"
  res@tmYLFormat = "@*+^sg"
  res@tmYRFormat = "@*+^sg"

end

;-- head of xylineplot -------------------------------------------------

undef("xylineplot")
function xylineplot(wks:graphic, x_array:numeric, y_array:numeric, \
                    s1:string, v1:integer, s2:string, v2:numeric,  \
                    s3:string, v3:integer)
local wks, x_array, y_array, s1, v1, s2, v2, s3, v3, \
      nx2, ny2, nl, xx, yy, vv1, vv2, vv3, res, obj

;-- body of xylineplot -------------------------------------------------

begin

  if ( s1 .ne. "D" .or. s2 .ne. "T" .or. s3 .ne. "C" ) then
    print("Set the order of variables in 'xylineplot' :  D-T-C")
    exit
  end if

  nx2 = dimsizes(x_array)
  ny2 = dimsizes(y_array)

  if ( dimsizes(nx2) .gt. 2 .or. dimsizes(ny2) .gt. 2 ) then
    print("Dimensions of input in 'xylineplot' should be less than 3.")
    exit
  end if

  if (dimsizes(nx2) .eq. 2) then
    nl = nx2(0)
  else
    nl = ny2(0)
  end if

  if ( dimsizes(nx2) .eq. 1 .and. nx2 .eq. 1 ) then
    xx = new(ny2,float)
    xx(:) = x_array
  else
    xx = x_array
  end if
  if ( dimsizes(ny2) .eq. 1 .and. ny2 .eq. 1 ) then
    yy = new(nx2,float)
    yy(:) = y_array
  else
    yy = y_array
  end if

  vv1 = new(max((/nl,dimsizes(v1)/)),integer)
  vv2 = new(max((/nl,dimsizes(v2)/)),float  )
  vv3 = new(max((/nl,dimsizes(v3)/)),integer)

  vv1(:dimsizes(v1)-1) = v1(:)
  vv2(:dimsizes(v2)-1) = v2(:)
  vv3(:dimsizes(v3)-1) = v3(:)

  if (dimsizes(v1) .eq. 1) then
    vv1(:) = v1
  end if
  if (dimsizes(v2) .eq. 1) then
    vv2(:) = v2
  end if
  if (dimsizes(v3) .eq. 1) then
    vv3(:) = v3
  end if

  res          = True
  res@gsnDraw  = False
  res@gsnFrame = False
  axisvalue_form(res)

  res@xyDashPatterns    = vv1
  res@xyLineThicknesses = vv2
  res@xyLineColors      = vv3
  res@tfPolyDrawOrder = "PreDraw"

  obj = gsn_xy(wks,xx,yy,res)

  return(obj)

end

;-- head of xyaddmarker ------------------------------------------------

undef("xyaddmarker")
procedure xyaddmarker(obj:graphic, i_mark:integer,                  \
                      s1:string, v1:numeric, s2:string, v2:numeric, \
                      s3:string, v3:integer)
local obj, i_mark, s1, v1, s2, v2, s3, v3, \
      spec, nl, i, mm, vv1, vv2, vv3, res, obj

;-- body of xyaddmarker ------------------------------------------------

begin

  if ( s1 .ne. "S" .or. s2 .ne. "T" .or. s3 .ne. "C" ) then
    print("Set the order of variables in 'xyaddmarker' :  S-T-C")
    exit
  end if

  getvalues obj
   "xyCoordDataSpec" : spec
  end getvalues

  nl = dimsizes(spec)

  mm  = new(max((/nl,dimsizes(i_mark)/)),integer)
  vv1 = new(max((/nl,dimsizes(v1    )/)),float  )
  vv2 = new(max((/nl,dimsizes(v2    )/)),float  )
  vv3 = new(max((/nl,dimsizes(v3    )/)),integer)

  mm (:dimsizes(i_mark)-1) = i_mark(:)
  vv1(:dimsizes(v1    )-1) = v1    (:)
  vv2(:dimsizes(v2    )-1) = v2    (:)
  vv3(:dimsizes(v3    )-1) = v3    (:)

  if (dimsizes(i_mark) .eq. 1) then
    mm(:) = i_mark
  end if
  if (dimsizes(v1) .eq. 1) then
    vv1(:) = v1
  end if
  if (dimsizes(v2) .eq. 1) then
    vv2(:) = v2
  end if
  if (dimsizes(v3) .eq. 1) then
    vv3(:) = v3
  end if

  do i=0, nl-1
    setvalues spec(i)
     "xyMarkLineMode"      : "MarkLines"
     "xyMarker"            : mm (i)
     "xyMarkerSizeF"       : vv1(i)*0.01
     "xyMarkerThicknessF"  : vv2(i)*0.1
     "xyMarkerColor"       : vv3(i)
    end setvalues
  end do

  return(obj)

end

;-- head of xyaddline --------------------------------------------------

undef("xyaddline")
procedure xyaddline(obj:graphic, x_array:numeric, y_array:numeric, \
                    s1:string, v1:integer, s2:string, v2:numeric,  \
                    s3:string, v3:integer)
local obj, x_array, y_array, s1, v1, s2, v2, s3, v3, \
      nx2, ny2, nl, xx, yy, vv1, vv2, vv3, spec

;-- body of xyaddline --------------------------------------------------

begin

  if ( s1 .ne. "D" .or. s2 .ne. "T" .or. s3 .ne. "C" ) then
    print("Set the order of variables in 'xyaddline' :  D-T-C")
    exit
  end if

  nx2 = dimsizes(x_array)
  ny2 = dimsizes(y_array)

  if ( dimsizes(nx2) .gt. 2 .or. dimsizes(ny2) .gt. 2 ) then
    print("Dimensions of input in 'xyaddline' should be less than 3.")
    exit
  end if

  nl = 1
  if (dimsizes(nx2) .eq. 2) then
    nl = nx2(0)
  end if
  if (dimsizes(ny2) .eq. 2) then
    nl = ny2(0)
  end if

  if ( dimsizes(nx2) .eq. 1 .and. nx2 .eq. 1 ) then
    xx = new(ny2,float)
    xx(:) = x_array
  else
    xx = x_array
  end if
  if ( dimsizes(ny2) .eq. 1 .and. ny2 .eq. 1 ) then
    yy = new(nx2,float)
    yy(:) = y_array
  else
    yy = y_array
  end if

  addarray = create "ADD_LINE" coordArraysClass defaultapp
   "caXArray" : xx
   "caYArray" : yy
   "caXMissingV" : 1.e32
   "caYMissingV" : 1.e32
  end create

  vv1 = new(max((/nl,dimsizes(v1)/)),integer)
  vv2 = new(max((/nl,dimsizes(v2)/)),float  )
  vv3 = new(max((/nl,dimsizes(v3)/)),integer)

  vv1(:dimsizes(v1)-1) = v1(:)
  vv2(:dimsizes(v2)-1) = v2(:)
  vv3(:dimsizes(v3)-1) = v3(:)

  if (dimsizes(v1) .eq. 1) then
    vv1(:) = v1
  end if
  if (dimsizes(v2) .eq. 1) then
    vv2(:) = v2
  end if
  if (dimsizes(v3) .eq. 1) then
    vv3(:) = v3
  end if

  spec = NhlAddData(obj,"xyCoordData",addarray)

  setvalues spec
   "xyDashPatterns"    : vv1
   "xyLineThicknesses" : vv2
   "xyLineColors"      : vv3
  end setvalues

end

;-- head of lgplot -----------------------------------------------------

undef("lgplot")
function lgplot(wks:graphic, label:string, s1:string, v1:integer, \
                s2:string, v2:numeric, s3:string, v3:integer,     \
                lgjust:integer)
local wks, label, s1, v1, s2, v2, s3, v3, lgjust, \
      nl, vv1, vv2, vv3, juststr, obj

;-- body of lgplot -----------------------------------------------------

begin

  if ( s1 .ne. "D" .or. s2 .ne. "T" .or. s3 .ne. "C" ) then
    print("Set the order of variables in 'lgplot' :  D-T-C")
    exit
  end if
  
  nl = dimsizes(label)

  vv1 = new(max((/nl,dimsizes(v1)/)),integer)
  vv2 = new(max((/nl,dimsizes(v2)/)),float  )
  vv3 = new(max((/nl,dimsizes(v3)/)),integer)

  vv1(:dimsizes(v1)-1) = v1(:)
  vv2(:dimsizes(v2)-1) = v2(:)
  vv3(:dimsizes(v3)-1) = v3(:)

  if (dimsizes(v1) .eq. 1) then
    vv1(:) = v1
  end if
  if (dimsizes(v2) .eq. 1) then
    vv2(:) = v2
  end if
  if (dimsizes(v3) .eq. 1) then
    vv3(:) = v3
  end if

  if (lgjust .le. 1) then
    juststr = "TopRight"
  end if
  if (lgjust .eq. 2) then
    juststr = "TopLeft"
  end if
  if (lgjust .eq. 3) then
    juststr = "BottomLeft"
  end if
  if (lgjust .eq. 4) then
    juststr = "BottomRight"
  end if

  obj = create "lg1" legendClass wks
   "lgAutoManage"        : False
   "lgItemCount"         : nl
   "lgItemOrder"         : ispan(nl-1,0,1)
   "lgLabelStrings"      : label
   "lgLabelJust"         : "CenterLeft"
   "lgLabelOffsetF"      : 0.15
   "lgMonoLineThickness" : False
   "lgDashIndexes"       : vv1
   "lgLineThicknesses"   : vv2
   "lgLineColors"        : vv3
   "lgPerimOn"           : False
   "lgJustification"     : juststr
   "lgTopMarginF"        : 0.15
   "lgBottomMarginF"     : 0.15
   "lgRightMarginF"      : 0.15
   "lgLeftMarginF"       : 0.15
  end create

  return(obj)

end

;-- head of cnlineplot -------------------------------------------------

undef("cnlineplot")
function cnlineplot(wks:graphic, x_array:numeric, y_array:numeric, \
                    field:numeric, lev:numeric, llab:integer)
local wks, x_array, y_array, field, lev, llab, \
      res, tag_lev, obj

;-- body of cnlineplot -------------------------------------------------

begin

  res          = True
  res@gsnDraw  = False
  res@gsnFrame = False
  axisvalue_form(res)

  res@sfXArray = x_array
  res@sfYArray = y_array

  tag_lev = True
  if ( dimsizes(dimsizes(lev)) .eq. 1 .and. dimsizes(lev) .eq. 1 ) then
    if (abs(lev) .eq. 999) then
      tag_lev = False
    end if
  end if

  if ( tag_lev ) then
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLevels             = lev
  end if

  res@cnInfoLabelOn = False

  if (llab .eq. 0) then
    res@cnLineLabelsOn = False
  else
    res@cnLineLabelsOn             = True
    res@cnLineLabelPlacementMode   = "Constant"
    res@cnLineLabelFont            = 21
    res@cnLineLabelFontHeightF     = p2n(8)
    res@cnLineLabelInterval        = 1
    res@cnLineDashSegLenF          = 0.2
    res@cnLineLabelFormat          = "@*+^sg"
    res@cnLineLabelBackgroundColor = -1
  end if

  obj = gsn_contour(wks,field,res)

  obj@xarrminmax = (/min(x_array),max(x_array)/)
  obj@yarrminmax = (/min(y_array),max(y_array)/)

  return(obj)

end

;-- head of cnshadeplot ------------------------------------------------

undef("cnshadeplot")
function cnshadeplot(wks:graphic, x_array:numeric, y_array:numeric, \
                     field:numeric, lev:numeric, lcol:integer, \
                     label)
local wks, x_array, y_array, field, lev, lcol, label, \
      res, tag_lev, tag_lcol, obj

;-- body of cnshadeplot ------------------------------------------------

begin

  res          = True
  res@gsnDraw  = False
  res@gsnFrame = False
  axisvalue_form(res)

  res@sfXArray = x_array
  res@sfYArray = y_array

  res@cnLinesOn = False
  res@cnFillOn  = True

  if ( label .ne. 0 .and. label .ne. "" ) then
    res@lbAutoManage             = False
    res@pmLabelBarDisplayMode    = "Always"
    res@lbOrientation            = "Horizontal"
    res@pmLabelBarSide           = "Bottom"
    res@pmLabelBarOrthogonalPosF = 0.04
    res@pmLabelBarKeepAspect     = True
    res@pmLabelBarHeightF        = 0.05
    res@pmLabelBarWidthF         = 0.5
    res@lbLabelAlignment         = "InteriorEdges"
    res@lbLabelJust              = "CenterCenter"
    res@lbLabelAutoStride        = True
    res@lbLabelFont              = 21
    res@lbLabelFontHeightF       = p2n(11)
    res@cnLineLabelFormat        = "@*+^sg"
    res@lbPerimOn                = False
  end if

  if ( typeof(label) .eq. "string" ) then
    res@lbTitleString      = label
    res@lbTitleFont        = 21
    res@lbTitleFontHeightF = p2n(13)
    res@lbTitlePosition    = "Bottom"
    res@lbTitleOffsetF     = 0.2
  end if

  tag_lev = True
  if ( dimsizes(dimsizes(lev)) .eq. 1 .and. dimsizes(lev) .eq. 1 ) then
    if (abs(lev) .eq. 999) then
      tag_lev = False
    end if
  end if
  tag_lcol = True
  if ( dimsizes(dimsizes(lcol)) .eq. 1 .and. dimsizes(lcol) .eq. 1 ) then
    if (abs(lcol) .eq. 999) then
      tag_lcol = False
    end if
  end if

  if ( tag_lev ) then
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLevels             = lev
  end if
  if ( tag_lcol ) then
    res@cnFillColors = lcol
  end if

  res@cnInfoLabelOn  = False
  res@cnLineLabelsOn = False

  obj = gsn_contour(wks,field,res)

  obj@xarrminmax = (/min(x_array),max(x_array)/)
  obj@yarrminmax = (/min(y_array),max(y_array)/)

  return(obj)

end

;-- head of cnline_att -------------------------------------------------

undef("cnline_att")
procedure cnline_att(obj:graphic, s1:string, v1:integer, s2:string, \
                     v2:numeric, s3:string, v3:integer)
local obj, s1, v1, s2, v2, s3, v3, \
      lev, nl, itag, ierr

;-- body of cnline_att -------------------------------------------------

begin

  if ( s1 .ne. "D" .or. s2 .ne. "T" .or. s3 .ne. "C" ) then
    print("Set the order of variables in 'cnline_att' :  D-T-C")
    exit
  end if

  setvalues obj
   "cnLinesOn" : True
  end setvalues

  getvalues obj
   "cnLevels" : lev
  end getvalues
  nl = dimsizes(lev)
  itag = new(nl,integer)
  do i=0, nl-1
    if (lev(i) .eq. 0.) then
      itag(i) = 1
    else
      if (lev(i) .gt. 0.) then
        itag(i) = 0
      else
        itag(i) = 2
      end if
    end if
  end do

  ierr = 0
  if (dimsizes(v1) .eq. 1) then
    setvalues obj
     "cnLineDashPattern" : v1
    end setvalues
  else
    if (dimsizes(v1) .eq. 3) then
      setvalues obj
       "cnMonoLineDashPattern" : False
       "cnLineDashPatterns"    : v1(itag(:))
      end setvalues
    else
      ierr = 1
    end if
  end if
  if (dimsizes(v2) .eq. 1) then
    setvalues obj
     "cnLineThicknessF" : v2
    end setvalues
  else
    if (dimsizes(v2) .eq. 3) then
      setvalues obj
       "cnMonoLineThickness" : False
       "cnLineThicknesses"   : v2(itag(:))
      end setvalues
    else
      ierr = 1
    end if
  end if
  if (dimsizes(v3) .eq. 1) then
    setvalues obj
     "cnLineColor"          : v3
     "cnLineLabelFontColor" : v3
    end setvalues
  else
    if (dimsizes(v3) .eq. 3) then
      setvalues obj
       "cnMonoLineColor"          : False
       "cnLineColors"             : v3(itag(:))
       "cnMonoLineLabelFontColor" : False
       "cnLineLabelFontColors"    : v3(itag(:))
      end setvalues
    else
      ierr = 1
    end if
  end if

  if (ierr .ne. 0) then
    print("The dimension of variables in 'cnline_att' must be 1 or 3.")
    exit
  end if

end

;-- head of cnaddline --------------------------------------------------

undef("cnaddline")
procedure cnaddline(obj:graphic, wks:graphic, x_array:numeric, \
                    y_array:numeric, field:numeric, lev:numeric, \
                    llab:integer, s1:string, v1:integer, s2:string, \
                    v2:numeric, s3:string, v3:integer)
local obj, wks, x_array, y_array, field, lev, s1, v1, s2, v2, s3, v3, \
      obj2

;-- body of cnaddline --------------------------------------------------

begin

  obj2 = cnlineplot(wks,x_array,y_array,field,lev,llab)
  cnline_att(obj2,s1,v1,s2,v2,s3,v3)
  NhlAddOverlay(obj,obj2,obj)

end

;-- head of cnfinalize -------------------------------------------------

undef("cnfinalize")
procedure cnfinalize(obj:graphic, wks:graphic)
local obj, wks, \
      obj_ol, xf, yf, wi, he

;-- body of cnfinalize -------------------------------------------------

begin

  obj_ol = obj

  obj = create "contour" logLinPlotClass wks
  end create
  axis_range(obj,obj_ol@xleft,obj_ol@xright,obj_ol@ybottom,obj_ol@ytop, \
             obj_ol@xlog,obj_ol@ylog)
  if ( obj_ol@xleft .eq. 999 .or. obj_ol@xright .eq. 999) then
    setvalues obj
     "trXMaxF" : obj@xarrminmax(1)
    end setvalues
  end if
  if ( obj_ol@xleft .eq. -999 .or. obj_ol@xright .eq. -999) then
    setvalues obj
     "trXMinF" : obj@xarrminmax(0)
    end setvalues
  end if
  if ( obj_ol@ybottom .eq. 999 .or. obj_ol@ytop .eq. 999) then
    setvalues obj
     "trYMaxF" : obj@yarrminmax(1)
    end setvalues
  end if
  if ( obj_ol@ybottom .eq. -999 .or. obj_ol@ytop .eq. -999) then
    setvalues obj
     "trYMinF" : obj@yarrminmax(0)
    end setvalues
  end if
  if ( obj_ol@xlog ) then
    setvalues obj
     "trXLog" : True
    end setvalues
    setvalues obj_ol
     "tmXBMinorPerMajor" : 8
     "tmXTMinorPerMajor" : 8
    end setvalues
  end if
  if ( obj_ol@ylog ) then
    setvalues obj
     "trYLog" : True
    end setvalues
    setvalues obj_ol
     "tmYLMinorPerMajor" : 8
     "tmYRMinorPerMajor" : 8
    end setvalues
  end if

  getvalues obj_ol
   "vpXF"      : xf
   "vpYF"      : yf
   "vpWidthF"  : wi
   "vpHeightF" : he
  end getvalues
  setvalues obj
   "vpXF"      : xf
   "vpYF"      : yf
   "vpWidthF"  : wi
   "vpHeightF" : he
  end setvalues

  NhlAddOverlay(obj,obj_ol,obj)

end

