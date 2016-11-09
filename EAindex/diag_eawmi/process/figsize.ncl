p_xpg = 810.  ; modified from 612. for letter (default)
p_ypg = 792.

;xspace = 170  ; 83, 170 - QJRMS  ; 80, 140, 165 - AMS
;yspace = 240  ; 240     - QJRMS  ; 230          - AMS
if (.not. isdefined("xspace")) then
  xspace = 165
end if
if (.not. isdefined("yspace")) then
  yspace = 230
end if

undef("m2p")
function m2p(mm:numeric)
local vpt
begin
  vpt = mm*612./215.9
  return(vpt)
end

undef("m2n")
function m2n(mm:numeric)
local vn
begin
  if (typeof(mm) .eq. "integer") then
    vn = int2flt(mm)/max((/xspace,yspace/))
  else
    vn = mm/max((/xspace,yspace/))
  end if
  return(vn)
end

undef("p2n")
function p2n(p:numeric)
local vn
begin
  vn = p*215.9/612./max((/xspace,yspace/))
  return(vn)
end

axisfontratio = 1.0
sz_mstr = p2n(6)
sz_xstr = p2n(5)
sz_ystr = sz_xstr
sz_xlb  = p2n(4)
sz_ylb  = sz_xlb

sz_majtck = -sz_xlb*0.6
sz_mintck = sz_majtck*0.7

sz_lbt = p2n(4.5)
sz_lb  = sz_xlb

;-- head of cr_wks -----------------------------------------------------

undef("cr_wks")
function cr_wks(out_fname:string, res:integer, i_outline:integer)
local out_fname, res, i_outline, \
      wkdlx, wkdux, wkdly, wkduy, wks, fname, fform, outline

;-- body of cr_wks -----------------------------------------------------

begin

  wkdlx = 0.5*(p_xpg-m2p(max((/xspace,yspace/))))
  wkdux = p_xpg - wkdlx
  wkdly = 0.5*(p_ypg-m2p(max((/xspace,yspace/))))
  wkduy = p_ypg - wkdly

  if ( isStrSubset(out_fname,".ps") ) then
    fname = out_fname
    fform = "ps"
  else
    fname = out_fname+".eps"
    fform = "eps"
  end if

  wks = create "ps" psWorkstationClass defaultapp
   "wkPSFileName"   : fname
   "wkPSFormat"     : fform
   "wkPSResolution" : res
   "wkDeviceLowerX" : wkdlx
   "wkDeviceUpperX" : wkdux
   "wkDeviceLowerY" : wkdly
   "wkDeviceUpperY" : wkduy
   "wkBackgroundColor" : (/1.,1.,1./)
  end create

  if (i_outline .ne. 0) then
    outline =  create "space" xyPlotClass wks
     "tmXBOn" : False
     "tmXTOn" : False
     "tmYLOn" : False
     "tmYROn" : False
     "vpXF" : 0.
     "vpYF" : 1.
     "vpWidthF"  : m2n(xspace)
     "vpHeightF" : m2n(yspace)
    end create
    draw(outline)
  end if

  return(wks)

end

;-- head of set_size0 --------------------------------------------------

undef("set_size0")
procedure set_size0(obj:graphic, xf:numeric, yf:numeric, \
                    widthf:numeric, heightf:numeric)
local obj, xf, yf, widthf, heightf

;-- body of set_size0 --------------------------------------------------

begin

  setvalues obj
   "vpXF"    : m2n(xf)
   "vpYF"    : 1.-m2n(yf)
   "vpWidthF"  : m2n(widthf)
   "vpHeightF" : m2n(heightf)
  end setvalues

end

;-- head of set_size ---------------------------------------------------

undef("set_size")
procedure set_size(obj:graphic, xf:numeric, yf:numeric, \
                   widthf:numeric, heightf:numeric)
local obj, xf, yf, widthf, heightf, \
      classname

;-- body of set_size ---------------------------------------------------

begin

  set_size0(obj,xf,yf,widthf,heightf)

  classname = NhlClassName(obj)

  if (classname .ne. "logLinPlotClass") then

    setvalues obj
     "tmBorderThicknessF"  : 1.0
     "tmXBMajorThicknessF" : 1.0
     "tmXTMajorThicknessF" : 1.0
     "tmYLMajorThicknessF" : 1.0
     "tmYRMajorThicknessF" : 1.0
     "tmXBMinorThicknessF" : 0.5
     "tmXTMinorThicknessF" : 0.5
     "tmYLMinorThicknessF" : 0.5
     "tmYRMinorThicknessF" : 0.5

     "tmXBMajorLengthF"  : sz_majtck*axisfontratio
     "tmXTMajorLengthF"  : sz_majtck*axisfontratio
     "tmYLMajorLengthF"  : sz_majtck*axisfontratio
     "tmYRMajorLengthF"  : sz_majtck*axisfontratio
     "tmXBMinorLengthF"  : sz_mintck*axisfontratio
     "tmXTMinorLengthF"  : sz_mintck*axisfontratio
     "tmYLMinorLengthF"  : sz_mintck*axisfontratio
     "tmYRMinorLengthF"  : sz_mintck*axisfontratio

     "tmXBLabelFont" : 21
     "tmXTLabelFont" : 21
     "tmYLLabelFont" : 21
     "tmYRLabelFont" : 21
     "tmXBLabelFontHeightF" : sz_xlb*axisfontratio
     "tmXTLabelFontHeightF" : sz_xlb*axisfontratio
     "tmYLLabelFontHeightF" : sz_ylb*axisfontratio
     "tmYRLabelFontHeightF" : sz_ylb*axisfontratio

     "tiMainFont"  : 22
     "tiXAxisFont" : 22
     "tiYAxisFont" : 22
     "tiMainFontHeightF"  : sz_mstr*axisfontratio
     "tiXAxisFontHeightF" : sz_xstr*axisfontratio
     "tiYAxisFontHeightF" : sz_ystr*axisfontratio
    end setvalues

  end if

end

;-- head of set_size_mp ------------------------------------------------

undef("set_size_mp")
procedure set_size_mp(obj:graphic, xf:numeric, yf:numeric, \
                      widthf:numeric, heightf:numeric)
local obj, xf, yf, widthf, heightf

;-- body of set_size_mp ------------------------------------------------

begin

  set_size0(obj,xf,yf,widthf,heightf)
  setvalues obj
   "mpPerimLineThicknessF" : 1.0
   "tmBorderThicknessF"  : 1.0
   "tmXBMajorThicknessF" : 1.0
   "tmXTMajorThicknessF" : 1.0
   "tmYLMajorThicknessF" : 1.0
   "tmYRMajorThicknessF" : 1.0
   "tmXBMinorThicknessF" : 0.5
   "tmXTMinorThicknessF" : 0.5
   "tmYLMinorThicknessF" : 0.5
   "tmYRMinorThicknessF" : 0.5

   "tmXBMajorLengthF"  : sz_majtck*axisfontratio
   "tmXTMajorLengthF"  : sz_majtck*axisfontratio
   "tmYLMajorLengthF"  : sz_majtck*axisfontratio
   "tmYRMajorLengthF"  : sz_majtck*axisfontratio
   "tmXBMinorLengthF"  : sz_mintck*axisfontratio
   "tmXTMinorLengthF"  : sz_mintck*axisfontratio
   "tmYLMinorLengthF"  : sz_mintck*axisfontratio
   "tmYRMinorLengthF"  : sz_mintck*axisfontratio

   "tmXBLabelFont" : 21
   "tmXTLabelFont" : 21
   "tmYLLabelFont" : 21
   "tmYRLabelFont" : 21
   "tmXBLabelFontHeightF" : sz_xlb*axisfontratio
   "tmXTLabelFontHeightF" : sz_xlb*axisfontratio
   "tmYLLabelFontHeightF" : sz_ylb*axisfontratio
   "tmYRLabelFontHeightF" : sz_ylb*axisfontratio
  end setvalues

end

;-- head of set_size_lb ------------------------------------------------

undef("set_size_lb")
procedure set_size_lb(obj:graphic, xf:numeric, yf:numeric, \
                      widthf:numeric, heightf:numeric)
local obj, xf, yf, widthf, heightf

;-- body of set_size_lb ------------------------------------------------

begin

  set_size0(obj,xf,yf,widthf,heightf)
  setvalues obj
   "lbAutoManage" : False

   "lbBoxLineThicknessF" : 0.6
   "lbPerimThicknessF"   : 0.6

   "lbLabelFont"        : 21
   "lbLabelFontHeightF" : sz_lb*axisfontratio
   "lbTitleFont"        : 22
   "lbTitleFontHeightF" : sz_lbt*axisfontratio
  end setvalues

end

;-- head of set_size_lg ------------------------------------------------

undef("set_size_lg")
procedure set_size_lg(obj:graphic, xf:numeric, yf:numeric, \
                      widthf:numeric, heightf:numeric)
local obj, xf, yf, widthf, heightf

;-- body of set_size_lg ------------------------------------------------

begin

  set_size0(obj,xf,yf,widthf,heightf)
  setvalues obj
   "lgAutoManage" : False

   "lgBoxLineThicknessF" : 0.6
   "lgPerimThicknessF"   : 0.6

   "lgLabelFont"        : 21
   "lgLabelFontHeightF" : sz_lbt*axisfontratio
   "lgTitleFont"        : 22
   "lgTitleFontHeightF" : sz_lbt*axisfontratio

   "lgLineDashSegLenF"    : 0.05
  end setvalues

end

;-- head of set_size_csmtck --------------------------------------------

undef("set_size_csmtck")
procedure set_size_csmtck(obj:graphic, xf:numeric, yf:numeric, \
                          widthf:numeric, heightf:numeric)
local obj, xf, yf, widthf, heightf

;-- body of set_size_csmtck --------------------------------------------

begin

  set_size0(obj,xf,yf,widthf,heightf)
  setvalues obj
   "tmBorderThicknessF"  : 1.0
   "tmXBMajorThicknessF" : 1.0
   "tmXTMajorThicknessF" : 1.0
   "tmYLMajorThicknessF" : 1.0
   "tmYRMajorThicknessF" : 1.0
   "tmXBMinorThicknessF" : 0.5
   "tmXTMinorThicknessF" : 0.5
   "tmYLMinorThicknessF" : 0.5
   "tmYRMinorThicknessF" : 0.5

   "tmXBLabelFont" : 21
   "tmXTLabelFont" : 21
   "tmYLLabelFont" : 21
   "tmYRLabelFont" : 21
   "tmXBLabelFontHeightF" : sz_xlb*axisfontratio
   "tmXTLabelFontHeightF" : sz_xlb*axisfontratio
   "tmYLLabelFontHeightF" : sz_ylb*axisfontratio
   "tmYRLabelFontHeightF" : sz_ylb*axisfontratio

   "tiMainFont"  : 22
   "tiXAxisFont" : 22
   "tiYAxisFont" : 22
   "tiMainFontHeightF"  : sz_mstr*axisfontratio
   "tiXAxisFontHeightF" : sz_xstr*axisfontratio
   "tiYAxisFontHeightF" : sz_ystr*axisfontratio
  end setvalues

end

;-- head of dr_txt -----------------------------------------------------

undef("dr_txt")
procedure dr_txt(wks:graphic, str:string, xf:numeric, yf:numeric, \
                 font:integer, fonthf:numeric, ang:numeric)
local obj, wks, str, xf, yf, font, fonthf, ang

;-- body of dr_txt -----------------------------------------------------

begin

  obj = create "Text" textItemClass wks
    "txString" : str
    "txPosXF"  : m2n(xf)
    "txPosYF"  : 1.-m2n(yf)
    "txFont"   : font
    "txFontHeightF" : p2n(fonthf)
    "txAngleF"      : ang
  end create
  draw(obj)

end

;-- head of dr_txt2 ----------------------------------------------------

undef("dr_txt2")
procedure dr_txt2(wks:graphic, str:string, xf:numeric, yf:numeric, \
                  font:integer, fonthf:numeric, ang:numeric,       \
                  just:string, col:numeric)
local obj, wks, str, xf, yf, font, fonthf, ang, just, col

;-- body of dr_txt2 ----------------------------------------------------

begin

  obj = create "Text" textItemClass wks
    "txString" : str
    "txPosXF"  : m2n(xf)
    "txPosYF"  : 1.-m2n(yf)
    "txFont"   : font
    "txFontHeightF" : p2n(fonthf)
    "txAngleF"      : ang
    "txJust"        : just
    "txFontColor"   : col
  end create
  draw(obj)

end


