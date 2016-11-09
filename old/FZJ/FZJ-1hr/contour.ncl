;***********************************************************************;
; Function : gsn_csm_contour_map_polar                                  ;
;                   wks: workstation object                             ;
;                  data: 2-dimensional data                             ;
;               resources: optional resources                           ;
;                                                                       ;
; This function creates and draws a contour plot over a polar           ;
; stereographic map plot to the workstation "wks" (the variable         ;
; returned from a previous call to "gsn_open_wks").  "data" is the      ;
; 2-dimensional data to be contoured, and "resources" is an optional    ;
; list of resources. The id of the map plot is returned.                ;
;                                                                       ;
; This function behaves differently from gsn_contour in that it will    ;
; create a special kind of plot if certain attributes are set.          ;
;                                                                       ;
;    1. By default, a cyclic point is added.  If gsnAddCyclic is set to ;
;       False, then the cyclic point is not added.                      ;
;    2. The longitude lines are labeled.                                ;
;    3. If any of the special GSUN resources "gsnLeftString,"           ;
;       "gsnCenterString," and/or "gsnRightString" are set, then they   ;
;       are used to title the top of the plot.                          ;
;       (in addition, the regular resource "tiMainString" can be set to ;
;       create a regular title).                                        ;
;    4. If the resource "cnFillOn" is set to True, then a labelbar is   ;
;       drawn and line labels are turned off.                           ;
;    5. If data has an attribute called "long_name," and "gsnLeftString";
;       hasn't explicitly been set, then the value of this attribute    ;
;       is used for the left string title.                              ;
;    6. If data has an attribute called "units," and "gsnRightString"   ;
;       hasn't explicitly been set, then the value of this attribute    ;
;       is used for the right string title.                             ;
;    7. If either of the resources "gsnPolarNH" or "gsnPolarSH" are set ;
;       to True then only the northern or southern hemisphere is        ;
;       displayed.                                                      ;
;***********************************************************************;

function contour_angle(wks:graphic,data:numeric,resources:logical)
local i, contour_object, map_object, res, lbar_pos, lbar_zone, \
lbar_on, calldraw, callframe, main_offset, datanew, res2, cnres, mpres, \
vpwf, vphf, font_height, levels, colors, zonal_zone, lbres, \
infolabel_on, infolabel_zone, lbar_zone, lbar_height, lbar_width
begin

; Initialize.
    lbar_on       = False     ; Labelbar flag
    res2          = get_resources(resources) ; Copy of resources
    lbar_zone     = 3         ; Zone for labelbar (may change later)
;    mpres         = True      ; Will hold map resources
    mpres         = False      ; Will hold map resources

; Check for cyclic point (default is True if we are not setting both
; lon2d/lat2d attributes).

;    set_cyclic = .not.(isatt(data,"lat2d").and.isatt(data,"lon2d"))
;    if(get_res_value(res2,"gsnAddCyclic",set_cyclic)) then
;      datanew = gsn_add_cyclic_point(data)
;    else
      datanew = data
;    end if

; Check for left/right titles at top. Use valid "long_name" type attributes
; and "units" if they exist.

    set_right_subtitle(datanew,res2,mpres)
    set_left_subtitle(datanew,res2,mpres)

; Check for draw and frame.

    calldraw  = get_res_value(res2,"gsnDraw", True)
    callframe = get_res_value(res2,"gsnFrame",True)
    maxbb     = get_bb_res(res2)

;
; Check for type of polar plot and polar labels desired.
;
    mpres@gsnPolar     = get_polar_type(res2)
    mpres@gsnPolarTime = get_res_value(res2,"gsnPolarTime",False)
    mpres@gsnPolarUT   = get_res_value(res2,"gsnPolarUT",0.)

; Check for subtitles at top and add to mpres if exist.

;    set_subtitles_res(res2,mpres)
;
; Tickmarks.
;
    mpres@gsnTickMarksOn = get_res_value(res2,"gsnTickMarksOn",True)

; Create data object and use coordinate variables if they exist.

    check_for_y_lat_coord(datanew,res2,"contour_map")
    check_for_lon_coord(datanew,res2,"contour_map")

; Set some more contour plot resources.

    set_attr(res2,"cnLineLabelBackgroundColor",-1)
    set_attr(res2,"cnInfoLabelOrthogonalPosF",0.06)
    set_attr(res2,"cnInfoLabelZone",2)

;
; By default, mpOutlineOn is False, unless cnFillOn is set to True
; or mpFillOn is set to False, then it is set back to True.
;
    if(check_attr(res2,"cnFillOn",True,False).or.\
       check_attr(res2,"mpFillOn",False,False))
      set_attr(res2,"mpOutlineOn",True)
    end if

; This section tests for regular resources.

    lbres = get_res_eq(res2,(/"lb","pm"/))
    mpres = get_res_eq(res2,(/"mp","ti","vp","tx","am","pmA","pmO","pmT",\
                              "gsnPolar"/))
    cnres = get_res_ne(res2,(/"mp","vp","lb","tx","am","gsnPolar"/))

    if(cnres)
      if(check_attr(cnres,"cnFillOn",True,False))
        set_attr(cnres,"cnInfoLabelOn",False)
        if(.not.isatt(lbres,"lbLabelBarOn").or.\
          check_attr(lbres,"lbLabelBarOn",True,False))
          set_attr(cnres,"cnLineLabelsOn",False)
          lbar_on = True   ; Turn on a labelbar
        end if
      end if
    end if
;
; Compute zone for labelbar if it is supposed to get drawn.
; Zone for labelbar changes only if an info label is not drawn.
;
    if(check_attr(cnres,"cnInfoLabelOn",False,False))
      lbar_zone = 2
    end if

; Before we create the objects, turn off draw and frame for them.
    cnres           = True 
    mpres           = True
    cnres@gsnDraw   = False
    cnres@gsnFrame  = False
    mpres@gsnDraw  = False
    mpres@gsnFrame = False

    contour_object = gsn_contour(wks,datanew,cnres) ; Create contours.
    map_object = gsn_csm_map(wks,mpres)             ; Create map.
    overlay(map_object,contour_object)              ; Overlay contour plot
                                                    ; on map.
;
; Retrieve the view port location of the map plot so we know where
; to put titles and labelbar (if there is one).
;
    getvalues map_object
      "vpWidthF"  : vpwf
      "vpHeightF" : vphf
    end getvalues

; Make sure axes labels are the same size.

    font_height = 0.02 * vphf     ; Make various label sizes a function
                                  ; of the height of the viewport.
;
; Check if user setting own font heights.
;
    main_font_height = get_res_value_keep(res2,"tiMainFontHeightF", \
                                          2.*font_height)
    main_font = get_res_value_keep(res2,"tiMainFont","helvetica-bold")

    setvalues map_object
        "tiMainFont"           : main_font        ; main title font
        "tiMainFontHeightF"    : main_font_height ; main title size
    end setvalues

; Create a labelbar.

    if(lbar_on)
      add_labelbar(wks,contour_object,lbar_zone,font_height,"polar",lbres)
    end if
  
    draw_and_frame(wks,map_object,calldraw,callframe,0,maxbb)

; Return plot object and data object (as attribute of plot object).

    map_object@data    = contour_object@data
    map_object@contour = contour_object
    return(map_object)
end
