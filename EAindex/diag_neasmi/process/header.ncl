load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

setvalues NhlGetWorkspaceObjectId()
 "wsMaximumSize" : 64000000
end setvalues

res0          = True
res0@gsnDraw  = False
res0@gsnFrame = False


