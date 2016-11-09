


f0 = addfile("u8_y.nc","r")
f1 = addfile("u8a_y.nc","w")

u = f0->u

um = dim_avg(u)

ua = u
do i=0, 191
  ua(:,:,:,i) = u(:,:,:,i) - um(:,:,:)
end do

f1->u = ua

do n=0, 147
  print( abs( 0.5*(um(n+1,:,:) - um(n,:,:)) ) - dim_max(abs(ua(n,:,:,:))) )
end do
