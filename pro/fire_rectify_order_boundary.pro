function fire_rectify_order_boundary,bstr,im

  npix = 2048
  x = findgen(npix)
  xx = findgen(npix)#replicate(1,npix)
  yy = replicate(1,npix)#findgen(npix)

  ;; Get subimage
  x1 = x[bstr.lo:bstr.hi]
  ylo = poly(x1,bstr.y0coef) > 0
  yhi = poly(x1,bstr.y1coef) < (npix-1)
  y0 = round(min(ylo)) > 0
  y1 = round(max(yhi)) < (npix-1)
  im2 = im[bstr.lo:bstr.hi,y0:y1]
  xx2 = xx[bstr.lo:bstr.hi,y0:y1]
  yy2 = yy[bstr.lo:bstr.hi,y0:y1]    
  mask = (yy2 ge (poly(xx2,bstr.y0coef)>0) and $
          yy2 le (poly(xx2,bstr.y1coef)<(npix-1)))

  ;; Rectify
  nx = bstr.hi-bstr.lo+1
  ny = max((round(yhi)<(npix-1)) - (round(ylo)>0))+1
  recim2 = fltarr(nx,ny)

  for j=0,nx-1 do begin
    ysrt = round(ylo[j])-y0
    yend = round(yhi[j])-y0
    recim2[j,0] = (im2*mask)[j,ysrt:yend]
  endfor

  return, recim
  
end
