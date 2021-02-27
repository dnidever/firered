function fire_getorderimage,tstr,im

  ;; maybe use fire_rectify_order.pro here

  npix = 2048

  ;; Get aperture subimage using boundary
  x = findgen(npix)
  xx = findgen(npix)#replicate(1,npix)
  yy = replicate(1,npix)#findgen(npix)

  x1 = x[tstr.bndx0:tstr.bndx1]
  ylo = poly(x1,tstr.bndy0coef) > 0
  yhi = poly(x1,tstr.bndy1coef) < (npix-1)
  y0 = round(min(ylo)) > 0
  y1 = round(max(yhi)) < (npix-1)
  subim = im[tstr.bndx0:tstr.bndx1,y0:y1]
  subxx = xx[tstr.bndx0:tstr.bndx1,y0:y1]
  subyy = yy[tstr.bndx0:tstr.bndx1,y0:y1]    
  submask = (subyy ge (poly(subxx,tstr.bndy0coef)>0) and $
             subyy le (poly(subxx,tstr.bndy1coef)<(npix-1)))
  subsz = size(subim)
  subnx = subsz[1]
  subny = subsz[2]
  ysub = findgen(subny)

  ;; Masked "aperture" image
  apim = subim*submask

  out = {im:apim,mask:submask,x0:tstr.bndx0,x1:tstr.bndx1,y0:y0,y1:y1,ylo:ylo,yhi:yhi}
  
  return,out

end
