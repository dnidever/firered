function fire_getorderimage,tstr,im

  ;; maybe use fire_rectify_order.pro here

  npix = 2048

  if im.subimage eq 1 then begin
    print,'Input must be original image'
    return,-1
  endif
  
  ;; Get aperture subimage using boundary
  x = findgen(npix)
  xx = findgen(npix)#replicate(1,npix)
  yy = replicate(1,npix)#findgen(npix)

  x1 = x[tstr.bndx0:tstr.bndx1]
  ylo = poly(x1,tstr.bndy0coef) > 0
  yhi = poly(x1,tstr.bndy1coef) < (npix-1)
  y0 = round(min(ylo)) > 0
  y1 = round(max(yhi)) < (npix-1)
  subflux = im.flux[tstr.bndx0:tstr.bndx1,y0:y1]
  suberr = im.err[tstr.bndx0:tstr.bndx1,y0:y1]
  submask = im.mask[tstr.bndx0:tstr.bndx1,y0:y1]  
  subxx = xx[tstr.bndx0:tstr.bndx1,y0:y1]
  subyy = yy[tstr.bndx0:tstr.bndx1,y0:y1]    
  mask = (subyy ge (poly(subxx,tstr.bndy0coef)>0) and $
          subyy le (poly(subxx,tstr.bndy1coef)<(npix-1)))
  subsz = size(subim)
  subnx = subsz[1]
  subny = subsz[2]

  ;; Masked "aperture" image
  apflux = subflux*mask
  aperr = suberr*mask + (1-mask)*1e30
  apmask = submask*mask
  apx = subxx[*,0]
  apy = reform(subyy[0,*])
  
  
  out = {file:im.file,flux:apflux,err:aperr,mask:apmask,$
         x:apx,y:apy,nx:subnx,ny:subny,head:im.head,$
         exptype:im.exptype,subimage:1,$
         tstr:tstr}
  ;x0:tstr.bndx0,x1:tstr.bndx1,y0:y0,y1:y1,ylo:ylo,yhi:yhi}
  
  return,out

end
