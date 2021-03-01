function fire_getorderimage,tstr,im

  ;; maybe use fire_rectify_order.pro here

  npix = 2048

  if im.subimage eq 1 then begin
    print,'Input must be original image'
    return,-1
  endif

  ;; Trace structure
  if tag_exist(tstr,'bndy0coef') then begin
    bndy0coef = tstr.bndy0coef
    bndy1coef = tstr.bndy1coef
    bndx0 = tstr.bndx0
    bndx1 = tstr.bndx1    
     
  ;; Boundary structure
  endif else begin
    bndy0coef = tstr.y0coef
    bndy1coef = tstr.y1coef
    bndx0 = tstr.lo
    bndx1 = tstr.hi
  endelse
  
  ;; Get aperture subimage using boundary
  x = lindgen(npix)
  xx = lindgen(npix)#replicate(1,npix)
  yy = replicate(1,npix)#lindgen(npix)

  x1 = x[bndx0:bndx1]
  ylo = poly(x1,bndy0coef) > 0
  yhi = poly(x1,bndy1coef) < (npix-1)
  y0 = round(min(ylo)) > 0
  y1 = round(max(yhi)) < (npix-1)
  subflux = im.flux[bndx0:bndx1,y0:y1]
  suberr = im.err[bndx0:bndx1,y0:y1]
  submask = im.mask[bndx0:bndx1,y0:y1]  
  subxx = xx[bndx0:bndx1,y0:y1]
  subyy = yy[bndx0:bndx1,y0:y1]    
  mask = (subyy ge (poly(subxx,bndy0coef)>0) and $
          subyy le (poly(subxx,bndy1coef)<(npix-1)))
  subsz = size(subflux)
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
         exptype:im.exptype,subimage:1,rectified:0,$
         tstr:tstr}
  
  return,out

end
