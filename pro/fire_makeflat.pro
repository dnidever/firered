pro fire_makeflat,files,outfile

  ;;  --- making flat ---
  
  ;files = 'ut131222/fire_0011.fits'
  files = 'ut131222/fire_00'+['11','12','13','14']+'.fits'    
  outfile = 'flat.fits'
  
  nfiles = n_elements(files)
  npix = 2048
  cube = fltarr(npix,npix,nfiles)
  exptime = fltarr(nfiles)
  for i=0,nfiles-1 do begin
    fits_read,files[i],im,head
    exptime[i] = sxpar(head,'exptime')
    print,files[i],exptime[i],median(im),mad(im)
    newim = FIRE_REFCORRECT(im)
    ;; scale by exptime
    cube[*,*,i] = float(newim)/exptime[i]
    displayc,float(newim)/exptime[i],/z,tit=files[i]
    wait,0.3
  endfor
  ;; Median
  if nfiles gt 1 then begin
    med = median(cube,dim=3) > 0
  endif else med = cube
  im = FIRE_MAKEIMAGE(med)
  bpmfile = 'bpm3.fits'
  bpm = FIRE_READIMAGE(bpmfile)
  im = FIRE_BPMCORRECT(im,bpm)
  
  ;; Work on the orders
  if n_elements(bndfile) eq 0 then bndfile = 'fire_boundary_0011.fits'
  bstr = mrdfits(bndfile,1)
  norders = n_elements(bstr)

  flat = im
  flat.flux[*,*] = 1
  for i=0,norders-1 do begin
    apim = FIRE_GETORDERIMAGE(bstr[i],im)

    ;; Get order mask
    bndy0coef = bstr[i].y0coef
    bndy1coef = bstr[i].y1coef
    bndx0 = bstr[i].lo
    bndx1 = bstr[i].hi
  
    ;; Get aperture subimage using boundary
    x = lindgen(npix)
    xx = lindgen(npix)#replicate(1,npix)
    yy = replicate(1,npix)#lindgen(npix)
    x1 = x[bndx0:bndx1]
    ylo = poly(x1,bndy0coef) > 0
    yhi = poly(x1,bndy1coef) < (npix-1)
    y0 = round(min(ylo)) > 0
    y1 = round(max(yhi)) < (npix-1)
    subxx = xx[bndx0:bndx1,y0:y1]
    subyy = yy[bndx0:bndx1,y0:y1]    
    mask1 = (subyy ge (poly(subxx,bndy0coef)>0) and $
             subyy le (poly(subxx,bndy1coef)<(npix-1)))
    ;; Grow "bad" region by 1-2 pixels
    kernel = fltarr(5,5)+1
    mask2 = convol(mask1,kernel)
    mask = (mask2 gt total(kernel)-1)
    
    ;; Median Smooth for bad pixels
    temp = apim.flux
    bdpix = where(apim.mask eq 0,nbdpix)
    temp[bdpix] = !values.f_nan
    medy = median(temp,dim=2)
    medy2d = medy#replicate(1,apim.ny)
    temp2 = apim.flux*mask + (1-mask)*medy2d

    ;; Smooth in vertical direction
    smy = medfilt2d(temp2,11,dim=2,/even,/edge_copy)
    apflat = temp2 / smy

    ;; Stuff it back in
    x0 = apim.x[0]
    x1 = apim.x[-1]
    y0 = apim.y[0]
    y1 = apim.y[-1]
    flatflux = flat.flux[x0:x1,y0:y1]
    flatflux = flatflux*(1-mask) + mask*apflat
    flat.flux[x0:x1,y0:y1] = flatflux

    displayc,apflat,/z,tit=i
    wait,0.3
    ;stop
  endfor

  stop
  
  print,'Writing to ',outfile
  head = flat.head
  sxaddpar,head,'exptype','Flat'
  flat = FIRE_UPDATEHEADER(flat,head)  
  FIRE_WRITEIMAGE,flat,outfile

  stop

end
