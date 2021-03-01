function fire_extract_order,tstr,im,arc=arc,recenter=recenter,yrecenter=yrecenter

  ;; maybe use fire_rectify_order.pro here

  npix = 2048

  if n_elements(yrecenter) eq 0 then yrecenter=0.0
  
  ;; Get aperture subimage using boundary
  apim = fire_getorderimage(tstr,im)

  ;;; Get trace information
  ;ytrace = poly(apim.x,tstr.tycoef)-apim.y[0]
  
  x = findgen(npix)
  xx = findgen(npix)#replicate(1,npix)
  yy = replicate(1,npix)#findgen(npix)

  ;x1 = x[tstr.bndx0:tstr.bndx1]
  ylo = poly(apim.x,tstr.bndy0coef) > 0
  yhi = poly(apim.x,tstr.bndy1coef) < (npix-1)
  y0 = round(min(ylo)) > 0
  ;y1 = round(max(yhi)) < (npix-1)
  ;subim = im[tstr.bndx0:tstr.bndx1,y0:y1]
  ;suberr = err[tstr.bndx0:tstr.bndx1,y0:y1]  
  ;subxx = xx[tstr.bndx0:tstr.bndx1,y0:y1]
  ;subyy = yy[tstr.bndx0:tstr.bndx1,y0:y1]    
  ;submask = (subyy ge (poly(subxx,tstr.bndy0coef)>0) and $
  ;           subyy le (poly(subxx,tstr.bndy1coef)<(npix-1)))
  ;subsz = size(subim)
  ;subnx = subsz[1]
  ;subny = subsz[2]
  ;ysub = findgen(subny)

  ;;; Masked "aperture" image
  ;apim = subim*submask
  ;aperr = suberr*submask
  
  ;; Trace and sigma information
  ;xsub = findgen(subnx)
  ytrace = poly(apim.x,tstr.tycoef)-y0 + yrecenter
  sigtrace = poly(apim.x,tstr.tsigcoef)

  ;; Median for each column
  ;temp = apim*0 + !values.f_nan
  ;temp[where(submask eq 1)] = apim[where(submask eq 1)]
  ;medy = median(temp,dim=1)
  ;bd = where(finite(medy) eq 0,nbd,comp=gd)
  ;medy[bd] = median(medy[gd])
  
  ;; Extract each column
  extstr = replicate({num:0L,pars:fltarr(4),perror:fltarr(4),flux:0.0,err:0.0,status:0,boxflux:0.0},apim.nx)
  extstr.num = lindgen(apim.nx)+1
  for i=0,apim.nx-1 do begin
    parinfo = replicate({limited:[0,0],limits:[0.0,0.0],fixed:0},4)
    parinfo[0].limited[0] = 1  & parinfo[0].limits[0] = 0.0
    parinfo[1].fixed = 1
    parinfo[2].fixed = 1
    parinfo[3].limited = 1     & parinfo[3].limits = [-100,500]
    estimates = [apim.flux[i,round(ytrace[i])]>1, ytrace[i], sigtrace[i], 0.0]
    ysrt = round(ylo[i])-y0 + 3
    yend = round(yhi[i])-y0 - 3

    pars = mpfitfun('gaussian',apim.y[ysrt:yend],reform(apim.flux[i,ysrt:yend]),reform(apim.err[i,ysrt:yend])>1,estimates,$
                     parinfo=parinfo,perror=perror,status=status,/quiet)
    extstr[i].status = status
    if status gt 0 then begin
      extstr[i].pars = pars
      extstr[i].perror = perror
      extstr[i].flux = pars[0]*pars[2]*sqrt(2*!dpi)
      extstr[i].err = extstr[i].flux * sqrt( (perror[0]/pars[0])^2 + (perror[2]/pars[2])^2 )
    endif
  endfor

  return, extstr

end
