function fire_extract_order,tstr,im,err,arc=arc,recenter=recenter,yrecenter=yrecenter

  ;; maybe use fire_rectify_order.pro here

  npix = 2048

  ;; Recenter by default for object spectra
  if not keyword_set(arc) and n_elements(recenter) eq 0 then recenter=1

  ;; Recenter
  if keyword_set(recenter) then begin
    ;; Rectify
    recim = fire_rectify_order(tstr,im)
    recsz = size(recim) 
    psf = total(recim>0h,1)/recsz[1]
    psf -= median(psf)
    nyhalf = recsz[2]/2
    bestind = first_el(maxloc(psf))
    ;yrecenter = bestind-nyhalf
    ;; Fit gaussian to PSF
    estimates = [psf[bestind],bestind,2.0,0.0]
    parinfo = replicate({limited:[0,0],limits:[0.0,0.0],fixed:0},4)
    parinfo[0].limited[0] = 1  & parinfo[0].limits[0] = 0.0
    parinfo[1].limited = 1     & parinfo[1].limits = [-2,2]+estimates[1]
    parinfo[2].limited = 1     & parinfo[2].limits = [1,10]
    parinfo[3].limited = 1     & parinfo[3].limits = [-100,100]
    psfpars = mpfitfun('gaussian',findgen(recsz[2]),psf,psf*0+1,estimates,$
                       parinfo=parinfo,perror=psfperror,yfit=yfit,status=status,/quiet)
    yrecenter = psfpars[1]-nyhalf


    ;; What if I did something similar to cross-correlation in 2D
    ;; create a 2D PSF (unit height?) and x-correlated with observed image
    ;; then shift by +/-2 pixels (relative to bestind) in 0.1 pix step
    ;;   subtract median/background for each column
    ;;   only use pixels with values above some value (>1?)
    ;; then interpolate to get the best value
    ;; should I do a boxcar first and use that for the PSF heights?
    
    stop
    
    ;;; fit with Gaussian
    ;estimates = [psf[bestind]-median(psf),bestind,2.0,0.0]
    ;yfit1 = mpfitpeak(fingen(n_elements(psf)),psf-median(psf),psfpars,nterms=4,/gaussian,estimates=estimates)
    print,'Recenter = ',strtrim(yrecenter,2)
  endif

  if n_elements(yrecenter) eq 0 then yrecenter=0

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
  suberr = err[tstr.bndx0:tstr.bndx1,y0:y1]  
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
  aperr = suberr*submask
  
  ;; Trace and sigma information
  xsub = findgen(subnx)
  ytrace = poly(xsub,tstr.tycoef)-y0 + yrecenter
  sigtrace = poly(xsub,tstr.tsigcoef)

  ;; Median for each column
  ;temp = apim*0 + !values.f_nan
  ;temp[where(submask eq 1)] = apim[where(submask eq 1)]
  ;medy = median(temp,dim=1)
  ;bd = where(finite(medy) eq 0,nbd,comp=gd)
  ;medy[bd] = median(medy[gd])
  
  ;; Extract each column
  extstr = replicate({num:0L,pars:fltarr(4),perror:fltarr(4),flux:0.0,fluxerr:0.0,status:0,boxflux:0.0},subnx)
  extstr.num = lindgen(subnx)+1
  for i=0,subnx-1 do begin
    parinfo = replicate({limited:[0,0],limits:[0.0,0.0],fixed:0},4)
    parinfo[0].limited[0] = 1  & parinfo[0].limits[0] = 0.0
    parinfo[1].fixed = 1
    parinfo[2].fixed = 1
    parinfo[3].limited = 1     & parinfo[3].limits = [-100,500]
    estimates = [apim[i,round(ytrace[i])]>1, ytrace[i], sigtrace[i], 0.0]
    ysrt = round(ylo[i])-y0 + 3
    yend = round(yhi[i])-y0 - 3

    pars = mpfitfun('gaussian',ysub[ysrt:yend],reform(apim[i,ysrt:yend]),reform(aperr[i,ysrt:yend])>1,estimates,$
                     parinfo=parinfo,perror=perror,status=status,/quiet)
    extstr[i].status = status
    if status gt 0 then begin
      extstr[i].pars = pars
      extstr[i].perror = perror
      extstr[i].flux = pars[0]*pars[2]*sqrt(2*!dpi)
      extstr[i].fluxerr = extstr[i].flux * sqrt( (perror[0]/pars[0])^2 + (perror[2]/pars[2])^2 )
    endif
  endfor

  return, extstr

end
