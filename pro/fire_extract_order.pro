pro fire_extract_order,tstr,im,extstr,model,arc=arc,recenter=recenter,yrecenter=yrecenter

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
  ;ysub = findgen(apim.ny)

  ;;; Masked "aperture" image
  ;apim = subim*submask
  ;aperr = suberr*submask
  
  ;; Trace and sigma information
  ;xsub = findgen(subnx)
  ytrace = poly(apim.x,tstr.tycoef) + yrecenter
  sigtrace = poly(apim.x,tstr.tsigcoef)

  ;; Median for each column
  ;temp = apim*0 + !values.f_nan
  ;temp[where(submask eq 1)] = apim[where(submask eq 1)]
  ;medy = median(temp,dim=1)
  ;bd = where(finite(medy) eq 0,nbd,comp=gd)
  ;medy[bd] = median(medy[gd])

  ;; Rectified image for boxcar flux
  recim = fire_rectify_order(tstr,im,/exact)
  
  ;; Extract each column
  extstr = replicate({num:0L,x:0.0,pars:fltarr(4),perror:fltarr(4),chisq:0.0,rchisq:0.0,$
                      flux:0.0,err:0.0,mask:0,status:0,boxflux:0.0},apim.nx)
  extstr.num = lindgen(apim.nx)+1
  extstr.x = apim.x
  model = apim.flux*0
  for i=0,apim.nx-1 do begin
    parinfo = replicate({limited:[0,0],limits:[0.0,0.0],fixed:0},4)
    parinfo[0].limited[0] = 1  & parinfo[0].limits[0] = 0.0
    parinfo[1].fixed = 1
    parinfo[2].fixed = 1
    parinfo[3].limited = 1     & parinfo[3].limits = [-100,500]
    estimates = [apim.flux[i,round(ytrace[i]-y0)]>1, ytrace[i], sigtrace[i], 0.0]
    ;ysrt = round(ytrace[i]-y0-5*sigtrace[i])
    ;yend = round(ytrace[i]-y0+5*sigtrace[i])    
    ysrt = round(ylo[i])-y0 + 3
    yend = round(yhi[i])-y0 - 3

    pars = mpfitfun('gaussian',apim.y[ysrt:yend],reform(apim.flux[i,ysrt:yend]),reform(apim.err[i,ysrt:yend])>1,$
                    estimates,parinfo=parinfo,perror=perror,status=status,yfit=yfit,bestnorm=chisq,/quiet)
    rchisq = chisq/n_elements(apim.flux[i,ysrt:yend])
    extstr[i].chisq = chisq
    extstr[i].rchisq = rchisq
    flux = pars[0]*pars[2]*sqrt(2*!dpi)
    ;; Boxcar flux
    recflux = reform(recim.flux[i,*])
    recmask = reform(recim.mask[i,*])
    yblo = round(recim.ny/2-3*sigtrace[i])
    ybhi = round(recim.ny/2+3*sigtrace[i])
    back = [recflux[0:yblo-1],recflux[ybhi+1:*]]
    backmask = [recmask[0:yblo-1],recmask[ybhi+1:*]]
    gdback = where(backmask eq 1,ngdback)
    medback = 0.0
    if ngdback gt 0 then medback=median(back[gdback])
    boxflux = total((recflux[yblo:ybhi]-medback)>0)
    extstr[i].boxflux = boxflux
    
    ;; If the Gaussian fit flux is much lower than the boxcar flux
    ;; then refit with outlier rejection
    if status gt 0 and flux lt boxflux*0.8 then begin
      pars1 = pars
      perror1 = perror
      yfit1 = yfit
      estimates2 = pars1
      yin = apim.y[ysrt:yend]
      fluxin = reform(apim.flux[i,ysrt:yend])
      errin = reform(apim.err[i,ysrt:yend])>1
      sig = mad((fluxin-yfit1)/errin,/zero)
      gd = where(abs((fluxin-yfit1)/errin) lt 5*sig,ngd,comp=bd,ncomp=nbd)
      if nbd gt 0 then begin
        pars = mpfitfun('gaussian',yin[gd],fluxin[gd],errin[gd],estimates2,$
                        parinfo=parinfo,perror=perror,status=status,yfit=yfit,bestnorm=chisq,/quiet)
        rchisq = chisq/ngd
        extstr[i].chisq = chisq
        extstr[i].rchisq = rchisq
        flux = pars[0]*pars[2]*sqrt(2*!dpi)
        ;;plot,fluxin
        ;;oplot,yfit1,co=150
        ;;oplot,gd,yfit,co=250
      endif
    endif
    
    extstr[i].status = status
    if status gt 0 then begin
      extstr[i].pars = pars
      extstr[i].perror = perror
      extstr[i].flux = pars[0]*pars[2]*sqrt(2*!dpi)
      extstr[i].err = extstr[i].flux * sqrt( (perror[0]/pars[0])^2 + (perror[2]/pars[2])^2 )
      extstr[i].mask = 1
      ;; add the model
      model[i,*] = gaussian(apim.y,[pars[0],pars[1],pars[2],0.0])
    endif else begin
      extstr[i].flux = 1e30
      extstr[i].err = 1e30
      extstr[i].mask = 0      
    endelse
  endfor

end
