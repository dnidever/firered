function fire_scale_trace,tstr,im,pl=pl

  ;; Recenter and scale the aperture
  ;; tstr - trace structure for this order
  ;; im  - original image

  npix = 2048
  add_tag,tstr,'recenter',0.0,tstr
  add_tag,tstr,'rescale',0.0,tstr
  
  ;; Rectify
  recim = FIRE_RECTIFY_ORDER(tstr,im,/exact)
  flux = recim.flux
  flux -= median(flux)
  flux = medfilt2d(recim.flux,11,dim=1,/edge_copy)
  ;psf = total(recim.flux>0,1)/recim.nx
  ;psf = median(recim.flux>0,dim=1)
  ;psf1 -= median(psf1)
  psf = median(flux,dim=1)
  psf -= median(psf)
  nyhalf = recim.ny/2
  bestind = first_el(maxloc(psf))
  
  ;; Fit gaussian to PSF
  estimates = [psf[bestind],bestind,2.0,0.0]
  parinfo = replicate({limited:[0,0],limits:[0.0,0.0],fixed:0},4)
  parinfo[0].limited[0] = 1  & parinfo[0].limits[0] = 0.0
  parinfo[1].limited = 1     & parinfo[1].limits = [-5,5]+estimates[1]
  parinfo[2].limited = 1     & parinfo[2].limits = [1,10]
  parinfo[3].limited = 1     & parinfo[3].limits = [-10,10]
  x = lindgen(recim.ny)
  totmask = total(recim.mask,1)
  gd = where(totmask eq recim.nx,ngd)
  xin = x[min(gd):max(gd)]  ; mask first few and last lines
  psfin = psf[min(gd):max(gd)]
  errin = psfin*0+1
  psfpars = mpfitfun('gaussian',xin,psfin,errin,estimates,$
                     parinfo=parinfo,perror=psfperror,yfit=yfit,status=status,/quiet)
  yrecenter = psfpars[1]-nyhalf
  
  ;; Recenter aperture 
  tstr.tycoef[0] += yrecenter
  print,'Recenter = ',strtrim(yrecenter,2)
  tstr.recenter = yrecenter
  

  ;; Rescale sigma based on this fit
  sig = poly(recim.x,tstr.tsigcoef)
  medsig = median(sig)
  scale = psfpars[2]/medsig
  tstr.tsigcoef *= scale
  tstr.thwhmcoef *= scale
  tstr.tmoffcoef *= scale  
  print,'Rescale = ',strtrim(scale,2)
  tstr.rescale = scale

  if keyword_set(pl) then begin
    plot,xin,psfin,tit='Order '+strtrim(tstr.order,2)
    oplot,xin,yfit,co=250
    oplot,[0,0]+nyhalf,[-1e5,1e5],co=80
    oplot,[0,0]+nyhalf+yrecenter,[-1e5,1e5],co=150
    ;stop
  endif
  
  return, tstr

end
