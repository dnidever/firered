function fire_recenter,tstr,im

  ;; Recenter the aperture
  ;; tstr - trace structure for this order
  ;; im  - original image

  npix = 2048

  ;; Rectify
  recim = FIRE_RECTIFY_ORDER(tstr,im,/exact)
  psf = total(recim.flux>0,1)/recim.nx
  psf -= median(psf)
  nyhalf = recim.ny/2
  bestind = first_el(maxloc(psf))
  
  ;; Fit gaussian to PSF
  estimates = [psf[bestind],bestind,2.0,0.0]
  parinfo = replicate({limited:[0,0],limits:[0.0,0.0],fixed:0},4)
  parinfo[0].limited[0] = 1  & parinfo[0].limits[0] = 0.0
  parinfo[1].limited = 1     & parinfo[1].limits = [-2,2]+estimates[1]
  parinfo[2].limited = 1     & parinfo[2].limits = [1,10]
  parinfo[3].limited = 1     & parinfo[3].limits = [-100,100]
  psfpars = mpfitfun('gaussian',findgen(recim.ny),psf,psf*0+1,estimates,$
                     parinfo=parinfo,perror=psfperror,yfit=yfit,status=status,/quiet)
  yrecenter = psfpars[1]-nyhalf

  ;; Recenter aperture 
  tstr.tycoef[0] += yrecenter
  print,'Recenter = ',strtrim(yrecenter,2)
  

  ;; Rescale sigma based on this fit
  sig = poly(recim.x,tstr.tsigcoef)
  medsig = median(sig)
  scale = psfpars[2]/medsig
;  tstr.tsigcoef = XXX
  print,'Rescale = ',strtrim(scale,2)
  
  ;; MAYBE RESCALE/RECENTER THE TRACE STRUCTURE AND RETURN IT!!!
stop
  
  return, tstr

end
