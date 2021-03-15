function fire_scale_trace,tstr,im,norescale=norescale,pl=pl,psfile=psfile

  ;; Recenter and scale the aperture
  ;; tstr - trace structure for this order
  ;; im  - original image

  setdisp
  !p.font = 0
  
  npix = 2048
  add_tag,tstr,'recenter',0.0,tstr
  add_tag,tstr,'rescale',1.0,tstr
  yrecenter = 0.0
  rescale = 1.0
  
  ;; Rectify
  recim = FIRE_RECTIFY_ORDER(tstr,im,/exact)
  flux = recim.flux
  flux = medfilt2d(recim.flux,11,dim=1,/edge_copy)  
  gd = where(recim.mask eq 1,ngd,comp=bd,ncomp=nbd)
  flux -= median(flux[gd])
  if nbd gt 0 then flux[bd] = 0.0
  ;;psf = total(recim.flux>0,1)/recim.nx
  ;;psf = median(recim.flux>0,dim=1)
  ;;psf1 -= median(psf1)
  gdpsf = where(total(recim.mask,1) gt 0.95*recim.nx,ngdpsf,comp=bdpsf,ncomp=nbdpsf)
  x = findgen(recim.ny)
  psf = median(flux,dim=1)
  psf -= median(psf[gdpsf])
  xin = x[min(gdpsf):max(gdpsf)]
  psfin = psf[min(gdpsf):max(gdpsf)]
  nyhalf = recim.ny/2
  bestind = xin[first_el(maxloc(psfin))]

  if keyword_set(norescale) then begin
    xin = lindgen(recim.ny)
    psfin = psf
    goto,plotting
  endif
     
  ;; Fit gaussian to PSF
  estimates = [psf[bestind],bestind,2.0,0.0]
  parinfo = replicate({limited:[0,0],limits:[0.0,0.0],fixed:0},4)
  parinfo[0].limited[0] = 1  & parinfo[0].limits[0] = 0.0
  parinfo[1].limited = 1     & parinfo[1].limits = [-5,5]+estimates[1]
  parinfo[2].limited = 1     & parinfo[2].limits = [1,10]
  parinfo[3].limited = 1     & parinfo[3].limits = [-10,10]
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

  PLOTTING:

  green = fsc_color('forest green',1)
  if keyword_set(pl) then begin
    plot,xin,psfin,tit='Order '+strtrim(tstr.order,2)
    if not keyword_set(norescale) then oplot,xin,yfit,co=250
    oplot,[0,0]+nyhalf,[-1e5,1e5],co=80
    oplot,[0,0]+nyhalf+yrecenter,[-1e5,1e5],co=green
  endif

  if n_elements(psfile) gt 0 then begin
    ps_open,psfile,/color,thick=4,/encap
    device,/inches,xsize=8.5,ysize=6.5
    plot,xin,psfin,tit='Order '+strtrim(tstr.order,2),xtit='Y',ytit='Median flux'
    if not keyword_set(norescale) then oplot,xin,yfit,co=250
    oplot,[0,0]+nyhalf,[-1e5,1e5],co=80
    oplot,[0,0]+nyhalf+yrecenter,[-1e5,1e5],co=green
    legend_old,['Recenter = '+stringize(yrecenter,ndec=2)],textcolor=green,/top,/left,box=0
    ps_close
    ps2png,psfile+'.eps',/eps
    if file_test(psfile+'.pdf') then file_delete,psfile+'.pdf'
    spawn,['epstopdf',psfile+'.eps'],/noshell
  endif
  
  return, tstr

end
