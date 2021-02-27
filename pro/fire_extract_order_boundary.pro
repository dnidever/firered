function fire_extract_order_boundary,tstr,im,err,arc=arc

  ;; maybe use fire_rectify_order.pro here
  
  npix = 2048
  x = findgen(npix)
  xx = findgen(npix)#replicate(1,npix)
  yy = replicate(1,npix)#findgen(npix)

  x1 = x[tstr.lo:tstr.hi]
  ylo = poly(x1,tstr.y0coef) > 0
  yhi = poly(x1,tstr.y1coef) < (npix-1)
  y0 = round(min(ylo)) > 0
  y1 = round(max(yhi)) < (npix-1)
  im2 = im[tstr.lo:tstr.hi,y0:y1]
  err2 = err[tstr.lo:tstr.hi,y0:y1]  
  xx2 = xx[tstr.lo:tstr.hi,y0:y1]
  yy2 = yy[tstr.lo:tstr.hi,y0:y1]    
  mask = (yy2 ge (poly(xx2,tstr.y0coef)>0) and $
          yy2 le (poly(xx2,tstr.y1coef)<(npix-1)))

  ;; rectify
  nx = tstr.hi-tstr.lo+1
  ny = max((round(yhi)<(npix-1)) - (round(ylo)>0))+1
  recim2 = fltarr(nx,ny)
  recerr2 = fltarr(nx,ny)+999999.  

  for j=0,nx-1 do begin
    ysrt = round(ylo[j])-y0
    yend = round(yhi[j])-y0
    recim2[j,0] = (im2*mask)[j,ysrt:yend]
    recerr2[j,0] = (err2*mask)[j,ysrt:yend]    
  endfor
  bd = where(recerr2 eq 0.0,nbd)
  if nbd gt 0 then recerr2[bd] = 999999.

  ;; Remove some background
  ;med2 = median(recimobj2,5) ;; median smooth
  ;medy = median(med2,dim=2)

  ;; Normal object spectrum
  if not keyword_set(arc) then begin
  
    ;; collapse in wavelength to get PSF
    psf = total(recim2>0,1)/nx

    ;; Fit gaussian to PSF
    estimates = [max(psf),ny/2,2.0,median(psf)]
    parinfo = replicate({limited:[0,0],limits:[0.0,0.0],fixed:0},4)
    parinfo[0].limited[0] = 1  & parinfo[0].limits[0] = 0.0
    parinfo[1].limited = 1     & parinfo[1].limits = [0.8,1.2]*estimates[1]
    parinfo[2].limited = 1     & parinfo[2].limits = [0.8,1.2]*estimates[2]
    parinfo[3].limited = 1     & parinfo[3].limits = [0.5,1.5]*estimates[3]
    psfpars = mpfitfun('gaussian',findgen(ny),psf,psf*0+1,estimates,$
                       parinfo=parinfo,perror=perror1,status=status,/quiet)
    ;yfit = mpfitpeak(findgen(ny),psf,psfpars,nterms=4,estimates=estimates,/gaussian)
    y0 = round(psfpars[1]-3*psfpars[2]) > 0
    y1 = round(psfpars[1]+3*psfpars[2]) < (ny-1)

  endif else begin
    y0 = ny/2-5
    y1 = ny/2+5
    psfpars = [1.0, ny/2, 2.0, 10.0]
  endelse
     
  y = findgen(ny)

    
  ;; Extract each column
  extstr = replicate({num:0L,pars:fltarr(4),perror:fltarr(4),flux:0.0,status:0,boxflux:0.0},nx)
  extstr.num = lindgen(nx)+1
  for i=0,nx-1 do begin
    parinfo = replicate({limited:[0,0],limits:[0.0,0.0],fixed:0},4)
    parinfo[0].limited[0] = 1  & parinfo[0].limits[0] = 0.0
    parinfo[1].limited = 1     & parinfo[1].limits = [0.8,1.2]*psfpars[1]
    parinfo[2].limited = 1     & parinfo[2].limits = [0.8,1.2]*psfpars[2]
    parinfo[3].limited = 1     & parinfo[3].limits = [0.5,1.5]*psfpars[3]
    pars1 = mpfitfun('gaussian',y[y0:y1],reform(recim2[i,y0:y1]),reform(recerr2[i,y0:y1]),psfpars,$
                     parinfo=parinfo,perror=perror1,status=status,/quiet)
    if status gt 0 then begin
      extstr[i].pars = pars1
      extstr[i].perror = perror1    
      extstr[i].flux = pars1[0]*pars1[2]*sqrt(2*!dpi)
    endif
    extstr[i].status = status
  endfor

  ;; also do boxcar
  tot = total(recim2[*,y0:y1],2)  ; subtract median!!!
  extstr.boxflux = tot
  ;; use upper/lower background regions

  ;; can I just shift the pixels in x-direction to rectify things
  ;; either shift by integer pixel values or interpolate

  ;; shifts 12 pixels in x in 55 pixels in y, 0.2 slope
  ;; we could leave the observed spectrum as is, but fit the skylines
  ;; as a sloping line
  
  ;stop
  
  ;; ARC, find the lines and fit them
  if keyword_set(arc) then begin
    mid = median(recim2[*,ny/2-2:ny/2+1],dim=2)
    maxima,mid,minarr,maxarr
    sig = mad(mid)
    gd = where(mid[maxarr] gt 10*sig and mid[maxarr-1] gt 5*sig and mid[maxarr+1] gt 5*sig,ngd)
    maxarr = maxarr[gd]
    nlines = ngd
    ;; Fit the lines
    xmid = findgen(n_elements(mid))
    if nlines gt 0 then begin
      linestr = replicate({num:0L,pars:fltarr(4)},nlines)
      linestr.num = lindgen(nlines)+1
      for j=0,nlines-1 do begin
        estimates = [mid[maxarr[j]],maxarr[j],2.0,0.0]
        x0 = (maxarr[j]-10) > 0
        x1 = (maxarr[j]+10) < (n_elements(mid)-1)
        yfit1 = mpfitpeak(xmid[x0:x1],mid[x0:x1],pars1,nterms=4,/gaussian)
        linestr[j].pars = pars1
        ;; Fit these in Y
      endfor
    endif else begin
      linestr = -1
    endelse
  endif
  
  
  ;; output structure
  if not keyword_set(arc) then begin
    outstr = {x0:tstr.lo,x1:tstr.hi,y0:y0,y1:y1,im:im2,mask:mask,recim:recim2,recerr:recerr2,extstr:extstr,psf:psf,psfpars:psfpars}
  endif else begin
    outstr = {x0:tstr.lo,x1:tstr.hi,y0:y0,y1:y1,im:im2,mask:mask,recim:recim2,recerr:recerr2,extstr:extstr,linestr:linestr}
  endelse
    
  return, outstr

end
