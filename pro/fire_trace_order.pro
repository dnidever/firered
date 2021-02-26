function fire_trace_order,bstr,im,err

  ;; Trace a stellar aperture in one order

  sz = size(im)
  nx = sz[1]
  ny = sz[2]
  x = findgen(nx)
  y = findgen(ny)
  xx = x#replicate(1,ny)
  yy = replicate(1,nx)#y

  ;; Subimage
  x1 = x[bstr.lo:bstr.hi]
  ylo = poly(x1,bstr.y0coef) > 0
  yhi = poly(x1,bstr.y1coef) < (ny-1)
  bndy0 = round(min(ylo)) > 0
  bndy1 = round(max(yhi)) < (ny-1)
  subim = im[bstr.lo:bstr.hi,bndy0:bndy1]
  suberr = err[bstr.lo:bstr.hi,bndy0:bndy1]  
  subxx = xx[bstr.lo:bstr.hi,bndy0:bndy1]
  subyy = yy[bstr.lo:bstr.hi,bndy0:bndy1]    
  submask = (subyy ge (poly(subxx,bstr.y0coef)>0) and $
             subyy le (poly(subxx,bstr.y1coef)<(ny-1)))
  subsz = size(subim)
  subnx = subsz[1]
  subny = subsz[2]
  
  apim = subim*submask  ;; aperture image
  aperr = suberr*submask
  
  ;; Loop over columns and find peaks
  step = 20
  nbin = 40
  tstr = replicate({num:0,xmed:0.0,pars:fltarr(4),perror:fltarr(4),flux:0.0,status:0},subnx/step)
  For i=0,subnx/step-1 do begin
    xmed = i*step+step/2
    xlo = xmed-nbin/2 > 0
    xhi = xmed+nbin/2 < (subnx-1)
    ylo1 = ceil(ylo[xmed])+3 - bndy0
    yhi1 = floor(yhi[xmed])-3 - bndy0
    med = median(apim[xlo:xhi,ylo1:yhi1],dim=1)
    sm = smooth(med,3,/edge_truncate)
    mederr = sqrt( total(aperr[xlo:xhi,ylo1:yhi1]^2,1)/(yhi1-ylo1+1) )
    ysm = reform(subyy[xmed,ylo1:yhi1])
    
    ;; Find maximum
    maxind = first_el(maxloc(sm))
    ;; Fit Gaussian
    estimates = [max(sm)-median(sm),ysm[maxind],2.0,median(sm)]
    parinfo = replicate({limited:[0,0],limits:[0.0,0.0],fixed:0},4)
    parinfo[0].limited[0] = 1  & parinfo[0].limits[0] = 0.0
    parinfo[1].limited = 1     & parinfo[1].limits = [0.8,1.2]*estimates[1]
    parinfo[2].limited = 1     & parinfo[2].limits = [2,5]
    parinfo[3].limited = 1     & parinfo[3].limits = [0.5,1.5]*estimates[3]
    pars = mpfitfun('gaussian',ysm,sm,mederr,estimates,parinfo=parinfo,$
                    perror=perror,yfit=yfit,status=status,/quiet)
    tstr[i].num = i+1
    tstr[i].xmed = xmed
    tstr[i].status = status
    if status gt 0 then begin
      tstr[i].pars = pars
      tstr[i].perror = perror
      tstr[i].flux = pars[0]*pars[1]*sqrt(2*!dpi)
    endif
  Endfor

  ;; Fit polynomial to order
  tcoef = robust_poly_fit(tstr.xmed,tstr.pars[1],3)

  ;; Should the trace fit use X-values WITHIN the subregion or the
  ;; FULL array????
  
  
  stop

 
end
