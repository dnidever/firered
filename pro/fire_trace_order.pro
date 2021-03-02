pro fire_trace_order,bstr,im,tstr,coeffs

  ;; Trace a stellar aperture in one order

  sz = size(im.flux)
  nx = sz[1]
  ny = sz[2]
  x = findgen(nx)
  y = findgen(ny)
  xx = x#replicate(1,ny)
  yy = replicate(1,nx)#y

  ;; Subimage
  apim = fire_getorderimage(bstr,im)
  x0 = apim.x[0]
  
  ;x1 = x[bstr.lo:bstr.hi]
  ylo = poly(apim.x,bstr.y0coef) > 0
  yhi = poly(apim.x,bstr.y1coef) < (ny-1)
  bndy0 = round(min(ylo)) > 0
  ;bndy1 = round(max(yhi)) < (ny-1)
  ;subim = im[bstr.lo:bstr.hi,bndy0:bndy1]
  ;suberr = err[bstr.lo:bstr.hi,bndy0:bndy1]  
  ;subxx = xx[bstr.lo:bstr.hi,bndy0:bndy1]
  ;subyy = yy[bstr.lo:bstr.hi,bndy0:bndy1]    
  ;submask = (subyy ge (poly(subxx,bstr.y0coef)>0) and $
  ;           subyy le (poly(subxx,bstr.y1coef)<(ny-1)))
  ;subsz = size(subim)
  ;subnx = subsz[1]
  ;subny = subsz[2]
  
  ;apim = subim*submask  ;; aperture image
  ;aperr = suberr*submask
  
  ;; Loop over columns and find peaks
  step = 20
  nbin = 40
  tstr = replicate({num:0,xmed:0.0,gpars:fltarr(4),gperror:fltarr(4),gstatus:0,gflux:0.0,$
                    mpars:fltarr(5),mperror:fltarr(5),mstatus:0,mflux:0.0,$
                    bad:0},apim.nx/step)
  
  For i=0,apim.nx/step-1 do begin
    xmed = i*step+step/2
    xlo = xmed-nbin/2 > 0
    xhi = xmed+nbin/2 < (apim.nx-1)
    ylo1 = ceil(ylo[xmed])+3 - bndy0
    yhi1 = floor(yhi[xmed])-3 - bndy0
    med = median(apim.flux[xlo:xhi,ylo1:yhi1],dim=1)
    sm = smooth(med,3,/edge_truncate)
    mederr = sqrt( total(apim.err[xlo:xhi,ylo1:yhi1]^2,1)/(yhi1-ylo1+1) )
    ysm = apim.y[ylo1:yhi1]
    ;ysm = reform(subyy[xmed,ylo1:yhi1])
    
    ;; Find maximum
    maxind = first_el(maxloc(sm))
    ;; Fit Gaussian
    gestimates = [max(sm)-median(sm),ysm[maxind],2.0,median(sm)]
    gparinfo = replicate({limited:[0,0],limits:[0.0,0.0],fixed:0},4)
    gparinfo[0].limited[0] = 1  & gparinfo[0].limits[0] = 0.0
    gparinfo[1].limited = 1     & gparinfo[1].limits = [0.8,1.2]*gestimates[1]
    gparinfo[2].limited = 1     & gparinfo[2].limits = [2,10]
    gparinfo[3].limited = 1     & gparinfo[3].limits = [0.5,1.5]*gestimates[3]
    gpars = mpfitfun('gaussian',ysm,sm,mederr,gestimates,parinfo=gparinfo,$
                    perror=gperror,yfit=gyfit,status=gstatus,/quiet)
    tstr[i].num = i+1
    tstr[i].xmed = xmed+x0
    tstr[i].gstatus = gstatus

    if gstatus gt 0 then begin
      tstr[i].gpars = gpars
      tstr[i].gperror = gperror
      tstr[i].gflux = gpars[0]*gpars[1]*sqrt(2*!dpi)
      ;; bad solution at edges
      if gpars[0] lt 2*median(mederr) then tstr[i].bad = 1
    endif
    
    mestimates = [gpars[0:2],3.0,gpars[3]]
    mparinfo = replicate({limited:[0,0],limits:[0.0,0.0],fixed:0},5)
    mparinfo[0].limited[0] = 1  & mparinfo[0].limits[0] = 0.0
    mparinfo[1].limited = 1     & mparinfo[1].limits = [0.8,1.2]*mestimates[1]
    mparinfo[2].limited = 1     & mparinfo[2].limits = [2,20]
    mparinfo[3].limited = 1     & mparinfo[3].limits = [0,15]
    mparinfo[4].limited = 1     & mparinfo[4].limits = minmax(sm)
    
    mpars = mpfitfun('fire_moffat',ysm,sm,mederr,mestimates,parinfo=mparinfo,$
                     perror=mperror,yfit=myfit,status=mstatus,/quiet)
    tstr[i].mpars = mpars
    tstr[i].mstatus = mstatus
    
    if mstatus gt 0 then begin
      tstr[i].mpars = mpars
      tstr[i].mperror = mperror
      tstr[i].mflux = mpars[0]*mpars[1]*sqrt(2*!dpi)
      ;; bad solution at edges
      ;if mpars[0] lt 2*median(mederr) then tstr[i].bad = 1
    endif
  Endfor

  ;; Fit polynomial to Y position
  ;gd = where(tstr.mstatus gt 0 and tstr.bad eq 0 and tstr.mperror[1] lt 1 and tstr.mperror[1] gt 0,ngd)
  gd = where(tstr.mstatus gt 0 and tstr.mperror[1] lt 1 and tstr.mperror[1] gt 0,ngd)  
  tcoef = robust_poly_fit(tstr[gd].xmed,tstr[gd].mpars[1],2)

  ;; Fit Gaussian sigma and Moffat values
  sigcoef = robust_poly_fit(tstr[gd].xmed,tstr[gd].gpars[2],2)
  hwhmcoef = robust_poly_fit(tstr[gd].xmed,tstr[gd].mpars[2],2)  
  moffcoef = robust_poly_fit(tstr[gd].xmed,tstr[gd].mpars[3],2)  
  ;; the fit to the moffat values are not great at the edges, but good enough

  ;; Put coefficients together
  coeffs = {tcoef:tcoef,sigcoef:sigcoef,hwhmcoef:hwhmcoef,moffcoef:moffcoef}
  
  ;plot,tstr[gd].xmed,tstr[gd].gpars[2],ps=1
  ;oplot,tstr.xmed,poly(tstr.xmed,sigcoef),co=60

  ;plot,tstr[gd].xmed,tstr[gd].mpars[2],ps=1
  ;oplot,tstr.xmed,poly(tstr.xmed,hwhmcoef),co=250
  ;plot,tstr[gd].xmed,tstr[gd].mpars[3],ps=1
  ;oplot,tstr.xmed,poly(tstr.xmed,moffcoef),co=250  
  
  ;; X values are LOCAL and Y values are GLOBAL
  ;displayc,apim,findgen(subnx),reform(subyy[0,*]),/z
  ;oplot,tstr.xmed,tstr.pars[1],ps=1
  ;oplot,tstr.xmed,poly(tstr.xmed,tcoef),co=60
  
  ;; Should the trace fit use X-values WITHIN the subregion or the
  ;; FULL array????
  
  ;stop

 
end
