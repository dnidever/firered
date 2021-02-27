pro fire_linefit2d_boundary,im,err,linestr,model,arc=arc

  ;; Fit emissions lines in 2D
  ;; input rectified image

  ;; Find the peaks
  sz = size(im)
  nx = sz[1]
  ny = sz[2]
  x = findgen(nx)
  y = findgen(ny)

  ;; Find peaks at lower and upper half
  med1 = median(im[*,5:ny/2-5],dim=2)
  sm1 = medfilt1d(med1,101)
  med1 -= sm1
  maxima,med1,minarr1,maxarr1
  sig1 = mad(med1)
  gd1 = where(med1[maxarr1] gt 5*sig1 and med1[maxarr1-1] gt 3*sig1 and med1[maxarr1+1] gt 3*sig1,ngd1)
  maxarr1 = maxarr1[gd1]

  med2 = median(im[*,ny/2+5:ny-6],dim=2)  
  sm2 = medfilt1d(med2,101)
  med2 -= sm2
  maxima,med2,minarr2,maxarr2
  sig2 = mad(med2)
  gd2 = where(med2[maxarr2] gt 5*sig2 and med2[maxarr2-1] gt 3*sig2 and med2[maxarr2+1] gt 3*sig2,ngd2)
  maxarr2 = maxarr2[gd2]  

  ;; Find the shift between them
  xcorlb,med1,med2,15,xsh

  ;; Match up the lines, need to find them in both
  srcmatch,maxarr1,maxarr1*0,maxarr2+xsh,maxarr2*0,3.0,ind1,ind2
  maxarr1 = maxarr1[ind1]
  maxarr2 = maxarr2[ind2]
  nlines = n_elements(maxarr1)
  
  ;; Get initial estimate for slope
  slp = -xsh/(mean([ny/2+5,ny-6])-mean([5,ny/2-5]))

  ;; Fit the lines
  residim = im
  model = im*0
  if nlines gt 0 then begin
    linestr = replicate({num:0L,pars:fltarr(6),status:0},nlines)
    linestr.num = lindgen(nlines)+1
    for j=0,nlines-1 do begin
      ;; Get the subimage
      xlo = min([maxarr1[j],maxarr2[j]])-12 > 0
      xhi = max([maxarr1[j],maxarr2[j]])+12 < (nx-1)
      subim = residim[xlo:xhi,*]
      suberr = err[xlo:xhi,*]
      subsz = size(subim)
      nxsub = subsz[1]
      nysub = subsz[2]
      ;; parameters: Gaussian height, Gaussian width, central X position,
      ;;             central Y position, slope (deltaX/deltaY), offset
      ;; X and Y should be 2D arrays
      xx = x[xlo:xhi]#replicate(1,ny)
      yy = replicate(1,nxsub)#y

      estimates = [mean([med1[maxarr1[j]], med2[maxarr2[j]]]), 2.0, mean(x[xlo:xhi]), ny/2, slp, 0.0]      
      parinfo = replicate({limited:[0,0],limits:[0.0,0.0],fixed:0},6)
      parinfo[0].limited[0] = 1  & parinfo[0].limits[0] = 0.0                   ;; height
      parinfo[1].limited = 1     & parinfo[1].limits = [1,6]                    ;; sigma
      parinfo[2].fixed = 1                                                      ;; cenX
      parinfo[3].limited = 1     & parinfo[3].limits = [-2,2]+estimates[3]      ;; cenY
      parinfo[4].limited = 1     & parinfo[4].limits = [0.7,1.3]*estimates[4]   ;; slope
      parinfo[5].limited = 1     & parinfo[5].limits = minmax(subim)            ;; offset
      ;; Mask central region
      if not keyword_set(arc) then begin
        mask = subim*0+1
        mask[*,ny/2-5:ny/2+5] = 0
        subim = subim*mask
      endif else mask=subim*0+1
      ;; Mask edges
      fa = {y0:3,y1:nysub-3,mask:mask}
      subim[*,0:3] = 0
      subim[*,nysub-3:*] = 0
      pars1 = mpfit2dfun('fire_gauss2d',xx,yy,subim,suberr,estimates,parinfo=parinfo,status=status,functargs=fa,/quiet)
      linestr[j].status = status
      if status gt 0 then begin
        linestr[j].pars = pars1
        gmodel1 = fire_gauss2d(xx,yy,pars1,_extra=fa)  ;; model with offset and trimming/masking
        tpars1 = pars1    ; no offset
        tpars1[5] = 0.0 
        gmodel = fire_gauss2d(xx,yy,tpars1)
        model[xlo:xhi,*] += gmodel
        ;; subtract
        residim[xlo:xhi,*] -= gmodel
      endif
    endfor

  ;; no lines
  endif else begin
    linestr = -1
  endelse

  end
