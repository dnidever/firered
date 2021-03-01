pro fire_linefit2d,tstr,im,err,linestr,model,residim,yrecenter=yrecenter,arc=arc

  ;; Fit emission lines in 2D
  npix = 2048
  yrecenter = 0
  
  ;; Get order image
  outim = fire_getorderimage(tstr,im)
  outerr = fire_getorderimage(tstr,err)
  apim = outim.im
  aperr = outerr.im
  subsz = size(apim)
  subnx = subsz[1]
  subny = subsz[2]  
  xsub = findgen(subnx)
  ysub = (findgen(npix))[outim.y0:outim.y1]

  ;; Rectify to find peaks
  recim = fire_rectify_order(tstr,im)
  recerr = fire_rectify_order(tstr,err)  

  ;; Get trace information
  ytrace = poly(xsub,tstr.tycoef)-outim.y0
  
  ;; Find the peaks
  recsz = size(recim)
  recnx = recsz[1]
  recny = recsz[2]
  xrec = findgen(recnx)
  yrec = findgen(recny)

  ;; Find peaks at lower and upper half
  med1 = median(recim[*,5:recny/2-5],dim=2)
  sm1 = medfilt1d(med1,101)
  med1 -= sm1
  maxima,med1,minarr1,maxarr1
  sig1 = mad(med1)
  gd1 = where(med1[maxarr1] gt 5*sig1 and med1[maxarr1-1] gt 3*sig1 and med1[maxarr1+1] gt 3*sig1,ngd1)
  maxarr1 = maxarr1[gd1]

  med2 = median(recim[*,recny/2+5:recny-6],dim=2)  
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

  ;; No lines found
  if nlines eq 0 then begin
    print,'No lines found'
    linestr = -1
    return
  endif
  
  ;; Get initial estimate for slope
  slp = -xsh/(mean([recny/2+5,recny-6])-mean([5,recny/2-5]))

  
  ;; FIRST, fit the lines in rectified space
  ;;-----------------------------------------
  
  ;; Fit the lines
  residim1 = recim
  model1 = recim*0
  
  linestr1 = replicate({num:0L,pars:fltarr(6),status:0},nlines)
  linestr1.num = lindgen(nlines)+1
  for j=0,nlines-1 do begin
    ;; Get the subimage
    xlo = min([maxarr1[j],maxarr2[j]])-12 > 0
    xhi = max([maxarr1[j],maxarr2[j]])+12 < (recnx-1)
    subim1 = residim1[xlo:xhi,*]
    suberr1 = recerr[xlo:xhi,*]
    subsz1 = size(subim1)
    nxsub1 = subsz1[1]
    nysub1 = subsz1[2]
    ;; parameters: Gaussian height, Gaussian width, central X position,
    ;;             central Y position, slope (deltaX/deltaY), offset
    ;; X and Y should be 2D arrays
    xx = xrec[xlo:xhi]#replicate(1,nysub1)
    yy = replicate(1,nxsub1)#yrec

    ;; fix zero error
    bderr = where(suberr1 eq 0,nbderr)
    if nbderr gt 0 then suberr1[bderr]=999999.
    
    estimates = [mean([med1[maxarr1[j]], med2[maxarr2[j]]]), 2.0, mean(xrec[xlo:xhi]), recny/2, slp, 0.0]      
    parinfo = replicate({limited:[0,0],limits:[0.0,0.0],fixed:0},6)
    parinfo[0].limited[0] = 1  & parinfo[0].limits[0] = 0.0                   ;; height
    parinfo[1].limited = 1     & parinfo[1].limits = [1,6]                    ;; sigma
    parinfo[2].fixed = 1                                                      ;; cenX
    parinfo[3].limited = 1     & parinfo[3].limits = [-2,2]+estimates[3]      ;; cenY
    parinfo[4].limited = 1     & parinfo[4].limits = [0.7,1.3]*estimates[4]   ;; slope
    parinfo[5].limited = 1     & parinfo[5].limits = minmax(subim1)            ;; offset
    ;; Mask central region
    if not keyword_set(arc) then begin
      mask1 = subim1*0+1
      mask1[*,recny/2-5:recny/2+5] = 0
      subim1 = subim1*mask1
    endif else mask1=subim1*0+1
    ;; Mask edges
    fa = {y0:3,y1:nysub1-3,mask:mask1}
    subim1[*,0:3] = 0
    subim1[*,nysub1-3:*] = 0
    pars1 = mpfit2dfun('fire_gauss2d',xx,yy,subim1,suberr1,estimates,parinfo=parinfo,status=status,functargs=fa,/quiet)
    linestr1[j].status = status
    if status gt 0 then begin
      linestr1[j].pars = pars1
      gmodel1 = fire_gauss2d(xx,yy,pars1,_extra=fa)  ;; model with offset and trimming/masking
      tpars1 = pars1    ; no offset
      tpars1[5] = 0.0 
      gmodel = fire_gauss2d(xx,yy,tpars1)
      model1[xlo:xhi,*] += gmodel
      ;; subtract
      residim1[xlo:xhi,*] -= gmodel
    endif
  endfor

  ;; SECOND, fit lines in original (curved) space
  ;;---------------------------------------------

  residim = apim
  model = apim*0
  
  linestr = replicate({num:0L,pars:fltarr(6),xtrace:0.0,ytrace:0.0,status:0},nlines)
  linestr.num = lindgen(nlines)+1
  lmodel = fltarr(subnx)
  for j=0,nlines-1 do begin
    ;; Get the subimage
    xlo = min([maxarr1[j],maxarr2[j]])-12 > 0
    xhi = max([maxarr1[j],maxarr2[j]])+12 < (recnx-1)
    
    ysrt = round(min(outim.ylo[xlo:xhi]))-outim.y0 + 3
    yend = round(max(outim.yhi[xlo:xhi]))-outim.y0 - 3

    subim1 = residim[xlo:xhi,ysrt:yend]
    suberr1 = aperr[xlo:xhi,ysrt:yend]
    subsz1 = size(subim1)
    nxsub1 = subsz1[1]
    nysub1 = subsz1[2]
    ;; parameters: Gaussian height, Gaussian width, central X position,
    ;;             central Y position, slope (deltaX/deltaY), offset
    ;; X and Y should be 2D arrays
    xx = xrec[xlo:xhi]#replicate(1,nysub1)
    yy = replicate(1,nxsub1)#ysub[ysrt:yend]-outim.y0

    ;; fix zero error
    bderr = where(suberr1 eq 0,nbderr)
    if nbderr gt 0 then suberr1[bderr]=999999.
    
    estimates = linestr1[j].pars
    ;; force ycen to be on the trace line, and hold that fixed
    ycen = poly(estimates[2],tstr.tycoef)-outim.y0
    estimates[3] = ycen
    parinfo = replicate({limited:[0,0],limits:[0.0,0.0],fixed:0},6)
    parinfo[0].limited[0] = 1  & parinfo[0].limits[0] = 0.0                   ;; height
    parinfo[1].limited = 1     & parinfo[1].limits = [1,6]                    ;; sigma
    parinfo[2].fixed = 1                                                      ;; cenX
    parinfo[3].limited = 1     & parinfo[3].limits = [-5,5]+estimates[3]      ;; cenY
    parinfo[4].limited = 1     & parinfo[4].limits = [0.7,1.3]*estimates[4]   ;; slope
    parinfo[5].limited = 1     & parinfo[5].limits = [-10,10]+estimates[5]    ;; offset
    ;; Mask central region
    ;if not keyword_set(arc) then begin
    ;  mask1 = subim1*0+1
    ;  mask1[*,recny/2-5:recny/2+5] = 0
    ;  subim1 = subim1*mask1
    ;endif else mask1=subim1*0+1
    mask1=subim1*0+1
    ;; Mask edges
    ;;fa = {y0:3,y1:nysub1-3,mask:mask1}
    fa = {mask:mask1}    
    ;;subim1[*,0:3] = 0
    ;;subim1[*,nysub1-3:*] = 0
    pars1 = mpfit2dfun('fire_gauss2d',xx,yy,subim1,suberr1,estimates,parinfo=parinfo,status=status,functargs=fa,/quiet)
    linestr[j].status = status
    if status gt 0 then begin
      linestr[j].pars = pars1
      gmodel1 = fire_gauss2d(xx,yy,pars1,_extra=fa)  ;; model with offset and trimming/masking
      tpars1 = pars1    ; no offset
      tpars1[5] = 0.0 
      gmodel = fire_gauss2d(xx,yy,tpars1)
      model[xlo:xhi,ysrt:yend] += gmodel
      ;; subtract
      residim[xlo:xhi,ysrt:yend] -= gmodel
      ;; get X/Y position on the trace curve
      ;; where the two curves intersect
      ;; trace curve is pretty flat cross the ~10 pixels we are
      ;; workign with, assume totally flat to get this
      ;; where does the line y-value cross trace ycen (at what x)
      ;; Y = delta_X / slope + Ycen
      ;; X = xcen+(y[i]-ycen)*slope
      ycen1 = ycen
      xcen1 = pars1[2]+(ycen1-pars1[3])*pars1[4]
      ycen2 = poly(xcen1,tstr.tycoef)-outim.y0  ;; iterate
      xcen2 = pars1[2]+(ycen2-pars1[3])*pars1[4]
      linestr[j].xtrace = xcen2
      linestr[j].ytrace = ycen2
      ;; make 1D model of line
      ;;lmodel[xlo:xhi] = gaussian(lindgen(xlo:xhi)...,[linestr[j].pars[0],xcen2,linestr[j].pars[1]])
    endif
  endfor
  
  ;stop

  end
