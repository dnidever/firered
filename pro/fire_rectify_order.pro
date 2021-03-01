function fire_rectify_order,tstr,im,exact=exact

  ;; tstr - trace structure for one order
  ;; im - original image structure
  
  if n_elements(exact) eq 0 then exact=1  ; default
  
  npix = 2048
  x = findgen(npix)
  xx = findgen(npix)#replicate(1,npix)
  yy = replicate(1,npix)#findgen(npix)

  ;; Get subimage
  apim = fire_getorderimage(tstr,im)
  
  ;x1 = x[tstr.bndx0:tstr.bndx1]
  ylo = poly(apim.x,tstr.bndy0coef) > 0
  yhi = poly(apim.x,tstr.bndy1coef) < (npix-1)
  y0 = apim.y[0]
  ;y0 = round(min(ylo)) > 0
  ;y1 = round(max(yhi)) < (npix-1)
  ;subim = im[tstr.bndx0:tstr.bndx1,y0:y1]
  ;subxx = xx[tstr.bndx0:tstr.bndx1,y0:y1]
  ;subyy = yy[tstr.bndx0:tstr.bndx1,y0:y1]    
  ;submask = (subyy ge (poly(subxx,tstr.bndy0coef)>0) and $
  ;           subyy le (poly(subxx,tstr.bndy1coef)<(npix-1)))
  ;subsz = size(subim)
  ;subnx = subsz[1]
  ;subny = subsz[2]  
  
  ;x1 = findgen(subnx)
  ytrace = poly(apim.x,tstr.tycoef)-apim.y[0]
  
  ;; Rectify
  nx = apim.nx
  ;nx = tstr.bndx1-tstr.bndx0+1
  ;ny = max((round(yhi)<(npix-1)) - (round(ylo)>0))+1
  ny = apim.ny
  ;; add extra space in case the trace is not in the middle of the boundary
  ny += 20
  ;; ny must be odd, center in the middle
  if ny mod 2 eq 0 then ny+=1
  nyhalf = ny/2
  recim = fltarr(nx,ny)
  recerr = fltarr(nx,ny)+1e30
  recmask = bytarr(nx,ny)
    
  ;; Interpolate
  if keyword_set(exact) then begin
    used = lonarr(ny)
    for j=0,nx-1 do begin
      ysrt = round(ylo[j])-y0
      yend = round(yhi[j])-y0
      ny1 = yend-ysrt+1
      ;yoff = nyhalf-(round(ytrace[j])-ysrt)
      yoff = ytrace[j]-nyhalf
      ;; Only use good pixels
      gdin = where(reform(apim.mask[j,ysrt:yend]) eq 1,ngdin,comp=bdin,ncomp=nbdin)
      ;; input and output y-arrays
      yin_orig = apim.y[ysrt:yend]-y0
      yin = yin_orig-yoff
      yout0 = ceil(min(yin[gdin]))
      yout1 = floor(max(yin[gdin]))
      yout = lindgen(yout1-yout0+1)+yout0
      ;; Interpolate flux
      fluxout = spline(yin[gdin],reform(apim.flux[j,yin_orig[gdin]]),yout)
      recim[j,yout] = fluxout
      ;; Interpolate the mask
      interp,yin,float(reform(apim.mask[j,yin_orig])),yout,maskoutint
      maskout = bytarr(n_elements(maskoutint))+1
      bdmask = where(maskoutint lt 0.5,nbdmask)
      if nbdmask gt 0 then maskout[bdmask] = 0
      recmask[j,yout] = maskout
      ;; Interpolate err
      errin = reform(apim.err[j,ysrt:yend])
      if nbdin gt 0 then begin
        interp,gdin,errin[gdin],bdin,bderr
        errin[bdin] = bderr
      endif
      interp,yin[gdin],errin[gdin],yout,errout
      ;errout = spline(yin[gdin],errin[gdin],yout)
      bderrout = where(errout gt 1e20,nbderrout)
      if nbderrout gt 0 then errout[bderrout] = 1e30
      if nbdmask gt 0 then errout[bdmask] = 1e30
      recerr[j,yout] = errout
      ;; Rows that we used
      used[yout] OR= 1
    endfor
    ;; Fix bad pixels
    bderr = where(recerr le 0.0,nbderr)
    if nbderr gt 0 then begin
      recerr[bderr] = 1e30
      recmask[bderr] = 0
    endif
    
  ;; Pixel shifts
  endif else begin
    
    used = lonarr(ny)
    for j=0,nx-1 do begin
      ysrt = round(ylo[j])-y0
      yend = round(yhi[j])-y0
      ny1 = yend-ysrt+1
      yoff = nyhalf-(round(ytrace[j])-ysrt)
      recim[j,yoff:yoff+ny1-1] = reform(apim.flux[j,ysrt:yend])
      recerr[j,yoff:yoff+ny1-1] = reform(apim.err[j,ysrt:yend])
      recmask[j,yoff:yoff+ny1-1] = reform(apim.mask[j,ysrt:yend])      
      used[yoff:yoff+ny1-1] OR= 1 
   endfor
  endelse
    
  ;; Trim unused rows, but make sure the "center" is still along the
  ;; central row at the end, and that there are an odd number of rows
  recim0 = recim
  recerr0 = recerr
  recmask0 = recmask  
  used0 = used
  ;; Trim at beginning
  ycenter = nyhalf
  while (used[0] eq 0) do begin
    recim = recim[*,1:*]
    recerr = recerr[*,1:*]
    recmask = recmask[*,1:*]
    used = used[1:*]
    ycenter -= 1
  endwhile
  ;; Trim at end
  while (used[-1] eq 0) do begin
    recim = recim[*,0:-2]
    recerr = recerr[*,0:-2]
    recmask = recmask[*,0:-2]    
    used = used[0:-2]
    ;; ycenter unchanged
  endwhile

  ;; rows above and below the "center" row
  nabove = n_elements(recim[0,*])-ycenter-1
  nbelow = ycenter
  
  ;; Center row NOT in the middle row of the array
  if nabove ne nbelow then begin

    ;; Add rows at the bottom
    if nbelow lt nabove then begin
      nadd = nabove-nbelow
      tempim = recim
      temperr = recerr
      tempmask = recmask      
      recim = fltarr(nx,n_elements(tempim[0,*])+nadd)
      recerr = fltarr(nx,n_elements(temperr[0,*])+nadd)+1e30
      recmask = fltarr(nx,n_elements(tempmask[0,*])+nadd)      
      recim[*,nadd:*] = tempim
      recerr[*,nadd:*] = temperr
      recmask[*,nadd:*] = tempmask         
    
    ;; Add rows at the top
    endif else begin
      nadd = nbelow-nabove
      tempim = recim
      temperr = recerr
      tempmask = recmask      
      recim = fltarr(nx,n_elements(tempim[0,*])+nadd)
      recerr = fltarr(nx,n_elements(temperr[0,*])+nadd)+1e30
      recmask = fltarr(nx,n_elements(tempmask[0,*])+nadd)      
      recim[*,0:n_elements(tempim[0,*])-1] = tempim
      recerr[*,0:n_elements(temperr[0,*])-1] = temperr
      recmask[*,0:n_elements(tempmask[0,*])-1] = tempmask
    endelse
  endif
  
  ;; Create image object
  sz = size(recim)
  out = {file:im.file,flux:recim,err:recerr,mask:recmask,$
         x:apim.x,y:lindgen(sz[2])-ycenter,nx:apim.nx,ny:sz[2],head:im.head,$
         exptype:im.exptype,subimage:1,rectified:1,$
         tstr:tstr}
  
  return, out
  
end
