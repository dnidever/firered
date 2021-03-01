function fire_rectify_order,tstr,im,exact=exact

  npix = 2048
  x = findgen(npix)
  xx = findgen(npix)#replicate(1,npix)
  yy = replicate(1,npix)#findgen(npix)

  ;; Get subimage
  x1 = x[tstr.bndx0:tstr.bndx1]
  ylo = poly(x1,tstr.bndy0coef) > 0
  yhi = poly(x1,tstr.bndy1coef) < (npix-1)
  y0 = round(min(ylo)) > 0
  y1 = round(max(yhi)) < (npix-1)
  subim = im[tstr.bndx0:tstr.bndx1,y0:y1]
  subxx = xx[tstr.bndx0:tstr.bndx1,y0:y1]
  subyy = yy[tstr.bndx0:tstr.bndx1,y0:y1]    
  submask = (subyy ge (poly(subxx,tstr.bndy0coef)>0) and $
             subyy le (poly(subxx,tstr.bndy1coef)<(npix-1)))
  subsz = size(subim)
  subnx = subsz[1]
  subny = subsz[2]  
  
  x1 = findgen(subnx)
  ytrace = poly(x1,tstr.tycoef)-y0
  
  ;; Rectify
  nx = tstr.bndx1-tstr.bndx0+1
  ny = max((round(yhi)<(npix-1)) - (round(ylo)>0))+1
  ;; add extra space in case the trace is not in the middle of the boundary
  ny += 20
  ;; ny must be odd, center in the middle
  if ny mod 2 eq 0 then ny+=1
  nyhalf = ny/2
  recim = fltarr(nx,ny)

  exact = 1
  
  ;; Interpolate
  if keyword_set(exact) then begin
    used = lonarr(ny)
    for j=0,nx-1 do begin
      ysrt = round(ylo[j])-y0
      yend = round(yhi[j])-y0
      ny1 = yend-ysrt+1
      yoff = nyhalf-(round(ytrace[j])-ysrt)
stop
      recim[j,yoff:yoff+ny1-1] = reform( (subim*submask)[j,ysrt:yend] )
      used[yoff:yoff+ny1-1] OR= 1 
    endfor
     
  ;; Pixel shifts
  endif else begin
    used = lonarr(ny)
    for j=0,nx-1 do begin
      ysrt = round(ylo[j])-y0
      yend = round(yhi[j])-y0
      ny1 = yend-ysrt+1
      yoff = nyhalf-(round(ytrace[j])-ysrt)
      recim[j,yoff:yoff+ny1-1] = reform( (subim*submask)[j,ysrt:yend] )
      used[yoff:yoff+ny1-1] OR= 1 
   endfor
  endelse
    
  ;; Trim unused rows, but make sure the "center" is still along the
  ;; central row at the end, and that there are an odd number of rows
  recim0 = recim
  used0 = used
  ;; Trim at beginning
  ycenter = nyhalf
  while (used[0] eq 0) do begin
    recim = recim[*,1:*]
    used = used[1:*]
    ycenter -= 1
  endwhile
  ;; Trim at end
  while (used[-1] eq 0) do begin
    recim = recim[*,0:-2]
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
      temp = recim
      recim = fltarr(nx,n_elements(temp[0,*])+nadd)
      recim[*,nadd:*] = temp
    
    ;; Add rows at the top
    endif else begin
      nadd = nbelow-nabove
      temp = recim
      recim = fltarr(nx,n_elements(temp[0,*])+nadd)
      recim[*,0:n_elements(temp[0,*])-1] = temp
    endelse
  endif
  
  return, recim
  
end
