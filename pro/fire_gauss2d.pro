function fire_gauss2d,xx,yy,pars,mask=mask,y0=y0,y1=y1

  ;; Created an elongated 2D Gaussian to simulate a line

  ;; parameters: Gaussian height, Gaussian width, central X position,
  ;;             central Y position, slope (deltaX/deltaY), offset
  ;; X and Y should be 2D arrays, but on a regular cartesian grid

  ;; need to know if we want to mask anything
  ;; how about upper/lower limits in Y?

  sz = size(xx)
  nx = sz[1]
  ny = sz[2]
  x = xx[*,0]
  y = reform(yy[0,*])

  if n_elements(mask) eq 0 then mask=xx*0+1
  if n_elements(y0) eq 0 then y0=0
  if n_elements(y1) eq 0 then y1=ny-1
  
  ;; Loop over
  im = xx*0
  height = pars[0]
  sigma = pars[1]
  xcen = pars[2]
  ycen = pars[3]
  slope = pars[4]   ;; deltaX/deltaY
  offset = pars[5]
  for i=y0,y1 do begin
    ;; Get x-position
    cen = xcen+(y[i]-ycen)*slope
    pars1 = [height, cen, sigma, offset]
    im[*,i] = gaussian(x,pars1)
  endfor

  return, im*mask
  
  end
