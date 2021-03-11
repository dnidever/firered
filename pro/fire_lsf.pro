pro fire_lsf


  ;; Get the FIRE LSF for each order from the sky/arc lines
  norders = 22

  lines = mrdfits('/Users/nidever/observing/fire/dec2013/red/fire_lines.fits',1)

  x = findgen(2048)
  for i=1,norders-1 do begin
    ind = where(lines.order eq i,nind)
    lines1 = lines[ind]
    gd = where(lines1.perror[1] gt 1e-10,ngd)
    lines1 = lines1[gd]
    if ngd gt 1 then begin
      coef1 = poly_fit(lines1.xtrace,lines1.pars[1],1,measure_errors=lines1.perror[1]>1e-5)
      coef2 = poly_fit(lines1.xtrace,lines1.pars[1],2,measure_errors=lines1.perror[1]>1e-5)
      coef1b = robust_poly_fit(lines1.xtrace,lines1.pars[1],1)
      coef2b = robust_poly_fit(lines1.xtrace,lines1.pars[1],2)        
      plot,lines1.xtrace,lines1.pars[1],ps=8,tit=i
      oplot,x,poly(x,coef1),co=150
      oplot,x,poly(x,coef2),co=250
      print,i,coef2
    endif else print,'only 1 line'
    stop
  endfor

  ;; Fitting Gaussian Sigma in PIXEL values

  ;; 5 orders only have 1 line!
  
  stop

  end
