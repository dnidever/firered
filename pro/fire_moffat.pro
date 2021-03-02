function fire_moffat,x,pars

  ;; Moffat function
  ;; copied from mpfitpeak.pro
  wid = abs(pars[2]) > 1e-20
  u = ((x-pars[1])/wid)^2
  if n_elements(pars) GE 5 then f = pars[4] else f = 0
  if n_elements(pars) GE 6 then f = f + pars[5]*x
  return, f + pars[0] / (u + 1)^pars[3]

end
