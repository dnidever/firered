function fire_refcorrect,im

  ;; Correct image with reference pixels

  newim = im
  
  ;; Use reference pixel around edge to correct
  for j=0,3 do begin
    refim1 = im[j*512:(j+1)*512-1,0:3]
    refim2 = im[j*512:(j+1)*512-1,2044:2047]
    back = median([refim1,refim2])
    newim[j*512:(j+1)*512-1,*] -= back
  endfor

  ;; use a ramp?

  return, newim

end
