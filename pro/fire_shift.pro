function fire_shift,im,sh,dim

  ;;  Shift in Y a different amount for each column
  ;;  use array indexing

  if n_elements(im) eq 0 or n_elements(sh) eq 0 or n_elements(dim) eq 0 then begin
    print,'Syntax - imsh = fire_shift(im,shift,dim)'
    return,-1
  endif
  
  sz = size(im)
  nx = sz[1]
  ny = sz[2]
  
  ;; Shift in X
  if dim eq 1 then begin

    ;; Double-check sh array
    if n_elements(sh) ne ny then begin
      print,'SHIFT array must be equal to NY'
      return,-1
    endif
    
    ;; Original indices
    xx1 = lindgen(nx)#replicate(1,ny)
    yy1 = replicate(1,nx)#lindgen(ny)

    ;; same for right side, but y values shift by sh
    xx2 = xx1 + replicate(1,nx)#sh
    yy2 = yy1
    bd = where(xx2 lt 0,nbd)
    if nbd gt 0 then xx2[bd] += nx
    bd = where(xx2 gt (nx-1),nbd)
    if nbd gt 0 then xx2[bd] -= nx  

    imsh = im*0
    imsh[xx2,yy2] = im[xx1,yy1]
  endif
  
  ;; Shift in Y
  if dim eq 2 then begin

    ;; Double-check sh array
    if n_elements(sh) ne nx then begin
      print,'SHIFT array must be equal to NX'
      return,-1
    endif
     
    ;; Original indices
    xx1 = lindgen(nx)#replicate(1,ny)
    yy1 = replicate(1,nx)#lindgen(ny)

    ;; same for right side, but y values shift by sh
    xx2 = xx1
    yy2 = yy1 + sh#replicate(1,ny)
    bd = where(yy2 lt 0,nbd)
    if nbd gt 0 then yy2[bd] += ny
    bd = where(yy2 gt (ny-1),nbd)
    if nbd gt 0 then yy2[bd] -= ny  

    imsh = im*0
    imsh[xx2,yy2] = im[xx1,yy1]
  endif
    
  return,imsh
  
  end
