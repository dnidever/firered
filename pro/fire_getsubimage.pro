function fire_getsubimage,im,xr,yr

  ;; Get sub-image of image object

  ;;   flux, err, mask (good mask), nx, ny, x (global), y (global),
  ;;   type (arc, object, quartz)

  if n_elements(im) eq 0 or n_elements(xr) eq 0 or n_elements(yr) eq 0 then begin
    print,'Syntax - subim = fire_getsubimage(im,xr,yr)'
    return,-1
  endif
  
  ;; Check XR and YR
  if n_elements(xr) ne 2 or n_elements(yr) ne 2 then begin
    print,'XR and YR must be [LOWER_INDEX, UPPER_INDEX]'
    return,-1
  endif
  if min(xr) lt 0 or max(xr) gt im.nx-1 or min(yr) lt 0 or max(yr) gt im.ny-1 then begin
    print,'XR and YR must be between 0 and NX/NY-1'
    return,-1
  endif

  newim = {file:im.file,flux:im.flux[xr[0]:xr[1],yr[0]:yr[1]],$
           err:im.err[xr[0]:xr[1],yr[0]:yr[1]],$
           mask:im.mask[xr[0]:xr[1],yr[0]:yr[1]],$
           x:im.x[xr[0]:xr[1]],y:im.y[yr[0]:yr[1]],$
           nx:xr[1]-xr[0]+1,ny:yr[1]-yr[0]+1,$
           head:im.head,exptype:im.exptype,subimage:0}

  return,newim
  
  end
