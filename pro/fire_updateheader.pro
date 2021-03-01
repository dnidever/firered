function fire_updateheader,im,head

  ;; Update header, have to construct new structure

  newim = {file:im.file,flux:im.flux,err:im.err,mask:im.mask,$
           x:im.x,y:im.y,nx:im.nx,ny:im.ny,$
           head:head,exptype:im.exptype,subimage:im.subimage}

  return,newim

  end
