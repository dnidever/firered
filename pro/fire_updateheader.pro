function fire_updateheader,im,head

  ;; Update header, have to construct new structure

  ntags = n_tags(im)
  tags = tag_names(im)
  newim = create_struct(tags[0],im.(0))
  for i=1,ntags-1 do begin
    if tags[i] ne 'HEAD' then begin

       newim = create_struct(newim,tags[i],im.(i))
    endif else begin
      newim = create_struct(newim,'head',head)
    endelse
  endfor
  
  return,newim

  end
