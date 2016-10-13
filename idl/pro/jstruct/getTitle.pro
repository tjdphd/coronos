FUNCTION getTitle, qty, action, desc_label, step, threshold, width, prefix, global_j_minmax

  str_jmin              = STRTRIM(global_j_minmax[0],2)
  str_jmax              = STRTRIM(global_j_minmax[1],2)

  ds                    = 'ts'
  time                  = scan_parameters('tstart',step, ds, prefix)
  str_time              = STRTRIM(time,2)

  IF (qty    EQ 'j')     THEN BEGIN
    vol_qty_plotted     = 'Volumetric Current Density'
    iso_surfs_plotted   = 'Current Density Isosurfaces'
    qty_label           = 'J'
  ENDIF ELSE BEGIN
    IF (qty    EQ 'o')   THEN BEGIN
      vol_qty_plotted   = 'Volumetric Vorticity '
      iso_surfs_plotted = 'Vorticity Isosurfaces'
    qty_label           = '!9W!X'
    ENDIF
  ENDELSE

  IF (action EQ 'vol' ) THEN  BEGIN
    str_thr             = STRTRIM(threshold,  2)
    str_thr             = STRMID(str_thr, 0,  4)
    str_title           = vol_qty_plotted + ' at t = ' + str_time + $
    '!C threshold = ' + str_thr + '; ' + qty_label + '!Dmax!N = ' + str_jmax + '; ' + qty_label + '!Dmin!N = ' + str_jmin
  ENDIF ELSE BEGIN
    IF (action EQ 'isf' ) THEN  BEGIN
      str_width         = STRTRIM(width,  2)
      str_width         = STRMID(str_width, 0,  6)
      PRINT, 'getTitle: str_width = ', str_width
      str_title         = iso_surfs_plotted + ' at t = ' + str_time + $
      '!C width = ' + str_width + '; ' + qty_label + '!Dmax!N = ' + str_jmax + '; ' + qty_label + '!Dmin!N = ' + str_jmin
    ENDIF
  ENDELSE

  RETURN, str_title
END
