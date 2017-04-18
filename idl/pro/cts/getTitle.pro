FUNCTION getTitle, qty, desc_label, n_step, n_slice, prefix, global_q_minmax

  str_qmin              = STRTRIM(global_q_minmax[0],2)
  str_qmax              = STRTRIM(global_q_minmax[1],2)

  time                  = scan_parameters('tstart', n_step, desc_label)

  str_time              = STRMID(STRTRIM(time,2),0,5)

  mp                    = scan_parameters('np', 0, desc_label)
  n3                    = scan_parameters('p3', 0, desc_label)
  zl                    = scan_parameters('zl', 0, desc_label)

  dz                    = zl / (mp * n3)

  z_loc                 = STRMID(STRTRIM(n_slice*dz, 2),0,5)
  
  CASE qty OF
  'p':BEGIN
        qty_plotted     = 'Stream Function '
        qty_label       = '!9F!X'
      END
  'a':BEGIN
        qty_plotted     = 'Flux Function '
        qty_label       = 'A'
      END
  'o':BEGIN
        qty_plotted     = 'Vorticity '
        qty_label       = '!9W!X'
      END
  'j':BEGIN
        qty_plotted     = 'Current Density '
        qty_label       = 'J'
      END
  ELSE:
  ENDCASE

    str_title = 'Contours of ' + qty_plotted + qty_label + ' in layer z = ' + z_loc + ': t = ' + str_time 

  RETURN, str_title
END
