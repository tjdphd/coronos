FUNCTION init_jvox, first_slice, last_slice

  n3               = scan_parameters('n3',  0, desc_label )
  ip1              = scan_parameters('ip1', 0, desc_label )
  ip2              = scan_parameters('ip2', 0, desc_label )

  i_x_res          = LONG64(2^ip1)
  i_y_res          = LONG64(2^ip2)  
  i_z_res          = LONG64(last_slice - first_slice + 1)

  JVOX             = BYTARR(i_z_res, i_x_res, i_y_res)

  RETURN, JVOX

END
