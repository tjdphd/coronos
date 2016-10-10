FUNCTION getResString, desc_label

  n3               = scan_parameters('n3',  0, desc_label )
  mp               = scan_parameters('mp',  0, desc_label )
  ip1              = scan_parameters('ip1', 0, desc_label )
  ip2              = scan_parameters('ip2', 0, desc_label )

  i_x_res          = LONG64(2^ip1)
  i_y_res          = LONG64(2^ip2)
  i_z_res          = LONG64(mp*n3)

  IF (i_x_res EQ i_y_res) THEN BEGIN
    res_str        = STRTRIM(i_x_res,2) + '_' + STRTRIM(i_z_res,2)
  ENDIF ELSE BEGIN
    res_str        = STRTRIM(i_x_res,2) + 'X' + STRTRIM(i_y_res,2) + '_' + STRTRIM(i_z_res,2)
  ENDELSE

  RETURN, res_str

END
