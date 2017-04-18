FUNCTION getSpectrumTitle, sfld, layer, str_dt

  n3                          = scan_parameters('n3' , 0, dsc_lab )                     ; number of slices per data file
  mp                          = scan_parameters('mp' , 0, dsc_lab )                     ; number of processors used in run
  zl                          = scan_parameters('zl' , 0, dsc_lab )                     ; total height along z of integration volume

  str_zl                      = STRTRIM(UINT(zl),2)

  z_res                       = n3 * mp                                                 ; resolution in z

  str_tot_layers              = STRTRIM(UINT(z_res), 2)
  dz                          = zl / z_res
  zcoord                      =    dz * layer
  str_zpos                    = STRMID(STRTRIM(zcoord,2),0,5)

; str_ave                     = ' Average for '            + str_dt
  str_layer                   = STRTRIM(layer,2)
  WHILE (STRLEN(str_layer) LT 3) DO BEGIN
    str_layer                 = '0' + str_layer
  ENDWHILE


; str_ave                     = ' Averaged over '            + str_dt

  str_ps                      = ' Power Spectrum for layer ' + str_layer     + ' of ' + str_tot_layers
; str_line_two                = 'at z = '   + str_zpos       + ' (L = '      + str_zl + '),' + str_ave
  str_line_two                = 'at z = '   + str_zpos       + ' (L = '      + str_zl + '),' + str_dt
  aumlaut                     = STRING(228B)
  str_els                     = ' Els' + aumlaut + 'sser '
  str_elsp                    = 'Z!E+!N Els' + aumlaut + 'sser' 
  str_elsm                    = 'Z!E-!N Els' + aumlaut + 'sser'

  log_of_average              = " (Logarithm of the Average) "
  average_of_log              = " (Average of the Logarithm) "

  i_sfld                      = getSfldIndex(sfld)
  CASE i_sfld OF
  '1'  : str_title            = 'Kinetic Energy'             + str_ps + '!C' + str_line_two + '!C' + log_of_average
  '2'  : str_title            = 'Magnetic Energy'            + str_ps + '!C' + str_line_two + '!C' + log_of_average
  '3'  : str_title            = 'Total Energy'               + str_ps + '!C' + str_line_two + '!C' + log_of_average
  '4'  : str_title            = str_elsp                     + str_ps + '!C' + str_line_two + '!C' + log_of_average
  '5'  : str_title            = str_elsm                     + str_ps + '!C' + str_line_two + '!C' + log_of_average
  '6'  : str_title            = 'Total' + str_els + 'Energy' + str_ps + '!C' + str_line_two + '!C' + log_of_average
  '7'  : str_title            = 'Kinetic Energy'             + str_ps + '!C' + str_line_two + '!C' + average_of_log
  '8'  : str_title            = 'Magnetic Energy'            + str_ps + '!C' + str_line_two + '!C' + average_of_log
  '9'  : str_title            = str_elsp                     + str_ps + '!C' + str_line_two + '!C' + average_of_log
  '10' : str_title            = str_elsm                     + str_ps + '!C' + str_line_two + '!C' + average_of_log
  ELSE : str_title            = 'Something is terribly wrong'
  ENDCASE

  RETURN, str_title

END
