FUNCTION plotSpVsZ,    desc_label,                   $
                       sfld,                         $ ; needed
                       str_dt,                       $ ; needed
                       SPIDX,                        $ 
                       first_step, last_step,        $
                       plot_mode

; IF ( (STRCMP(plot_mode, "idxvsz",   6, /FOLD_CASE) EQ 1)                                                                  $
; OR   (STRCMP(plot_mode, "idxvsnch", 6, /FOLD_CASE) EQ 1)                                                                  $
;    )                                                                                                                      $
; THEN BEGIN

  data_out_dir                = getDataOutDir(plot_mode)

  aumlaut                     = STRING(228B)
  str_els                     = ' Els' + aumlaut + 'sser '
  str_elsp                    = 'Z!E+!N Els' + aumlaut + 'sser' 
  str_elsm                    = 'Z!E-!N Els' + aumlaut + 'sser'
  log_of_average              = " (Logarithm of the Average) "
  average_of_log              = " (Average of the Logarithm) "

  i_sfld                      = getSfldIndex(sfld)

  CASE i_sfld OF
  '1'  : str_for_case         = ' for Kinetic Energy'             + '!C' + str_dt + log_of_average
  '2'  : str_for_case         = ' for Magnetic Energy'            + '!C' + str_dt + log_of_average
  '3'  : str_for_case         = ' for Total Energy'               + '!C' + str_dt + log_of_average
  '4'  : str_for_case         = ' for '      + str_elsp           + '!C' + str_dt + log_of_average
  '5'  : str_for_case         = ' for '      + str_elsm           + '!C' + str_dt + log_of_average
  '6'  : str_for_case         = ' for Total' + str_els + 'Energy' + '!C' + str_dt + log_of_average
  '7'  : str_for_case         = ' for Kinetic Energy'             + '!C' + str_dt + average_of_log
  '8'  : str_for_case         = ' for Magnetic Energy'            + '!C' + str_dt + average_of_log
  '9'  : str_for_case         = ' for ' + str_elsp                + '!C' + str_dt + average_of_log
  '10' : str_for_case         = ' for ' + str_elsm                + '!C' + str_dt + average_of_log
  ELSE : str_for_case         = ' Call a doctor.'
  ENDCASE
  str_srs_ave                 = getSubrunRangeString(first_step,last_step)
  i_y                         = 1
  IF (STRCMP(plot_mode, "idxvsz",     6, /FOLD_CASE) EQ 1) THEN BEGIN
    sp_vs_z_eps_out           = data_out_dir +  sfld + '/' + 'spidx_vs_z_' +  sfld + '_' + str_srs_ave + '.eps'
    PRINT, "plotting spectral index vs z to file ", sp_vs_z_eps_out
    i_x                       = 0
    str_x_title               = 'z'
    str_title                 = 'Spectral Index (a) vs. z' + str_for_case
    
  ENDIF ELSE BEGIN


    IF (STRCMP(plot_mode, "idxvsnch", 6, /FOLD_CASE) EQ 1) THEN BEGIN
      sp_vs_z_eps_out         = data_out_dir +  sfld + '/' + 'spidx_vs_nch_' +  sfld + '_' + str_srs_ave + '.eps'
      PRINT, "plotting spectral index vs normalized cross helicity to file ", sp_vs_z_eps_out
      i_x                     = 4
      str_x_title             = '!7r!X'
      str_title               = 'Spectral Index (a) vs. Normalized Cross-Helicity (!7r!X)' + '!C' + str_for_case
    ENDIF
  ENDELSE

  str_y_title                 = 'a'

  x_max                       = MAX(SPIDX[*,i_x])
  x_min                       = MIN(SPIDX[*,i_x])

  y_max                       = MAX(SPIDX[*,i_y])
  y_min                       = MIN(SPIDX[*,i_y])

  x_rng                       = [x_min, x_max]
  y_rng                       = [y_min, y_max]

  SET_PLOT, 'PS'
  DEVICE, /COLOR, /ENCAPSULATED
  DEVICE, FILENAME            = sp_vs_z_eps_out

  PLOT, SPIDX[*,i_x], SPIDX[*,i_y],                                                                                       $
        CHARSIZE              = 1.2,                                                                                      $
        LINESTYLE             = 0,                                                                                        $
        YTICKFORMAT           = '(E8.1)',                                                                                 $
        XMARGIN               = [12,4],                                                                                   $
        YMARGIN               = [4,6],                                                                                    $
        XRANGE                = x_rng,                                                                                    $
        YRANGE                = y_rng,                                                                                    $
        XTITLE                = str_x_title,                                                                              $
        YTITLE                = str_y_title,                                                                              $
        THICK                 = 2,                                                                                        $
        TITLE                 = str_title


  ave_spec                    = MEAN(SPIDX[*,i_y])

  str_ave_spec                = STRMID(STRTRIM(ave_spec,2),0,6)
PRINT, "ave_spec = ", ave_spec

  AVEGX                       = [x_rng[0], x_rng[1]  ]
; AVEGY                       = [ ave_spec, ave_spec ]
  AVEGY                       = [ -2.0, -2.0 ]

  OPLOT, AVEGX, AVEGY,                                                                                                    $ 
    LINESTYLE                 = 2,                                                                                        $ 
    THICK                     = 2
 
  XYOUTS, 13.0, -0.25, 'Average Spectral Index = ' + str_ave_spec
  DEVICE, /CLOSE
          
  SET_PLOT, 'X'

  RETURN, 0
END
