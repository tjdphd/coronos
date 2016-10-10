PRO run_ave_spec, dsc_lab, sfld, layer, first_step, last_step, minmax_mode

COMMON step, time

  i_k             = 0 ; k-value
  i_pe            = 1 ; kinetic        energy spectrum in layer ? of time step ?
  i_ae            = 2 ; magnetic       energy spectrum in layer ? of time step ?
  i_ts            = 3 ; total          energy spectrum in layer ? of time step ?
  i_zp            = 4 ; elsasser+      energy spectrum in layer ? of time step ?
  i_zm            = 5 ; elsasser-      energy spectrum in layer ? of time step ?
  i_tz            = 6 ; total elsasser energy spectrum in layer ? of time step ?

  CASE sfld OF
  'pe' : i_sfld   = i_pe
  'ae' : i_sfld   = i_ae
  'ts' : i_sfld   = i_ts
  'zp' : i_sfld   = i_zp
  'zm' : i_sfld   = i_zm
  'tz' : i_sfld   = i_tz
  ELSE:  i_sfld   = -1

  ENDCASE
  CASE i_sfld OF
  '0' :  str_sfld = 'kk_'
  '1' :  str_sfld = 'pe_'
  '2' :  str_sfld = 'ae_'
  '3' :  str_sfld = 'ts_'
  '4' :  str_sfld = 'zp_'
  '5' :  str_sfld = 'zm_'
  '6' :  str_sfld = 'tz_'
  ELSE:  str_sfld = 'xx_'
  ENDCASE

  ip1             = scan_parameters('ip1', 0, dsc_lab )                     ; power of 2 giving x-resolution
  ip2             = scan_parameters('ip2', 0, dsc_lab )                     ; power of 2 giving y-resolution
  n3              = scan_parameters('n3' , 0, dsc_lab )                     ; number of slices per data file
  mp              = scan_parameters('mp' , 0, dsc_lab )                     ; number of processors used in run
  zl              = scan_parameters('zl' , 0, dsc_lab )                     ; total height along z of integration volume
   
  x_res           = 2^ip1                                                   ; resolution in x
  y_res           = 2^ip2                                                   ; resolution in y
  z_res           = n3 * mp                                                 ; resolution in z


  n_proc          = UINT(mp / (z_res / layer  )) - 1
  str_zpos        = STRTRIM(zl / (z_res/layer),2)
  loc_layer       = UINT(layer / ( n_proc+1 ))
  str_proc        = STRTRIM(n_proc,   2)
  len_str_proc    = STRLEN(str_proc)

  WHILE (STRLEN(str_proc) LT 3) DO BEGIN
    str_proc      = '0' + str_proc
  ENDWHILE
  
  str_layer       = STRTRIM(loc_layer,2)
  
  WHILE (STRLEN(str_layer) LT 3) DO BEGIN
    str_layer     = '0' + str_layer
  ENDWHILE
  
  i_x_res         = UINT(x_res)
  i_y_res         = UINT(y_res)
  i_z_res         = UINT(z_res)
   
  str_x_res       = STRTRIM(i_x_res, 2)
  str_y_res       = STRTRIM(i_y_res, 2)
  str_z_res       = STRTRIM(i_z_res, 2)
  
  str_zl          = STRTRIM(UINT(zl),2)

  cur_dir         = GETENV('PWD')
  eps_out_dir     = cur_dir + '/ra_spec/' + sfld + '/eps'

  res_str         = str_x_res + '!9X!X'   + str_y_res +'!9X!X'    + str_z_res

  case_res_str    = + ' (L = ' + str_zl    + ')'

  IF (minmax_mode EQ  'global') THEN BEGIN
    glb_minmax_y  = DBLARR(2)
    glb_minmax_y  = global_sfld_minmax(i_sfld, dsc_lab, layer, first_step, last_step)
  ENDIF

  E               = open_spec_data_file( dsc_lab, layer, first_step)
  size_E          = SIZE(E)

  n_lines         = size_E[1]
  n_cols          = size_E[2]

  tfirst          = time
  FOR I = first_step + 1, last_step DO BEGIN

       E_next     = open_spec_data_file( dsc_lab, layer, I)

       FOR J = 0, 6 DO BEGIN
         E[*,J]   = E[*, J] + E_next[*, J]
       ENDFOR

  ENDFOR
  tlast           = time

  str_steps       = '_' + STRTRIM(first_step,2) + '-' + STRTRIM(last_step,2)
  str_dt          = STRTRIM(tfirst,2) + ' < t < ' + STRTRIM(tlast,2)
  tot_steps       = FLOAT(last_step - first_step + 1)

  E[*,*]          = E[*,*] / tot_steps

  data_out_file   = cur_dir + '/ra_spec/' + 'ra_' + sfld + '_L-' + str_zl + '.dat'



  OpenW, data_unit, data_out_file, /GET_LUN

  FOR K = 0, n_lines - 1 DO BEGIN

     PRINTF, data_unit,  E[K, 0:6], FORMAT='(7(E16.8,1x),:/)'

  ENDFOR

  FREE_LUN, data_unit

  case_file_str   = str_sfld + "ra_spec_" + str_x_res + '_' + str_z_res + '_zl_' + str_zl $
                    + '_proc-' + str_proc + '_layer-' + str_layer + '_srs' + str_steps 
  eps_out         = eps_out_dir + '/' + case_file_str + '.eps'

  str_ave         = ' Averaged over ' + str_dt
  str_ps          = ' Power Spectrum for layer ' + str_layer + ' of process '
  str_line_two    = 'at z = ' + str_zpos + ' (L = ' + str_zl + '),' + str_ave
  aumlaut         = STRING(228B)
  str_els         = ' Els' + aumlaut + 'sser '
  str_elsp        = 'Z!E+!N Els' + aumlaut + 'sser' 
  str_elsm        = 'Z!E-!N Els' + aumlaut + 'sser'
  CASE i_sfld OF
  '1' : str_title = 'Kinetic Energy'             + str_ps + str_proc + '!C' + str_line_two
  '2' : str_title = 'Magnetic Energy'            + str_ps + str_proc + '!C' + str_line_two
  '3' : str_title = 'Total Energy'               + str_ps + str_proc + '!C' + str_line_two
  '4' : str_title = str_elsp                     + str_ps + str_proc + '!C' + str_line_two
  '5' : str_title = str_elsm                     + str_ps + str_proc + '!C' + str_line_two
  '6' : str_title = 'Total' + str_els + 'Energy' + str_ps + str_proc + '!C' + str_line_two
  ELSE: str_title = 'Something is terribly wrong'
  ENDCASE

  nz_idx          = WHERE(E[*, i_sfld] GT 0., nz_count)
  PRINT, 'nz_count = ', nz_count
  PRINT, 'i_sfld   = ', i_sfld



  min_x           = MIN(E[nz_idx, i_k   ])
  max_x           = MAX(E[*,      i_k   ])

  x_rng           = [min_x, max_x]

  min_y           = MIN(E[nz_idx, i_sfld])
  max_y           = MAX(E[*,      i_sfld])

  PRINT, 'min_y = ', min_y
  PRINT, 'max_y = ', max_y

  y_rng           = [min_y, max_y]

  SET_PLOT, 'PS' 
  DEVICE, /ENCAPSULATED
  DEVICE, FILENAME = eps_out

  PLOT, E[*,i_k], E[*, i_sfld],    $
        CHARSIZE    = 0.8,         $
        LINESTYLE   = 0,           $
        YTICKFORMAT = '(E8.1)',    $
        XMARGIN     =[12,4],       $
        YMARGIN     =[4,4],        $
        XRANGE      = x_rng,       $
        YRANGE      = y_rng,       $
        XTITLE      = 'k',         $
        YTITLE      = str_y_title, $
        THICK       = 2,           $
        /XLOG,                     $
        /YLOG,                     $
        TITLE       = str_title

 
  DEVICE, /CLOSE
        
  SET_PLOT, 'X'

END
