FUNCTION scan_parameters, param, n_step, desc_label

  FORWARD_FUNCTION scan_parameters

  cur_dir        = GETENV('PWD')

  IF (n_step EQ 0) THEN BEGIN

    par_file     = '/coronos.in'

  ENDIF ELSE BEGIN

    prefix       = scan_parameters('prefix',   0, desc_label )
    ip1          = scan_parameters('p1',       0, desc_label )
    n3           = scan_parameters('p3',       0, desc_label )
    mp           = scan_parameters('np',       0, desc_label )
    data_dir     = scan_parameters('data_dir', 0, desc_label )

    x_res        = 2^ip1
    z_res        = n3 * mp

    i_x_res      = UINT(x_res)
    i_z_res      = UINT(z_res)

    str_x_res    = STRTRIM(i_x_res, 2)
    str_z_res    = STRTRIM(i_z_res, 2)
    str_n_step   = STRTRIM(n_step,  2)

      str_res    = '_' + str_x_res + '_' + str_z_res

    par_file     = '/' + data_dir + '/' + prefix + str_res + '.00.' + 'o' + desc_label + str_n_step

  ENDELSE

  par_dat        = cur_dir + par_file

  IF (NOT FILE_TEST(par_dat)) THEN par_dat = cur_dir + '/coronos.in'

  n_lines        = FILE_LINES(par_dat)

  OPENR, par_unit, par_dat, /GET_LUN

  record         = 'dummy'

  str_type       = ''

  FOR I = 1, n_lines DO BEGIN

    READF, par_unit, FORMAT = '(A)', record
    split_name   = STRSPLIT(record)
    str_name     = STRTRIM(STRMID(record,split_name[0],split_name[1]),2)
  
    IF (param EQ str_name) THEN BEGIN

      sub_record = STRTRIM(STRMID(record,split_name[1],split_name[2]),2)
      split_val  = STRSPLIT(sub_record)
      str_val    = STRTRIM(STRMID(sub_record, split_val[0],split_val[1]),2)


      s_sub_rec  = STRTRIM(STRMID(sub_record,split_val[1],split_val[2]),2)
      split_type = STRSPLIT(s_sub_rec)
      str_type   = STRTRIM(STRMID(s_sub_rec, split_type[0],split_type[1]),2)

      
    ENDIF
  ENDFOR

  FREE_LUN, par_unit

  CASE str_type of
    'int': value = UINT(str_val)
    'dbl': value = DOUBLE(str_val)
    'str': value = str_val
    'log': value = UINT(str_val)
    ELSE : PRINT, 'scan_parameters: WARNING - parameter ', str_name, ' is of unknown type.'
  ENDCASE
  

  RETURN, value

END
