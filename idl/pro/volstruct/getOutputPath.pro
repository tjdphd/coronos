FUNCTION getOutputPath, n_step, str_pfx, str_ifx, last_step

; would like this or getInfix to know about:

; - pm
; - n_cntrs
; - inc_cntrs 
; - action
; - rotation angle

  str_last_stp          = STRTRIM(last_step,2)
  len_steps             = STRLEN(str_last_stp)
  str_n_step            = STRTRIM(n_step,2)
  len_n_step            = STRLEN(str_n_step)
  IF (len_steps LT 3) THEN len_steps = 3
  dif_len               = len_steps - len_n_step

  FOR I_Len             = 0, dif_len - 1 DO str_n_step = '0' + str_n_step

  out_dir               = GETENV('PWD') + '/' + str_pfx + '/eps'

  output_path           = out_dir + '/' + str_pfx + '_' + str_ifx + str_n_step + '.eps'
  
  RETURN, output_path

END
