FUNCTION construct_eps_outfile_name, prefix, suffix, n_slice, n_step, tot_steps
COMMON loc, eps_out_dir

str_n_slice     = STRTRIM(n_slice, 2)

mp              = scan_parameters('mp', 0, desc_label)
n3              = scan_parameters('n3', 0, desc_label)
z_res           = n3 * mp
i_z_res         = UINT(z_res)
str_z_res       = STRTRIM(i_z_res, 2)

max_z_digits    = STRLEN(str_z_res)
z_digits        = STRLEN(str_n_slice)
zero_pad        = max_z_digits - z_digits
zero_str        = ''

IF (zero_pad NE 0) THEN BEGIN
  FOR I = 1, zero_pad DO BEGIN
    zero_str    = zero_str + '0'
  ENDFOR
ENDIF
str_n_slice     = zero_str + str_n_slice

str_tot_steps   = STRTRIM(tot_steps,2)
str_n_step      = STRTRIM(n_step,2)
max_z_digits    = STRLEN(str_tot_steps)
z_digits        = STRLEN(str_n_step)
zero_pad        = max_z_digits - z_digits
zero_str        = ''

IF (zero_pad NE 0) THEN BEGIN
  FOR I = 1, zero_pad DO BEGIN
    zero_str    = zero_str + '0'
  ENDFOR
ENDIF
str_n_step      = zero_str + str_n_step

IF (STRLEN(eps_out_dir) EQ 0) THEN BEGIN
   out_dir      = GETENV('PWD')
ENDIF ELSE BEGIN
   out_dir      = eps_out_dir
      ENDELSE
PRINT, 'eps output will be place in the directory: ', out_dir
eps_file_out    = out_dir + '/' + prefix + '_slc-' + str_n_slice + '_stp-' + str_n_step + suffix + '.eps'

RETURN, eps_file_out
END
