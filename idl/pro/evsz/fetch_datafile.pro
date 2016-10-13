FUNCTION fetch_datafile, n_slice, step, desc_label

   ip1          = scan_parameters('ip1',0, desc_label)
   ip2          = scan_parameters('ip2',0, desc_label)
   n3           = scan_parameters('n3' ,0, desc_label)
   mp           = scan_parameters('mp' ,0, desc_label)
   zl           = scan_parameters('zl' ,0, desc_label)
    
   x_res        = 2^ip1
   y_res        = 2^ip2   
   z_res        = (n3 * mp)
   
   z_res        = (n3 * mp)
   
   IF (n_slice EQ 1) THEN BEGIN
      n_slc_file = 0
   ENDIF ELSE BEGIN
       IF ( n3 GT 1) THEN BEGIN
         s           = n_slice - 1
         n_slc_file  = (s - ( (s + n3) mod n3 )) / n3
       ENDIF ELSE BEGIN
         n_slc_file     = n_slice - 1 
             ENDELSE
         ENDELSE
   i_x_res      = UINT(x_res)
   i_y_res      = UINT(y_res)
   i_z_res      = UINT(z_res)
   
   str_x_res    = STRTRIM(i_x_res, 2)
   str_z_res    = STRTRIM(i_z_res, 2)
   

;IF ( inc_res EQ 'y' ) THEN BEGIN

   str_res      = '_' + str_x_res + '_' + str_z_res
 
;ENDIF ELSE BEGIN

;         str_res = ''

;     ENDELSE

str_n_slice    = STRTRIM(n_slice, 2)
str_n_slc_file = STRTRIM(n_slc_file,2)
str_n_step     = STRTRIM(step,  2)

max_z_digits = STRLEN(str_z_res)
IF (max_z_digits LT 3) THEN max_z_digits=3
z_digits     = STRLEN(str_n_slc_file)
zero_pad     = max_z_digits - z_digits
zero_str     = ''
IF (zero_pad NE 0) THEN BEGIN
  FOR I = 1, zero_pad DO BEGIN
    zero_str = zero_str + '0'
  ENDFOR
ENDIF

str_n_slc_file  = zero_str + str_n_slc_file
cur_dir      = GETENV('PWD')
datafile     = cur_dir + '/' + 'rmct2' + str_res + '.' + str_n_slc_file + '.o' + desc_label + str_n_step + '.gz'

IF ( NOT FILE_TEST(datafile) ) THEN BEGIN
  datafile   = cur_dir + '/' + 'rmct2' + str_res + '.' + str_n_slc_file + '.o' + desc_label + str_n_step
  IF ( NOT FILE_TEST(datafile) ) THEN BEGIN
    datafile = 'file_not_found'
  ENDIF
ENDIF


;print, "fetch_datafile: datafile = ", FILE_BASENAME(datafile)

RETURN, datafile
END
