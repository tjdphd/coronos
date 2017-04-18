FUNCTION open_spec_data_file, desc_label, layer, n_step

COMMON step, time

ip1          = scan_parameters('p1',       0, desc_label )                ; power of 2 giving x-resolution
ip2          = scan_parameters('p2',       0, desc_label )                ; power of 2 giving y-resolution
n3           = scan_parameters('p3',       0, desc_label )                ; number of slices per data file
mp           = scan_parameters('np',       0, desc_label )                ; number of processors used in run
zl           = scan_parameters('zl',       0, desc_label )                ; total height along z of integration volume
data_dir     = scan_parameters('data_dir', 0, desc_label )
 
x_res        = 2^ip1                                                      ; resolution in x
y_res        = 2^ip2                                                      ; resolution in y
z_res        = n3 * mp                                                    ; resolution in z

IF (layer GT n3) THEN BEGIN
  is_layer_even = (layer MOD 2) EQ 0
  IF (is_layer_even) THEN BEGIN
    n_proc     = CEIL(UINT(layer/n3)-1)
  ENDIF ELSE BEGIN
    n_proc     = FLOOR(UINT(layer/n3))
  ENDELSE
ENDIF ELSE BEGIN
  n_proc     = UINT(0)
ENDELSE

IF ( ((n3 + layer) MOD n3 ) NE 0 ) THEN BEGIN
  loc_layer  =  ( (n3 + layer) MOD n3 )
ENDIF ELSE BEGIN
  loc_layer  = n3
ENDELSE

str_zpos     = STRTRIM(zl / (z_res/layer),2)
str_proc     = STRTRIM(n_proc,   2)
len_str_proc = STRLEN(str_proc)
WHILE (STRLEN(str_proc) LT 3) DO BEGIN
  str_proc   = '0' + str_proc
ENDWHILE


str_layer    = STRTRIM(loc_layer,2)

WHILE (STRLEN(str_layer) LT 3) DO BEGIN
  str_layer  = '0' + str_layer
ENDWHILE

i_x_res      = UINT(x_res)
i_y_res      = UINT(y_res)
i_z_res      = UINT(z_res)
 
str_x_res    = STRTRIM(i_x_res, 2)
str_z_res    = STRTRIM(i_z_res, 2)

str_res_lab  = str_x_res + '_' + str_z_res

str_n_step   = STRTRIM(n_step,2)

; e.g. 'spectra_128_128.015_004.ots1'

data_file    = 'spectra_' + str_res_lab + '.' + str_proc + '_' + str_layer + '.o' + desc_label + str_n_step 

cur_dir      = GETENV('PWD')
data_file    = cur_dir + '/' + data_dir + '/' + data_file

OPENR, data_unit, data_file, /GET_LUN

nlines       = LONG64(0)
nlines       = FILE_LINES(data_file) - 2

time         = 0.0D
z_coord      = 0.0D
READF, data_unit, time
READF, data_unit, z_coord

line         = ""
READF, data_unit, line
cols         = N_ELEMENTS(StrSplit(line))
str_cols     = STRTRIM(cols,2)
str_fmt      = '(' + str_cols + '(e24.12,1x),:/)'

Point_Lun, data_unit, 0
SKIP_LUN,  data_unit, 2, /LINES

Eline        = DBLARR(cols, 1)
E            = DBLARR(cols, nlines)

FOR I = 0, nlines - 1 DO BEGIN

  READF, data_unit, FORMAT = str_fmt, Eline
   E[*, I]   = Eline
  
ENDFOR

E            = TRANSPOSE(E)  

FREE_LUN, data_unit

RETURN, E

END
