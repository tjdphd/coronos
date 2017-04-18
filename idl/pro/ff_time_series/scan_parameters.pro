<<<<<<< HEAD
FUNCTION scan_parameters, param,          $
                          n_step,         $
                          desc_label,     $
                          PREFIX = prefix

FORWARD_FUNCTION scan_parameters

cur_dir = GETENV('PWD')
IF (n_step EQ 0) THEN BEGIN
   par_file = '/rmct2.00.tmp-bak'
ENDIF ELSE BEGIN
         ip1        = scan_parameters('ip1', 0, desc_label )
         n3         = scan_parameters('n3',  0, desc_label )
         mp         = scan_parameters('mp',  0, desc_label )

         x_res      = 2^ip1
         z_res      = n3 * mp

         i_x_res    = UINT(x_res)
         i_z_res    = UINT(z_res)

         str_x_res  = STRTRIM(i_x_res, 2)
         str_z_res  = STRTRIM(i_z_res, 2)
         str_n_step = STRTRIM(n_step,  2)

         str_res    = '_' + str_x_res + '_' + str_z_res
          
         par_file   = '/' + prefix + str_res + '.00.' + 'o' + desc_label + str_n_step
      ENDELSE
par_dat = cur_dir + par_file

n_lines   = FILE_LINES(par_dat)

OPENR, par_unit, par_dat, /GET_LUN
var_name  = "var"
record    = 'dummy'
param_len=STRLEN(param)
FOR I = 1, n_lines DO BEGIN
  READF, par_unit, FORMAT = '(A)', record
  IF ( I LE 2) THEN CONTINUE
  rec_len = STRLEN(record)
  FOR J = 0, rec_len - 1 DO BEGIN 
     sub_str    = STRMID(record, J, param_len)
     match_name = STRMATCH(sub_str, param)
     IF (match_name[0] GT 0) THEN BREAK
  ENDFOR
  IF (match_name[0] GT 0) THEN BEGIN
        J_Start = J + param_len + 1  
        str_val = ''
        str_inc = 'x'
        WHILE(str_inc NE '') DO BEGIN
          str_inc = STRMID(record, J_Start, 1)
          J_Start = J_Start + 1
          IF (str_inc EQ ' ') THEN BEGIN
             IF (STRLEN(str_val) EQ 0 ) THEN BEGIN           
               CONTINUE
             ENDIF ELSE BEGIN
                BREAK
             ENDELSE
          ENDIF ELSE BEGIN
             CASE str_inc OF
              '.':  str_val = str_val + str_inc
              '+':  str_val = str_val + str_inc
              '-':  str_val = str_val + str_inc
              'e':  str_val = str_val + str_inc
              'E':  str_val = str_val + str_inc
              '0':  str_val = str_val + str_inc
              '1':  str_val = str_val + str_inc
              '2':  str_val = str_val + str_inc
              '3':  str_val = str_val + str_inc
              '4':  str_val = str_val + str_inc
              '5':  str_val = str_val + str_inc
              '6':  str_val = str_val + str_inc
              '7':  str_val = str_val + str_inc
              '8':  str_val = str_val + str_inc
              '9':  str_val = str_val + str_inc
             ELSE:  
             ENDCASE
                ENDELSE
        ENDWHILE
     BREAK
  ENDIF
ENDFOR
FREE_LUN, par_unit

CASE param OF
    'tstart' : value = FLOAT(str_val)
    'ip1'    : value = FLOAT( str_val)
    'ip2'    : value = FLOAT( str_val)
    'n3'     : value = UINT( str_val)
    'mp'     : value = UINT( str_val)
    'nw'     : value = UINT( str_val)
    'cnw'    : value = FLOAT(str_val)   
    'zl'     : value = FLOAT(str_val)      
    'nu'     : value = FLOAT(str_val)
    'eta'    : value = FLOAT(str_val)    
    'esp'    : value = FLOAT(str_val)    
    'dt'     : value = FLOAT(str_val)     
    'ndt'    : value = UINT( str_val)    
    's1'     : value = FLOAT(str_val)
    's2'     : value = FLOAT(str_val)     
    'tta'    : value = FLOAT(str_val)    
    'ilnr'   : value = UINT( str_val)       
    'iptest' : value = UINT( str_val)     
    'ilt'    : value = UINT( str_val)        
    'q1'     : value = FLOAT(str_val)     
    'q2'     : value = FLOAT(str_val)     
    'dtr'    : value = FLOAT(str_val)    
    'imain'  : value = UINT( str_val)      
    'ss1'    : value = FLOAT(str_val)    
    'qp'     : value = FLOAT(str_val)     
    'ftpow'  : value = FLOAT(str_val)  
    'ediss1' : value = FLOAT(str_val) 
    'kc'     : value = FLOAT(str_val)     
    'numold' : value = UINT( str_val)     
    'tauE'   : value = FLOAT(str_val)   
    'f'      : value = FLOAT(str_val)      
    'dtau'   : value = FLOAT(str_val)   
    'oldnum' : value = UINT( str_val)    
    'deng'   : value = FLOAT(str_val)   
    'tauC'   : value = FLOAT(str_val)   
    'aveKE'  : value = FLOAT(str_val)  
    'aveME'  : value = FLOAT(str_val)  
    'AVEz'   : value = FLOAT(str_val)   
    'AVEpn'  : value = FLOAT(str_val)  
    'AVEpv'  : value = FLOAT(str_val)  
    'Beta'   : value = FLOAT(str_val)   
    'igen'   : value = UINT( str_val)       
    'bdrys'  : value = UINT( str_val)      
    'rcount' : value = UINT( str_val)     
    'aveK2'  : value = FLOAT(str_val)  
    'nf'     : value = UINT( str_val)          
    'AVEpe'  : value = FLOAT(str_val)  
 ENDCASE
 
 
 
RETURN, value
=======
FUNCTION scan_parameters, param, n_step, desc_label

  FORWARD_FUNCTION scan_parameters

  this_step        = n_step

  IF ( NOT ( ( STRCMP(param, "srun",4)) AND (n_step = -1) )) THEN BEGIN

    next_step      = scan_parameters("srun", -1, desc_label)
    last_step      = next_step - 1
    n_step         = this_step

  ENDIF ELSE BEGIN
    last_step      = -1
  ENDELSE

  cur_dir          = GETENV('PWD')

  IF (n_step EQ 0) THEN BEGIN

    par_file       = '/crs_init.in'

  ENDIF ELSE BEGIN

    IF (n_step GT 0 AND n_step LT last_step)  THEN BEGIN

      prefix       = scan_parameters('prefix', 0, desc_label )
      ip1          = scan_parameters('p1',     0, desc_label )
      n3           = scan_parameters('p3',     0, desc_label )
      mp           = scan_parameters('np',     0, desc_label )

      x_res        = 2^ip1
      z_res        = n3 * mp

      i_x_res      = UINT(x_res)
      i_z_res      = UINT(z_res)

      str_x_res    = STRTRIM(i_x_res, 2)
      str_z_res    = STRTRIM(i_z_res, 2)
      str_n_step   = STRTRIM(n_step,  2)

      str_res      = '_' + str_x_res + '_' + str_z_res

      par_file     = '/' + prefix + str_res + '.00.' + 'o' + desc_label + str_n_step

    ENDIF ELSE BEGIN
      IF (STRCMP(param,"srun",4) AND (n_step EQ -1)) THEN BEGIN
        par_file   = '/coronos.in'
      ENDIF ELSE BEGIN
        IF (n_step EQ last_step) THEN BEGIN
          par_file = '/coronos.in'
          n_step   = this_step
        ENDIF ELSE BEGIN
          PRINT, "scan_parameters: ERROR - Something is wrong, this should never be output"
        ENDELSE
      ENDELSE
    ENDELSE
  ENDELSE

  par_dat        = cur_dir + par_file

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
>>>>>>> 2828a96c7de58aa9dbd3162460334aa15b677d54
END
