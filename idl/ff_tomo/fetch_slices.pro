FUNCTION fetch_slices, first_slice, last_slice, n_step, desc_label

COMMON fnc, prefix, inc_res, nfld
COMMON loc, eps_out_dir

ip1          = scan_parameters('ip1', 0, desc_label )
ip2          = scan_parameters('ip2', 0, desc_label )
n3           = scan_parameters('n3',  0, desc_label )
mp           = scan_parameters('mp',  0, desc_label )
zl           = scan_parameters('zl',  0, desc_label )
 
current_file = 'empty'

FOR kth_slc  = first_slice, last_slice DO BEGIN

   datafile  = fetch_datafile(kth_slc, n_step, desc_label)

   IF (kth_slc EQ 1) THEN BEGIN

         n_slc_file                         = 0
         n_prior_files                      = 0
         line_set                           = 1

   ENDIF ELSE BEGIN
             IF ( n3 GT 1) THEN BEGIN
                   s                        = kth_slc - 1
                   n_slc_file               = (s - ( (s + n3) mod n3 )) / n3
                   IF (n_slc_file GT 0) THEN BEGIN
                      n_prior_files         = n_slc_file
                      line_set              = kth_slc - (n_prior_files * n3)
                   ENDIF ELSE BEGIN
                             n_prior_files  = 0
                             line_set       = kth_slc - (n_prior_files * n3)
                         ENDELSE
             ENDIF ELSE BEGIN
                       n_slc_file           = kth_slc - 1
                       n_prior_files        = n_slc_file
                       line_set             = 1
                   ENDELSE
         ENDELSE

   IF (kth_slc EQ first_slice) THEN BEGIN

      x_res           = 2^ip1
      y_res           = 2^ip2   
      z_res           = (n3 * mp)

      i_x_res         = LONG64(x_res)
      i_y_res         = LONG64(y_res)
      i_z_res         = LONG64(z_res)

      dx              = 1  / ( x_res )
      dy              = 1  / ( y_res )
      dz              = zl / ( z_res )
        
      OPENR, data_unit, datafile, /COMPRESS, /GET_LUN, ERROR = op_err

      IF (op_err NE 0) THEN BEGIN

         df_len       = STRLEN(datafile)
         datafile_unc = STRMID(datafile, 0, df_len - 3)

         OPENR, data_unit, datafile_unc, /GET_LUN, ERROR = op_err

         IF ( op_err NE 0 ) THEN BEGIN

            PRINT, FORMAT = '(a57,1x,a28)',                                              $
                            'fetch_slices: ERROR - cannot return any requested slices.', $
                            'Cannot open associated file.'
            
            Slices  = FLTARR(4, 1, n3)
            Slices  = TRANSPOSE(Slices)

            BREAK

         ENDIF ELSE BEGIN

                    datafile   = datafile_unc
                    PRINT,  FORMAT = '(a50,1x,a16,1x,a50,1x,a8)',                $
                           'fetch_slices: WARNING - Compressed file not found.', $
                           'Opening the file', datafile, 'instead.'

                    file_lines = LONG64(0)
                    file_lines = FILE_LINES(datafile)

               ENDELSE

      ENDIF    ELSE BEGIN

                  file_lines   = LONG64(0)
                  file_lines   = FILE_LINES(datafile, /COMPRESS)

               ENDELSE
      
      line                     = ""               
      READF, data_unit, line

      cols                     = N_ELEMENTS(StrSplit(line))

      Point_Lun, data_unit, 0

      current_file             = datafile

      slc_lines                = LONG64(0)
      slc_lines                = file_lines / n3

      expected_filelines       = i_x_res * i_y_res * n3

      IF (expected_filelines NE file_lines) THEN BEGIN

         expected_slc_lines    = i_x_res * i_y_res
         n_slices              = UINT(file_lines / expected_slc_lines)
         extra_lines           = (file_lines - (n_slices * expected_slc_lines))

         IF (file_lines NE (extra_lines + (n_slices * expected_slc_lines)) ) THEN BEGIN
            PRINT,'fetch_slices: WARNING - Something is wrong with the math'
            PRINT, ' '
            PRINT, 'file_lines         = ', file_lines
            PRINT, 'expected_slc_lines = ', expected_slc_lines
            PRINT, 'n_slices           = ', n_slices
            PRINT, 'extra_lines        = ', extra_lines
         ENDIF
         
         PRINT, 'fetch_slices: WARNING - Found fewer lines than expected in the file', $
                 datafile

         IF (n_slices LT line_set ) THEN BEGIN

            PRINT, FORMAT = '(a58,1x,a18)',                                              $
                            'fetch_slices: ERROR - cannot return any requested slices.', $
                            'Insufficent data.'
            Slices     = FLTARR(4, 1, n3)
            Slices     = TRANSPOSE(Slices)
            BREAK
         ENDIF

      ENDIF ELSE BEGIN
              n_slices = (last_slice - first_slice + 1) 
            ENDELSE
      
      A                = FLTARR(cols, slc_lines)
      X                = FLTARR(1,    slc_lines)
      Y                = FLTARR(1,    slc_lines)  

      X                = TRANSPOSE(X)
      Y                = TRANSPOSE(Y)

      cur_idx          = LONG64(0)

      FOR I = 0, x_res - 1 DO BEGIN
         FOR J = 0, y_res - 1 DO BEGIN

            cur_idx    = LONG64((I * x_res) + J)
            X[cur_idx] = LONG64(I) * dx
            Y[cur_idx] = LONG64(J) * dx

         ENDFOR
      ENDFOR

      Slices           = FLTARR(cols + 2, slc_lines, n3)
      Slices           = TRANSPOSE(Slices)

   ENDIF ELSE BEGIN

             IF ( current_file NE datafile) THEN BEGIN          ; This should be modified
               FREE_LUN, (data_unit)                            ; as above to account for the
               OPENR, data_unit, datafile, /COMPRESS, /GET_LUN  ; possibility of incomplete       
               current_file = datafile                          ; datafiles              
             ENDIF   
         ENDELSE

   skip_lines          = slc_lines * (line_set - 1)
   
   PRINT, FORMAT='(a28,1x,i2,1x,a2,1x,i2,1x,a10,1x,i8,1x,a2,1x,i8,1x,a12,1x,a150,a3)',               $
                 'fetch_slices: fetching slice', line_set, 'of', n3, 'from lines', skip_lines + 1,   $
                 'to', skip_lines + slc_lines, 'of data file', datafile, '...'
   
   col_str             = STRTRIM(cols,2)
   fmt_str             = '(' + col_str + '(e24.16,1x),:)'

   print, "fetch_slices: fmt_str = ", fmt_str

   READF, data_unit, FORMAT = fmt_str , A

   A                   = TRANSPOSE(A)
   
   IF ( ( kth_slc - first_slice + 1) LE n3 ) THEN BEGIN
      
      req_slices       = (last_slice - first_slice + 1)

      IF (n3 GT 1) THEN BEGIN

         slc_idx       = kth_slc MOD req_slices

         Slices[slc_idx, *, 0] = X[*,0]
         Slices[slc_idx, *, 1] = Y[*,0]

         FOR ith_col = 0, cols - 1 DO BEGIN
           Slices[slc_idx, *, ith_col + 2] = A[*, ith_col]
         ENDFOR

      ENDIF ELSE BEGIN

                Slices[*, 0]   = X[*,0]
                Slices[*, 1]   = Y[*,0]

         FOR ith_col = 0, cols - 1 DO BEGIN
           Slices[*, ith_col + 2] = A[*, ith_col]
         ENDFOR

            ENDELSE
      A                        = TRANSPOSE(A)
   ENDIF ELSE BEGIN
             PRINT, 'I = ', I
             PRINT, 'fetch_slices: WARNING - Slice array full. Returning only slices ', first_slice, ' to ', kth_slc, '.'
         ENDELSE

   n_slices                    = n_slices - 1

   IF ( n_slices EQ 0 ) THEN BREAK

ENDFOR

IF (N_ELEMENTS(data_unit) GE 1) THEN FREE_LUN, (data_unit)

RETURN, Slices
END
