FUNCTION fetch_layer, layer, step, label, KK

  null_check         = !null

  n_lines            = LONG64(0)

  datafile           = fetch_datafile(layer, step, label)

  IF (STRCMP(datafile, 'file_not_found')) THEN BEGIN
    PRINT, "fetch_layer: ERROR - cannot open file for layer ", layer, " of step ", step
    LYR              = DBLARR(1,1)
  ENDIF ELSE BEGIN
    IF (STRMATCH(datafile, '*.gz', /FOLD_CASE) EQ 1) THEN BEGIN
      OPENR, data_unit, datafile, /COMPRESS, /GET_LUN, ERROR = op_err
      n_lines        = FILE_LINES(datafile, /COMPRESS)
    ENDIF ELSE BEGIN
      OPENR, data_unit, datafile, /GET_LUN, ERROR = op_err
      n_lines        = FILE_LINES(datafile)
    ENDELSE
  ENDELSE

  ip1                = scan_parameters( 'p1', 0, label)
  ip2                = scan_parameters( 'p2', 0, label)
  n3                 = scan_parameters( 'p3', 0, label)
  mp                 = scan_parameters( 'np', 0, label)

  x_res              = 2^ip1
  y_res              = 2^ip2   

  dx                 = 1.0D / DOUBLE( x_res )
  dy                 = 1.0D / DOUBLE( y_res )

  local_layer        =  ((layer - 1) mod n3 ) + 1

  lines_per_layer    = n_lines / n3
  expected_lines     = LONG64(x_res * y_res)

  IF (expected_lines NE lines_per_layer)  THEN BEGIN
    PRINT, 'fetch_layer: WARNING - lines_per layer does not match expected number of lines.'
    PRINT, ' '
    PRINT, 'fetch_layer: lines_per_layer = ', lines_per_layer
    PRINT, 'fetch_layer: expected_lines  = ', expected_lines
    PRINT, 'fetch_layer: x_res           = ', x_res
    PRINT, 'fetch_layer: y_res           = ', y_res
  ENDIF

  X                  = FLTARR(lines_per_layer)
  Y                  = FLTARR(lines_per_layer)

  idx                = LONG64(0)

  FOR I = 0, LONG64(x_res) - 1 DO BEGIN
    FOR J = 0, LONG64(y_res) - 1 DO BEGIN

      idx            = LONG64((I*x_res) + J)
      X[idx]         = LONG64(I) * dx
      Y[idx]         = LONG64(J) * dy

    ENDFOR
  ENDFOR

  l_start            = (local_layer - 1) * lines_per_layer + 1
  l_stop             = l_start + lines_per_layer - 1

  IF (op_err NE 0) THEN BEGIN
    PRINT, "fetch_layer: ERROR - file exists but cannot be opened"
    LYR              = DBLARR(1,1)
  ENDIF ELSE BEGIN

    line             = ""
    READF, data_unit, line
    cols             = N_ELEMENTS(StrSplit(line))
    Point_Lun, data_unit, 0
    LYR              = DBLARR(cols, lines_per_layer)
    col_str          = STRTRIM(cols,2)
    fmt_str          = '(' + col_str + '(e26.18,1x),:)'

    SKIP_LUN, data_unit, l_start - 1, /LINES

    READF, data_unit, FORMAT = fmt_str, LYR

    LYR              = TRANSPOSE(LYR)

    FREE_LUN, data_unit

    IF (cols EQ 4 ) THEN BEGIN

      Slice          = FLTARR(cols + 2, lines_per_layer)
      Slice          = TRANSPOSE(Slice)

      Slice[*,0]     = X[  *    ]
      Slice[*,1]     = Y[  *    ]
      Slice[*,2:5]   = LYR[*,0:3]

    ENDIF ELSE BEGIN

      IF (cols EQ 2 AND KK NE null_check) THEN BEGIN ; for backward compatibility

        Slice        = FLTARR(cols + 4, lines_per_layer)
        Slice        = TRANSPOSE(Slice)

        P            = REFORM(LYR[*,0], x_res, y_res)
        A            = REFORM(LYR[*,1], x_res, y_res)

        p_inf        = WHERE(FINITE(P), count_p_inf)
        a_inf        = WHERE(FINITE(A), count_a_inf)

        IF (count_p_inf NE 0) THEN  PRINT, 'fetch_layer: WARNING - found nans in P), first nan at ', p_inf[0] & STOP
        IF (count_a_inf NE 0) THEN  PRINT, 'fetch_layer: WARNING - found nans in A), first_nan at ', a_inf[0] & STOP

        FTP          = FFT(P, /DOUBLE)
        FTA          = FFT(A, /DOUBLE)

        FTO          = FTP * ( KK^2 )
        FTJ          = FTA * ( KK^2 )

        O_Sqr        = REAL_PART(FFT(FTO, 1, /DOUBLE))
        J_Sqr        = REAL_PART(FFT(FTJ, 1, /DOUBLE))


        O            = REFORM(O_Sqr, x_res*y_res)
        J            = REFORM(J_Sqr, x_res*y_res)

        Slice[*,0]   = X[*]
        Slice[*,1]   = Y[*]
        Slice[*,2] = LYR[*,0]
        Slice[*,3] = LYR[*,1]
        Slice[*,4]   = O[*]
        Slice[*,5]   = J[*]
      
      ENDIF ELSE BEGIN
        PRINT, 'fetch_layer: ERROR - two-column data set with k-vector array undefined'
        Slice        = null_check
      ENDELSE

    ENDELSE
  ENDELSE

  RETURN, Slice

END
