FUNCTION fetch_layer, layer, step, label

  ip1             = scan_parameters('p1', 0, label)
  ip2             = scan_parameters('p2', 0, label)
   n3             = scan_parameters( 'p3', 0, label)

  n_lines         = LONG64(0)
  
  datafile        = fetch_datafile(layer, step, label)

  IF (STRCMP(datafile, 'file_not_found')) THEN BEGIN
    PRINT, "fetch_layer: ERROR - cannot open file ", datafile, " for layer ", layer, " of step ", step
    A             = DBLARR(1,1)
  ENDIF ELSE BEGIN

    PRINT, "fetch_layer: opening file ", datafile, "..."

    IF (STRMATCH(datafile, '*.gz', /FOLD_CASE) EQ 1) THEN BEGIN
      OPENR, data_unit, datafile, /COMPRESS, /GET_LUN, ERROR = op_err
      n_lines     = FILE_LINES(datafile, /COMPRESS)
    ENDIF ELSE BEGIN
      OPENR, data_unit, datafile, /GET_LUN, ERROR = op_err
      n_lines     = FILE_LINES(datafile)
    ENDELSE
  ENDELSE

  local_layer     =  ((layer - 1) mod n3 ) + 1

  lines_per_layer = n_lines / n3

  l_start         = (local_layer - 1) * lines_per_layer + 1
  l_stop          = l_start + lines_per_layer - 1
 
  IF (op_err NE 0) THEN BEGIN
    PRINT, "fetch_layer: ERROR - file exists but cannot be opened"
    A             = DBLARR(1,1)
  ENDIF ELSE BEGIN

    line          = ""
    READF, data_unit, line
    cols          = N_ELEMENTS(StrSplit(line))
    Point_Lun, data_unit, 0
    A             = DBLARR(cols, lines_per_layer)
    col_str       = STRTRIM(cols,2)
    fmt_str       = '(' + col_str + '(e24.16,1x),:)'

    SKIP_LUN, data_unit, l_start - 1, /LINES

    READF, data_unit, FORMAT = fmt_str, A

    A             = TRANSPOSE(A)

    FREE_LUN, data_unit

  ENDELSE


  RETURN, A

END
