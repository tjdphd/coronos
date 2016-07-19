FUNCTION global_q_minmax_by_slice, n_step, i_qty, first_slice, last_slice, desc_label, KK

  i_x                             = 6 ; 4
  i_j                             = 5 ; 3
  i_o                             = 4 ; 2
  i_a                             = 3 ; 1
  i_p                             = 2 ; 0

  IF (i_qty EQ -1) THEN BEGIN
  
        PRINT, 'max_q: WARNING - value ', qty, ' for the variable qty not recognized.'
        PRINT, 'Setting indext i_qty to default value for contours of the current density j.'
  
        i_qty                     = 5

  ENDIF
  
  n3                              = scan_parameters('n3',  0, desc_label )
  ip1                             = scan_parameters('ip1', 0, desc_label )
  ip2                             = scan_parameters('ip2', 0, desc_label )
  
  x_res                           = 2^ip1
  y_res                           = 2^ip2   
  
  i_x_res                         = LONG64(x_res)
  i_y_res                         = LONG64(y_res)


  IF (i_qty NE i_x) THEN BEGIN

    max_q_out                     = -1.0e-100
    min_q_out                     =  1.0e+100
    min_qp_out                    =  1.0e+100
    max_qm_out                    = -1.0e-100

    max_q_in                      = -1.0e-100
    min_q_in                      =  1.0e+100
    min_qp_in                     =  1.0e+100
    max_qm_in                     = -1.0e-100

  ENDIF ELSE BEGIN

    max_q_out                     = FLTARR(4)
    min_q_out                     = FLTARR(4)
    min_qp_out                    = FLTARR(4)
    max_qm_out                    = FLTARR(4)

    max_q_in                      = FLTARR(4)
    min_q_in                      = FLTARR(4)
    min_qp_in                     = FLTARR(4)
    max_qm_in                     = FLTARR(4)

    max_q_out[*]                  = -1.0e-100
    min_q_out[*]                  =  1.0e+100
    min_qp_out[*]                 =  1.0e+100
    max_qm_out[*]                 = -1.0e-100


    max_q_in[*]                   = -1.0e-100
    min_q_in[*]                   =  1.0e+100
    min_qp_in[*]                  =  1.0e+100
    max_qm_in[*]                  = -1.0e-100

  ENDELSE
  
  FOR I = first_slice, last_slice DO BEGIN

    Slice                         = fetch_layer(I, n_step, desc_label, KK)

    size_slice                    = SIZE(Slice)
    nlines                        = size_slice[1]
  
    IF (nlines NE (i_x_res * i_y_res)) THEN BEGIN
   
       PRINT, 'max_q: size_slices   = ', size_slices
       PRINT, FORMAT   = '(a40,1x,i4,1x,a18)', 'max_q_by_slice: WARNING - skipping slice', I, 'data not available'
 
       CONTINUE
 
    ENDIF ELSE BEGIN
   
      IF (i_qty NE i_x) THEN BEGIN

        PRINT, 'global_q_minmax_by_slice: i_qty = ', i_qty

        QTY                       = Slice[   *, i_qty]
 
        max_q_in                  = MAX(QTY)
        min_q_in                  = MIN(QTY)
 
        gtz_idx                   = WHERE(QTY GT 0., gtz_cnt)
        ltz_idx                   = WHERE(QTY LT 0., ltz_cnt)
 
        IF ( gtz_cnt GT 0) THEN BEGIN
           min_qp_in              = MIN(QTY[WHERE(QTY GT 0., /NULL)])
        ENDIF ELSE BEGIN
           min_qp_in              = +1.0e-100
        ENDELSE
        IF ( ltz_cnt GT 0) THEN BEGIN
           max_qm_in              = MAX(QTY[WHERE(QTY LT 0., /NULL)])
        ENDIF ELSE BEGIN
           max_qm_in              = -1.0e-100
        ENDELSE

        IF max_q_in  GT max_q_out  THEN max_q_out    = max_q_in
        IF min_q_in  LT min_q_out  THEN min_q_out    = min_q_in
 
        IF max_qm_in GT max_qm_out THEN max_qm_out   = max_qm_in
        IF min_qp_in LT min_qp_out THEN min_qp_out   = min_qp_in
 
      ENDIF ELSE BEGIN

        FOR I_q = 2, 5 DO BEGIN

          QTY                     = Slice[   *, I_q]
 
          max_q_in[i_q - 2]       = MAX(QTY)
          min_q_in[i_q - 2]       = MIN(QTY)
 
          gtz_idx                 = WHERE(QTY GT 0., gtz_cnt)
          ltz_idx                 = WHERE(QTY LT 0., ltz_cnt)
 
          IF ( gtz_cnt GT 0) THEN BEGIN
             min_qp_in[I_q - 2]   = MIN(QTY[WHERE(QTY GT 0., /NULL)])
          ENDIF ELSE BEGIN
             min_qp_in[I_q - 2]   = +1.0e-100
          ENDELSE
          IF ( ltz_cnt GT 0) THEN BEGIN
             max_qm_in[I_q - 2]   = MAX(QTY[WHERE(QTY LT 0., /NULL)])
          ENDIF ELSE BEGIN
             max_qm_in[I_q - 2]   = -1.0e-100
          ENDELSE

          IF max_q_in[I_q  - 2] GT max_q_out[I_q - 2]  THEN max_q_out[I_q  - 2] = max_q_in[I_q  - 2]
          IF min_q_in[I_q  - 2] LT min_q_out[I_q - 2]  THEN min_q_out[I_q  - 2] = min_q_in[I_q  - 2]
 
          IF max_qm_in[I_q - 2] GT max_qm_out[I_q - 2] THEN max_qm_out[I_q - 2] = max_qm_in[I_q - 2]
          IF min_qp_in[I_q - 2] LT min_qp_out[I_q - 2] THEN min_qp_out[I_q - 2] = min_qp_in[I_q - 2]

        ENDFOR

      ENDELSE

    ENDELSE

  ENDFOR

  IF (i_qty NE i_x) THEN BEGIN
   
    global_q_minmax_by_slice      = FLTARR(4)
   
    global_q_minmax_by_slice[0]   = min_q_out
    global_q_minmax_by_slice[1]   = max_q_out
    global_q_minmax_by_slice[2]   = min_qp_out
    global_q_minmax_by_slice[3]   = max_qm_out

  ENDIF ELSE BEGIN

    global_q_minmax_by_slice      = FLTARR(4,4)

    global_q_minmax_by_slice[0,*] = min_q_out[ *] ; second index is i_qty
    global_q_minmax_by_slice[1,*] = max_q_out[ *]
    global_q_minmax_by_slice[2,*] = min_qp_out[*]
    global_q_minmax_by_slice[3,*] = max_qm_out[*]

  ENDELSE
  
  RETURN, global_q_minmax_by_slice

END
