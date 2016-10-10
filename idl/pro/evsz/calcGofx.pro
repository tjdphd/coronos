FUNCTION jZeroOfX, x

  jzofx       = BESELJ(x, 0)

  RETURN, jzofx

END

FUNCTION calcGofx, KK, desc_label

  ZERO                                             = 0.0D                      ; Define some constants for convenience
  HALF                                             = 0.5D                      ; everywhere
  ONE                                              = 1.0D                      ; everywhere
  TWO                                              = 2.0D
  THREE                                            = 3.0D
  PI                                               = 3.141592653589793D
  TWO_THIRDS                                       = TWO / THREE
  TWO_PI                                           = TWO*PI

  CZERO                                            = DCOMPLEX(ZERO,ZERO)       ; and to make sure that DOUBLE's are used
  CI                                               = DCOMPLEX(ZERO,ONE)        ; and to make sure that DOUBLE's are used

  size_kk                                          = SIZE(KK, /DIMENSIONS)


  n1                                               = size_kk[0]
  n2                                               = size_kk[1]

  GOFX                                             = FLTARR(n1,n2)
  ZZLL                                             = FLTARR(n1,n2)
  ZZLL[*,*]                                        = ZERO

  multithread                                      = 1

  IF (multithread) THEN BEGIN

    col_rng                                        = getSliceRange(0, n1-1)
    size_cr                                        = SIZE(col_rng, /DIMENSIONS)
    n_threads                                      = size_cr[0]

    oBridge                                        = OBJARR(n_threads - 1)

    FOR J = 0, n_threads - 1 DO BEGIN

      first_col_sec                                = col_rng[J,0]
      final_col_sec                                = col_rng[J,1]
      n_secs_this_thread                           = final_col_sec - first_col_sec + 1

      IF (J EQ n_threads - 1) THEN BEGIN

        thread_gofx                                = FLTARR(n_secs_this_thread, n2)
        local_gofx                                 = FLTARR(n_secs_this_thread, n2)

        thread_gofx[0:n_secs_this_thread-1,0:n2-1] = QSIMP( 'jZeroOfx',                           $
                                                     ZZLL[first_col_sec:final_col_sec,   0:n2-1], $
                                                     HALF*KK[first_col_sec:final_col_sec,0:n2-1]  $
                                                   )

        local_gofx[*,*]                            = thread_gofx[*,*]

      ENDIF ELSE BEGIN

        func_name                                  = 'jZeroOfx'

        oBridge[J]                                 = OBJ_NEW('IDL_IDLBRIDGE')
        oBridge[J]                                 -> setVar, 'func_name',     func_name
        oBridge[J]                                 -> setVar, 'ZZLL',          ZZLL
        oBridge[J]                                 -> setVar, 'KK',            KK 
        oBridge[J]                                 -> setVar, 'first_col_sec', col_rng[J,0]
        oBridge[J]                                 -> setVar, 'final_col_sec', col_rng[J,1]
        oBridge[J]                                 -> setVar, 'n_secs_this_thread', n_secs_this_thread
        oBridge[J]                                 -> setVar, 'n2', n2
        oBridge[J]                                 -> Execute, "thread_gofx = FLTARR(n_secs_this_thread,n2)"
        oBridge[J]                                 -> Execute,                                              $
        "thread_gofx[0:n_secs_this_thread-1,0:n2-1] = QSIMP( func_name,"                                    $
                                                      +     "ZZLL[first_col_sec:final_col_sec,0:n2-1],"     $
                                                      +     "HALF*KK[first_col_sec:final_col_sec,0:n2-1])", $
        /nowait
      ENDELSE

    ENDFOR

    IF (n_threads GT 1) THEN BEGIN
      notdone                                      = 1
      WHILE notdone DO BEGIN
        done                                       = 0
        FOR J = 0, n_threads - 2 DO done           = done + oBridge[J] -> Status()
        IF (done EQ 0 ) THEN notdone               = done
      ENDWHILE
    ENDIF

    FOR J = 0, n_threads - 1 DO BEGIN

      first_col_sec                                = col_rng[J,0]
      final_col_sec                                = col_rng[J,1]
      n_secs_this_thread                           = final_col_sec - first_col_sec + 1

      IF (J EQ n_threads - 1) THEN BEGIN

        GOFX[first_col_sec:final_col_sec, *]       = local_gofx[*, *]

      ENDIF ELSE BEGIN

        thread_gofx                                = FLTARR(n_secs_this_thread, n2)
        thread_gofx                                = oBridge[J] -> getVar('thread_gofx')
        GOFX[first_col_sec:final_col_sec,0:n2-1]   = thread_gofx[0:n_secs_this_thread-1,0:n2-1] 
        obj_destroy, oBridge[J]

      ENDELSE

    ENDFOR

  ENDIF ELSE BEGIN

    GOFX[0:n1-1,0:n2-1]                            = QSIMP('jZeroOfx', ZZLL[0:n1-1,0:n2-1], HALF*KK[0:n1-1,0:n2-1])

  ENDELSE

  PRINT, 'calcGofx: completed...'

  GOFX_LIN = REFORM(GOFX,n1*n2)

  cur_dir  = GETENV('pwd')
  res_str  = getResString(desc_label)
  gofx_dat = cur_dir + '/ra_evsz/' + 'gofx_' + res_str + '.dat'

  OpenW, gofx_unit, gofx_dat, /GET_LUN

  FOR I = 0, n1*n2 -1 DO BEGIN
    PRINTF, gofx_unit, FORMAT = '(E16.8),:/)', GOFX_LIN[I]
  ENDFOR

  FREE_LUN, gofx_unit
  
  RETURN, GOFX
END
