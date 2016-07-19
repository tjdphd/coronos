FUNCTION jZeroOfX, x

  jzofx     = BESELJ(x, 0)

  RETURN, jzofx

END

PRO plotGofx, nx, x_f

  ZERO      = 0.0D
  TWO       = 2.0D
  THREE     = 3.0D

    X       = FLTARR(nx)
    GOFX    = FLTARR(nx)

    x_inc   = x_f / FLOAT(nx - 1)
  
    X[0]    = 0

    FOR I   = 1, nx - 1 DO X[I] = X[I - 1] + x_inc 
   
    GOFX[*] = QSIMP('jZeroOfx', ZERO, X[*])


    x_min   = MIN(X)
    x_max   = MAX(X)

    y_min   = MIN(GOFX)
    y_max   = MAX(GOFX)

    x_rng   = [x_min, x_max]
    y_rng   = [y_min, (THREE / TWO)*y_max]

    plt     = PLOT(X, GOFX,            $
                   THICK      = 2,     $
                   AXIS_STYLE = 2,     $
                   LINESTYLE  = 0,     $
                   XRANGE     = x_rng, $
                   YRANGE     = y_rng  $
                  )

    
    txt         = TEXT(10.0, 1.75, '$G(x) = \int_0^x dx\prime J_0(x\prime)$',/DATA, FONT_SIZE=18,FONT_NAME='Helvetica')
    ax          = plt.AXES
    ax[0].TITLE = '$x$'
    ax[1].TITLE = '$G(x)$'

END
