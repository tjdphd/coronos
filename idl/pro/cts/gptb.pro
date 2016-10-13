FUNCTION TRVSR, r

  COMMON share, mu, npmult

  PI    = DOUBLE(!PI)

  r     = DOUBLE(r)
  mu    = DOUBLE(mu)

  lhs   = (r*mu - ATAN((2.0D*mu*r*(r^2-1.0D)) / ((2.0D*mu*(3.0D*r^2-1.0D)) - (r^2-1.0D)) )) / PI

  rhs   = DOUBLE(npmult)

  RETURN, DOUBLE(lhs - rhs)

END

PRO gptb, sigscale, pimult

  COMMON share, mu, npmult

  r_tol       = DOUBLE(1.0E-09)
  npmult      = pimult

  str_n       = STRTRIM(npmult,2)

  PI          = DOUBLE(!PI)

  ZERO        =  0.0D
  FOURTH      =  0.25D
  HALF        =  0.5D
  ONE         =  1.0D
  TWO         =  2.0D
  THREE       =  3.0D
  FOUR        =  4.0D
  FIVE        =  5.0D
  SIX         =  6.0D
  SEVEN       =  7.0D
  EIGHT       =  8.0D
  NINE        =  9.0D
  TEN         = 10.0D
  TWO_THIRDS  = TWO / THREE
  TWO_PI      = TWO*PI
  CI          = DCOMPLEX(ZERO, ONE )
  CZERO       = DCOMPLEX(ZERO, ZERO)

  sigscale    = DOUBLE(sigscale)
  npts        = LONG64(128)

  l_perp      = ONE
  l_perp_sqd  = l_perp * l_perp

  sigma_x     = sigscale*l_perp/SQRT(TWO)
  sigma_y     = sigscale*l_perp/SQRT(TWO)

  str_sig_x   = STRTRIM(sigma_x,2)
  str_sig_y   = STRTRIM(sigma_y,2)

  PRINT, 'str_sig_x = ', str_sig_x
  PRINT, 'str_sig_y = ', str_sig_y

  sigma_x_sqd = sigma_x * sigma_x
  sigma_y_sqd = sigma_y * sigma_y

  kappa_x     = ONE / (TWO*sigma_x_sqd)
  kappa_y     = ONE / (TWO*sigma_y_sqd)

  mu_x        = kappa_x * l_perp_sqd / FOUR
  mu_y        = kappa_y * l_perp_sqd / FOUR

  r_guess     = [-0.1000,0.0000, 0.1000]
  
  mu          = mu_x
  r_x         = FX_ROOT(r_guess,'TRVSR',/DOUBLE)
  type_rx     = TYPENAME(r_x)

  IF (type_rx NE 'DOUBLE') THEN BEGIN
    PRINT, 'root is type ', type_rx
    STOP
  ENDIF
  
  r_x_check   = TRVSR(r_x)

  IF (ABS(r_x_check) GE r_tol) THEN BEGIN
    PRINT, 'gptb: ERROR - r_x_check exceeds tolerance'
    STOP
  ENDIF

  mu          = mu_y
  r_y         = FX_ROOT(r_guess,'TRVSR',/DOUBLE)
  type_ry     = TYPENAME(r_y)

  IF (type_ry NE 'DOUBLE') THEN BEGIN
    PRINT, 'root is type ', type_ry
    STOP
  ENDIF

  r_y_check   = TRVSR(r_y)

  IF (ABS(r_y_check) GE r_tol) THEN BEGIN
    PRINT, 'gptb: ERROR - r_y_check exceeds tolerance'
    STOP
  ENDIF

  varsigma_x  = r_x * kappa_x
  varsigma_y  = r_y * kappa_y

  PRINT, ''
  PRINT, 'gptb: r_x_check   = ', r_x_check
  PRINT, 'gptb: r_x         = ', r_x
  PRINT, 'gptb: mu_x        = ', mu_x
  PRINT, 'gptb: kappa_x     = ', kappa_x
  PRINT, 'gptb: varsigma_x  = ', varsigma_x
  PRINT, ''
  PRINT, 'gptb: r_y_check   = ', r_y_check
  PRINT, 'gptb: r_y         = ', r_y
  PRINT, 'gptb: mu_y        = ', mu_y
  PRINT, 'gptb: kappa_y     = ', kappa_y
  PRINT, 'gptb: varsigma_y  = ', varsigma_y
  PRINT, ''

  ; r_x  plot

  RATX        = DINDGEN(npts)
  RATX[*]     = TWO*r_x*RATX[*] / RATX[npts-1]

  rx_x_min    = MIN(RATX)
  rx_x_max    = MAX(RATX)

; PRINT, 'rx_x_min = ', rx_x_min
; PRINT, 'r_x      = ', r_x
; PRINT, 'rx_x_max = ', rx_x_max

  rx_x_range  = [rx_x_min, rx_x_max]

  RXRHS       = DBLARR(npts)
  RXLHS       = DBLARR(npts)

  RXRHS[*]    = DOUBLE(npmult)
  RXLHS[*]    = (   RATX[*]*mu_x                                $
                  - ATAN( (TWO*mu_x*RATX[*]*(RATX[*]^2 - ONE) ) $
                            /                                   $
                          ( TWO*mu_x*(THREE*RATX[*]^2 - ONE) -(RATX[*]^2 - ONE)) )) / PI

  rx_y_min    = MIN(RXLHS)
  rx_y_max    = MAX(RXLHS)

; PRINT, 'rx_y_min = ', rx_y_min
; PRINT, 'rx_y_max = ', rx_y_max

  rx_y_range  = [rx_y_min, rx_y_max]

  ; r_y  plot

  RATY        = DINDGEN(npts)
  RATY[*]     = TWO*r_y*RATY[*] / RATY[npts-1]

  ry_x_min    = MIN(RATY)
  ry_x_max    = MAX(RATY)

; PRINT, 'ry_x_min = ', ry_x_min
; PRINT, 'r_y      = ', r_y
; PRINT, 'ry_x_max = ', ry_x_max

  ry_x_range  = [ry_x_min,ry_x_max]

  RYRHS       = DBLARR(npts)
  RYLHS       = DBLARR(npts)

  RYRHS[*]    = DOUBLE(npmult)
  RYLHS[*]    = (   RATY[*]*mu_y                                                                          $
                  - ATAN( ( TWO*mu_y*RATY[*]*(RATY[*]^2 - ONE) )                                          $
                            /                                                                             $
                          ( TWO*mu_y*(THREE*RATY[*]^2 - ONE) -(RATY[*]^2 - ONE)) )) / PI

  ry_y_min    = MIN(RYLHS)
  ry_y_max    = MAX(RYLHS)

; PRINT, 'ry_y_min = ', ry_y_min
; PRINT, 'ry_y_max = ', ry_y_max

  ry_y_range  = [ry_y_min, ry_y_max]

  plt_title   = 'Graphical Solution for $r_x$; $\sigma_x = $' + str_sig_x + ', $m=$' + str_n

  pltlx       = PLOT(RATX, RXLHS,                                                                         $
                     THICK     = 2,                                                                       $
                     XRANGE    = rx_x_range,                                                              $
                     YRANGE    = rx_y_range,                                                              $
                     LINESTYLE = 0,                                                                       $
                     TITLE     = plt_title                                                                $
                    )

  pltrx       = PLOT(RATX, RXRHS,                                                                         $
                     /OVERPLOT,                                                                           $
                     COLOR     = 'red',                                                                   $
                     THICK     = 4,                                                                       $
                     LINESTYLE = 2                                                                        $
                    )

  linerx      = POLYLINE( [r_x,r_x], [ZERO,0.8*rx_y_max], /DATA,                                          $
                          TARGET    = pltlx,                                                              $
                          COLOR     = !COLOR.RED,                                                         $
                          THICK     = 3,                                                                  $
                          LINESTYLE = 2                                                                   $
                        )

  ax          = pltlx.AXES
  ax[0].TITLE = '$r_x$'
  rx_txt_ytitle  = TEXT( 0.025,4.0,                                                                       $
                        '$(1/\pi)\{r_x\mu_x-arctan[2\mu_x r_x(r_x^2-1)/[2\mu_x(3r_x^2-1)-(r_x^2-1)]]\}$', $
                        /DATA, FONT_SIZE=16, TARGET=pltlx)

  plt_title   = 'Graphical Solution for $r_y$; $\sigma_y = $' + str_sig_y + ', $n=$' + str_n

  pltly       = PLOT(RATY, RYLHS,                                                                         $
                     THICK     = 2,                                                                       $
                     XRANGE    = ry_x_range,                                                              $
                     YRANGE    = ry_y_range,                                                              $
                     LINESTYLE = 0,                                                                       $
                     TITLE     = plt_title                                                                $
                    )

  pltry       = PLOT(RATY, RYRHS,                                                                         $
                     /OVERPLOT,                                                                           $
                     COLOR     = 'red',                                                                   $
                     THICK     = 4,                                                                       $
                     LINESTYLE = 2                                                                        $
                    )

  ax          = pltly.AXES
  ax[0].TITLE = '$r_y$'

  ry_txt_ytitle  = TEXT( 0.025,4.0,                                                                       $
                        '$(1/\pi)\{r_y\mu_y-arctan[2\mu_y r_y(r_y^2-1)/[2\mu_y(3r_y^2-1)-(r_y^2-1)]]\}$', $
                        /DATA, FONT_SIZE=16, TARGET=pltly)

  linery      = POLYLINE([r_y, r_y],[ZERO, 0.8*ry_y_max], /DATA,                                          $
                          TARGET    = pltly,                                                              $
                          COLOR     = !COLOR.RED,                                                         $
                          THICK     = 3,                                                                  $
                          LINESTYLE = 2                                                                   $
                        )

  phi_x       = ACOS( kappa_x / (SQRT(varsigma_x^2+kappa_x^2)))
  phi_y       = ACOS( kappa_y / (SQRT(varsigma_y^2+kappa_y^2)))

  PRINT, 'type of phi_x = ', TYPENAME(phi_x)
  PRINT, 'type of phi_y = ', TYPENAME(phi_y)

  gamma_1_x   = EXP(-mu_x) * SIN((varsigma_x * l_perp_sqd / FOUR) + phi_x )
  gamma_1_y   = EXP(-mu_y) * SIN((varsigma_y * l_perp_sqd / FOUR) + phi_y )

   PRINT, ''
   PRINT, 'gptb: phi_x       = ', phi_x
   PRINT, 'gptb: phi_y       = ', phi_y
   PRINT, ''
 
;  PRINT, ''
;  PRINT, 'gptb: gamma_1_x   = ', gamma_1_x
;  PRINT, 'gptb: gamma_1_y   = ', gamma_1_y
;  PRINT, ''
 
   X          = DINDGEN(npts)
   Y          = DINDGEN(npts)
   X[*]       = X[*]/X[npts-1]
   Y[*]       = Y[*]/Y[npts-1]
 
   pnt_inc    = l_perp / (DOUBLE(npts))
 
   n1         = npts
   n2         = npts

   Xf               = DINDGEN((n1-1)/2) + ONE
   is_n1_even       = (n1 MOD 2) EQ 0
   IF (is_n1_even) THEN KX = TWO_PI * [ ZERO, Xf,   DOUBLE(n1)/TWO, DOUBLE(-n1)/TWO  + Xf] /(pnt_inc * n1) $
   ELSE                 KX = TWO_PI * [ ZERO, Xf, -(DOUBLE(n1)/TWO + ONE) + Xf ] /(pnt_inc*n1)
   Yf               = DINDGEN((n2-1)/2) + ONE
   is_n2_even       = (n2 MOD 2) EQ 0
   IF (is_n2_even) THEN KY = TWO_PI * [ ZERO, Yf,   DOUBLE(n2)/TWO, DOUBLE(-n2)/TWO  + Yf] /(pnt_inc * n2) $
   ELSE                 KY = TWO_PI * [ ZERO, Yf, -(DOUBLE(n2)/TWO + ONE) + Yf ]  /(pnt_inc * n2)
 
   XARG        = DBLARR(npts)
   YARG        = DBLARR(npts)
 
   XARG[*]     = (varsigma_x * ( X[*]  - (HALF*l_perp) )^2) + phi_x
   YARG[*]     = (varsigma_y * ( Y[*]  - (HALF*l_perp) )^2) + phi_y
 
   GAMX        = DBLARR(npts)
   GAMY        = DBLARR(npts)
 
   GAMX[*]     = EXP( -kappa_x * (X[*] - (HALF*l_perp))^2 ) * SIN( XARG[*] ) 
   GAMY[*]     = EXP( -kappa_y * (Y[*] - (HALF*l_perp))^2 ) * SIN( YARG[*] )
 
;  PRINT, ''
;  PRINT, 'gptb: GAMY at 0   = ', GAMY[0]
;  PRINT, 'gptb: GAMY at L   = ', GAMY[npts-1]
;  PRINT, 'gptb: gamma_1_y   = ', gamma_1_y
;  PRINT, ''
 
   A_X         = DBLARR(npts)
   B_X         = DBLARR(npts)
 
   A_X         = varsigma_x * ( ONE - ( FOUR * kappa_x *( X[*] - (HALF*l_perp) )^2) )
   B_X         = kappa_x    * ( ONE - ( TWO * (kappa_x + ( (varsigma_x*varsigma_x) / kappa_x ) )*(X[*]-HALF*l_perp)^2) )
 
   A_Y         = DBLARR(npts)
   B_Y         = DBLARR(npts)
 
   A_Y         = varsigma_y * ( ONE - ( FOUR * kappa_y *( Y[*] - (HALF*l_perp) )^2) )
   B_Y         = kappa_y    * ( ONE - ( TWO * (kappa_y + ( (varsigma_y*varsigma_y) / kappa_y ) )*(Y[*]-HALF*l_perp)^2) )
 
   by_should   = kappa_y * (ONE - (TWO*(ONE + (varsigma_y / kappa_y)^2) * (kappa_y*l_perp_sqd / FOUR)))
   ay_should   = varsigma_y* (ONE - kappa_y*l_perp_sqd)
 
;  PRINT, ''
;  PRINT, 'gptb: A_Y at 0   = ', A_Y[0]
;  PRINT, 'gptb: A_Y at L   = ', A_Y[npts-1]
;  PRINT, 'gptb: ay_should  = ', ay_should  
;  PRINT, ''
;  PRINT, 'gptb: B_Y at 0   = ', B_Y[0]
;  PRINT, 'gptb: B_Y at L   = ', B_Y[npts-1]
;  PRINT, 'gptb: by_should  = ', by_should  
;  PRINT, ''
 
   GAMDDX      = DBLARR(npts)
   GAMDDY      = DBLARR(npts)
 
   GAMDDX [*]  = TWO * EXP( -kappa_x * ( X[*]-(HALF*l_perp) )^2 ) * ( A_X * COS(XARG[*]) - B_X * SIN(XARG[*]) ) 
   GAMDDY [*]  = TWO * EXP( -kappa_y * ( Y[*]-(HALF*l_perp) )^2 ) * ( A_Y * COS(YARG[*]) - B_Y * SIN(YARG[*]) ) 
 
;  PRINT, ''
;  PRINT, 'gptb: GAMDDY at 0   = ', GAMDDY[0]
;  PRINT, 'gptb: GAMDDY at L   = ', GAMDDY[npts-1]
;  PRINT, ''
 
   v_norm      =  ONE /SQRT( ((varsigma_x * varsigma_x) + (kappa_x*kappa_x)) * ((varsigma_y * varsigma_y) + (kappa_y*kappa_y)) )
   v_norm      = v_norm * DOUBLE(1.0E04)
   PRINT, 'gptb: v_norm        = ', v_norm
 
   VOR         = DBLARR(npts, npts)
   PHI         = DBLARR(npts, npts)
 
   FOR IX = 0, npts - 1 DO BEGIN
     FOR JY = 0, npts - 1 DO BEGIN
 
       VOR[IX, JY] = -( ( (GAMY[JY] - gamma_1_y) * GAMDDX[IX] ) + ( ( GAMX[IX] - gamma_1_x ) *  GAMDDY[JY] ) )
       PHI[IX, JY] =   (GAMX[IX] - gamma_1_x) * (GAMY[JY]-gamma_1_y)
 
     ENDFOR
   ENDFOR

   PHI[*,*]    = PHI[*,*] * v_norm
   VOR[*,*]    = VOR[*,*] * v_norm

   FFTPHI      = FFT(PHI, /DOUBLE)
   FFTVOR      = FFT(VOR, /DOUBLE)

   FFTOME      = DCOMPLEXARR(n1, n2)
   FFTPSI      = DCOMPLEXARR(n1, n2)

   FOR IX = 0, n1 - 1 DO BEGIN
     FOR JY = 0, n2 - 1 DO BEGIN

       k_sqd           = KX[IX]^2 + KY[JY]^2
       FFTOME[IX,JY]   = FFTPHI[IX,JY] * k_sqd

       IF (k_sqd NE ZERO) THEN BEGIN
         FFTPSI[IX,JY] = FFTVOR[IX,JY] / k_sqd
       ENDIF ELSE BEGIN
        ;FFTPSI[IX,JY] = CZERO
         FFTPSI[IX,JY] = FFTPHI[IX,JY]
       ENDELSE

     ENDFOR
   ENDFOR

  PSI          = REAL_PART(FFT(FFTPSI, /DOUBLE, /INVERSE))
  OME          = REAL_PART(FFT(FFTOME, /DOUBLE, /INVERSE))

  vor_abs_max  = MAX(ABS(VOR))
  ome_abs_max  = MAX(ABS(OME))

  vor_abs_min  = MIN(ABS(VOR))
  ome_abs_min  = MIN(ABS(OME))

  PRINT, ""
  PRINT, "vor_abs_max = ", vor_abs_max
  PRINT, "ome_abs_max = ", ome_abs_max
  PRINT, ""
  PRINT, "vor_abs_min = ", vor_abs_min
  PRINT, "ome_abs_min = ", ome_abs_min
  PRINT, ""

  vor_nz          = WHERE(ABS(VOR) NE ZERO)

  vor_nz_min      = MIN(ABS(VOR[vor_nz]))
  vor_nz_max      = MAX(ABS(VOR[vor_nz]))
  delta_vor_nz    = vor_nz_max - vor_nz_min

  ome_nz_min      = MIN(ABS(OME[vor_nz]))
  ome_nz_max      = MAX(ABS(OME[vor_nz]))
  delta_ome_nz    = ome_nz_max - ome_nz_min

  PRINT, "vor_nz_max   = ", vor_nz_max
  PRINT, "vor_nz_min   = ", vor_nz_min
  PRINT, "delta_vor_nz = ", delta_vor_nz
  PRINT, ""
  PRINT, "ome_nz_max   = ", ome_nz_max
  PRINT, "ome_nz_min   = ", ome_nz_min
  PRINT, "delta_ome_nz = ", delta_ome_nz
  PRINT, ""

  delta_omega_max = ZERO
   big_delta_tol  = DOUBLE(1.0E-07)
   big_delta      = LONG64(0)
   

   FOR IX = 0, npts - 1 DO BEGIN
     FOR JY = 0, npts - 1 DO BEGIN

        IF (ABS(VOR[IX,JY]) GT big_delta_tol AND ABS(OME[IX,JY]) GT big_delta_tol) THEN BEGIN
          delta_omega = ABS(( VOR[IX,JY] - OME[IX,JY] ) / VOR[IX, JY])
          IF (delta_omega GT ONE) THEN BEGIN
            big_delta = big_delta+1
            PRINT, "VOR[",IX,",",JY,"] = ", VOR[IX,JY], "     ", $ 
            "OME[",IX,",",JY,"] = ", OME[IX,JY]
          ENDIF
        ENDIF ELSE BEGIN
          delta_omega = ZERO
        ENDELSE
        IF (delta_omega GT delta_omega_max) THEN delta_omega_max = delta_omega

     ENDFOR
   ENDFOR
   PRINT, 'delta_omega = ', delta_omega
   PRINT, 'big_delta   = ', big_delta
   PRINT, 'npts_sqd    = ', npts * npts
   PRINT, 'npts        = ', npts
   PRINT, FORMAT = '(A14,E24.18)', 'r_x         = ', r_x
   PRINT, FORMAT = '(A14,E24.18)', 'r_y         = ', r_y
   PRINT, FORMAT = '(A14,E24.18)', 'sigma_x     = ', sigma_x
   PRINT, FORMAT = '(A14,E24.18)', 'sigma_y     = ', sigma_y
   PRINT, 'varsigma_x  = ', varsigma_x
   PRINT, 'varsigma_y  = ', varsigma_y

  x_min        = MIN(X)
  x_max        = MAX(X)
               
  y_min        = MIN(Y)
  y_max        = MAX(Y)
               
  x_range      = [x_min, x_max]
  y_range      = [y_min, y_max]
               
  phi_max      = MAX(PHI)
  phi_min      = MIN(PHI)
               
  phi_range    = [phi_min, phi_max]
               
  vor_max      = MAX(VOR)
  vor_min      = MIN(VOR)
               
  vor_range    = [vor_min, vor_max]

  ome_max      = MAX(OME)
  ome_min      = MIN(OME)

  ome_range    = [ome_min, ome_max]

  psi_max      = MAX(PSI)
  psi_min      = MIN(PSI)

  psi_range    = [psi_min, psi_max]

  sfc_title    = "$\phi(x,y,z)$; $E(z)=1$, $\sigma_x=\sigma_y$ = " + str_sig_x
               
  phi_surf     = SURFACE(PHI, X, Y,                                                        $
                         STYLE        = 2,                                                 $
                         PERSPECTIVE  = 1,                                                 $
                         TITLE        = sfc_title,                                         $
                         DEPTH_CUE    = [0,2],                                             $
                         AXIS_STYLE   = 1,                                                 $
                         FONT_SIZE    = 16,                                                $
                         COLOR        = 'blue',                                            $
                         XRANGE       = x_range,                                           $
                         YRANGE       = y_range,                                           $
                         ZRANGE       = phi_range,                                         $
                         TRANSPARENCY = 00,                                                $
                         CLIP         = 0                                                  $
                       )
               
  ax           = phi_surf.AXES
  ax[0].TITLE  = '$X$'
  ax[1].TITLE  = '$Y$'

  sfc_title    = '$\phi(x,y,z)$ from Fourier transform of $\Omega$'

  psi_surf     = SURFACE(PSI, X, Y,                                                    $
                     STYLE        = 2,                                                 $
                     PERSPECTIVE  = 1,                                                 $
                     TITLE        = sfc_title,                                         $
                     DEPTH_CUE    = [0,2],                                             $
                     AXIS_STYLE   = 1,                                                 $
                     FONT_SIZE    = 16,                                                $
                     COLOR        = 'blue',                                            $
                     XRANGE       = x_range,                                           $
                     YRANGE       = y_range,                                           $
                     ZRANGE       = phi_range,                                         $
                     TRANSPARENCY = 00,                                                $
                     CLIP         = 0                                                  $
                   )

  ax           = psi_surf.AXES
  ax[0].TITLE  = '$X$'
  ax[1].TITLE  = '$Y$'
               
  sfc_title    = "$\Omega(x,y,z)$; $E(z)=1$, $\sigma_x=\sigma_y$ = " + str_sig_x

  vor_surf     = SURFACE(VOR, X, Y, $
                         STYLE        = 2,                                                 $
                         PERSPECTIVE  = 1,                                                 $
                         TITLE        = sfc_title,                                         $
                         DEPTH_CUE    = [0,2],                                             $
                         AXIS_STYLE   = 1,                                                 $
                         FONT_SIZE    = 16,                                                $
                         COLOR        = 'blue',                                            $
                         XRANGE       = x_range,                                           $
                         YRANGE       = y_range,                                           $
                         ZRANGE       = vor_range,                                         $
                         TRANSPARENCY = 00,                                                $
                         CLIP         = 0                                                  $
                       )
                  
  ax           = vor_surf.AXES
  ax[0].TITLE  = '$X$'
  ax[1].TITLE  = '$Y$'

  sfc_title    = '$\Omega(x,y,z)$ from fourier transform of $\phi$'

  ome_surf     = SURFACE(OME, X, Y, $
                         STYLE        = 2,                                                 $
                         PERSPECTIVE  = 1,                                                 $
                         TITLE        = sfc_title,                                         $
                         DEPTH_CUE    = [0,2],                                             $
                         AXIS_STYLE   = 1,                                                 $
                         FONT_SIZE    = 16,                                                $
                         COLOR        = 'blue',                                            $
                         XRANGE       = x_range,                                           $
                         YRANGE       = y_range,                                           $
                         ZRANGE       = vor_range,                                         $
                         TRANSPARENCY = 00,                                                $
                         CLIP         = 0                                                  $
                       )

  ax           = ome_surf.AXES
  ax[0].TITLE  = '$X$'
  ax[1].TITLE  = '$Y$'

  sfc_title    = '$\phi(x,y,z)$ from Fourier transform of $\Omega$'

; psi_surf     = SURFACE(PSI, X, Y,                                                    $
;                    STYLE        = 2,                                                 $
;                    PERSPECTIVE  = 1,                                                 $
;                    TITLE        = sfc_title,                                         $
;                    DEPTH_CUE    = [0,2],                                             $
;                    AXIS_STYLE   = 1,                                                 $
;                    FONT_SIZE    = 16,                                                $
;                    COLOR        = 'blue',                                            $
;                    XRANGE       = x_range,                                           $
;                    YRANGE       = y_range,                                           $
;                    ZRANGE       = psi_range,                                         $
;                    TRANSPARENCY = 00,                                                $
;                    CLIP         = 1                                                  $
;                  )

  ax           = psi_surf.AXES
  ax[0].TITLE  = '$X$'
  ax[1].TITLE  = '$Y$'

END
