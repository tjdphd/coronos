PRO FFT

ZERO             = 0.0D                            ; Define some constants for convenience
CZERO            = DCOMPLEX(ZERO,ZERO)             ; and to make sure that DOUBLE's are used
ONE              = 1.0D                            ; everywhere
TWO              = 2.0D
THREE            = 3.0D
PI               = 3.141592653589793D
TWO_THIRDS       = TWO / THREE

TWO_PI           = TWO*PI

; Defining the Gaussian in Cartesian space

n1               = 64                             ; number of data samples along x
n2               = 64                             ; number of data samples along y

X                = DBLARR(n1)                      ; x-coordinates
Y                = DBLARR(n2)                      ; y-coordinates

dx               = ONE/DOUBLE(n1)                  ; x-increment
dy               = ONE/DOUBLE(n2)                  ; y-increment

X                = DINDGEN(n1) * dx                ; range in x
Y                = DINDGEN(n2) * dy                ; range in y

U                = DBLARR(n1,n2)                   ; rank 2 array for Gaussian in Cartesian space

sigma_x          = 0.0325D                         ; Gaussian width in x
sigma_y          = 0.0325D                         ; Gaussian width in y

 x_0             = 0.0D                            ; x-coordinate of maximum U
 y_0             = 0.0D                            ; y0coordinate of mininum U

;sigma_x       = ONE / SQRT(TWO_PI)                ; Alternative witdth
;sigma_y       = ONE / SQRT(TWO_PI)                ; choices

sigma_x_sqd      = sigma_x * sigma_x
sigma_y_sqd      = sigma_y * sigma_y

FOR I = 0, n1 - 1 DO BEGIN                         ; Initialize the Gaussian in
  FOR J = 0, n2 - 1 DO BEGIN                       ; Cartesian space

    del_x_sqd    = (X[I] - x_0) * (X[I]-x_0)
    del_y_sqd    = (Y[J] - y_0) * (Y[J]-y_0)

    U[I,J]       = 1.0D*EXP(-del_x_sqd/(TWO*sigma_x_sqd))*EXP(-del_y_sqd/(TWO*sigma_y_sqd))

  ENDFOR
ENDFOR

; Now define the wave numbers kx and ky
; Notice that X and Y range in [0,1]
; this means that the "L" defined on the
; matlab webpage is L=0.5 and that
; we require the multiplicative value of
; TWO_PI when defining out K vectors

Xf               = DINDGEN((n1-1)/2) + ONE
is_n1_even       = (n1 MOD 2) EQ 0
IF (is_n1_even) THEN KX = TWO_PI * [ ZERO, Xf,   DOUBLE(n1)/TWO, DOUBLE(-n1)/TWO  + Xf]/ ( DOUBLE(n1)*dx) $
ELSE                 KX = TWO_PI * [ ZERO, Xf, -(DOUBLE(n1)/TWO + ONE) + Xf ] / (DOUBLE(n1)*dx)
Yf               = DINDGEN((n2-1)/2) + ONE
is_n2_even       = (n2 MOD 2) EQ 0
IF (is_n2_even) THEN KY = TWO_PI * [ ZERO, Yf,   DOUBLE(n2)/TWO, DOUBLE(-n2)/TWO  + Yf]/ ( DOUBLE(n2)*dx) $
ELSE                 KY = TWO_PI * [ ZERO, Yf, -(DOUBLE(n2)/TWO + ONE) + Yf ] / (DOUBLE(n2)*dx)

FTU              = FFT(U)                          ; FTU now contains the Fourier components
                                                   ; of U. It is a DCOMPLEX array of rank two
                                                   ; and dimensional sizes of n1 and n2 - i.e.
                                                   ; FTU contains the same number of complex
                                                   ; numbers as U contains reals. 

x_max            = MAX(X)                          ; If X and Y were defined to range over
y_max            = MAX(Y)                          ; [0,x_max] and [0,y_max] for example I 
                                                   ; would have used x_max and y_max to normalize
x_min            = MIN(X)                          ; the X and Y vectors to [0,1].
y_min            = MIN(Y)

x_range          = [x_min, x_max]                  ; Then if, for example x_max = y_max = 10R_J
y_range          = [y_min, y_max]                  ; where R_J is a Jovian radius, my K vectors
                                                   ; as defined above would be in units of (10R_J)^(-1)
kx_max           = MAX(KX)           
ky_max           = MAX(KY)                         ; The other thing to be understood about what's
                                                   ; happening here is that I'm establishing the plotting
kx_min           = MIN(KX)                         ; ranges in Cartesian space and Fourier space to
ky_min           = MIN(KY)                         ; be used later in the contour plots

kx_range         = [kx_min, kx_max]
ky_range         = [ky_min, ky_max]

IF ABS(kx_min) LT kx_max THEN kx_max = ABS(kx_min) ; This ensures a symmetric range in KX and KY
IF ABS(ky_min) LT ky_max THEN ky_max = ABS(ky_min) ; and helps to enforce the de-aliasing that comes next

k_sqd_max        = kx_max*kx_max + ky_max*ky_max   ; Now I'm getting ready to de-alias the DFT using
;threshold        = 1.0D*TWO_THIRDS * k_sqd_max     ; the "2/3" method which is similar to but different from
                                                   ; the discussion on the matlab website


IF ABS(kx_max) LE ABS(ky_max) THEN k_max = kx_max
IF ABS(ky_max) LE ABS(kx_max) THEN k_max = ky_max

threshold        = TWO_THIRDS * k_max

PRINT, 'threshold = ', threshold
PRINT, 'k_max     = ', k_max

 FOR L = 0, n1 - 1 DO BEGIN
   FOR M = 0, n2 - 1 DO BEGIN

;   k_sqd         = KX[L]*KX[L] + KY[M]*KY[M]       ; the DFT components corresponding to wave numbers larger
;   IF (k_sqd GT threshold ) THEN FTU[L,M] = CZERO  ; than 2/3 of the maximum are set to complex zero

    IF ( ABS(KY[M]) EQ ZERO      AND ABS(KX[L])  EQ k_max   ) THEN FTU[L,M] = CZERO
    IF ( ABS(KY[M]) GT threshold OR  ABS(KX[L]) GT threshold) THEN FTU[L,M] = CZERO   

   ENDFOR
 ENDFOR

;FOR L = 0, n1 - 1 DO BEGIN                         ; This is for diagnostic purposes and can
;   FOR M = 0, n2 - 1 DO BEGIN                      ; be commented out
;
;     IF ( KX[L] GE ZERO ) THEN BEGIN
;
;     rpart       =  REAL_PART(FTU[L,M])
;     ipart       =  IMAGINARY(FTU[L,M])
;
;     PRINT, FORMAT = '(A9,F7.1,A7,F7.1,A5,F8.5,A1,F8.5,A1,2x,A4, I4,A6,I4)', $
;      'FTU[kx = ', KX[L]/TWO_PI, ', ky = ', KY[M]/TWO_PI, '] = (',rpart, ',', ipart,')', 'L = ', L,', M = ', M
;
;     ENDIF
;  ENDFOR
;ENDFOR

kx_sort          = SORT(KX)                        ; Recall that the ordering of KX and KY is not
ky_sort          = SORT(KY)                        ; monotonically ascending. To make sure the CONTOUR
                                                   ; function behaves properly we have to re-order
KX[*]            = KX[kx_sort]                     ; KX, KY, and ...
KY[*]            = KY[ky_sort]

size_ftu         = SIZE(FTU)
type_ftu         = TYPENAME(FTU)

FTU_SORT         = DCOMPLEXARR(n1,n2)              ; ...FTU

FOR L = 0, n1 - 1 DO BEGIN                         ; This loop ensures that the correct value of
   FOR M = 0, n2 - 1 DO BEGIN                      ; FTU is associated with the correct values
                                                   ; of KX and KY
   FTU_SORT[L,M] = FTU[kx_sort[L],ky_sort[M]]

  ENDFOR
ENDFOR

FTU              = FTU_SORT                        ; Now we can set FTU to its sorted value



MFTU             = SQRT(REAL_PART(FTU)^2 + IMAGINARY(FTU)^2)    ; Calculating the moduli of the FTU[L,M]

u_max            = MAX(U)                                       ; I need max's and min's of the various
u_min            = MIN(U)                                       ; arrays I want to plot. I also need
 
rftu_max         = MAX(REAL_PART(FTU))                          ; to know what range of values my contours
rftu_min         = MIN(REAL_PART(FTU))                          ; to lie in. I'm working those things out now.

iftu_max         = MAX(IMAGINARY(FTU))
iftu_min         = MIN(IMAGINARY(FTU))

mftu_max         = MAX(MFTU)
mftu_min         = MIN(MFTU)

n_ctrs           = 31                                           ; The nominal number of contours to be used
inc_uctrs        = (u_max    - u_min)    / DOUBLE(n_ctrs-1)     ; Contour increments for U
inc_rftuctrs     = (rftu_max - rftu_min) / DOUBLE(n_ctrs-1)     ; Contour increments for real part of FTU
inc_iftuctrs     = (iftu_max - iftu_min) / DOUBLE(n_ctrs-1)     ; Contour increments for imaginary part of FTU
inc_mftuctrs     = (mftu_max - mftu_min) / DOUBLE(n_ctrs-1)     ; Contour increments for modulus of FTU

UCTRS            = DBLARR(n_ctrs)                               ; array for U contour values
RFTUCTRS         = DBLARR(n_ctrs)                               ; array for real part of FTU contour values
IFTUCTRS         = DBLARR(n_ctrs)                               ; array for imaginary part of FTU contour values
MFTUCTRS         = DBLARR(n_ctrs)                               ; array for modulus of FTU contour values

FOR K = 0, n_ctrs - 1 DO BEGIN 

    UCTRS[K]     = u_min    + (K * inc_uctrs   )                ; Fill the above arrays
 RFTUCTRS[K]     = rftu_min + (K * inc_rftuctrs)                ; appropriately
 IFTUCTRS[K]     = iftu_min + (K * inc_iftuctrs) 
 MFTUCTRS[K]     = Mftu_min + (K * inc_Mftuctrs) 

ENDFOR

 UCNTR            = CONTOUR( U, X, Y,                                  $ ; Contour plot of U on X Y plane
                             RGB_TABLE = 6,                            $
                             C_VALUE   = UCTRS,                        $
                             XRANGE    = x_range,                      $
                             YRANGE    = y_range,                      $
                             TITLE     = 'U(X,Y)',                     $
                             XTITLE    = 'X',                          $
                             YTITLE    = 'Y'                           $
                           )
 
 
RFTUCNTR         = CONTOUR( REAL_PART(FTU), KX, KY,                   $ ; Contour plot of real part of FTU on Kx Ky plane
                            RGB_TABLE = 6,                            $
                            C_VALUE   = RFTUCTRS[2:n_ctrs-1],         $
                            XRANGE    = kx_range,                     $
                            YRANGE    = ky_range,                     $
                            TITLE     = 'Re[FTU(K!Dx!N,K!Dy!N)]',     $
                            XTITLE    = 'K!Dx!N',                     $
                            YTITLE    = 'K!Dy!N'                      $
                          )

IFTUCNTR         = CONTOUR( IMAGINARY(FTU), KX, KY,                   $ ; Contour plot of imaginary part of FTU on Kx Ky plane
                            RGB_TABLE = 6,                            $
                            C_VALUE   = IFTUCTRS[2:n_ctrs-1],         $
                            XRANGE    = kx_range,                     $
                            YRANGE    = ky_range,                     $
                            TITLE     = 'Im[FTU(K!Dx!N,K!Dy!N)]',     $
                            XTITLE    = 'K!Dx!N',                     $
                            YTITLE    = 'K!Dy!N'                      $
                          )
MFTUCNTR         = CONTOUR( MFTU, KX, KY,                             $ ; Contour plot of modulus of FTU on Kx Ky plane
                            RGB_TABLE = 6,                            $
                            C_VALUE   = MFTUCTRS[2:n_ctrs-1],         $
                            XRANGE    = kx_range,                     $
                            YRANGE    = ky_range,                     $
                            TITLE     = '|FTU(K!Dx!N,K!Dy!N)|',       $
                            XTITLE    = 'K!Dx!N',                     $
                            YTITLE    = 'K!Dy!N'                      $
                          )

;USF              = SURFACE(U, X, Y,                           $ ; alternatively, or additionally one
;                           AXIS_STYLE = 2,                    $ ; could plot these as surfaces
;                           STYLE      = 2,                    $
;                           DEPTH_CUE  = [0,2],                $
;                           COLOR      = 'green',              $
;                           XRANGE     = x_range,              $
;                           YRANGE     = y_range               $
;                         )

;RFTUSF           = SURFACE(REAL_PART(FTU), KX, KY,            $
;                           AXIS_STYLE = 2,                    $
;                           STYLE      = 2,                    $
;                           DEPTH_CUE  = [0,2],                $
;                           COLOR      = 'red',                $
;                           XRANGE     = kx_range,             $
;                           YRANGE     = ky_range              $
;                          )

;IFTUSF           = SURFACE(IMAGINARY(FTU), KX, KY,            $
;                           AXIS_STYLE = 2,                    $
;                           STYLE      = 2,                    $
;                           DEPTH_CUE  = [0,2],                $
;                           COLOR      = 'red',                $
;                           XRANGE     = kx_range,             $
;                           YRANGE     = ky_range              $
;                          )

;MFTUSF           = SURFACE(MFTU,           KX, KY,            $
;                           AXIS_STYLE = 2,                    $
;                           STYLE      = 2,                    $
;                           DEPTH_CUE  = [0,2],                $
;                           COLOR      = 'red',                $
;                           XRANGE     = kx_range,             $
;                           YRANGE     = ky_range              $
;                          )

END
