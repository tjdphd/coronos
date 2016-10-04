FUNCTION min_q, n_slice, qty, first_step, last_step, desc_label
COMMON fnc, prefix, inc_res, nfld

CASE qty OF
'vz': i_qty      =  5
'bz': i_qty      =  4
'a' : i_qty      =  3
'p' : i_qty      =  2
'j' : i_qty      =  0
ELSE: i_qty      = -1
ENDCASE

IF (i_qty EQ -1) THEN BEGIN

      PRINT, 'min_q: WARNING - value ', qty, ' for the variable qty not recognized.'
      PRINT, 'Setting indext i_qty to default value for contours of the flux function a.'

      i_qty      = 3
ENDIF

min_q_in         = 1.0e+100
min_q_out        = 1.0e+100

n3               = scan_parameters('n3',  0, desc_label )
ip1              = scan_parameters('ip1', 0, desc_label )
ip2              = scan_parameters('ip2', 0, desc_label )

x_res            = 2^ip1
y_res            = 2^ip2   

i_x_res          = LONG64(x_res)
i_y_res          = LONG64(y_res)

FOR I = first_step, last_step DO BEGIN

   Slices        = fetch_slices(n_slice, n_slice, I, desc_label)
   size_slices   = SIZE(Slices)   

   IF (size_slices[0] EQ 3) THEN nlines = size_slices[2]
   IF (size_slices[0] EQ 2) THEN nlines = size_slices[1]

   IF (nlines NE (i_x_res * i_y_res)) THEN BEGIN

      PRINT, FORMAT = '(a31,1x,i4,1x,a11,1x,i4,1x,a18)', 'min_a: WARNING - skipping slice', $
                                                          n_slice, 'for sub-run', I,        $
                                                         'data not available'
      CONTINUE
   ENDIF ELSE BEGIN

             IF (qty NE 'j') THEN BEGIN
                IF (n3 GT 1) THEN BEGIN
                        min_q_in = MIN(Slices[0, *, i_qty] )
                ENDIF ELSE BEGIN
                        min_q_in = MIN(Slices[   *, i_qty] )
                      ENDELSE
             ENDIF ELSE BEGIN
                       J_z       = calc_jz(Slices, n_slice, I, desc_label )
                       min_q_in  = MIN(J_z)
                   ENDELSE
            IF min_q_in LT min_q_out THEN min_q_out = min_q_in
         ENDELSE
  ENDFOR
  RETURN, min_q_out
END 
