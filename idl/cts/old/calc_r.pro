FUNCTION calc_r, Slices, n_slice, n_step, desc_label

COMMON Save_r, Pre
COMMON First, first_call
COMMON fnc, prefix, inc_res, nfld
COMMON loc, eps_out_dir

; Pre        = FLTARR(1,7)
; Pre        = TRANSPOSE(Pre)
  
frac         = 0.07
symmetric    = 0

ip1          = scan_parameters('ip1', 0, desc_label )      ; power of 2 giving x-resolution
ip2          = scan_parameters('ip2', 0, desc_label )      ; power of 2 giving y-resolution
n3           = scan_parameters('n3' , 0, desc_label )      ; number of slices per data file
mp           = scan_parameters('mp' , 0, desc_label )      ; number of processors used in run
zl           = scan_parameters('zl' , 0, desc_label )      ; total height along z of integration volume
 
x_res        = 2^ip1                                       ; resolution in x
y_res        = 2^ip2                                       ; resolution in y
z_res        = n3 * mp                                     ; resolution in z
 
nx           = LONG64(x_res - 1)                           ; integer resolutions for looping
ny           = LONG64(y_res - 1)
nz           = LONG64(z_res - 1)  

fst_step = n_step                                          ; for calls requiring these quantities
lst_step = n_step                                          ; to insure just the n_step'th slice is processed

                                                           ; Slices is an array containing n3 X x_res X y_res
slices_size = SIZE(Slices)                                 ; records each of which has data for x, y, a, and p
rank = slices_size[0]                                      ; where x and y are the point coordinates within one
Slc  = FLTARR(slices_size[rank], slices_size[rank - 1], 2) ; of the n3 slices. 'slc_size' is a 1-D array whose
Slc  = TRANSPOSE(Slc)                                      ; elements 0, 1, 2 and 3 respectively are, the
IF (n3 GT 1) THEN BEGIN                                    ; rank of Slices, the number of slices, the number
   Slc[0,*,*] = Slices[0,*,*]                              ; records per slice and the number of fields per
ENDIF ELSE BEGIN                                           ; record. ( unless n3 = 1 in which case elements
   Slc[0,*,*] = Slices[*,*]                                ; 0, 1, and 2, are the rank, the number of records
      ENDELSE                                              ; and the fields per record.)
                                                           
slc_size      = SIZE(Slc)                                  
n_lins        = slices_size[rank-1]                        ; for ease of use we name this element of slc_size

hlf_lin        = n_lins/2                                  ; Used to section the array Slc and search for maxima
                                                           ; of 'a' separately in the left- and right-hand halves of
                                                           ; the domain.
AUX            = FLTARR(4, n_lins)
AUX            = TRANSPOSE(AUX)
 
i_skip         = LONG64(0)
i_aux          = LONG64(-1)

IF (first_call EQ 0) THEN BEGIN
      Pre[*]   = 0  
    r_last     = Pre[0]
 
    min_a_l    = MIN(Slc[0,0:hlf_lin,3],SUBSCRIPT_MAX=i_max_a) ; left-hand min and array index in Slc of left-hand max_a
    max_a_l    = MAX(Slc[0,0:hlf_lin,3],SUBSCRIPT_MIN=i_min_a) ; left-hand max and array index in Slc of right_hand min_a
    
    x_a_max_l  = Slc[0,i_max_a, 0]                             ; x and y coordinates of left-hand a_max
    y_a_max_l  = Slc[0,i_max_a, 1]

    x_a_min_l  = Slc[0,i_min_a, 0]                             ; x and y coordinates of left-hand a_min
    y_a_min_l  = Slc[0,i_min_a, 1]

    min_a_r    = MIN(Slc[0,hlf_lin + 1:n_lins - 1,3], SUBSCRIPT_MAX=i_max_a)   ; Now we start over, this time 
    max_a_r    = MAX(Slc[0,hlf_lin + 1:n_lins - 1,3], SUBSCRIPT_MIN=i_min_a)   ; in the right-hand plane


    i_max_a    = i_max_a + hlf_lin + 1                         ; the indexes i_max_and i_min_a  are returned by
    i_min_a    = i_min_a + hlf_lin + 1                         ; MIN and MAX (go figure) relative to hlf_lin. So
                                                               ; we have to add hlf_lin back in before using the
                                                               ; index values in slice

    x_a_max_r  = Slc[0,i_max_a, 0]                             ; Now we can use the indexes to fetch the x and y
    y_a_max_r  = Slc[0,i_max_a, 1]                             ; coordinates of a_max and a_min from Slc.


    x_a_min_r  = Slc[0,i_min_a, 0]                                           
    y_a_min_r  = Slc[0,i_min_a, 1]
    
    B_T_Slc    = FLTARR(4, n_lins)                             ; new array for rearranged slice
    B_T_Slc    = TRANSPOSE(B_T_Slc)
    
    cur_idx    = LONG64(0)

    FOR J = 0, nx DO BEGIN
       FOR I = 0L, LONG64(n_lins - 1), (nx + 1) DO BEGIN       ; re-order data in Slc such that the
         B_T_Slc[cur_idx, *] = Slc[0,I+J, *]                   ; data changes most rapidly along x rather
         cur_idx             = cur_idx + 1                     ; than along y. 
      ENDFOR
    ENDFOR 

    min_a_b   = MIN(B_T_Slc[0:hlf_lin,3],SUBSCRIPT_MAX=i_max_a) ; bottom-half min and array index in 
                                                                ; Slc of bottom-half max_a
    max_a_b   = MAX(B_T_Slc[0:hlf_lin,3],SUBSCRIPT_MIN=i_min_a) ; bottom-half max and array index in 
                                                                ; Slc of right_hand min_a
           
    x_a_max_b = B_T_Slc[i_max_a, 0]                             ; x and y coordinates of left-hand a_max
    y_a_max_b = B_T_Slc[i_max_a, 1]
    
    x_a_min_b = B_T_Slc[i_min_a, 0]                             ; x and y coordinates of left-hand a_min
    y_a_min_b = B_T_Slc[i_min_a, 1]
    
    min_a_t   = MIN(B_T_Slc[hlf_lin + 1:(n_lins-1),3], SUBSCRIPT_MAX=i_max_a) ; Top-half min and array index for max
    max_a_t   = MAX(B_T_Slc[hlf_lin + 1:(n_lins-1),3], SUBSCRIPT_MIN=i_min_a) ; Top-half max and array index for min
    
    i_max_a   = i_max_a + hlf_lin + 1                           ; the indexes i_max_and i_min_a  are returned by
    i_min_a   = i_min_a + hlf_lin + 1                           ; MIN and MAX (go figure) relative to hlf_lin. So
                                                                ; we have to add hlf_lin back in before using the
                                                                ; index values in slice

    x_a_max_t = B_T_Slc[i_max_a, 0]                             ; Now we can use the indexes to fetch the x and y
    y_a_max_t = B_T_Slc[i_max_a, 1]                             ; coordinates of a_max and a_min from Slc.
    
    x_a_min_b = B_T_Slc[i_min_a, 0]                                           
    y_a_min_b = B_T_Slc[i_min_a, 1]

ENDIF ELSE BEGIN
         IF (symmetric EQ 0) THEN BEGIN
   
         r_last = Pre[0]
         x_1    = Pre[2]
         y_1    = Pre[3]
         x_2    = Pre[5]
         y_2    = Pre[6]

         
         FOR I=0L, LONG64(n_lins - 1) DO BEGIN        
            x      = Slc[0,I,0]
            y      = Slc[0,I,1]  
            r_test_1 = sqrt( (x - x_1)^2 + (y - y_1)^2 )
            r_test_2 = sqrt( (x - x_2)^2 + (y - y_2)^2 )
            x_test_1 = abs(x - x_1)
            x_test_2 = abs(x - x_2)
            y_test_1 = abs(y - y_1)
            y_test_2 = abs(y - y_2)
            reject = 0
            IF ( x_test_1 LE (frac * r_last) ) THEN BEGIN
               IF (y_test_1 LE (frac * r_last)) THEN BEGIN
                  IF (r_test_1 LE (frac * r_last)) THEN BEGIN
                    i_aux = i_aux + 1
                    AUX[i_aux, *] = Slc[0,I, *]
                  ENDIF ELSE BEGIN
                             reject = reject + 1
                        ENDELSE
               ENDIF ELSE BEGIN
                             reject = reject + 1
                     ENDELSE
            ENDIF ELSE BEGIN
                             reject = reject + 1
                  ENDELSE
            IF (x_test_2 LE (frac * r_last) ) THEN BEGIN
               IF (y_test_2 LE (frac * r_last) ) THEN BEGIN  
                  IF (r_test_2 LE (frac * r_last)) THEN BEGIN
                             IF (reject EQ 0) THEN BEGIN
                                PRINT, 'WARNING test_1: search regions overlapping'
                             ENDIF
                                i_aux = i_aux + 1
                                AUX[i_aux, *] = Slc[0,I,*]
                  ENDIF ELSE BEGIN
                             reject = reject + 1
                        ENDELSE
               ENDIF ELSE BEGIN
                             reject = reject + 1
                     ENDELSE
            ENDIF ELSE BEGIN    
                             reject = reject + 1
                  ENDELSE
            IF (reject EQ 2) THEN BEGIN
               i_skip = i_skip + 1
            ENDIF
         ENDFOR
       
         AUX[i_aux + 1: n_lins-1, 3 ] = -1.0e-20
         aux_order      = SORT(AUX[*,3], /L64)
         size_aux_order = SIZE(aux_order)

         rev_aux_order  = FLTARR(1, LONG64(size_aux_order[1]))
         rev_aux_order  = TRANSPOSE(rev_aux_order)
        
         rev_idx = LONG64(size_aux_order[1])
         FOR I=0L, LONG64(size_aux_order[1] - 1) DO BEGIN
             rev_idx = LONG64(rev_idx - 1)
             rev_aux_order[rev_idx] = aux_order[I]
         ENDFOR
         aux_order = rev_aux_order
             
         n_lins    = LONG64(n_lins - i_skip)
         WHILE ( n_lins LE 3) DO BEGIN
            n_lins = n_lins + 1
         ENDWHILE
      
         hlf_lin   = LONG64(n_lins/2)
        
         L_R_AUX   = FLTARR(4, n_lins)                   ; new array for rearranged slice
         L_R_AUX   = TRANSPOSE(L_R_AUX)
         
         O_L_R_AUX = FLTARR(4, n_lins)                   ; new array for rearranged slice
         O_L_R_AUX = TRANSPOSE(O_L_R_AUX)
         
         B_T_AUX   = FLTARR(4, n_lins)
         B_T_AUX   = TRANSPOSE(B_T_AUX)
         
         O_B_T_AUX = FLTARR(4, n_lins)
         O_B_T_AUX = TRANSPOSE(O_B_T_AUX)
         
         
         FOR I = 0L, LONG64(n_lins - 1) DO BEGIN
            J  = aux_order[I]
            L_R_AUX[I, *] =  AUX[J,*]
         ENDFOR
        
         x_order      = SORT(L_R_AUX[*,0], /L64)
        
         FOR I = 0L, LONG64(n_lins - 1) DO BEGIN
            J               = x_order[I]
            O_L_R_AUX[I, *] = L_R_AUX[J,*]
         ENDFOR
       
         min_a_l    = MIN(O_L_R_AUX[0:hlf_lin,3],SUBSCRIPT_MAX=i_max_a)  ; left-hand min and array index in Slc of left-hand max_a
         max_a_l    = MAX(O_L_R_AUX[0:hlf_lin,3],SUBSCRIPT_MIN=i_min_a)  ; left-hand max and array index in Slc of right_hand min_a
            
         x_a_max_l  = O_L_R_AUX[i_max_a, 0]                              ; x and y coordinates of left-hand a_max
         y_a_max_l  = O_L_R_AUX[i_max_a, 1]

         x_a_min_l  = O_L_R_AUX[i_min_a, 0]                              ; x and y coordinates of left-hand a_min
         y_a_min_l  = O_L_R_AUX[i_min_a, 1]

         line_span  = i_max_a - i_min_a + 1

         min_a_r    = MIN(O_L_R_AUX[(hlf_lin + 1):n_lins-1, 3], SUBSCRIPT_MAX=i_max_a)
         max_a_r    = MAX(O_L_R_AUX[(hlf_lin + 1):n_lins-1, 3], SUBSCRIPT_MIN=i_min_a)

         i_max_a    = i_max_a + (hlf_lin + 1)                     ; the indexes i_max_and i_min_a  are returned by
         i_min_a    = i_min_a + (hlf_lin + 1)                     ; MIN and MAX (go figure) relative to hlf_lin. So
                                                                  ; we have to add hlf_lin back in before using the
                                                                  ; index values in slice

         x_a_max_r  = O_L_R_AUX[i_max_a, 0]                       ; Now we can use the indexes to fetch the x and y
         y_a_max_r  = O_L_R_AUX[i_max_a, 1]                       ; coordinates of a_max and a_min from Slc.


         x_a_min_r  = O_L_R_AUX[i_min_a, 0]                                           
         y_a_min_r  = O_L_R_AUX[i_min_a, 1]
         
         O_B_T_AUX   = FLTARR(4, n_lins)                          ; new array for rearranged slice
         O_B_T_AUX   = TRANSPOSE(O_B_T_AUX)
         
         x_cnt     = LONG64(0)
         cur_idx   = LONG64(-1)
         
         y_order = SORT(L_R_AUX[*, 1], /L64)
            
         FOR I = 0L, LONG64(n_lins - 1) DO BEGIN
            J  = y_order[I]
            O_B_T_AUX[I, *] = L_R_AUX[J, *]
         ENDFOR

        min_a_b   = MIN(O_B_T_AUX[0:hlf_lin,3],SUBSCRIPT_MAX=i_max_a) ; bottom-half min and array index in 
                                                                      ; Slc of bottom-half max_a
        max_a_b   = MAX(O_B_T_AUX[0:hlf_lin,3],SUBSCRIPT_MIN=i_min_a) ; bottom-half max and array index in 
                                                                      ; Slc of right_hand min_a
    
        x_a_max_b = O_B_T_AUX[i_max_a, 0]                             ; x and y coordinates of left-hand a_max
        y_a_max_b = O_B_T_AUX[i_max_a, 1]
   
        x_a_min_b = O_B_T_AUX[i_min_a, 0]                             ; x and y coordinates of left-hand a_min
        y_a_min_b = O_B_T_AUX[i_min_a, 1]
   
        min_a_t   = MIN(O_B_T_AUX[hlf_lin + 1:(n_lins-1),3], SUBSCRIPT_MAX=i_max_a) ; Top-half min and array index for max
        max_a_t   = MAX(O_B_T_AUX[hlf_lin + 1:(n_lins-1),3], SUBSCRIPT_MIN=i_min_a) ; Top-half max and array index for min
   
        i_max_a    = i_max_a + hlf_lin + 1                            ; the indexes i_max_and i_min_a  are returned by
        i_min_a    = i_min_a + hlf_lin + 1                            ; MIN and MAX (go figure) relative to hlf_lin. So
                                                                      ; we have to add hlf_lin back in before using the
                                                                      ; index values in slice

        x_a_max_t  = O_B_T_AUX[i_max_a, 0]                            ; Now we can use the indexes to fetch the x and y
        y_a_max_t  = O_B_T_AUX[i_max_a, 1]                            ; coordinates of a_max and a_min from Slc.
   
        x_a_min_b  = O_B_T_AUX[i_min_a, 0]                                           
        y_a_min_b  = O_B_T_AUX[i_min_a, 1]
         
         ENDIF ELSE BEGIN ; symmetric = 0 if
      
                    x_a_max_r  = 1.0 - x_a_max_l
                    y_a_max_r  = 1.0 - y_a_max_l

                    x_a_max_b  = x_a_max_l
                    y_a_max_b  = y_a_max_l

                    x_a_max_t  = x_a_max_r
                    y_a_max_t  = y_a_max_r

                ENDELSE

         ENDELSE    

      r_lr        = SQRT( (x_a_max_l - x_a_max_r)^2 + (y_a_max_l - y_a_max_r)^2 ) ; r as determined from left-right split
      r_bt        = SQRT( (x_a_max_t - x_a_max_b)^2 + (y_a_max_t - y_a_max_b)^2 ) ; r as determined from bottom-top split
      

IF (r_lr GE r_bt) THEN BEGIN

   r              = r_lr
   max_a_1        = max_a_l
   x_a_max_1      = x_a_max_l
   y_a_max_1      = y_a_max_l
   max_a_2        = max_a_r
   x_a_max_2      = x_a_max_r
   y_a_max_2      = y_a_max_r

ENDIF ELSE BEGIN

           r         = r_bt
           max_a_1   = max_a_b
           x_a_max_1 = x_a_max_b
           y_a_max_1 = y_a_max_b
           max_a_2   = max_a_t
           x_a_max_2 = x_a_max_t
           y_a_max_2 = y_a_max_t

      ENDELSE

RAXYAXY    = FLTARR(1,7)
RAXYAXY    = TRANSPOSE(RAXYAXY)

RAXYAXY[0] = r
RAXYAXY[1] = max_a_1
RAXYAXY[2] = x_a_max_1
RAXYAXY[3] = y_a_max_1
RAXYAXY[4] = max_a_2
RAXYAXY[5] = x_a_max_2
RAXYAXY[6] = y_a_max_2

Pre[*]     = RAXYAXY[*]

IF (symmetric EQ 0 ) THEN first_call = 1

RETURN, RAXYAXY

END
