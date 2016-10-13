FUNCTION global_q_minmax, qty, first_step, last_step, first_slice, last_slice, desc_label, KK, $
                          GLOBAL_IN_MINMAX = global_in_minmax

  i_x                                          = 6
  i_j                                          = 5
  i_o                                          = 4
  i_a                                          = 3
  i_p                                          = 2

  CASE qty OF
  'x' : i_qty                                  =  i_x
  'j' : i_qty                                  =  i_j
  'o' : i_qty                                  =  i_o
  'a' : i_qty                                  =  i_a
  'p' : i_qty                                  =  i_p
  ELSE: i_qty                                  = -1
  ENDCASE

  IF (i_qty EQ -1) THEN BEGIN
  
    PRINT, 'max_q: WARNING - value ', qty, ' for the variable qty not recognized.'
    PRINT, 'Setting indext i_qty to default value for contours of the current density j.'
  
    i_qty                                      = 5

  ENDIF
  
  ; Prepare to multi-thread here

  total_steps                                  = last_step - first_step + 1
  
  IF (i_qty LT i_x) THEN BEGIN
    glbl_q_minmax_by_slice                     = FLTARR(4,   total_steps)
  ENDIF ELSE BEGIN
    glbl_q_minmax_by_slice                     = FLTARR(4,4, total_steps)
  ENDELSE

  slice_range                                  = getSliceRange(first_slice, last_slice)
  size_sr                                      = SIZE(slice_range, /DIMENSIONS)
  n_threads                                    = size_sr[0]

  IF (n_threads GT 1) THEN oBridge = OBJARR(n_threads - 1)

  IF (i_qty LT i_x) THEN BEGIN

    thread_q_minmax_by_slice                   = FLTARR(4)
    all_threads_q_minmax_by_slice              = FLTARR(4, n_threads)

  ENDIF ELSE BEGIN

    thread_q_minmax_by_slice                   = FLTARR(4,4)
    all_threads_q_minmax_by_slice              = FLTARR(4,4,n_threads)

  ENDELSE

  FOR I = first_step, last_step  DO BEGIN
  
    ndx                                        = I - first_step

    ; start multi-threading here.
    
    FOR J = 0, n_threads - 1 DO BEGIN
   
      IF (J EQ  n_threads - 1) THEN BEGIN

        thread_q_minmax_by_slice               = global_q_minmax_by_slice(I, i_qty, slice_range[J,0], slice_range[J,1], desc_label)

      ENDIF ELSE BEGIN

        oBridge[J]                             = OBJ_NEW('IDL_IDLBRIDGE')
        oBridge[J]                             -> SetVar, 'step_idx', I
        oBridge[J]                             -> SetVar, 'i_qty', i_qty
        oBridge[J]                             -> SetVar, 'start_slice', slice_range[J,0] 
        oBridge[J]                             -> SetVar, 'stop_slice',  slice_range[J,1] 
        oBridge[J]                             -> SetVar, 'desc_label',  desc_label

        IF (i_qty LT i_x) THEN BEGIN
          oBridge[J]                           -> Execute, "thread_q_minmax_by_slice = FLTARR(4)"
        ENDIF ELSE BEGIN
          oBridge[J]                           -> Execute, "thread_q_minmax_by_slice = FLTARR(4,4)"
        ENDELSE
        oBridge[J]                             -> Execute, $
        "thread_q_minmax_by_slice = global_q_minmax_by_slice(step_idx, i_qty, start_slice, stop_slice, desc_label)", /nowait

      ENDELSE
    ENDFOR
    notdone                                    = 1
    IF (n_threads GT 1) THEN BEGIN
      WHILE notdone DO BEGIN
        done                                   = 0
        FOR J = 0, n_threads - 2 DO done       = done + oBridge[J] -> Status()
        IF (done EQ 0) THEN notdone            = done
      ENDWHILE
    ENDIF

    FOR J = 0, n_threads - 1 DO BEGIN
      IF (J EQ n_threads - 1) THEN BEGIN
        IF (i_qty LT i_x) THEN BEGIN
          all_threads_q_minmax_by_slice[*,J]   = thread_q_minmax_by_slice[*]
        ENDIF ELSE BEGIN
          all_threads_q_minmax_by_slice[*,*,J] = thread_q_minmax_by_slice[*,*] ; second index is i_qty
        ENDELSE
      ENDIF ELSE BEGIN
        IF (i_qty LT i_x) THEN BEGIN
          all_threads_q_minmax_by_slice[*,J]   = oBridge[J] -> GetVar('thread_q_minmax_by_slice[*]')
        ENDIF ELSE BEGIN
          all_threads_q_minmax_by_slice[*,*,J] = oBridge[J] -> GetVar('thread_q_minmax_by_slice[*,*]')
        ENDELSE
        obj_destroy, oBridge[J]
      ENDELSE
    ENDFOR

    ; end multi-threading here.

    IF (i_qty LT i_x) THEN BEGIN

      glbl_q_minmax_by_slice[0, ndx]           = MIN(all_threads_q_minmax_by_slice[0,*])
      glbl_q_minmax_by_slice[1, ndx]           = MAX(all_threads_q_minmax_by_slice[1,*])
      glbl_q_minmax_by_slice[2, ndx]           = MIN(all_threads_q_minmax_by_slice[2,*])
      glbl_q_minmax_by_slice[3, ndx]           = MAX(all_threads_q_minmax_by_slice[3,*])

    ENDIF ELSE BEGIN

      glbl_q_minmax_by_slice[0,i_p-2, ndx]     = MIN(all_threads_q_minmax_by_slice[0,i_p-2,*])  ; min_q_out, qty, step
      glbl_q_minmax_by_slice[1,i_p-2, ndx]     = MAX(all_threads_q_minmax_by_slice[1,i_p-2,*])  ; max_q_out, qty, step
      glbl_q_minmax_by_slice[2,i_p-2, ndx]     = MIN(all_threads_q_minmax_by_slice[2,i_p-2,*])  ; min_qp_out, qty, step
      glbl_q_minmax_by_slice[3,i_p-2, ndx]     = MAX(all_threads_q_minmax_by_slice[3,i_p-2,*])  ; max_qm_out, qty, step
                                    
      glbl_q_minmax_by_slice[0,i_a-2, ndx]     = MIN(all_threads_q_minmax_by_slice[0,i_a-2,*])  ; min_q_out, qty, step
      glbl_q_minmax_by_slice[1,i_a-2, ndx]     = MAX(all_threads_q_minmax_by_slice[1,i_a-2,*])  ; max_q_out, qty, step
      glbl_q_minmax_by_slice[2,i_a-2, ndx]     = MIN(all_threads_q_minmax_by_slice[2,i_a-2,*])  ; min_qp_out, qty, step
      glbl_q_minmax_by_slice[3,i_a-2, ndx]     = MAX(all_threads_q_minmax_by_slice[3,i_a-2,*])  ; max_qm_out, qty, step
                                    
      glbl_q_minmax_by_slice[0,i_o-2, ndx]     = MIN(all_threads_q_minmax_by_slice[0,i_o-2,*])  ; min_q_out, qty, step
      glbl_q_minmax_by_slice[1,i_o-2, ndx]     = MAX(all_threads_q_minmax_by_slice[1,i_o-2,*])  ; max_q_out, qty, step
      glbl_q_minmax_by_slice[2,i_o-2, ndx]     = MIN(all_threads_q_minmax_by_slice[2,i_o-2,*])  ; min_qp_out, qty, step
      glbl_q_minmax_by_slice[3,i_o-2, ndx]     = MAX(all_threads_q_minmax_by_slice[3,i_o-2,*])  ; max_qm_out, qty, step
                                    
      glbl_q_minmax_by_slice[0,i_j-2, ndx]     = MIN(all_threads_q_minmax_by_slice[0,i_j-2,*])  ; min_q_out, qty, step
      glbl_q_minmax_by_slice[1,i_j-2, ndx]     = MAX(all_threads_q_minmax_by_slice[1,i_j-2,*])  ; max_q_out, qty, step
      glbl_q_minmax_by_slice[2,i_j-2, ndx]     = MIN(all_threads_q_minmax_by_slice[2,i_j-2,*])  ; min_qp_out, qty, step
      glbl_q_minmax_by_slice[3,i_j-2, ndx]     = MAX(all_threads_q_minmax_by_slice[3,i_j-2,*])  ; max_qm_out, qty, step

    ENDELSE

  ENDFOR
  
  IF (i_qty LT i_x) THEN BEGIN
 
    glbl_q_min                                 = MIN(glbl_q_minmax_by_slice[0, *] )
    glbl_q_max                                 = MAX(glbl_q_minmax_by_slice[1, *] )

    glbl_qp_min                                = MIN(glbl_q_minmax_by_slice[2, *] )
    glbl_qm_max                                = MAX(glbl_q_minmax_by_slice[3, *] )
    
    glbl_q_minmax                              = FLTARR(4)
    glbl_q_minmax                              = [glbl_q_min, glbl_q_max, glbl_qp_min, glbl_qm_max]

  ENDIF ELSE BEGIN

    glbl_q_minmax                              = FLTARR(4,4)

    glbl_q_min                                 = MIN(glbl_q_minmax_by_slice[0,0,*] )
    glbl_q_max                                 = MAX(glbl_q_minmax_by_slice[1,0,*] )

    glbl_qp_min                                = MIN(glbl_q_minmax_by_slice[2,0,*] )
    glbl_qm_max                                = MAX(glbl_q_minmax_by_slice[3,0,*] )

    glbl_q_minmax[*,i_p - 2]                   = [glbl_q_min, glbl_q_max, glbl_qp_min, glbl_qm_max]

    glbl_q_min                                 = MIN(glbl_q_minmax_by_slice[0,1,*] )
    glbl_q_max                                 = MAX(glbl_q_minmax_by_slice[1,1,*] )

    glbl_qp_min                                = MIN(glbl_q_minmax_by_slice[2,1,*] )
    glbl_qm_max                                = MAX(glbl_q_minmax_by_slice[3,1,*] )

    glbl_q_minmax[*,i_a - 2]                   = [glbl_q_min, glbl_q_max, glbl_qp_min, glbl_qm_max]

    glbl_q_min                                 = MIN(glbl_q_minmax_by_slice[0,2,*] )
    glbl_q_max                                 = MAX(glbl_q_minmax_by_slice[1,2,*] )

    glbl_qp_min                                = MIN(glbl_q_minmax_by_slice[2,2,*] )
    glbl_qm_max                                = MAX(glbl_q_minmax_by_slice[3,2,*] )

    glbl_q_minmax[*,i_o - 2]                   = [glbl_q_min, glbl_q_max, glbl_qp_min, glbl_qm_max]

    glbl_q_min                                 = MIN(glbl_q_minmax_by_slice[0,3,*] )
    glbl_q_max                                 = MAX(glbl_q_minmax_by_slice[1,3,*] )

    glbl_qp_min                                = MIN(glbl_q_minmax_by_slice[2,3,*] )
    glbl_qm_max                                = MAX(glbl_q_minmax_by_slice[3,3,*] )
    
    glbl_q_minmax[*,i_j - 2]                   = [glbl_q_min, glbl_q_max, glbl_qp_min, glbl_qm_max]

    IF (KEYWORD_SET(global_in_minmax)) THEN BEGIN
      PRINT, 'global_q_minmax: received previous results'

      FOR i_q = 0, 3 DO BEGIN

        IF (global_in_minmax[0,i_q] LT glbl_q_minmax[0,i_q]) THEN glbl_q_minmax[0,i_q] = global_in_minmax[0,i_q]
        IF (global_in_minmax[1,i_q] GT glbl_q_minmax[1,i_q]) THEN glbl_q_minmax[1,i_q] = global_in_minmax[1,i_q]
        IF (global_in_minmax[2,i_q] LT glbl_q_minmax[2,i_q]) THEN glbl_q_minmax[2,i_q] = global_in_minmax[2,i_q]
        IF (global_in_minmax[3,i_q] GT glbl_q_minmax[3,i_q]) THEN glbl_q_minmax[3,i_q] = global_in_minmax[3,i_q]

      ENDFOR

    ENDIF

  ENDELSE
  
  RETURN, glbl_q_minmax

END
