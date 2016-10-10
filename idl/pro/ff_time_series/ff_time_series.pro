PRO ff_time_series, dsc_lab, fld_one, fld_two

i_time            = 0  ; current time
i_pe              = 1  ; perpendicular kinetic energy
i_ae              = 2  ; perpendicular magnetic energy
i_mo              = 3  ; maximum vorticity
i_imo             = 4  ; location in z for maximum vorticity
i_jmo             = 5  ; location in plain for maximum vorticity
i_oe              = 6  ; square magnitude vorticity Laplacian
i_mj              = 7  ; maximum current
i_imj             = 8  ; location in z for maximum current
i_jmj             = 9  ; location in plain for maximum current
i_ce              = 10 ; square magnitude current Laplacian
i_nuoe            = 11 ; perpendicular viscous dissipation
i_etace           = 12 ; perpendicular resistive dissipation
i_fe              = 13 ; Poynting flux
i_ftpow           = 14 ; average injected foot-point energy per unit time
i_ediss           = 15 ; average dissipated energy per unit time
i_deng            = 16 ; average rage-of-change of total energy
i_dedt            = 17 ; instantaneous rate-of-change of total energy
i_wtf             = 18 ; average rate of energy gain minus energy loss
i_consE           = 19 ; consE
i_tntauc          = 20 ; current time in units of correlation time
i_dt              = 21 ; time increment
i_dtvb            = 22 ; dtvb
i_ceocee          = 23 ; sqrt(ce/cee) whatever that is
i_aveK2           = 24 ; aveK2/t
i_aveME           = 25 ; time averaged perpendicular magnetic field
i_avepe           = 26 ; sqrt(2.*AVEpe/t)
i_kze	          = 27 ; parallel kinetic energy
i_aze             = 28 ; parallel magnetic energy
i_ne              = 29 ; internal energy
i_bze             = 30 ; Z Laplacian
i_vze             = 31 ; V_z Laplacian
i_kapbze          = 32 ; parallel conductive heat loss
i_nuvze           = 33 ; parallel viscous dissipation
 
E_case            = open_energy_data_file(dsc_lab, 'ff')

out_dev           = 'X'

ip1               = scan_parameters('ip1', 0, desc_label )                     ; power of 2 giving x-resolution
ip2               = scan_parameters('ip2', 0, desc_label )                     ; power of 2 giving y-resolution
n3                = scan_parameters('n3' , 0, desc_label )                     ; number of slices per data file
mp                = scan_parameters('mp' , 0, desc_label )                     ; number of processors used in run
zl                = scan_parameters('zl' , 0, desc_label )                     ; total height along z of integration volume
 
x_res             = 2^ip1                                                      ; resolution in x
y_res             = 2^ip2                                                      ; resolution in y
z_res             = n3 * mp                                                    ; resolution in z

i_x_res           = UINT(x_res)
i_y_res           = UINT(y_res)
i_z_res           = UINT(z_res)
 
str_x_res         = STRTRIM(i_x_res, 2)
str_y_res         = STRTRIM(i_y_res, 2)
str_z_res         = STRTRIM(i_z_res, 2)

str_zl            = STRTRIM(UINT(zl),2)

i_field           = fld_one

 CASE i_field OF 
 '1': str_field   = 'pe'
 '2': str_field   = 'ae'
 '3': str_field   = 'mo'
 '4': str_field   = 'imo'
 '5': str_field   = 'jmo'
 '6': str_field   = 'oe'
 '7': str_field   = 'mj'
 '8': str_field   = 'imj'
 '9': str_field   = 'jmj'
'10': str_field   = 'ce'
'11': str_field   = 'nuoe'
'12': str_field   = 'etace'
'13': str_field   = 'fe'
'14': str_field   = 'ftpow'
'15': str_field   = 'ediss'
'16': str_field   = 'deng'
'17': str_field   = 'dedt'
'18': str_field   = 'wtf'
'19': str_field   = 'consE'
'20': str_field   = 'tntauc'
'21': str_field   = 'dt'
'22': str_field   = 'dtvb'
'23': str_field   = 'ceocee'
'24': str_field   = 'aveK2'
'25': str_field   = 'aveME'
'26': str_field   = 'avepe'
'27': str_field   = 'kze'
'28': str_field   = 'aze'
'29': str_field   = 'ne'
'30': str_field   = 'bze'
'31': str_field   = 'vze'
'32': str_field   = 'kapbze'
'33': str_field   = 'nuvze'
 ELSE: str_field  = 'Something_is_wrong'
 ENDCASE

case_file_str     = '_'       + str_field + '_'       + str_x_res + '_'        + str_z_res + '_zl_' + str_zl 
res_str           = str_x_res + '!9X!X'   + str_y_res +'!9X!X'    + str_z_res
case_res_str      = res_str   + ' (zl = ' + str_zl    + ')'

   
E_size            = SIZE(E_case)
n_lines           = E_size[1]
                  
i_field           = fld_one
                 
max_y             = MAX(E_case[*, i_field])
min_y             = MIN(E_case[*, i_field])
                 
y_range           = [min_y, max_y]
                  
min_t             = MIN(E_case[*, i_time])
max_t             = MAX(E_case[*, i_time])
                 
t_range           = [min_t, max_t]
                 
cur_dir           = GETENV('PWD')


 CASE i_field OF 
 '1': str_title   = 'Total Perpendicular Kinetic Energy vs.t: '
 '2': str_title   = 'Total Perpendicular Magnetic Energy vs.t: '
 '3': str_title   = 'Maximum Vorticity vs.t: '
 '4': str_title   = 'Should not be looking at this vs.t: '
 '5': str_title   = 'Should not be looking at this vs.t: '
 '6': str_title   = 'Square Magnitude of Vorticity Laplacian vs.t: '
 '7': str_title   = 'Maximum Current vs.t: '
 '8': str_title   = 'Should not be looking at this vs.t: '
 '9': str_title   = 'Should not be looking at this vs.t: '
'10': str_title   = 'Square Magnitude of Current Laplacian vs.t: '
'11': str_title   = 'Perpendicular Viscous Dissipation vs.t: '
'12': str_title   = 'Perpendicular Resistive Dissipation vs.t: '
'13': str_title   = 'Poynting Flux vs.t: '
'14': str_title   = 'Average Injected Foot-Point Energy Per Unit Time vs.t: '
'15': str_title   = 'Average Dissipated Energy Per Unit Time vs.t: '
'16': str_title   = 'Average Rate-of-Change of Total Energy vs.t: '
'17': str_title   = 'Instantaneous Rate-of-Change of Total Energy vs.t: '
'18': str_title   = 'Average Rage of Energy Gains Minus Energy Losses vs.t: '
'19': str_title   = 'consE vs.t: '
'20': str_title   = 'Current Time in Units of Correlation Time vs.t: '
'21': str_title   = 'dt vs.t: '
'22': str_title   = 'dtvb vs.t: '
'23': str_title   = 'I really do not know what this is vs.t: '
'24': str_title   = 'aveK2/t vs.t: '
'25': str_title   = 'Time-Averaged Perpendicular Field Strength vs.t: '
'26': str_title   = 'Mean Photospheric Flow Velocity vs.t: '
'27': str_title   = 'Total Parallel Kinetic Energy vs.t: '
'28': str_title   = 'Total Parallel Magnetic Energy vs.t: '
'29': str_title   = 'Total Internal Energy vs.t: '
'30': str_title   = 'Square Magnitude of Z vs.t: '
'31': str_title   = 'Square Magnitude of V_z vs.t: '
'32': str_title   = 'Parallel Conductive Heat Loss vs.t: '
'33': str_title   = 'Parallel Viscous Dissipation vs.t: '
 ELSE: str_title  = 'Something_is_wrong'
 ENDCASE


       max_ediss1 = MAX(E_case[*,  i_ediss], i_max_ediss1)
     t_max_ediss1 =     E_case[i_max_ediss1, i_time]

     PRINT, '  max_ediss1 = ',   max_ediss1
     PRINT, 't_max_ediss1 = ', t_max_ediss1
     PRINT, 'i_max_ediss1 = ', i_max_ediss1

     AVL          = FLTARR(n_lines,1)
     AVL          = E_case[*, i_field]

;    i_beg        = 0
;    i_beg        = n_lines - 1 - (i_max_ediss1)

     i_beg        = i_max_ediss1
     i_end        = n_lines - 1

     ave          = MEAN(  AVL[i_beg:i_end])
     sigma        = STDDEV(AVL[i_beg:i_end])

     PRINT, FORMAT = '(A31,i2,A4,E24.16)', 'average value for field:       ', i_field, ' is ', ave
     PRINT, FORMAT = '(A31,i2,A4,E24.16)', 'standard deviation for field:  ', i_field, ' is ', sigma

     AVL[*]       = ave

   str_title      = str_title + ' For ' + case_res_str

IF (out_dev EQ 'X') THEN BEGIN
   SET_PLOT, 'X'
ENDIF ELSE BEGIN

          eps_out = cur_dir + '/ffts' + case_file_str + '.eps'

          SET_PLOT, 'PS'
          DEVICE  , /ENCAPSULATED
          DEVICE  ,  FILENAME       = eps_out
      ENDELSE

PLOT, E_case[0:(n_lines-1),i_time], E_case[0:(n_lines-1),i_field], $
      CHARSIZE    = 1.1,                                           $
      LINESTYLE   = 0,                                             $
      YTICKFORMAT = '(E7.1)',                                      $
      XMARGIN     = [12,4],                                        $
      YMARGIN     = [4,4],                                         $
      TITLE       = str_title,                                     $
      XTITLE      = 't',                                           $
      YTITLE      = ' ',                                           $
      XRANGE      = t_range,                                       $
      YRANGE      = y_range,                                       $
      THICK       = 2;,                                            $
;     /YLOG

  OPLOT, E_case[0:(n_lines-1),i_time], AVL[0:(n_lines-1)],         $
  LINESTYLE = 1,                                                   $
  THICK     = 5

IF (fld_two NE '0') THEN BEGIN
  i_field = fld_two
  OPLOT, E_case[0:(n_lines-1),i_time], E_case[0:(n_lines-1),i_field], $
         LINESTYLE   = 2
  i_field = fld_one
ENDIF

Q = ''
READ, Q, PROMPT = 'Save to postscript? [y/n]:'

IF ( Q EQ 'y') THEN BEGIN

   eps_out        = cur_dir + '/ffts' + case_file_str + '.eps'

   SET_PLOT, 'PS'
   DEVICE  , /ENCAPSULATED
   DEVICE  ,  FILENAME = eps_out, /COLOR, BITS_PER_PIXEL=8
   LOADCT, 13

   PLOT, E_case[0:(n_lines-1),i_time], E_case[0:(n_lines-1),i_field], $
         CHARSIZE    = 0.9,                                           $
         LINESTYLE   = 0,                                             $
         YTICKFORMAT = '(E7.1)',                                      $
         XMARGIN     = [12,4],                                        $
         YMARGIN     = [4,4],                                         $
         TITLE       = str_title,                                     $
         XTITLE      = 't',                                           $
         YTITLE      = ' ',                                           $
         XRANGE      = t_range,                                       $
         YRANGE      = y_range,                                       $
         THICK       = 2,                                             $
         /NODATA       ;,                                             $
;        /YLOG,                                                       $

   OPLOT, E_case[0:(n_lines-1),i_time], E_case[0:(n_lines-1),i_field],   $
         COLOR       = 255
   
   IF (fld_two NE '0') THEN BEGIN

     i_field = fld_two
     OPLOT, E_case[0:(n_lines-1),i_time], E_case[0:(n_lines-1),i_field], $
         COLOR       = 219
     i_field = fld_one

   ENDIF

   PRINT, 'Plotted to file ', eps_out

   DEVICE, /CLOSE

   SET_PLOT, 'X'

ENDIF

READ, Q, PROMPT = 'enter <Return> to quit: '
              
;OPLOT, E_wn4[*, i_time], E_wn4[*, i_pe], LINESTYLE = 2, THICK = 2

;NNLGD             = FLTARR(2,2)
;NNLGD[0,0]        = 0.20
;NNLGD[0,1]        = 0.30
;NNLGD[1,*]        = 8.0e-03
;
;WNLGD             = FLTARR(2,2)
;WNLGD[0,0]        = 0.20
;WNLGD[0,1]        = 0.30
;WNLGD[1,*]        = 7.00e-03
;
;OPLOT, NNLGD[0,*], NNLGD[1,*], LINESTYLE=0, THICK=2
;OPLOT, WNLGD[0,*], WNLGD[1,*], LINESTYLE=2, THICK=2
;
;XYOUTS, CHARSIZE  = 1.5, 0.31, 8.0e-03, 'No Noise'
;XYOUTS, CHARSIZE  = 1.5, 0.31, 7.0e-03, case_legend_str

IF (out_dev EQ 'PS') THEN DEVICE, /CLOSE

END
