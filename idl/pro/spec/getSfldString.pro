FUNCTION getSfldString, sfld

  i_sfld          = getSfldIndex(sfld)

  CASE i_sfld OF
  '0'  :  str_sfld = 'kk_'
  '1'  :  str_sfld = 'pe_'
  '2'  :  str_sfld = 'ae_'
  '3'  :  str_sfld = 'ts_'
  '4'  :  str_sfld = 'zp_'
  '5'  :  str_sfld = 'zm_'
  '6'  :  str_sfld = 'tz_'
  '7'  :  str_sfld = 'pf_'
  '8'  :  str_sfld = 'af_'
  '9'  :  str_sfld = 'zpf_'
  '10' :  str_sfld = 'zmf_'
  ELSE :  str_sfld = 'xx_'
  ENDCASE

  RETURN, str_sfld

END
