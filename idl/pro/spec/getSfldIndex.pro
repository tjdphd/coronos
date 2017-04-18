FUNCTION getSfldIndex, sfld

  i_k                           =  0 ; k-value
  i_pe                          =  1 ; kinetic        energy spectrum in layer ? of time step ?
  i_ae                          =  2 ; magnetic       energy spectrum in layer ? of time step ?
  i_ts                          =  3 ; total          energy spectrum in layer ? of time step ?
  i_zp                          =  4 ; elsasser+      energy spectrum in layer ? of time step ?
  i_zm                          =  5 ; elsasser-      energy spectrum in layer ? of time step ?
  i_tz                          =  6 ; total elsasser energy spectrum in layer ? of time step ?
  i_pef                         =  7 ;
  i_aef                         =  8 ;
  i_zpf                         =  9 ;
  i_zmf                         = 10 ;


  CASE sfld OF
  'pe'  : i_sfld                = i_pe
  'ae'  : i_sfld                = i_ae
  'ts'  : i_sfld                = i_ts
  'zp'  : i_sfld                = i_zp
  'zm'  : i_sfld                = i_zm
  'tz'  : i_sfld                = i_tz
  'pf'  : i_sfld                = i_pef
  'af'  : i_sfld                = i_aef
  'zpf' : i_sfld                = i_zpf
  'zmf' : i_sfld                = i_zmf
  ELSE:  i_sfld                 = -1
  ENDCASE

  RETURN, i_sfld
  
END
