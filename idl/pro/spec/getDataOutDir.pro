FUNCTION getDataOutDir, plot_mode

  cur_dir          = GETENV('PWD')
  data_out_dir     = cur_dir                + '/ra_spec/'

  RETURN, data_out_dir

END
