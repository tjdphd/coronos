FUNCTION show_contour, n_slice, qty, n_step, tot_steps, out_dev, desc_label

COMMON save_r, Pre
COMMON First, last_call
COMMON loc, eps_out_dir
COMMON CC, Cntrs
COMMON JC, J_Cntrs
COMMON fnc, prefix, inc_res, nfld

          Pic   = construct_jcontour(n_slice, qty, n_step, tot_steps, desc_label)

          RETURN, 0
END
