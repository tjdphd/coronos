%.cu.o:	%.cu
	$(NVCC) -arch=sm_20 $(CRS_CPPFLAGS) -o $@ -c $<
bin_PROGRAMS     = coronos
coronos_SOURCES  = coronos.cpp cls_run_instance.cpp  cls_run_instance.hpp \
			             cls_canvas.cpp        cls_canvas.hpp \
				      cls_stack.cpp         cls_stack.hpp \
				    cls_lcsolve.cpp       cls_lcsolve.hpp \
			         cls_redhallmhd.cpp    cls_redhallmhd.hpp \
			      cls_parameter_map.cpp cls_parameter_map.hpp \
			          cls_parameter.cpp     cls_parameter.hpp \
                                        cls_fft.cpp           cls_fft.hpp \
                                                        nsp_constants.hpp

EXTRA_coronos_SOURCES = cls_run_instance_cuda_ext.cu   cls_run_instance_cuda_ext.hpp \
                        cls_canvas_cuda_ext.cu         cls_canvas_cuda_ext.hpp       \
                        cls_stack_cuda_ext.cu          cls_stack_cuda_ext.hpp        \
                        cls_lcsolve_cuda_ext.cu        cls_lcsolve_cuda_ext.hpp      \
                        cls_redhallmhd_cuda_ext.cu     cls_redhallmhd_cuda_ext.hpp   \
                        cls_fft_cuda_ext.cu            cls_fft_cuda_ext.hpp
if CUDA
CORONOS_SYSTEM        = cls_run_instance_cuda_ext.cu.$(OBJEXT)  \
                        cls_canvas_cuda_ext.cu.$(OBJEXT)        \
                        cls_stack_cuda_ext.cu.$(OBJEXT)         \
                        cls_lcsolve_cuda_ext.cu.$(OBJEXT)       \
                        cls_redhallmhd_cuda_ext.cu.$(OBJEXT)    \
                        cls_fft_cuda_ext.cu.$(OBJEXT)
else
endif

coronos_LDADD         = $(CORONOS_SYSTEM) $(CRS_LIBS)
coronos_DEPENDENCIES  = $(CORONOS_SYSTEM)

coronos_CPPFLAGS      = $(CRS_CPPFLAGS)
coronos_LDFLAGS       = $(CRS_LDFLAGS)
