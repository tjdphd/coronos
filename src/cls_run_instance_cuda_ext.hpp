#include <cuda.h>

class run_instance_cuda_ext
{
  friend class run_instance;

  private:

  public:

  run_instance_cuda_ext();

  static int getDeviceCount();

  ~run_instance_cuda_ext();
};
