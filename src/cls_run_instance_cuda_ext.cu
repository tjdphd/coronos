#include "cls_run_instance_cuda_ext.hpp"

run_instance_cuda_ext::run_instance_cuda_ext() {

}

int run_instance_cuda_ext::getDeviceCount(void) {

  int devices     = 0;

  cudaError_t err = cudaGetDeviceCount(&devices);

  return devices;

}

run_instance_cuda_ext::~run_instance_cuda_ext() {

}
