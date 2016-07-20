project(e-scatter/cuda_common)
list(APPEND SOURCE_FILES cuda_safe_call.cu)
cuda_add_library(common ${SOURCE_FILES})

