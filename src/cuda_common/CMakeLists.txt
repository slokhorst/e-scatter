project(e-scatter/cuda_common)
list(APPEND SOURCE_FILES cuda_safe_call.cu)
add_library(cuda_common ${SOURCE_FILES})
