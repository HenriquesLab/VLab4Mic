__constant sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

// 1D convolution along x-axis (columns) for buffer (handles 3D data)
__kernel void
conv1d_x(__global float *image, __global float *image_out, __global float *kernel_array, int kernel_size, int nrows, int ncols){

  int row = get_global_id(0);
  int col = get_global_id(1);
  int z = get_global_id(2);

  // Return if out of bounds
  if (row >= nrows || col >= ncols) return;

  int kernel_center = (kernel_size-1)/2;
  float acc = 0.0f;

  // Calculate the linear index for this z-slice
  int slice_offset = z * nrows * ncols;

  for (int k = 0; k < kernel_size; k++) {
    int localcol = col + (k - kernel_center);
    // Zero-padding at boundaries (mode='same')
    if (localcol >= 0 && localcol < ncols) {
      acc += kernel_array[k] * image[slice_offset + row * ncols + localcol];
    }
  }

  image_out[slice_offset + row * ncols + col] = acc;
}

// 1D convolution along y-axis (rows) for buffer (handles 3D data)
__kernel void
conv1d_y(__global float *image, __global float *image_out, __global float *kernel_array, int kernel_size, int nrows, int ncols){

  int row = get_global_id(0);
  int col = get_global_id(1);
  int z = get_global_id(2);

  // Return if out of bounds
  if (row >= nrows || col >= ncols) return;

  int kernel_center = (kernel_size-1)/2;
  float acc = 0.0f;

  // Calculate the linear index for this z-slice
  int slice_offset = z * nrows * ncols;

  for (int k = 0; k < kernel_size; k++) {
    int localrow = row + (k - kernel_center);
    // Zero-padding at boundaries (mode='same')
    if (localrow >= 0 && localrow < nrows) {
      acc += kernel_array[k] * image[slice_offset + localrow * ncols + col];
    }
  }

  image_out[slice_offset + row * ncols + col] = acc;
}

// 1D convolution along z-axis (depth) for buffer (handles 3D data)
__kernel void
conv1d_z(__global float *image, __global float *image_out, __global float *kernel_array, int kernel_size, int nrows, int ncols, int ndepth){

  int row = get_global_id(0);
  int col = get_global_id(1);

  // Return if out of bounds
  if (row >= nrows || col >= ncols) return;

  int kernel_center = (kernel_size-1)/2;

  // For each depth slice at this row,col position
  for (int z = 0; z < ndepth; z++) {
    float acc = 0.0f;

    for (int k = 0; k < kernel_size; k++) {
      int localz = z + (k - kernel_center);
      // Zero-padding at boundaries (mode='same')
      if (localz >= 0 && localz < ndepth) {
        acc += kernel_array[k] * image[localz * nrows * ncols + row * ncols + col];
      }
    }

    image_out[z * nrows * ncols + row * ncols + col] = acc;
  }
}