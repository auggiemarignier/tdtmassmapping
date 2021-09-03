#include "darkmapper/types.h"
#include "darkmapper/ssht_transforms.h"

#include <sopt/chained_operators.h>
#include <sopt/linear_transform.h>
#include <sopt/wavelets.h>

#include <cassert>
#include <iostream>
#include <string>
#include <fftw3.h>

namespace darkmapper {

//! Kernels related to functions on the sphere.
namespace spherical_kernels {
//! Covariance weighting applied only to the masked coefficients.
std::tuple<sopt::OperatorFunction<Vector<t_complex>>, sopt::OperatorFunction<Vector<t_complex>>> 
covariance_weighting(const t_int& L, const Vector<t_real>& covariance, const t_int& Size_Mask);
//! Masking operator 
std::tuple<sopt::OperatorFunction<Vector<t_complex>>, sopt::OperatorFunction<Vector<t_complex>>> 
masking(const t_int& L, const Vector<t_int>& mask);
//! Applying mask to complex vectors.
Vector<t_complex> apply_mask(const Vector<t_complex>& x, const Vector<t_int>& mask);
//! Applying mask to real vectors.
Vector<t_real> apply_mask_r(const Vector<t_real>& x, const Vector<t_int>& mask);
//! Computes harmonic space mapping kernel
Vector<t_complex> harmonic_space_kernel(const t_int& L);
//! Constructs lambda functions for applying kernel in harmonic space.
std::tuple<sopt::OperatorFunction<Vector<t_complex>>, sopt::OperatorFunction<Vector<t_complex>>> 
harmonic_mapping(const t_int& L);
}   // namespace spherical_kernels

//! operator that will perform spherical mass-mapping inversion.
namespace spherical_model {
//! Construct linear operator for masked analysis darkmapper_m mass-mapping  algorithm.
std::shared_ptr<sopt::LinearTransform<Vector<t_complex>>> 
darkmapperSphere(const t_int& L, const Vector<t_real>& covariance, const Vector<t_int>& mask, 
            const t_int& ssht_verbose, const ssht_dl_method_t dl_method);
//! Constructs linear operator for Spherical Kaiser Squires mass-mapping.
std::shared_ptr<sopt::LinearTransform<Vector<t_complex>>>
spherical_kaiser_squires(const t_int& L, const t_int& ssht_verbose, const ssht_dl_method_t dl_method);
}   // namespace spherical_model

//! operator that will perform planar mass-mapping inversion.
namespace planar_model {
//! 2D FFT operator
std::tuple<sopt::OperatorFunction<Vector<t_complex>>, sopt::OperatorFunction<Vector<t_complex>>>
init_fft_2d(const t_uint& imsizey_, const t_uint& imsizex_, const t_real& oversample_ratio = 1);
//! KS operator
std::tuple<sopt::OperatorFunction<Vector<t_complex>>, sopt::OperatorFunction<Vector<t_complex>>, sopt::OperatorFunction<Vector<t_complex>>>
fourier_kernels(const t_uint& imsizey_, const t_uint& imsizex_);
//! Standard masking operator
std::tuple<sopt::OperatorFunction<Vector<t_complex>>, sopt::OperatorFunction<Vector<t_complex>>> 
masking(const t_uint& imsizey_, const t_uint& imsizex_, const Vector<t_int>& mask);
//! combined operators of fft and KS.
std::shared_ptr<sopt::LinearTransform<Vector<t_complex>>> 
darkmapperPlane(const t_uint& imsizey_, const t_uint& imsizex_, const Vector<t_int>& mask, const Vector<t_real>& covariance);
//! combined operators of fft and KS.
std::shared_ptr<sopt::LinearTransform<Vector<t_complex>>> 
kaiser_squires(const t_uint& imsizey_, const t_uint& imsizex_);
}   // namespace planar_model
}   // namespace darkmapper
