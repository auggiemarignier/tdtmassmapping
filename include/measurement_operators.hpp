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
