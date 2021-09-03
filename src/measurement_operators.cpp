#include "darkmapper/measurement_operators.h"
#include <iostream>
#include <sopt/chained_operators.h>
#include <sopt/linear_transform.h>
#include <sopt/wavelets.h>

namespace darkmapper {
namespace planar_model {

std::tuple<sopt::OperatorFunction<Vector<t_complex>>, sopt::OperatorFunction<Vector<t_complex>>>
init_fft_2d(const t_uint& imsizey_, const t_uint& imsizex_, const t_real& oversample_ratio)
{
  t_uint const ftsizeu_ = std::floor(imsizex_ * oversample_ratio);
  t_uint const ftsizev_ = std::floor(imsizey_ * oversample_ratio);
  const t_int plan_flag = (FFTW_MEASURE | FFTW_PRESERVE_INPUT);
#ifdef darkmapper_OPENMP_FFTW
  DARKMAPPER_DEBUG("Using OpenMP threading with FFTW.");
  fftw_init_threads();
  fftw_plan_with_nthreads(omp_get_max_threads());
#endif
  Vector<t_complex> src = Vector<t_complex>::Zero(ftsizev_ * ftsizeu_);
  Vector<t_complex> dst = Vector<t_complex>::Zero(ftsizev_ * ftsizeu_);
  auto del = [](fftw_plan_s* plan) { fftw_destroy_plan(plan); };
  const std::shared_ptr<fftw_plan_s> m_plan_forward(
      fftw_plan_dft_2d(ftsizev_, ftsizeu_, reinterpret_cast<fftw_complex*>(src.data()),
          reinterpret_cast<fftw_complex*>(dst.data()), FFTW_FORWARD, plan_flag),
      del);
  const std::shared_ptr<fftw_plan_s> m_plan_inverse(
      fftw_plan_dft_2d(ftsizev_, ftsizeu_, reinterpret_cast<fftw_complex*>(src.data()),
          reinterpret_cast<fftw_complex*>(dst.data()), FFTW_BACKWARD, plan_flag),
      del);
  auto forward = [=](Vector<t_complex>& output, const Vector<t_complex>& input) {
    assert(input.size() == ftsizev_ * ftsizeu_);
    output = Vector<t_complex>::Zero(input.size());
    fftw_execute_dft(
        m_plan_forward.get(),
        const_cast<fftw_complex*>(reinterpret_cast<const fftw_complex*>(input.data())),
        reinterpret_cast<fftw_complex*>(output.data()));
    output /= std::sqrt(output.size());
  };
  auto backward = [=](Vector<t_complex>& output, const Vector<t_complex>& input) {
    assert(input.size() == ftsizev_ * ftsizeu_);
    output = Vector<t_complex>::Zero(input.size());
    fftw_execute_dft(
        m_plan_inverse.get(),
        const_cast<fftw_complex*>(reinterpret_cast<const fftw_complex*>(input.data())),
        reinterpret_cast<fftw_complex*>(output.data()));
    output /= std::sqrt(output.size());
  };
  return std::make_tuple(forward, backward);
}

//! Construct Forward model fourier space kernels
std::tuple<sopt::OperatorFunction<Vector<t_complex>>, sopt::OperatorFunction<Vector<t_complex>>, sopt::OperatorFunction<Vector<t_complex>>>
fourier_kernels(const t_uint& imsizey_, const t_uint& imsizex_)
{
  t_uint const N = imsizex_ * imsizey_;
  t_uint const n = std::sqrt(N);

  //! Precompute the KS forward/adjoint inversion kernels to reduce computation cost
  Vector<t_complex> D_ij_f = Vector<t_complex>::Zero(N);
  Vector<t_complex> D_ij_a = Vector<t_complex>::Zero(N);
    for ( int i = 0; i < n; i++){   
      for( int j = 0; j < n; j++){
          t_real kx = 0.0;
          t_real ky = 0.0;
          if(2*i-2<n){ ky = (static_cast<t_real>(i)); }
          else{ ky = (static_cast<t_real>(i-n)); }
          if(2*j-2<n){ kx = (static_cast<t_real>(j)); }
          else{ kx = (static_cast<t_real>(j-n)); }
          if( (kx != 0.0) || (ky != 0.0) ){
            D_ij_f(i*n+j) = t_complex( (ky*ky - kx*kx)/(kx*kx+ky*ky),(2.0*kx*ky)/(kx*kx+ky*ky) );
            D_ij_a(i*n+j) = t_complex( (ky*ky - kx*kx)/(kx*kx+ky*ky),(-2.0*kx*ky)/(kx*kx+ky*ky) ); } } }

  //import the ks_operators:
  auto forward = [=](Vector<t_complex>& output, const Vector<t_complex>& input) {
    assert(input.size() == N);
    output = Vector<t_complex>::Zero(input.size());
    for ( int i = 0; i < N; i++){ output(i) = input(i)*D_ij_f(i); }
  };
  auto inverse = [=](Vector<t_complex>& output, const Vector<t_complex>& input) {
    assert(input.size() == N);
    output = Vector<t_complex>::Zero(input.size());
    for ( int i = 0; i < N; i++){ output(i) = input(i)/D_ij_f(i); }
  };
  auto adjoint = [=](Vector<t_complex>& output, const Vector<t_complex>& input) {
    assert(input.size() == N);
    output = Vector<t_complex>::Zero(input.size());
    for ( int i = 0; i < N; i++){ output(i) = input(i)*D_ij_a(i); }
  };
  return std::make_tuple(forward, adjoint, inverse);
}

//! Standard masking operator
std::tuple<sopt::OperatorFunction<Vector<t_complex>>, sopt::OperatorFunction<Vector<t_complex>>> 
masking(const t_uint& imsizey_, const t_uint& imsizex_, const Vector<t_int>& mask) {
  //! Define realspace size.
  t_int const N = imsizey_ * imsizex_;
  t_int const n = std::sqrt(N);
  assert(mask.size() == N);

  //! Compute mask entry dimension
  t_int Size_Mask = 0;
  for(t_int i = 0; i < N; i++){ if(mask(i) != 0) { ++Size_Mask; } }

  std::vector<t_int> sampling;
  for(t_int i = 0; i < N; i++) { 
    if(mask(i) != 0) { sampling.push_back(i); } }

  const auto forward_mask = [=](Vector<t_complex> &out, const Vector<t_complex> &x) {
        out = Vector<t_complex>::Zero(Size_Mask);
        for (t_int i = 0; i < sampling.size(); i++) out(i) = x(sampling.at(i));
  };
  const auto inverse_mask = [=](Vector<t_complex> &out, const Vector<t_complex> &x) {
        out = Vector<t_complex>::Zero(N);
        for (t_int i = 0; i < sampling.size(); i++) out(sampling.at(i)) = x(i);
  };
  return std::make_tuple(forward_mask, inverse_mask);
}

//! combined operators of fft and KS.
std::shared_ptr<sopt::LinearTransform<Vector<t_complex>>> 
darkmapperPlane(const t_uint& imsizey_, const t_uint& imsizex_, 
  const Vector<t_int>& mask, const Vector<t_real>& covariance)
{
  t_int const N = imsizey_ * imsizex_;
  t_int const n = std::sqrt(N);

  const auto fft_ops = planar_model::init_fft_2d(imsizey_, imsizex_);
  const auto ks_ops = planar_model::fourier_kernels(imsizey_, imsizex_);

  //! Compute mask entry dimension
  t_int Size_Mask = 0;
  for(t_int i = 0; i < N; i++){ if(mask(i) != 0) { ++Size_Mask; } }

  //! Masking operators
  const auto masking_operator = planar_model::masking(imsizex_, imsizey_, mask);

  //! Covariance weighting to whiten noise for wavelet thresholding.
  const auto covariance_weighting = [=](Vector<t_complex>& output, const Vector<t_complex>& input) {
    assert(input.size() == Size_Mask);
    output = Vector<t_complex>::Zero(Size_Mask);
    for(t_int i = 0; i < Size_Mask; i++) { output(i) = input(i) / t_complex(covariance(i), covariance(i)); }
  };

  //! Chain into a tuple of operators
  const auto operator_tuple = std::make_tuple(
  sopt::chained_operators<Vector<t_complex>>(
    covariance_weighting,
    std::get<0>(masking_operator),
    std::get<1>(fft_ops), 
    std::get<0>(ks_ops), 
    std::get<0>(fft_ops) ),
  sopt::chained_operators<Vector<t_complex>>(
    std::get<1>(fft_ops), 
    std::get<1>(ks_ops), 
    std::get<0>(fft_ops),
    std::get<1>(masking_operator),
    covariance_weighting 
    )
  );

  //! Return a share pointer to this linear operator.
  return std::make_shared<sopt::LinearTransform<Vector<t_complex>>>(
    std::get<0>(operator_tuple), std::array<t_int, 3>{0, 1, static_cast<t_int>(Size_Mask)}, 
    std::get<1>(operator_tuple), std::array<t_int, 3>{0, 1, static_cast<t_int>(N)});
}

//! combined operators of fft and KS.
std::shared_ptr<sopt::LinearTransform<Vector<t_complex>>> 
kaiser_squires(const t_uint& imsizey_, const t_uint& imsizex_)
{
  t_int const N = imsizey_ * imsizex_;
  t_int const n = std::sqrt(N);

  const auto fft_ops = planar_model::init_fft_2d(imsizey_, imsizex_);
  const auto ks_ops = planar_model::fourier_kernels(imsizey_, imsizex_);

  //! Chain into a tuple of operators
  const auto operator_tuple = std::make_tuple(
  sopt::chained_operators<Vector<t_complex>>(
    std::get<1>(fft_ops), 
    std::get<0>(ks_ops), 
    std::get<0>(fft_ops) ),
  sopt::chained_operators<Vector<t_complex>>(
    std::get<1>(fft_ops), 
    std::get<2>(ks_ops), 
    std::get<0>(fft_ops)
    )
  );

  //! Return a share pointer to this linear operator.
  return std::make_shared<sopt::LinearTransform<Vector<t_complex>>>(
    std::get<0>(operator_tuple), std::array<t_int, 3>{0, 1, static_cast<t_int>(N)}, 
    std::get<1>(operator_tuple), std::array<t_int, 3>{0, 1, static_cast<t_int>(N)});
}
}   // namespace planar_model
}   // namespace darkmapper
