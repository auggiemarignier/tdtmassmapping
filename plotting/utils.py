import numpy as np
from numpy.fft import fft2, ifft2, ifftshift, fftshift


def meanvar_from_submeanvar(mu_x, mu_y, var_x, var_y, n_x, n_y):
    """
    Calculates the mean and varaince of a set of samples Z given
    the mean, variance and number of samples of two subsets X and Y
    """
    mu_z = (n_x * mu_x + n_y * mu_y) / (n_x + n_y)
    var_z = (1 / (n_x + n_y)) * (
        n_x * (mu_x ** 2 + var_x) + n_y * (mu_y ** 2 + var_y)
    ) - mu_z ** 2
    return mu_z, var_z


def snr(signal, noise):
    return 20 * np.log10(
        np.linalg.norm(signal) / np.linalg.norm(noise)
    )


class PlanarForwardModel:
    """
    Weak Gravitational Lensing planar Forward model

    Supports additional complexity which should simply
    be appended to the dir_op and adj_op objects
    appropriately (e.g. psf degridding step, psf
    deconvolution etc.)
    """

    def __init__(self, n, mask=None, ngal=None, supersample=1):
        """Construct class to hold the spherical forward and
        forward adjoint operators.

        Args:

                n (int): Pixel count along each axis (square images)
                mask (int array): Map of realspace masking.
                ngal (int array): Map of galaxy observation count.
                supersample (float): Degree of supersampling.

        Raises:

                ValueError: Raised if map is of size 0.
                ValueError: Raised if mask is the wrong shape.
                ValueError: Raised if ngal is the wrong shape.
                WarningLog: Raised if L is very large.
        """
        if n < 1:
            raise ValueError("Input map dimensions incorrect (null dimensions)!")

        # General class members
        self.n = n
        self.shape = (n, n)
        self.super = supersample
        self.ns = int((supersample - 1) * self.n / 2)

        self.get_resampling()

        # Define fourier transforms and lensing kernel
        self.fourier_kernels = self.compute_fourier_kernels()

        # Intrinsic ellipticity dispersion
        self.var_e = 0.37 ** 2

        # Define realspace masking
        if mask is None:
            self.mask = np.ones(self.shape, dtype=bool)
        else:
            self.mask = mask.astype(bool)

        # Define observational covariance
        if ngal is None:
            self.inv_cov = self.mask_forward(np.ones(self.shape))
        else:
            self.inv_cov = self.ngal_to_inv_cov(ngal)

        if self.mask.shape != self.shape:
            raise ValueError("Shape of mask map is incorrect!")

    def dir_op(self, kappa):
        """Planar weak lensing measurement operator

        Args:

                kappa (complex array): Convergence signal

        """
        # 1) Compute convergence fourier coefficients
        klm = fft2(kappa, norm="ortho")
        # 1b) Downsample map for superresolution
        if self.super > 1:
            klm = self.downsample(klm)
        # 2) Map to shear fourier coefficients
        ylm = self.fourier_kernels[0] * klm
        # 3) Compute shear realspace map
        y = ifft2(ylm, norm="ortho")
        # 4) Apply observational mask
        y_obs = self.mask_forward(y)
        # 5) Covariance weight the shear observations
        return self.cov_weight(y_obs)

    def adj_op(self, gamma):
        """Planar weak lensing adjoint measurement operator

        Args:

                gamma (complex array): Shear Observations (cov weighted)

        """
        # 1) Covariance weight the shear observations
        y_obs = self.cov_weight(gamma)
        # 2) Grid masked shear onto a full-sky
        y = self.mask_adjoint(y_obs)
        # 3) Compute shear fourier coefficients
        ylm = fft2(y, norm="ortho")
        # 4) Map to convergence fourier coefficients
        klm = self.fourier_kernels[1] * ylm
        # 4b) Pad mask for superresolution
        if self.super > 1:
            klm = self.upsample(klm)
        # 5) Compute convergence realspace map
        return ifft2(klm, norm="ortho")

    def get_resampling(self):
        N = self.n * self.super
        bounds = [self.n / 2, N - self.n // 2]
        resampling = []
        for i in range(N):
            for j in range(N):
                if i < bounds[0] or i >= bounds[1]:
                    if j < bounds[0] or j >= bounds[1]:
                        resampling.append(i * N + j)
        self.resampling = resampling

    def downsample(self, klm):
        return (
            ifftshift(
                fftshift(klm)[self.ns : self.n + self.ns, self.ns : self.n + self.ns]
            )
            / 1
        )

    def upsample(self, klm):
        return 1 * ifftshift(
            np.pad(fftshift(klm), ((self.ns, self.ns), (self.ns, self.ns)), "constant")
        )

    def ks_estimate(self, gamma):
        """Computes Kaiser-Squires estimator (for first estimate)

        Args:

                gamma (complex array): Shear Observations (patch)
        """
        ylm = fft2(gamma)
        klm = np.zeros_like(ylm)
        klm[self.fourier_kernels[0] != 0] = (
            ylm[self.fourier_kernels[0] != 0]
            / self.fourier_kernels[0][self.fourier_kernels[0] != 0]
        )
        return ifft2(klm)

    def compute_fourier_kernels(self):
        """Compuptes fourier space kernel mappings.

        Returns as a tuple {forward, inverse}.
        """
        D_f = np.zeros(self.shape, dtype=complex)

        for i in range(self.n):
            for j in range(self.n):
                kx = float(i) - float(self.n) / 2.0
                ky = float(j) - float(self.n) / 2.0
                k = kx ** 2.0 + ky ** 2.0
                if k > 0:
                    D_f[i, j] = (kx ** 2.0 - ky ** 2.0) + 1j * (2.0 * kx * ky)
                    D_f[i, j] /= k
        D_f = ifftshift(D_f)
        D_i = np.conjugate(D_f)

        return (D_f, D_i)

    def mask_forward(self, f):
        """Applies given mask to a field.

        Args:

                f (complex array): Realspace Signal

        Raises:

                ValueError: Raised if signal is nan
                ValueError: Raised if signal is of incorrect shape.

        """
        if f is not f:
            raise ValueError("Signal is NaN.")

        if f.shape != self.shape:
            raise ValueError("Signal shape is incorrect!")

        return f[self.mask]

    def mask_adjoint(self, x):
        """Applies given mask adjoint to observations

        Args:

                x (complex array): Set of observations.

        Raises:

                ValueError: Raised if signal is nan

        """
        if x is not x:
            raise ValueError("Signal is NaN.")

        f = np.zeros(self.shape, dtype=complex)
        f[self.mask] = x
        return f

    def ngal_to_inv_cov(self, ngal):
        """Converts galaxy number density map to
        data covariance.

        Assumes no intrinsic correlation between pixels.

        Args:
                ngal (real array): pixel space map of observation counts per pixel.

        """
        ngal_m = self.mask_forward(ngal)
        return np.sqrt((2.0 * ngal_m) / (self.var_e))

    def cov_weight(self, x):
        """Applies covariance weighting to observations.

        Assumes no intrinsic correlation between pixels.

        Args:
                x (array): pixel space map to be inverse covariance weighted.

        """
        return x * self.inv_cov
