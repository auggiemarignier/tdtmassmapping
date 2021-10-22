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
