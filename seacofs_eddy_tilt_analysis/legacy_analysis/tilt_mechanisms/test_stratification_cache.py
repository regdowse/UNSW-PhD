import numpy as np

import mechanism_tools as mech


def test_constant_potential_density_has_zero_n2():
    z = np.array([[-500.0, -200.0, -50.0, -2.0]])
    sigma0 = np.full_like(z, 1025.0)

    n2, z_mid = mech._n2_from_density(sigma0, z)

    np.testing.assert_allclose(n2, 0.0)
    np.testing.assert_allclose(z_mid, [[-350.0, -125.0, -26.0]])


def test_stable_linear_potential_density_gives_positive_n2():
    z = np.array([[-500.0, -200.0, -50.0, -2.0]])
    sigma0 = 1025.0 - 0.01 * z

    n2, _ = mech._n2_from_density(sigma0, z)

    np.testing.assert_allclose(n2, 9.81e-2 / 1025.0)


def test_profile_metrics_find_pycnocline_and_density_mld():
    z_rho = np.array([[-500.0, -300.0, -180.0, -100.0, -40.0, -2.0]])
    sigma0 = np.array([[1027.0, 1026.2, 1025.6, 1025.15, 1025.01, 1025.0]])
    n2, z_mid = mech._n2_from_density(sigma0, z_rho)

    metrics = mech._profile_stratification_metrics(
        n2, z_mid, sigma0, z_rho, depths=(200, 500),
        pycnocline_limit_m=500, pycnocline_half_width_m=50,
        mld_density_threshold=0.03,
    )

    assert metrics["N2_200m_mean_s2"][0] > 0
    assert metrics["N2_500m_integral_m_s2"][0] > 0
    assert metrics["N2_pycnocline_max_s2"][0] == np.nanmax(n2)
    assert 40.0 <= metrics["MLD_density_m"][0] <= 100.0


def test_cache_signature_changes_with_scientific_configuration():
    base = mech.N2CacheConfig(depths=(200, 500))
    changed_depth = mech.N2CacheConfig(depths=(300, 500))
    changed_window = mech.N2CacheConfig(depths=(200, 500), pycnocline_half_width_m=75)

    assert base.signature != changed_depth.signature
    assert base.signature != changed_window.signature
    assert base.partition_root != changed_depth.partition_root
