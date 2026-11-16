from types import SimpleNamespace

import numpy as np

import case_section_tools as cst


def test_horizontal_relative_vorticity_for_solid_body_rotation():
    x = np.linspace(-4, 4, 9)
    y = np.linspace(-3, 3, 7)
    X, Y = np.meshgrid(x, y, indexing="ij")
    omega = 2.5e-5
    # Coordinates are km; convert to metres in the analytic velocity field.
    u = -0.5 * omega * Y * 1e3
    v = 0.5 * omega * X * 1e3
    result = cst.horizontal_relative_vorticity(u, v, x, y)
    np.testing.assert_allclose(result, omega, rtol=1e-13, atol=1e-18)


def test_analytical_esp_vorticity_matches_velocity_curl():
    x = np.linspace(-15, 15, 101)
    y = np.linspace(-12, 12, 91)
    X, Y = np.meshgrid(x, y, indexing="ij")
    row = SimpleNamespace(
        xc=1.2, yc=-0.7, Rc=18.0, Omega=-3.2e-5,
        q11=1.5, q12=0.2, q22=0.75,
    )
    dx = (X-row.xc)*1e3
    dy = (Y-row.yc)*1e3
    rc = row.Rc*1e3
    rho2 = row.q11*dx**2 + 2*row.q12*dx*dy + row.q22*dy**2
    factor = row.Omega*np.exp(-rho2/rc**2)
    u = -factor*(row.q12*dx + row.q22*dy)
    v = factor*(row.q11*dx + row.q12*dy)
    numerical = cst.horizontal_relative_vorticity(u, v, x, y)
    analytical = cst.esp_relative_vorticity(X, Y, row)
    # Exclude the finite-difference boundary and allow second-order truncation.
    np.testing.assert_allclose(
        numerical[2:-2, 2:-2], analytical[2:-2, 2:-2],
        rtol=2e-3, atol=1.2e-8,
    )
