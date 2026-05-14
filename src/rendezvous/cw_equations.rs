/// Clohessy-Wiltshire (Hill's) equations for spacecraft rendezvous.
///
/// The CW equations describe relative motion of a chaser w.r.t. a target
/// in a circular reference orbit of mean motion `n`:
///
///   x'' - 2n*y' - 3n²x = fx
///   y'' + 2n*x'        = fy
///   z'' + n²z          = fz
///
/// State vector: [x, y, z, vx, vy, vz] (relative, in LVLH frame, km / km/s)

/// Relative state vector in LVLH frame.
pub type RelState = [f64; 6];

/// Propagate the CW state analytically over time `tau` (seconds)
/// given mean motion `n` (rad/s) and no thrust (free-drift).
///
/// Reference: Clohessy & Wiltshire (1960), "Terminal Guidance System for Satellite Rendezvous".
pub fn propagate_free(state: &RelState, n: f64, tau: f64) -> RelState {
    let (x0, y0, z0, vx0, vy0, vz0) = (
        state[0], state[1], state[2], state[3], state[4], state[5],
    );

    let nt = n * tau;
    let (s, c) = (nt.sin(), nt.cos());

    // CW state-transition matrix (STM)
    let x = (4.0 - 3.0 * c) * x0 + s / n * vx0 + 2.0 / n * (1.0 - c) * vy0;
    let y = 6.0 * (s - nt) * x0 + y0 - 2.0 / n * (1.0 - c) * vx0 + (4.0 * s - 3.0 * nt) / n * vy0;
    let z = z0 * c + vz0 / n * s;

    let vx = 3.0 * n * s * x0 + c * vx0 + 2.0 * s * vy0;
    let vy = -6.0 * n * (1.0 - c) * x0 - 2.0 * s * vx0 + (4.0 * c - 3.0) * vy0;
    let vz = -z0 * n * s + vz0 * c;

    [x, y, z, vx, vy, vz]
}

/// Compute the impulsive ΔV for a two-impulse CW rendezvous.
/// Returns (dv1, dv2) in km/s.
pub fn two_impulse_rendezvous(
    state: &RelState,
    target: &RelState,
    n: f64,
    tof: f64,
) -> ([f64; 3], [f64; 3]) {
    // Solve for initial ΔV using the CW STM inverse
    // (simplified: position match at tof; velocity boundary conditions from STM)
    let nt = n * tof;
    let (s, c) = (nt.sin(), nt.cos());

    let (x0, y0, z0) = (state[0], state[1], state[2]);
    let (xf, yf, zf) = (target[0], target[1], target[2]);

    // Position components of STM
    let s11 = 4.0 - 3.0 * c;
    let s12 = s / n;
    let s13 = 2.0 / n * (1.0 - c);
    let s21 = 6.0 * (s - nt);
    let s22 = -2.0 / n * (1.0 - c);
    let s23 = (4.0 * s - 3.0 * nt) / n;
    let s31 = c;
    let s32 = s / n;

    // Required velocities at t=0 (numerically from STM – simplified 3DOF)
    // In-plane: solve [s12 s13; s22 s23] * [vx0; vy0] = [xf - s11*x0; yf - s21*x0 - y0]
    let det = s12 * s23 - s13 * s22;
    let rhs_x = xf - s11 * x0;
    let rhs_y = yf - s21 * x0 - y0;
    let vx0_req = (rhs_x * s23 - rhs_y * s13) / det;
    let vy0_req = (s12 * rhs_y - s22 * rhs_x) / det;
    let vz0_req = (zf - s31 * z0) / s32;

    let dv1 = [
        vx0_req - state[3],
        vy0_req - state[4],
        vz0_req - state[5],
    ];

    // Propagate with required velocities to tof
    let intermediate = propagate_free(
        &[x0, y0, z0, vx0_req, vy0_req, vz0_req],
        n,
        tof,
    );
    let dv2 = [
        target[3] - intermediate[3],
        target[4] - intermediate[4],
        target[5] - intermediate[5],
    ];

    (dv1, dv2)
}
