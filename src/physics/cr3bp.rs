/// CR3BP (Circular Restricted Three-Body Problem) equations of motion.
///
/// The non-dimensional CR3BP equations are:
///   x'' - 2y' - x = -(1-μ)(x+μ)/r1³ - μ(x-1+μ)/r2³
///   y'' + 2x' - y = -(1-μ)y/r1³ - μy/r2³
///   z''           = -(1-μ)z/r1³ - μz/r2³
///
/// where μ = M2/(M1+M2), r1 = dist to primary, r2 = dist to secondary.

/// State vector: [x, y, z, vx, vy, vz]
pub type State = [f64; 6];

/// CR3BP mass parameter μ for the Earth–Moon system (non-dimensional).
pub const MU_EARTH_MOON: f64 = 0.012150585609624;

/// Compute the CR3BP accelerations for a given state and mass parameter.
pub fn derivatives(state: &State, mu: f64) -> State {
    let (x, y, z, vx, vy, vz) = (state[0], state[1], state[2], state[3], state[4], state[5]);

    let r1_sq = (x + mu).powi(2) + y * y + z * z;
    let r2_sq = (x - 1.0 + mu).powi(2) + y * y + z * z;
    let r1_cb = r1_sq.sqrt().powi(3);
    let r2_cb = r2_sq.sqrt().powi(3);

    let c1 = (1.0 - mu) / r1_cb;
    let c2 = mu / r2_cb;

    let ax = 2.0 * vy + x - c1 * (x + mu) - c2 * (x - 1.0 + mu);
    let ay = -2.0 * vx + y - c1 * y - c2 * y;
    let az = -c1 * z - c2 * z;

    [vx, vy, vz, ax, ay, az]
}

/// Jacobi constant (energy-like conserved quantity in CR3BP).
pub fn jacobi_constant(state: &State, mu: f64) -> f64 {
    let (x, y, z, vx, vy, vz) = (state[0], state[1], state[2], state[3], state[4], state[5]);
    let v2 = vx * vx + vy * vy + vz * vz;
    let r1 = ((x + mu).powi(2) + y * y + z * z).sqrt();
    let r2 = ((x - 1.0 + mu).powi(2) + y * y + z * z).sqrt();
    let omega = 0.5 * (x * x + y * y) + (1.0 - mu) / r1 + mu / r2;
    2.0 * omega - v2
}
