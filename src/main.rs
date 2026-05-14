mod physics;
mod rendezvous;
mod edl;
mod telemetry;

use physics::cr3bp::{jacobi_constant, MU_EARTH_MOON, State};
use physics::rk4::propagate;
use rendezvous::cw_equations::{propagate_free, two_impulse_rendezvous};
use edl::solver::{simulate, VehicleParams};
use telemetry::state::{StateVector, stats};

fn main() {
    println!("=== NRHO Engine ===\n");

    // ── 1. CR3BP / NRHO propagation ──────────────────────────────────────────
    // Approximate L2 southern NRHO initial conditions (non-dimensional, EM system)
    // Reference: Zimovan et al. (2017) – 9:2 NRHO
    let nrho_ic: State = [
        1.0211368985965216,   // x
        0.0,                  // y
        -0.1819855765419538,  // z
        0.0,                  // vx
        -0.1033633818602695,  // vy
        0.0,                  // vz
    ];

    let mu = MU_EARTH_MOON;
    let c0 = jacobi_constant(&nrho_ic, mu);
    println!("[CR3BP] Initial Jacobi constant: {:.10}", c0);

    // Propagate one approximate NRHO period (~1.5 non-dim time units)
    let dt = 0.0001;
    let period = 1.5;
    let traj = propagate(&nrho_ic, mu, dt, period);

    let c_final = jacobi_constant(&traj.last().unwrap().1, mu);
    println!("[CR3BP] Final   Jacobi constant: {:.10}", c_final);
    println!(
        "[CR3BP] Drift:  {:.2e}  (lower is better)",
        (c_final - c0).abs()
    );
    println!("[CR3BP] Propagated {} steps over {} non-dim time\n", traj.len(), period);

    // Collect telemetry
    let telem: telemetry::state::Telemetry = traj
        .iter()
        .step_by(100) // downsample
        .map(|(t, s)| {
            let j = jacobi_constant(s, mu);
            StateVector::from_cr3bp(*t, s, Some(j), "NRHO")
        })
        .collect();

    let s = stats(&telem);
    println!(
        "[Telemetry] {} points | r_mean={:.4}  r_min={:.4}  r_max={:.4}  duration={:.3}\n",
        s.n_points, s.mean_r, s.min_r, s.max_r, s.duration
    );

    // ── 2. Clohessy-Wiltshire rendezvous ────────────────────────────────────
    // Target: Gateway at rest, chaser 1 km ahead in radial (x), 0.5 km below (z)
    let n_moon = 2.665e-6_f64; // rad/s (lunar mean motion)
    let chaser: rendezvous::cw_equations::RelState = [1.0, 0.0, 0.5, 0.0, 0.0, 0.0]; // km, km/s
    let target: rendezvous::cw_equations::RelState = [0.0; 6];
    let tof = 3600.0 * 6.0; // 6 hours

    let (dv1, dv2) = two_impulse_rendezvous(&chaser, &target, n_moon, tof);
    let dv1_mag = (dv1[0]*dv1[0] + dv1[1]*dv1[1] + dv1[2]*dv1[2]).sqrt();
    let dv2_mag = (dv2[0]*dv2[0] + dv2[1]*dv2[1] + dv2[2]*dv2[2]).sqrt();
    println!("[CW Rendezvous] TOF: {} hours", tof / 3600.0);
    println!("[CW Rendezvous] ΔV₁ = {:.4} km/s  |  ΔV₂ = {:.4} km/s", dv1_mag, dv2_mag);
    println!("[CW Rendezvous] Total ΔV = {:.4} km/s\n", dv1_mag + dv2_mag);

    // Verify free-drift: propagate chaser with required initial ΔV
    let post_maneuver = [
        chaser[0], chaser[1], chaser[2],
        chaser[3] + dv1[0],
        chaser[4] + dv1[1],
        chaser[5] + dv1[2],
    ];
    let final_state = propagate_free(&post_maneuver, n_moon, tof);
    let miss = ((final_state[0]).powi(2) + (final_state[1]).powi(2) + (final_state[2]).powi(2)).sqrt();
    println!("[CW Rendezvous] Miss distance after free drift: {:.6} km\n", miss);

    // ── 3. EDL simulation ────────────────────────────────────────────────────
    let params = VehicleParams::default();
    let edl_traj = simulate(&params, 125_000.0, 5_800.0, 0.5);

    let landed = edl_traj.last().unwrap();
    println!("[EDL] Entry: alt=125 km, vel=5800 m/s");
    println!(
        "[EDL] Touchdown: t={:.1}s  v={:.1} m/s  phase={:?}",
        landed.time_s, landed.velocity_mps, landed.phase
    );
    let max_q = edl_traj.iter().map(|s| s.dynamic_pressure_pa).fold(0.0_f64, f64::max);
    let max_g = edl_traj.iter().map(|s| s.accel_mps2.abs()).fold(0.0_f64, f64::max);
    println!("[EDL] Max-Q: {:.1} Pa  |  Peak decel: {:.2} m/s²\n", max_q, max_g);

    println!("=== Done ===");
}
