/// 4th-order Runge-Kutta integrator for the CR3BP.

use super::cr3bp::{derivatives, State};

/// Advance `state` by time step `dt` using RK4.
pub fn rk4_step(state: &State, mu: f64, dt: f64) -> State {
    let k1 = scale(&derivatives(state, mu), dt);
    let k2 = scale(&derivatives(&add(state, &scale(&k1, 0.5)), mu), dt);
    let k3 = scale(&derivatives(&add(state, &scale(&k2, 0.5)), mu), dt);
    let k4 = scale(&derivatives(&add(state, &k3), mu), dt);

    let mut next = *state;
    for i in 0..6 {
        next[i] += (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
    }
    next
}

/// Propagate from `state` for `duration` time units using step `dt`.
/// Returns a vector of (time, state) snapshots.
pub fn propagate(state: &State, mu: f64, dt: f64, duration: f64) -> Vec<(f64, State)> {
    let steps = (duration / dt).ceil() as usize;
    let mut history = Vec::with_capacity(steps + 1);
    let mut s = *state;
    let mut t = 0.0_f64;

    history.push((t, s));
    for _ in 0..steps {
        s = rk4_step(&s, mu, dt);
        t += dt;
        history.push((t, s));
    }
    history
}

fn scale(s: &State, c: f64) -> State {
    [s[0] * c, s[1] * c, s[2] * c, s[3] * c, s[4] * c, s[5] * c]
}

fn add(a: &State, b: &State) -> State {
    [a[0]+b[0], a[1]+b[1], a[2]+b[2], a[3]+b[3], a[4]+b[4], a[5]+b[5]]
}
