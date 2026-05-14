use serde::{Deserialize, Serialize};

/// Full spacecraft state vector with metadata.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StateVector {
    /// Mission elapsed time (non-dimensional CR3BP units or seconds – context-dependent).
    pub time: f64,
    /// Position [x, y, z].
    pub position: [f64; 3],
    /// Velocity [vx, vy, vz].
    pub velocity: [f64; 3],
    /// Jacobi constant (conserved in CR3BP – monitors numerical drift).
    pub jacobi: Option<f64>,
    /// Human-readable mission phase label.
    pub phase: String,
}

impl StateVector {
    pub fn from_cr3bp(
        time: f64,
        state: &[f64; 6],
        jacobi: Option<f64>,
        phase: impl Into<String>,
    ) -> Self {
        Self {
            time,
            position: [state[0], state[1], state[2]],
            velocity: [state[3], state[4], state[5]],
            jacobi,
            phase: phase.into(),
        }
    }

    /// Euclidean distance from origin.
    pub fn radius(&self) -> f64 {
        self.position.iter().map(|x| x * x).sum::<f64>().sqrt()
    }
}

/// Telemetry stream: ordered collection of state vectors.
pub type Telemetry = Vec<StateVector>;

/// Compute basic telemetry statistics.
pub fn stats(telem: &Telemetry) -> TelemetryStats {
    if telem.is_empty() {
        return TelemetryStats::default();
    }
    let n = telem.len() as f64;
    let mean_r = telem.iter().map(|s| s.radius()).sum::<f64>() / n;
    let max_r = telem.iter().map(|s| s.radius()).fold(f64::NEG_INFINITY, f64::max);
    let min_r = telem.iter().map(|s| s.radius()).fold(f64::INFINITY, f64::min);
    let duration = telem.last().map(|s| s.time).unwrap_or(0.0)
        - telem.first().map(|s| s.time).unwrap_or(0.0);

    TelemetryStats { mean_r, max_r, min_r, duration, n_points: telem.len() }
}

#[derive(Debug, Default, Serialize, Deserialize)]
pub struct TelemetryStats {
    pub mean_r: f64,
    pub max_r: f64,
    pub min_r: f64,
    pub duration: f64,
    pub n_points: usize,
}
