/// Entry-Descent-Landing (EDL) solver using a simple ballistic drag model.
///
/// Models an aeroshell descending through a planetary atmosphere using:
///   - Ballistic coefficient β = m / (Cd * A)  [kg/m²]
///   - Exponential atmosphere: ρ(h) = ρ0 * exp(-h / H)
///   - Equations of motion in the vertical plane

use serde::{Deserialize, Serialize};

/// EDL vehicle parameters.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VehicleParams {
    /// Mass in kg.
    pub mass_kg: f64,
    /// Drag coefficient (dimensionless).
    pub cd: f64,
    /// Reference area m².
    pub ref_area_m2: f64,
    /// Surface gravity m/s².
    pub gravity_mps2: f64,
    /// Atmospheric scale height m.
    pub scale_height_m: f64,
    /// Sea-level density kg/m³.
    pub rho0_kg_m3: f64,
    /// Parachute deployment altitude m.
    pub chute_deploy_alt_m: f64,
    /// Parachute drag-area product m² (Cd * A after deploy).
    pub chute_cda_m2: f64,
}

impl Default for VehicleParams {
    /// Default: generic Mars EDL profile.
    fn default() -> Self {
        Self {
            mass_kg: 900.0,
            cd: 1.7,
            ref_area_m2: 15.9,
            gravity_mps2: 3.72,
            scale_height_m: 11_100.0,
            rho0_kg_m3: 0.020,
            chute_deploy_alt_m: 10_000.0,
            chute_cda_m2: 78.5,
        }
    }
}

/// EDL state snapshot.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EdlState {
    pub time_s: f64,
    pub altitude_m: f64,
    pub velocity_mps: f64,
    pub accel_mps2: f64,
    pub dynamic_pressure_pa: f64,
    pub phase: EdlPhase,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum EdlPhase {
    Hypersonic,
    Parachute,
    Landed,
}

/// Run the EDL simulation starting from `entry_alt_m` and `entry_vel_mps`.
/// Returns a time-series of EDL states.
pub fn simulate(
    params: &VehicleParams,
    entry_alt_m: f64,
    entry_vel_mps: f64,
    dt_s: f64,
) -> Vec<EdlState> {
    let mut alt = entry_alt_m;
    let mut vel = entry_vel_mps; // m/s downward (positive)
    let mut t = 0.0_f64;
    let mut history = Vec::new();

    loop {
        if alt <= 0.0 {
            history.push(EdlState {
                time_s: t,
                altitude_m: 0.0,
                velocity_mps: vel,
                accel_mps2: 0.0,
                dynamic_pressure_pa: 0.0,
                phase: EdlPhase::Landed,
            });
            break;
        }

        let rho = params.rho0_kg_m3 * (-alt / params.scale_height_m).exp();
        let cda = if alt < params.chute_deploy_alt_m {
            params.chute_cda_m2
        } else {
            params.cd * params.ref_area_m2
        };
        let phase = if alt < params.chute_deploy_alt_m {
            EdlPhase::Parachute
        } else {
            EdlPhase::Hypersonic
        };

        let q = 0.5 * rho * vel * vel;
        let drag = q * cda;
        let accel = drag / params.mass_kg - params.gravity_mps2;

        history.push(EdlState {
            time_s: t,
            altitude_m: alt,
            velocity_mps: vel,
            accel_mps2: accel,
            dynamic_pressure_pa: q,
            phase,
        });

        // Euler integration (sufficient for planning-level fidelity)
        vel -= accel * dt_s; // deceleration
        alt -= vel * dt_s;
        t += dt_s;

        if t > 3600.0 {
            break; // safety cap
        }
    }

    history
}
