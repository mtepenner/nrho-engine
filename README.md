# NRHO Rendezvous & Trajectory Engine

A fast, memory-safe astrodynamics simulation engine written in Rust, specifically designed for modeling Near-Rectilinear Halo Orbits (NRHO). This engine provides the critical math and physics solvers required for complex cislunar mission planning, orbital rendezvous, and atmospheric entry.

## Table of Contents
- [Features](#features)
- [Architecture](#architecture)
- [Technologies](#technologies)
- [Installation](#installation)
- [License](#license)

## 🚀 Features
- **CR3BP Physics Engine**: Implements the Circular Restricted Three-Body Problem to accurately model spacecraft behavior in cislunar space.
- **Runge-Kutta Integrator (RK4)**: High-precision 4th-order numerical integration for calculating complex orbital trajectories.
- **Orbital Rendezvous Math**: Solves Clohessy-Wiltshire (CW) equations for modeling proximity operations and docking.
- **EDL Solver**: Calculates Entry, Descent, and Landing dynamics for atmospheric return or surface operations.
- **Telemetry State Management**: Tracks and broadcasts real-time spacecraft state vectors.

## 🏗️ Architecture
The system is built as a highly modular Rust library/binary:
1.  **`physics/`**: The core computational engine handling the heavy astrodynamics (CR3BP, RK4).
2.  **`rendezvous/`**: Proximity operations logic utilizing CW equations.
3.  **`edl/`**: Aerodynamics and landing trajectory solvers.
4.  **`telemetry/`**: State vector serialization and management.

## 🛠️ Technologies
- **Language**: Rust
- **Build System**: Cargo
- **Domain**: Astrodynamics / Orbital Mechanics

## 📥 Installation
1. Clone the repository: `git clone https://github.com/mtepenner/nrho-engine.git`
2. Build the project using Cargo: `cargo build --release`
3. Run the engine: `cargo run`
*(Note: Requires Rust and Cargo to be installed on your system).*

## ⚖️ License
This project is licensed under the MIT License. See the `LICENSE` file for details.
