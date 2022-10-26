use eframe::egui::{
    self,
    plot::{Line, Plot, PlotPoints, Polygon},
    CentralPanel, SidePanel, TopBottomPanel, Slider,
};

use std::collections::VecDeque;

const DRIVE_STRENGTH: f64 = 1.5; // gamma
const DRIVING_FORCE_FREQ: f64 = 2.0 * core::f64::consts::PI; // omega
const NATURAL_FREQ: f64 = 1.5 * DRIVING_FORCE_FREQ; // omega_0
const DAMPING_COEFF: f64 = NATURAL_FREQ / 4.0; // beta

const LENGTH: f64 = 1.0;

const PHI_0: f64 = -core::f64::consts::FRAC_PI_2;
const EPSILON: f64 = 0.1;

#[derive(Default, Clone, Copy)]
struct Params {
    drive_strength: f64,
    driving_force_freq: f64,
    natural_freq: f64,
    damping_coeff: f64,
}

// Diffeq is \ddot{\phi} + 2\beta\dot{\phi} + \omega_0^2\sin\phi = \gamma\omega_0^2\sin\omega t
// We can rewrite this as a system of two first order ODEs:

#[derive(Default, Copy, Clone)]
struct State {
    phi: f64,
    phidot: f64,
}

#[derive(Default)]
struct Application {
    state: [(State, Params); 3],
    history: [VecDeque<(f64, State)>; 3],
    time: f64,

    show_pendulum: [bool; 3],
    show_diff: [bool; 3],
    show_log_diff: [bool; 3],
    steps_per_second: usize,
    dt: f64,
}

impl Application {
    fn diff(state: State, params: Params, time: f64) -> State {
        let State { phi, phidot } = state;
        let Params {
            drive_strength,
            driving_force_freq,
            natural_freq,
            damping_coeff,
        } = params;
        let cos_omega_t = (driving_force_freq * time).cos();
        let sin_phi = phi.sin();
        let phiddot = drive_strength * natural_freq.powi(2) * cos_omega_t
            - 2.0 * damping_coeff * phidot
            - natural_freq.powi(2) * sin_phi;
        State {
            phi: phidot,
            phidot: phiddot,
        }
    }

    fn update(&mut self, dt: f64) {
        for (state, params) in self.state.iter_mut() {
            let params = *params;

            let k1 = Self::diff(*state, params, self.time);
            let k2 = Self::diff(*state + k1 * (dt / 2.0), params, self.time + dt / 2.0);
            let k3 = Self::diff(*state + k2 * (dt / 2.0), params, self.time + dt / 2.0);
            let k4 = Self::diff(*state + k3 * dt, params, self.time + dt);
            *state += (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);
        }
        self.time += dt;
    }

    fn polygonize(&self) -> Vec<Polygon> {
        self.state
            .iter()
            .map(|(state, _)| {
                Polygon::new(vec![
                    [0.0, 0.0],
                    [LENGTH * state.phi.sin(), -LENGTH * state.phi.cos()],
                ])
            })
            .collect()
    }
}

fn main() -> anyhow::Result<()> {
    let options = eframe::NativeOptions {
        initial_window_size: Some(egui::Vec2::new(800.0, 600.0)),
        ..Default::default()
    };

    let app = Application {
        state: [(
            State {
                phi: 0.,
                phidot: 0.,
            },
            Params {
                drive_strength: DRIVE_STRENGTH,
                driving_force_freq: DRIVING_FORCE_FREQ,
                natural_freq: NATURAL_FREQ,
                damping_coeff: DAMPING_COEFF,
            },
        ), (
            State {
                phi: PHI_0,
                phidot: 0.,
            },
            Params {
                drive_strength: DRIVE_STRENGTH,
                driving_force_freq: DRIVING_FORCE_FREQ,
                natural_freq: NATURAL_FREQ,
                damping_coeff: DAMPING_COEFF,
            },
        ), (
            State {
                phi: PHI_0 + EPSILON,
                phidot: 0.,
            },
            Params {
                drive_strength: DRIVE_STRENGTH,
                driving_force_freq: DRIVING_FORCE_FREQ,
                natural_freq: NATURAL_FREQ,
                damping_coeff: DAMPING_COEFF,
            },
        )],
        show_pendulum: [false; 3],
        show_diff: [false; 3],
        show_log_diff: [false; 3],
        steps_per_second: 5000,
        dt: 1.0 / 600000.0,
        ..Default::default()
    };

    eframe::run_native("Duffing Pendulum", options, Box::new(|_cc| Box::new(app)));

    Ok(())
}

/********** PLOTTING **********/

impl eframe::App for Application {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        ctx.set_visuals(egui::Visuals::dark());

        SidePanel::left("Options").resizable(true).show(ctx, |ui| {
            ui.heading("Options");
            for (i, show) in self.show_pendulum.iter_mut().enumerate() {
                ui.checkbox(show, format!("Show {}", i));
            }

            for (i, show) in self.show_diff.iter_mut().enumerate() {
                ui.checkbox(show, format!("Show diff between {} and {}", i, (i + 1) % 3));
            }

            for (i, show) in self.show_log_diff.iter_mut().enumerate() {
                ui.checkbox(show, format!("Show log diff between {} and {}", i, (i + 1) % 3));
            }
            

            ui.add(Slider::new(&mut self.steps_per_second, 0..=100000).text("Steps per second").logarithmic(true));
            ui.add(Slider::new(&mut self.dt, 1e-5..=1e-1).text("Time Step").logarithmic(true));
        });

        TopBottomPanel::bottom("Base")
            .resizable(true)
            .show(ctx, |ui| {
                Plot::new("Graphs").show(ui, |plot_ui| {
                    for (history, enabled) in self.history.iter().zip(self.show_pendulum.iter()) {
                        if !enabled {
                            continue;
                        }

                        plot_ui.line(Line::new(PlotPoints::from_iter(
                            history.iter().map(|(t, s)| [*t, s.phi]),
                        )));
                    }

                    for (i, enabled) in self.show_diff.iter().enumerate() {
                        if !enabled {
                            continue;
                        }

                        let history = &self.history[i];
                        let history2 = &self.history[(i + 1) % 3];
                        plot_ui.line(Line::new(PlotPoints::from_iter(
                            history.iter().zip(history2.iter()).map(|((t, s), (_, s2))| [*t, s2.phi - s.phi]),
                        )));
                    }

                    for (i, enabled) in self.show_log_diff.iter().enumerate() {
                        if !enabled {
                            continue;
                        }

                        let history = &self.history[i];
                        let history2 = &self.history[(i + 1) % 3];
                        plot_ui.line(Line::new(PlotPoints::from_iter(
                            history.iter().zip(history2.iter()).map(|((t, s), (_, s2))| [*t, {
                                let logdiff = (s2.phi - s.phi).abs().log10();

                                if !logdiff.is_finite() {
                                    0.
                                } else {
                                    logdiff
                                }
                            }]),
                        )));
                    }
                });
            });

        CentralPanel::default().show(ctx, |ui| {
            Plot::new("Duffing Pendulum")
                .data_aspect(1.0)
                .show(ui, |plot_ui| {
                    if self.steps_per_second > 0 {
                        for _ in 0..self.steps_per_second {
                            self.update(self.dt);
                        }

                        for ((state, _), history) in self.state.iter_mut().zip(self.history.iter_mut()) {
                            history.push_back((self.time, *state));
                            // if history.len() > 10000 {
                            //     history.remove(0);
                            // }
                        }
                    }

                    let polygons = self.polygonize();

                    for ((history, enabled), polygon) in self
                        .history
                        .iter_mut()
                        .zip(self.show_pendulum.iter())
                        .zip(polygons)
                    {
                        if !enabled {
                            continue;
                        }

                        plot_ui.line(Line::new(PlotPoints::from_iter(
                            history
                                .iter()
                                .map(|(_, s)| [LENGTH * s.phi.sin(), -LENGTH * s.phi.cos()]),
                        )));

                        plot_ui.polygon(polygon);
                    }
                });
        });

        ctx.request_repaint();
    }
}

/********* MATH IMPLEMENTATION *********/

use std::ops::{Add, AddAssign, Mul, MulAssign};

impl Add for State {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        Self {
            phi: self.phi + other.phi,
            phidot: self.phidot + other.phidot,
        }
    }
}

impl Mul<f64> for State {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            phi: self.phi * rhs,
            phidot: self.phidot * rhs,
        }
    }
}

impl AddAssign for State {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

impl MulAssign<f64> for State {
    fn mul_assign(&mut self, rhs: f64) {
        *self = *self * rhs;
    }
}
