use eframe::egui::{
    self,
    plot::{Line, Plot, PlotPoints, Polygon},
    CentralPanel, SidePanel, TopBottomPanel, Slider,
};

const DRIVE_STRENGTH: f64 = 1.15; // gamma
const DRIVING_FORCE_FREQ: f64 = 0.667; // omega
const NATURAL_FREQ: f64 = 1.0; // omega_0
const DAMPING_COEFF: f64 = 0.5; // beta

const LENGTH: f64 = 1.0;

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
    steps_per_second: usize,
    dt: f64,
}

impl Application {
    fn update(&mut self, dt: f64) {
        for (state, params) in self.state.iter_mut() {
            let cos_omega_t = (params.driving_force_freq * self.time).cos();
            let sin_phi = state.phi.sin();

            let phiddot = params.drive_strength * params.natural_freq.powi(2) * cos_omega_t
                - 2.0 * params.damping_coeff * state.phidot
                - params.natural_freq.powi(2) * sin_phi;

            state.phi += state.phidot * dt + 0.5 * phiddot * dt.powi(2);

            let new_phiddot = params.drive_strength * params.natural_freq.powi(2) * cos_omega_t
                - 2.0 * params.damping_coeff * state.phidot
                - params.natural_freq.powi(2) * sin_phi;

            state.phidot += (phiddot + new_phiddot) / 2.0 * dt;
        }
        self.time += dt;
    }

    fn polygonize(&self) -> Vec<Polygon> {
        const WIDTH: f64 = 0.01;
        self.state
            .iter()
            .map(|(state, _)| {
                Polygon::new(vec![
                    [WIDTH / 2.0 * state.phi.cos(), WIDTH / 2.0 * state.phi.sin()],
                    [LENGTH * state.phi.sin() + WIDTH / 2.0 * state.phi.cos(), -LENGTH * state.phi.cos() + WIDTH / 2.0 * state.phi.sin()],
                    [LENGTH * state.phi.sin() - WIDTH / 2.0 * state.phi.cos(), -LENGTH * state.phi.cos() - WIDTH / 2.0 * state.phi.sin()],
                    [-WIDTH / 2.0 * state.phi.cos(), -WIDTH / 2.0 * state.phi.sin()],
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
        state: [
            (
                State {
                    phi: 0.,
                    phidot: 0.1,
                },
                Params {
                    drive_strength: DRIVE_STRENGTH,
                    driving_force_freq: DRIVING_FORCE_FREQ,
                    natural_freq: NATURAL_FREQ,
                    damping_coeff: DAMPING_COEFF,
                },
            ),
            (
                State {
                    phi: 0.5e-3,
                    phidot: 0.1,
                },
                Params {
                    drive_strength: DRIVE_STRENGTH + 0.01,
                    driving_force_freq: DRIVING_FORCE_FREQ + 0.01,
                    natural_freq: NATURAL_FREQ + 0.01,
                    damping_coeff: DAMPING_COEFF + 0.01,
                },
            ),
            (
                State {
                    phi: 1e-3,
                    phidot: 0.1,
                },
                Params {
                    drive_strength: DRIVE_STRENGTH + 0.02,
                    driving_force_freq: DRIVING_FORCE_FREQ + 0.02,
                    natural_freq: NATURAL_FREQ + 0.02,
                    damping_coeff: DAMPING_COEFF + 0.02,
                },
            ),
        ],
        show_pendulum: [true; 3],
        steps_per_second: 100000,
        dt: 1.0 / 6000000.0,
        ..Default::default()
    };

    eframe::run_native("Duffing Pendulum", options, Box::new(|_cc| Box::new(app)));

    Ok(())
}

/********** PLOTTING **********/

impl eframe::App for Application {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        SidePanel::left("Options").resizable(true).show(ctx, |ui| {
            ui.heading("Options");
            for (i, show) in self.show_pendulum.iter_mut().enumerate() {
                ui.checkbox(show, format!("Show {}", i));
            }

            ui.add(Slider::new(&mut self.steps_per_second, 1..=1000000).text("Steps per second").logarithmic(true));
            ui.add(Slider::new(&mut self.dt, 1e-6..=1e-1).text("Time Step").logarithmic(true));
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
                });
            });

        CentralPanel::default().show(ctx, |ui| {
            Plot::new("Duffing Pendulum")
                .data_aspect(1.0)
                .show(ui, |plot_ui| {
                    for _ in 0..self.steps_per_second {
                        self.update(self.dt);
                    }

                    for ((state, _), history) in self.state.iter_mut().zip(self.history.iter_mut()) {
                        history.push_back((self.time, *state));
                        // if history.len() > 10000 {
                        //     history.remove(0);
                        // }
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

use std::{
    collections::VecDeque,
    ops::{Add, AddAssign, Mul, MulAssign},
};

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
