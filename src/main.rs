#![feature(proc_macro_hygiene)]
extern crate inline_python;
extern crate peroxide;
use peroxide::*;
use inline_python::python;
use std::f64::consts::PI;

const x_f: f64 = 20f64;
const g_eff: f64 = 81f64;
const m_pl: f64 = 1.221e+19;
const m: f64 = 100f64;

fn main() {
    let mut ode_solver = ExplicitODE::new(|st| dydx(st, 1e-14f64));
    let mut ode_solver2 = ExplicitODE::new(|st| dydx(st, 1e-16f64));

    let init_state: State<f64> = State::new(1f64, c!(y_eq(1f64)), c!(0f64));

    ode_solver
        .set_method(ExMethod::RK4)
        .set_initial_condition(init_state.clone())
        .set_step_size(1e-5)
        .set_times(10000000);

    ode_solver2
        .set_method(ExMethod::RK4)
        .set_initial_condition(init_state.clone())
        .set_step_size(1e-5)
        .set_times(10000000);

    let result = ode_solver.integrate();
    let result2 = ode_solver2.integrate();

    let x = result.col(0);
    let y = result.col(1);
    let yeq = x.fmap(|t| y_eq(t));
    let y2 = result2.col(1);

    python! {
        import pylab as plt

        plt.rc("text", usetex=True)
        plt.rc("font", family="serif")

        plt.figure(figsize=(10,6), dpi=300)
        
        plt.loglog('x, 'y)
        plt.loglog('x, 'yeq)
        plt.loglog('x, 'y2)

        axes = plt.gca()
        axes.set_ylim([1e-20, 1])

        plt.grid()
        plt.savefig("plot.png")
    }
}

fn dydx(st: &mut State<f64>, sv: f64) {
    let x = st.param;
    let y = &st.value;
    let dy = &mut st.deriv;
    dy[0] = - (PI * g_eff / 45f64).sqrt() * m_pl * m / x.powi(2) * sv * (y[0].powi(2) - y_eq(x).powi(2));
    // dy[0] = - 1f64 / x.powi(2) * sv * (y[0].powi(2) - y_eq(x).powi(2));
}

fn y_eq(x: f64) -> f64 {
    45f64 / (2f64 * PI.powi(4)) * (PI / 8f64).sqrt() / g_eff * x.powf(1.5f64) * (-x).exp()
}