use std::{fmt, result::Result};
use grb::{
    Env, Var, Model, Status, attr, parameter, c, add_binvar,
    expr::{LinExpr, GurobiSum},
};
use rand::Rng;
use crate::{
    Error,
    model::{
        windows::REG_WINDOW_SHIFT,
        assgn::ReadAssignment,
    },
    ext::rand::XoshiroRng,
};

/// Gurobi ILP solver.
#[derive(Clone, Default)]
pub struct GurobiSolver;

/// Query read assignments from the ILP solution, and set them to `assignments`.
fn set_assignments(
    model: &Model,
    assignment_vars: &[Var],
    assignments: &mut ReadAssignment,
) -> Result<(), Error> {
    let vals = model.get_obj_attr_batch(attr::X, assignment_vars.iter().copied())?;
    let mut i = 0;
    assignments.try_init_assignments::<_, Error>(|locs| {
        let j = i + locs.len();
        let new_assgn = vals[i..j].iter().position(|&v| v >= 0.9999 && v <= 1.0001)
            .ok_or_else(|| Error::solver("Gurobi", "Cannot identify read assignment (output is not 0/1)"))?;
        i = j;
        Ok(new_assgn)
    })?;
    assert_eq!(i, vals.len(), "Numbers of total read assignments do not match");
    Ok(())
}

fn define_model(assignments: &ReadAssignment) -> Result<(Model, Vec<Var>), Error> {
    let mut env = Env::empty()?;
    env.set(parameter::IntParam::OutputFlag, 0)?;
    let mut model = Model::with_env("", &env.start()?)?;
    model.set_param(parameter::IntParam::Threads, 1)?;

    let total_windows = assignments.contig_windows().total_windows() as usize;
    let mut window_depth = vec![(0_u32, 0_u32); total_windows];
    let mut window_depth_constrs = vec![LinExpr::new(); total_windows];
    let mut objective = LinExpr::new();
    let mut assignment_vars = Vec::new();
    for rp in 0..assignments.total_reads() {
        let read_alns = assignments.possible_read_alns(rp);
        if read_alns.len() == 1 {
            let loc = &read_alns[0];
            for w in loc.windows().into_iter() {
                window_depth[w as usize].0 += 1;
            }
            objective.add_constant(loc.ln_prob());
            continue;
        }

        let prev_len = assignment_vars.len();
        for (j, loc) in read_alns.iter().enumerate() {
            let var = add_binvar!(model, name: &format!("R{:x}_{}", rp, j))?;
            assignment_vars.push(var);
            objective.add_term(loc.ln_prob(), var);
            let ws = loc.windows();
            let inc: u32 = if ws[0] == ws[1] { 2 } else { 1 };
            for w in ws.into_iter() {
                if w < REG_WINDOW_SHIFT {
                    window_depth[w as usize].0 += inc;
                } else {
                    window_depth[w as usize].1 += inc;
                    window_depth_constrs[w as usize].add_term(f64::from(inc), var);
                }
                if inc == 2 {
                    break;
                }
            }
        }
        model.add_constr(&format!("R{:x}", rp), c!( (&assignment_vars[prev_len..]).grb_sum() == 1 ))?;
    }

    let depth_contrib = assignments.depth_contrib();
    let mut depth_vars = Vec::new();
    for (w, mut depth_constr0) in window_depth_constrs.into_iter().enumerate() {
        let depth_distr = assignments.depth_distr(w);
        let (trivial_reads, non_trivial_reads) = window_depth[w];
        if non_trivial_reads == 0 {
            objective.add_constant(depth_contrib * depth_distr.ln_pmf(trivial_reads));
            continue;
        }

        depth_vars.clear();
        for depth_inc in 0..=non_trivial_reads {
            let depth = trivial_reads + depth_inc;
            let var = add_binvar!(model, name: &format!("D{:x}_{}", w, depth))?;
            depth_vars.push(var);
            if depth_inc > 0 {
                depth_constr0.add_term(-f64::from(depth_inc), var);
            }
            objective.add_term(depth_contrib * depth_distr.ln_pmf(depth), var);
        }
        model.add_constr(&format!("D{:x}_base", w), c!( (&depth_vars).grb_sum() == 1 ))?;
        model.add_constr(&format!("D{:x}_eq", w), c!( depth_constr0 == 0 ))?;
    }
    model.set_objective(objective, grb::ModelSense::Maximize)?;
    Ok((model, assignment_vars))
}

// GurobiSolver implements Solver and not MultiTrySolver, as the latter cannot store information between steps.

impl super::Solver for GurobiSolver {
    /// Distribute reads between several haplotypes in a best way.
    fn solve_nontrivial(&self, assignments: &mut ReadAssignment, rng: &mut XoshiroRng) -> Result<(), Error> {
        let (mut model, vars) = define_model(assignments)?;
        model.set_param(parameter::IntParam::Seed, rng.gen::<i32>().abs())?;
        model.optimize()?;
        let status = model.status()?;
        if status != Status::Optimal {
            log::error!("Gurobi achieved non-optimal status {:?}", status);
        }
        set_assignments(&model, &vars, assignments)?;
        let ilp_lik = model.get_attr(attr::ObjVal)?;
        let assgn_lik = assignments.likelihood();
        if (assgn_lik - ilp_lik).abs() > 1e-5 {
            log::error!("Gurobi likehood differs from the model likelihood: {} and {}", ilp_lik, assgn_lik);
        }
        Ok(())
    }
}

impl super::SetParams for GurobiSolver {
    fn set_params(&mut self, _: &json::JsonValue) -> Result<(), Error> {
        Ok(())
    }
}

impl fmt::Display for GurobiSolver {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Gurobi")
    }
}
