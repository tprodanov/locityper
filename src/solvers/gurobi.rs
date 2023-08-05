use std::{fmt, result::Result};
use grb::{
    Env, Var, Model, Status, attr, parameter, c, add_binvar,
    expr::{LinExpr, GurobiSum},
};
use rand::Rng;
use crate::{
    Error,
    model::assgn::{GenotypeAlignments, ReadAssignment},
    ext::rand::XoshiroRng,
};

fn define_model(gt_alns: &GenotypeAlignments) -> Result<(Model, Vec<Var>), Error> {
    let mut env = Env::empty()?;
    env.set(parameter::IntParam::OutputFlag, 0)?;
    let mut model = Model::with_env("", &env.start()?)?;
    model.set_param(parameter::IntParam::Threads, 1)?;

    let total_windows = gt_alns.total_windows();
    let mut window_depth = vec![(0_u32, 0_u32); total_windows];
    let mut window_depth_constrs = vec![LinExpr::new(); total_windows];
    let mut objective = LinExpr::new();
    let mut assignment_vars = Vec::new();
    for rp in 0..gt_alns.total_reads() {
        let read_alns = gt_alns.possible_read_alns(rp);
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
                window_depth[w as usize].1 += inc;
                window_depth_constrs[w as usize].add_term(f64::from(inc), var);
                if inc == 2 {
                    break;
                }
            }
        }
        model.add_constr(&format!("R{:x}", rp), c!( (&assignment_vars[prev_len..]).grb_sum() == 1 ))?;
    }

    let depth_contrib = gt_alns.depth_contrib();
    let mut depth_vars = Vec::new();
    for (w, mut depth_constr0) in window_depth_constrs.into_iter().enumerate() {
        if !gt_alns.use_window(w) {
            continue;
        }
        let depth_distr = gt_alns.depth_distr(w);
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

/// Query read assignments from the ILP solution, and create `ReadAssignment`.
fn get_assignments<'a>(
    gt_alns: &'a GenotypeAlignments,
    model: &Model,
    assignment_vars: &[Var],
) -> Result<ReadAssignment<'a>, Error> {
    let vals = model.get_obj_attr_batch(attr::X, assignment_vars.iter().copied())?;
    let mut i = 0;
    let assgns = ReadAssignment::try_new::<_, Error>(gt_alns, |locs| {
        let j = i + locs.len();
        let new_assgn = vals[i..j].iter().position(|&v| v >= 0.9999 && v <= 1.0001)
            .ok_or_else(|| Error::solver("Gurobi", "Cannot identify read assignment (output is not 0/1)"))?;
        i = j;
        Ok(new_assgn)
    })?;
    assert_eq!(i, vals.len(), "Numbers of total read assignments do not match");
    Ok(assgns)
}

/// Gurobi ILP solver.
#[derive(Clone, Default)]
pub struct GurobiSolver;

impl super::Solver for GurobiSolver {
    /// Distribute reads between several haplotypes in a best way.
    fn solve_nontrivial<'a>(
        &self,
        gt_alns: &'a GenotypeAlignments,
        rng: &mut XoshiroRng,
    ) -> Result<ReadAssignment<'a>, Error>
    {
        let (mut model, vars) = define_model(gt_alns)?;
        model.set_param(parameter::IntParam::Seed, rng.gen::<i32>().abs())?;
        model.optimize()?;
        let status = model.status()?;
        if status != Status::Optimal {
            log::error!("Gurobi achieved non-optimal status {:?}", status);
        }
        let assgns = get_assignments(gt_alns, &model, &vars)?;
        let ilp_lik = model.get_attr(attr::ObjVal)?;
        let assgn_lik = assgns.likelihood();
        if ((assgn_lik - ilp_lik) / ilp_lik).abs() > 1e-8 {
            log::error!("Gurobi likehood differs from the model likelihood: {} and {}", ilp_lik, assgn_lik);
        }
        Ok(assgns)
    }
}

impl super::SetParams for GurobiSolver {
    fn set_param(&mut self, key: &str, _val: &str) -> Result<(), Error> {
        log::error!("Gurobi solver: unknown parameter {:?}", key);
        Ok(())
    }
}

impl fmt::Display for GurobiSolver {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Gurobi")
    }
}
