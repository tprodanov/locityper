pub mod seq;
pub mod algo;
pub mod reconstr;

fn main() {
    use algo::math::{Phred, Ln};
    // for &q in &[1.0_f64, 10.0, 20.0, 30.0, 40.0, 50.0] {
    //     println!("Quality {}    Prob {:.10}   log(Prob) {:.10}   Log-prob {:.10}   Log10-prob {:.10}",
    //         q, Phred::to_prob(q), Phred::to_prob(q).ln(),
    //         Phred::to_lprob(q), algo::math::ln_to_log10(Phred::to_lprob(q)));
    // }

    // for &p in &[0.5_f64, 0.3, 0.1, 0.01, 0.001, 0.0001, 0.00001, 1e-8, 1e-16] {
    //     // println!("Prob {}    LProb {}   Phred {}   PhredL {}",
    //     //     p, p.ln(), Phred::from_prob(p), Phred::from_lprob(p.ln()));
    //     println!("Prob {}    log10(p) {}   ln(p) {}   log10->ln {}   ln->log10 {}",
    //         p, p.log10(), p.ln(), log10_to_ln(p.log10()), ln_to_log10(p.ln()));
    // }

    let a: f64 = 2.52143e-129;
    let b: f64 = 5.1235e-130;
    println!("-   {}   {}   {}", (a - b).ln(), Ln::sub(a.ln(), b.ln()),
        (a.ln().exp() - b.ln().exp()).ln());
    println!("+   {}   {}   {}", (a + b).ln(), Ln::add(a.ln(), b.ln()),
    (a.ln().exp() + b.ln().exp()).ln());

    // const LOG10: f64 = 2.302585092994046_f64;
    // println!("{:.40}", LOG10 - 10.0_f64.ln());
    // let x = vec![0.0, 1.0, 2.0, 3.0,  4.0,  5.0];
    // // let x = vec![5.0, 4.0, 3.0, 2.0, 1.0, 0.0];
    // let y = vec![0.0, 0.8, 0.9, 0.1, -0.8, -1.0];
    // let w = vec![1.0, 1.0, 1.0, 1.0,  1.0,  1.0];
    // let coefs = algo::loess::polyfit(&x, &y, &w, 3).unwrap();
    // println!("Coefs: {:?}", coefs);
    // for &z in &x {
    //     println!("   p({:.5}) = {:.5}", z, algo::loess::polyval(&coefs, z));
    // }

    // println!("Loess:");
    // let y2 = algo::loess::Loess::new()
    //     .set_degree(2)
    //     .calculate(&x, &y);
    // for (xval, yval) in x.iter().zip(y2) {
    //     println!("   {:.5} -> {:.5}", xval, yval);
    // }
}
