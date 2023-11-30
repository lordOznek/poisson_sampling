use poisson_sampling::*;

fn main() {
    // Tester la fonction avec un lambda de 10 et 100 points à générer
    let sample = poisson_sample_box2d(0.1, 30, 1.0, 1.0);

    let sample_disk = poisson_sample_disk(0.1, 30, 1.0);
    println!("Samples for disk shape {:?}", sample_disk);
}
