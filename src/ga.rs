/// This module contains the genetic algorithm implementation.
use super::constants::CROSSOVER_RATE;
use crate::chromosome::Chromosome;
use rand::rngs::StdRng;
use rand::Rng;

pub fn crossover(
    rng: &mut StdRng,
    individual1: &Chromosome,
    individual2: &Chromosome,
) -> (Chromosome, Chromosome) {
    let mut new_genes1 = Vec::new();
    let mut new_genes2 = Vec::new();

    let zipped_genes = individual1.genes.iter().zip(individual2.genes.iter());

    for (gene1, gene2) in zipped_genes {
        if rng.gen_range(0.0..1.0) <= CROSSOVER_RATE {
            new_genes1.push(*gene2);
            new_genes2.push(*gene1);
        } else {
            new_genes1.push(*gene1);
            new_genes2.push(*gene2);
        }
    }

    let mut new_individual1 = individual1.clone();
    let mut new_individual2 = individual2.clone();

    new_individual1.genes = new_genes1;
    new_individual2.genes = new_genes2;

    (new_individual1, new_individual2)
}
