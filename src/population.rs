use crate::chromosome::Chromosome;
use crate::constants;
use crate::ga;
use crate::restraints;
use crate::structure::Molecule;
use core::cmp::Ordering;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::Rng;

#[derive(Debug, Clone)]
pub struct Population {
    pub chromosomes: Vec<Chromosome>,
    pub receptor: Molecule,
    pub ligand: Molecule,
    pub generation: u64,
    pub restraints: Vec<restraints::Restraint>,
    pub weights: constants::EnergyWeights,
}

impl Population {
    pub fn new(
        individuals: Vec<Chromosome>,
        receptor: Molecule,
        ligand: Molecule,
        restraints: Vec<restraints::Restraint>,
        weights: constants::EnergyWeights,
    ) -> Population {
        Population {
            chromosomes: individuals,
            receptor,
            ligand,
            generation: 0,
            restraints,
            weights,
        }
    }

    pub fn eval_fitness(&mut self) {
        let weights = self.weights;
        self.chromosomes.iter_mut().for_each(|c| {
            c.fitness(
                &self.receptor,
                &self.ligand,
                &self.restraints,
                &weights,
            );
        });
    }

    pub fn evolve(&self, rng: &mut StdRng) -> Population {
        // ELITISM: Save best individuals before evolution
        let mut elite = self.chromosomes.clone();
        elite.sort_by(|a, b| a.fitness.partial_cmp(&b.fitness).unwrap_or(Ordering::Equal));
        let elite_individuals: Vec<Chromosome> = elite
            .iter()
            .take(constants::ELITISM_COUNT)
            .cloned()
            .collect();

        // Tournament selection
        let mut new_population = self.tournament_selection(rng);

        // Crossover population
        new_population.chromosomes.shuffle(rng);
        let mut new_chromosomes = Vec::new();
        for chunk in new_population.chromosomes.chunks(2) {
            if chunk.len() == 2 {
                let (offspring1, offspring2) = ga::crossover(rng, &chunk[0], &chunk[1]);
                new_chromosomes.push(offspring1);
                new_chromosomes.push(offspring2);
            }
        }
        new_population.chromosomes = new_chromosomes;

        // Mutation
        new_population
            .chromosomes
            .iter_mut()
            .for_each(|individual| {
                individual.mutate(rng, constants::MUTATION_RATE);
            });

        // ELITISM: Replace worst individuals with elite
        new_population
            .chromosomes
            .sort_by(|a, b| b.fitness.partial_cmp(&a.fitness).unwrap_or(Ordering::Equal));
        for (i, elite_individual) in elite_individuals.iter().enumerate() {
            if i < new_population.chromosomes.len() {
                new_population.chromosomes[i] = elite_individual.clone();
            }
        }

        new_population
    }

    fn tournament_selection(&self, rng: &mut StdRng) -> Population {
        let mut offspring = Population::new(
            Vec::with_capacity(self.size()),
            self.receptor.clone(),
            self.ligand.clone(),
            self.restraints.clone(),
            self.weights,
        );

        offspring.generation = self.generation + 1;

        for _ in 0..self.size() {
            let tournament_individuals = (0..constants::TOURNAMENT_SIZE)
                .map(|_| {
                    let random_index = rng.gen_range(0..self.size());
                    &self.chromosomes[random_index]
                })
                .collect::<Vec<_>>();

            let fittest = tournament_individuals
                .iter()
                .min_by(|a, b| {
                    a.fitness
                        .partial_cmp(&b.fitness)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .expect("No individuals in tournament");

            offspring.chromosomes.push((*fittest).clone());
        }

        offspring
    }

    pub fn get_min_fittest(&self) -> &Chromosome {
        self.chromosomes
            .iter()
            .min_by(|a, b| a.fitness.partial_cmp(&b.fitness).unwrap_or(Ordering::Equal))
            .expect("No chromosomes present")
    }

    pub fn get_mean_fitness(&self) -> f64 {
        if self.chromosomes.is_empty() {
            return f64::NAN;
        }
        let total_fitness: f64 = self.chromosomes.iter().map(|c| c.fitness).sum();
        total_fitness / self.chromosomes.len() as f64
    }

    pub fn size(&self) -> usize {
        self.chromosomes.len()
    }
}
