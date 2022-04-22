/*
 * Declarations for Deme class to evolve a genetic algorithm for the
 * travelling-salesperson problem.  A deme is a population of individuals.
 */

#include "chromosome.hh"
#include "deme.hh"
#include <cassert>

// Generate a Deme of the specified size with all-random chromosomes.
// Also receives a mutation rate in the range [0-1].
Deme::Deme(const Cities* cities_ptr, unsigned pop_size, double mut_rate)
    : mut_rate_(mut_rate)
{
    assert(pop_size>0);
    
    // As loop for the size of the population we want
    for(size_t u = 0; u < pop_size; ++ u){
	
	// Create a bunch of new chromosomes and add to our group
        Chromosome* ch1 = new Chromosome(cities_ptr);
        pop_.push_back(ch1);

    }

    // Fitness calculations //
    double total_fitness = 0;
    
    for(auto cur_chromosome: pop_){ //find total fitness
        total_fitness += cur_chromosome->get_fitness();
    }
    
    // Probability calculations //
    double cur_prob = 0;
    
    for(auto cur_chromosome: pop_){ // find proballity of each chromosome
        cur_prob += (cur_chromosome->get_fitness())/(total_fitness);
        fitness_intervals_.push_back(cur_prob);
    }
}

// Clean up as necessary
Deme::~Deme()
{
 for(size_t j = 0; j < pop_.size(); ++ j){
    delete pop_[j];
 }
}

double Deme::random_gen(double from, double to){ //random generator
  return(from + double(generator_()-generator_.min())/(generator_.max()-generator_.min())* to);
}

// Evolve a single generation of new chromosomes, as follows:
// We select pop_size/2 pairs of chromosomes (using the select() method below).
// Each chromosome in the pair can be randomly selected for mutation, with
// probability mut_rate, in which case it calls the chromosome mutate() method.
// Then, the pair is recombined once (using the recombine() method) to generate
// a new pair of chromosomes, which are stored in the Deme.
// After we've generated pop_size new chromosomes, we delete all the old ones.
void Deme::compute_next_generation()
{
  std::vector<Chromosome*> children_pop_;

  // Loop until we have half the population
  for(size_t n =0; n < (pop_.size()/2); ++n){

    // Assign parents
    Chromosome* parent_one = select_parent();
    Chromosome* parent_two = select_parent();
    
    // Randomly mutate
    if(random_gen(0,1)< mut_rate_){
        parent_one->mutate();
    }

    if(random_gen(0,1) < mut_rate_){
        parent_two->mutate();
    }
    
    // Take mutates parents and create children
    auto nu_children = parent_one->recombine(parent_two);
    
    // Add them to the cluster
    children_pop_.push_back(nu_children.first);
    children_pop_.push_back(nu_children.second);
  }

  // Remove the children
  for(size_t d = 0; d < pop_.size(); ++ d){
    delete pop_[d];
  }
  pop_ = children_pop_;
}

// Return a copy of the chromosome with the highest fitness.
const Chromosome* Deme::get_best() const
{
    double highest = std::numeric_limits<double>::min();
    
    Chromosome* highest_chromo = NULL;
    
    // As long as there are chromosomes to pop
    for(auto chromosome: pop_){
        // Pop and take the fitness of them
	double cur_fitness = chromosome->get_fitness();
        
	// if new record is reached
	if(cur_fitness > highest){
	    
	    // celebrate and make record higher
            highest = cur_fitness;
            highest_chromo = chromosome;
        }
    }
    // return the very best, like no one ever was...
    return highest_chromo;
}

// Randomly select a chromosome in the population based on fitness and
// return a pointer to that chromosome.

Chromosome* Deme::select_parent()
{
  // Generate some randome numbers
  double rand_numbers = random_gen(0,1);

  for(size_t p = 0; p < fitness_intervals_.size(); ++ p){
    
    // Grab some parents at random (pop them and then leave
    if(rand_numbers < fitness_intervals_[p]){
        return pop_[p];
    }
  }

  // Otherwise break and complain
  assert(false);
  return NULL;
}
