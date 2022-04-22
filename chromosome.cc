/*
 * Implementation for Chromosome class
 */

#include <algorithm>
#include <cassert>
#include <random>
#include <limits>

#include "chromosome.hh"

//////////////////////////////////////////////////////////////////////////////
// Generate a completely random permutation from a list of cities
Chromosome::Chromosome(const Cities* cities_ptr)
  : cities_ptr_(cities_ptr),
    order_(random_permutation(cities_ptr->size())),
    generator_(rand())
{
  assert(is_valid());
}

//////////////////////////////////////////////////////////////////////////////
// Clean up as necessary
Chromosome::~Chromosome()
{
  assert(is_valid());
}

double Chromosome::random_gen(double from, double to){ //random generator
  return(from + double(generator_()-generator_.min())/(generator_.max()-generator_.min())* to);
}

//////////////////////////////////////////////////////////////////////////////
// Perform a single mutation on this chromosome
void
Chromosome::mutate()
{
  int q = random_gen(0, order_.size()-1);
  int u = random_gen(0, order_.size()-1);
  auto thing = order_[q];
  order_[q] = order_[u];
  order_[u] = thing;
  assert(is_valid());
}

//////////////////////////////////////////////////////////////////////////////
// Return a pair of offsprings by recombining with another chromosome
// Note: this method allocates memory for the new offsprings
std::pair<Chromosome*, Chromosome*>
Chromosome::recombine(const Chromosome* other)
{
  assert(is_valid());
  assert(other->is_valid());

  unsigned a = random_gen(0, order_.size()-1);
  unsigned b = random_gen(0, order_.size()-1);
  unsigned beginning;
  unsigned ending;
  if(a<b){
    beginning = a;
    ending = b;
  }
  else{
    beginning = b;
    ending = a;
  }
  Chromosome* child_one = create_crossover_child(this, other, beginning, ending);
  Chromosome* child_two = create_crossover_child(other, this, beginning, ending);
  return std::pair<Chromosome*, Chromosome*>(child_one, child_two);
}

//////////////////////////////////////////////////////////////////////////////
// For an ordered set of parents, return a child using the ordered crossover.
// The child will have the same values as p1 in the range [b,e),
// and all the other values in the same order as in p2.
Chromosome*
Chromosome::create_crossover_child(const Chromosome* p1, const Chromosome* p2,
                                   unsigned b, unsigned e) const
{
  Chromosome* child = p1->clone(); // set pointer child to a clone of the first chromosome

  // We iterate over both parents separately, copying from parent1 if the
  // value is within [b,e) and from parent2 otherwise
  
  unsigned i = 0, j = 0; 					// create new unsigned ints

  for ( ; i < p1->order_.size() && j < p2->order_.size(); ++i) { // As long as both parents have values to check
    
    if (i >= b and i < e) { 		// If the counter value is within the intended range
     
	child->order_[i] = p1->order_[i]; // Set the child's value at the counter to the chromosome's
    
    }
    
    else { // Increment j as long as its value is in the [b,e) range of p1
      
	while (p1->is_in_range(p2->order_[j], b, e)) { // Loop for the second chromosome
      
	++j;
      
      	assert(j < p2->order_.size()); // Make sure we haven't gone out of bounds
      }

      child->order_[i] = p2->order_[j]; // set the child's order value to the second chromosome's
      j++;
    }
  }
  assert(child->is_valid()); // check we havent f***ed up
  return child; // return
}

// Return a positive fitness value, with higher numbers representing
// fitter solutions (shorter total-city traversal path).

//Helper function that builds non-random permutations for the fitness values

/*void all_permutations(Chromosome::perm_options_t options, std::vector<Cities::permutation_t>& every_perm, Cities::permutation_t& cur_permut){
  if(options.empty()){
    every_perm.push_back(cur_permut);
    return;
  }
  Chromosome::perm_options_t next_iter_options = options;
  for(auto element: options){
    cur_permut.push_back(element);
    next_iter_options.erase(element);
    all_permutations(next_iter_options, every_perm, cur_permut);
    next_iter_options.insert(element);
    cur_permut.pop_back();
  }
}*/

double
Chromosome::get_fitness() const 	// Function for determining the fitness of a chromosome
{
  double fitness = 1000000/(cities_ptr_->total_path_distance(order_)); // divide large number by total distance
  return fitness; 			// return it (higher is better)
}

// A chromsome is valid if it has no repeated values in its permutation,
// as well as no indices above the range (length) of the chromosome.
bool
Chromosome::is_valid() const
{
  if(cities_ptr_->size()!=order_.size()){ 	// If the chromosome isnt the same size as the list
    return false; 				// get mad (return false)
  }
  else{ 						// Assuming the chromosome is the same size
    
    for(long unsigned int i = 0; i<order_.size(); i++){ // Loop through the list
      int value_count = std::count(order_.begin(), order_.end(), i); // and count how many times each value appears
      if(value_count != 1) return false; // If it's more than one get mad (return false)
    }
    
    return true; 	// Otherwise we chilling
  
  }
}

// Find whether a certain value appears in a given range of the chromosome.
// Returns true if value is within the specified the range specified
// [begin, end) and false otherwise.
bool
Chromosome::is_in_range(unsigned value, unsigned begin, unsigned end) const
{
  
  for(unsigned g=begin; g < end; g++){ // Check in range for certain character
    
    if(value == order_[g]){ 	       // If that value is there then return true
      return true;
    }
  
  }
  return false;
}
