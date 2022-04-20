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
  Chromosome* child = p1->clone();

  // We iterate over both parents separately, copying from parent1 if the
  // value is within [b,e) and from parent2 otherwise
  unsigned i = 0, j = 0;

  for ( ; i < p1->order_.size() && j < p2->order_.size(); ++i) {
    if (i >= b and i < e) {
      child->order_[i] = p1->order_[i];
    }
    else { // Increment j as long as its value is in the [b,e) range of p1
      while (p1->is_in_range(p2->order_[j], b, e)) {
        ++j;
        assert(j < p2->order_.size());
      }
      child->order_[i] = p2->order_[j];
      j++;
    }
  }

  assert(child->is_valid());
  return child;
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
Chromosome::get_fitness() const
{
  double fitness = 1000000/(cities_ptr_->total_path_distance(order_));
  return fitness;
}

// A chromsome is valid if it has no repeated values in its permutation,
// as well as no indices above the range (length) of the chromosome.
bool
Chromosome::is_valid() const
{
  if(cities_ptr_->size()!=order_.size()){
    return false;
  }
  else{
    for(long unsigned int i = 0; i<order_.size(); i++){
      int value_count = std::count(order_.begin(), order_.end(), i);
      if(value_count != 1) return false;
    }
    return true;
  }
}

// Find whether a certain value appears in a given range of the chromosome.
// Returns true if value is within the specified the range specified
// [begin, end) and false otherwise.
bool
Chromosome::is_in_range(unsigned value, unsigned begin, unsigned end) const
{
  // Add your implementation here
  for(unsigned g=begin; g < end; g++){
    if(value == order_[g]){
      return true;
    }
  }
  return false;
}
