//Daniel Kryzhanovsky and Gabe Howland

TO COMPILE:

run the make file with command 'make'

Some errors may appear as unused variables, it depends which generator is active.

Change algorithms in tsp.cc, comment out one of the lines between lines 132-136.

Output can also be changed within that file, just below the algorithm choice 

///////////////////////////////////////////////////////////////////////////////

METHOD:

///////////////////////////////////////////////////////////////////////////////
Chromosome: For chromosome, we started off by making sure the is_valed returns 
an actual order list so that the chromosme's could be innitialized within the 
whole process. 

The fitness is found by taking a large number (one million in this case) and 
dividing it by the total distance of the ordered cities to get the chromosome's 
fitness level. 

Nothing had to be changed about chromosome's constructor and deconstuctor.
 
mutate() is able to mutate by taking a random part of the two chromosomes, 
and switching it given by the two chromosomes switching the city with each other. 

Recombine combines two city pointers 'a' and b by checking which of the two 
are bigger, whichever one is bigger, ends up as the beginning city, 
and the other as the end. 

It then uses the built in crossover child function by initalizing two new 
children with the parents, then returning those children as a pair. 

In order to check if the value it's combining is even in the chromosome, it 
uses the is_in_range() function, which checks if the currently checked value 
is anywhere in the chromosome, returning a true or fale depending on if it's 
there or not. 

In order to generate random nummbers for the various instances of end, we 
built a random helper function that generates a number by checking within the 
range of the stated minimum and the stated max, having the range check from the
min to the max.


///////////////////////////////////////////////////////////////////////////////
Deme: For Deme's constuctor and deconstructor, we built a loop in both that 
inserts a chromosome in the population, and deletes it from the deconstructor 
when done with it.

We also added a loop that gets the toal fitness of all the chromosomes 
and a loop that determines the probabillity of each chromosome. 
Both of those are used for select_parent(). 

We also intialized mutate_rate outside of {} to be a true class parameter that 
will be used later. 

With get_best(), we are able to achieve a selctor that allows to return the 
highest fitness level in the list, that is then used to compare chromosomes for. 

For select_parent(), we use the same random number generator as the one used 
in the chromosome class to generate the random intervals that will be compared 
to the generated chromosomes. 

We then essentially match the fitness to the probability, which allows us to 
select a random parent, biasdly. 

To create the next generation of chromosomes, we select two random random 
parents and then mutate them based on the given mutation rate. We then create 
the new children by recombining the two parents, being pushed back as either 
the first or second child. 

After the intial loop, we delete the current parents through a 'for' loop and 
end the program by adding the new created children to the population.
