# Quicksort-Pthread implementation

Quicksort using pthreads for concurrency.

To compile use Makefile
  $ make quicksort

To run
  $ ./quicksort number_of_elements, intervall, distribution (G for guassian, E for exponential, U for uniform), Number of threads

# Options
* number_of_elements = size of array
* intervall = largest number in array

# Output
The time for each algorithm 
* sequential quicksort
* threeway quicksort
* pthreads quicksort
