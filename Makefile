quicksort: quicksort.c
	gcc -o quicksort quicksort.c -O3 -pthread -lm
clean:
	rm -f quicksort
