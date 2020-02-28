#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>
#include <pthread.h>

//GLOBAL
int num_size, *index_low, *index_high, NUM_THREADS;
int * arr;
//STRUCTS
typedef struct data{
	long size;
	long start_pos;
	long end_pos;
	int *start_arr;
	int *out_arr;
}data_t;

typedef struct mergeData{
	int left_start;
	int right_start;
	int part_size;
	int *start_arr;
	int *out_arr;
}mergeData_t;

//PTHREAD VARS
void *status;
int NUM_THREADS, active_worker;
pthread_mutex_t lock;


// TIMERS
double total_run_time = 0;
double total_alg_time = 0;
double start_run_time;
double start_alg_time;
double start_rng_time;
double end_rng_time;
double total_rng_time;
double end_time;
double threeway_start_time;
double threeway_end_time;
double threeway_tot_time;
static double get_wall_seconds(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
	return seconds;
}

//DECLARATIONS 
void *qs_parallel(void *);
void *find_index(void *);
void print_arr(int *, int);
void *mergeWorkers(void *);


int *gen_rand_arr(int size, char dist){
	srand(time(0));	
	int *arr = malloc(size*sizeof(int));

	//generating uniform distribution
	if (dist == 'U'){
		for (int i = 0; i<size; i++){
			arr[i] = rand() % num_size + 1;
		}
	}

	// generating normal distibution using the box muller method
	else if (dist == 'G'){
		double gaussian_one, u, v, gaussian_two;
		for (int i = 0; i<size; i+=2){
			u = (double)rand() / (double)RAND_MAX;
			v = (double)rand() / (double)RAND_MAX;
			gaussian_one = sqrt(-2*log(u))*cos(2*M_PI*v);
			gaussian_two = sqrt(-2*log(u))*sin(2*M_PI*v);
			arr[i] = (int)(num_size*gaussian_one);
			arr[i+1] = (int)(num_size*gaussian_two);
		}
		if (size % 2 != 0){
			int divisable = size % 2;
			for (int i = 0; i<divisable; i++){
				v = (double)rand() / (double)RAND_MAX;
				gaussian_one = sqrt(-2*log(u))*cos(2*M_PI*v);
				arr[i] = (int)(num_size*gaussian_one);
			}
		}

	}

	// generating exponential distribution
	else if (dist == 'E'){
		double temp, expon;
		double lambda = 0.05;
		for (int i = 0; i<size; i++){
			temp = (double)rand() / (double)RAND_MAX;
			arr[i] = -log(1 - temp) / lambda;
		}
	}
	return arr;
}



//swaps numbers
void change_num(int *a, int *b){
	int temp = *a;
	*a = *b;
	*b = temp;
}
//partitions the array
int partition (int arr[], int low, int high){
	int pivot = arr[high];
	int i = (low - 1);

	for (int j = low; j <= high - 1; j++){
		if (arr[j] <= pivot){
			i++;
			change_num(&arr[i], &arr[j]);
		}
	}
	change_num(&arr[i + 1], &arr[high]);
	return (i+1);
}
//support method that aranges the different sub arrays into the final sorted array
int * qs_parallel_manager(int *arr, int low, int high, int size, int num_threads){

	if (num_threads/2 >= 1){
		pthread_t threads[num_threads/2];
		pthread_attr_t attr;
		pthread_attr_init(&attr);

		int remainder = size%num_threads;
		mergeData_t arg_arr[num_threads/2];
		int interval = size/num_threads;
		int count=0;
		for(int t=0; t<num_threads/2; t++){
			arg_arr[t].start_arr = arr;
			arg_arr[t].left_start = interval*count;
			count++;
			arg_arr[t].right_start = interval*count;
			count++;
			arg_arr[t].part_size = interval;

			if(t == (num_threads/2)-1)
				arg_arr[t].left_start+=remainder;

			pthread_create(&threads[t], &attr, mergeWorkers, (void *)&arg_arr[t]);
		}

		for(int t=0; t<num_threads/2; t++)
			pthread_join(threads[t], &status);
		arr = qs_parallel_manager(arr, low, high, size, num_threads/2);


	}else{
		return arr;
	}
}
//parallel section for merging into the final array
void * mergeWorkers(void * arg){
	mergeData_t args = *((mergeData_t *)arg);
	int * left = args.start_arr +args.left_start;
	int * right = args.start_arr + args.right_start;
	int * arr = args.start_arr;
	int * temp = malloc(args.part_size*2*sizeof(int));
	int size = args.part_size;	//possibly times 2
	int i=0 ,j=0 ,k=0;

	while(i < size && j < size){
		if(left[i] <= right[j]){
			temp[k++] = left[i++];
		}
		else{
			temp[k++] = right[j++];
		}
	}
	while(i < size)
		temp[k++] = left[i++];
	while(j < size)
		temp[k++] = right[j++];
	i=0;
	int s = args.left_start;
	for(; i < args.part_size*2; i++)
		args.start_arr[s++] = temp[i];
	pthread_exit(0);
}

//threeway partitioning
void partition_three_way(int arr[], int left, int right, int *i, int *j){
	*index_low = left, *index_high = right;
	int part = left-1, q = right;
	int last_element = arr[right];

	while(1) {
		while (arr[*index_low] < last_element)
			*index_low++;

		while (last_element < arr[*index_high]){
			*index_high--;
			if (*index_high == left)
				break;
		}
		if (*index_low >= *index_high)
			break;

		change_num(&arr[*index_low], &arr[*index_high]);

		if (arr[*index_low] == last_element){
			part++;
			change_num(&arr[part], &arr[*index_low]);
		}

		if (arr[*index_low] == last_element){
			q--;
			change_num(&arr[*index_high], &arr[q]);
		}

	}

	change_num(&arr[*index_low], &arr[right]);

	index_high = index_low - 1;
	for (int t = left; t > q; t++, *index_high--)
		change_num(&arr[t], &arr[*index_high]);
	index_low = index_high + 1;
	for (int t = left; t>q; t--, *index_low++)
		change_num(&arr[*index_low], &arr[t]);
}
//quicksort threeway
void qs_threeway(int arr[], int bound_low, int bound_high){
	if(bound_high <= bound_low){
		return;
	}

	int index_low = bound_low;
	int i = bound_low + 1;
	int index_high = bound_high;
	int pivot = arr[bound_low];

	while(i <= index_high){
		if (arr[i] < pivot){
			change_num(&arr[index_low++], &arr[i++]);
		}
		else if (arr[i] >= pivot){
			change_num(&arr[i], &arr[index_high--]);
		} else {
			i++;
		}
	}
	qs_threeway(arr, bound_low, index_low - 1);
	qs_threeway(arr, index_high + 1, bound_high);
}

//original quicksort
void quick_sort(int arr[], int low, int high){
	while (low < high){
		int part = partition(arr, low, high);
		if (part - low < high - part){
			quick_sort(arr, low, part - 1);
			low = part + 1;
		} else {
			quick_sort(arr, part + 1, high);
			high = part - 1;
		}
	}
}
//this method creates new threads and gives eachone equal amount of work
int *quick_sort_parallel(int N, int input_arr[]){
	pthread_t threads[NUM_THREADS];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_mutex_init(&lock, NULL);

	int remainder, *res, interval;
	long index, size=(long)N;

	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	data_t argarr[NUM_THREADS];

	remainder = size%NUM_THREADS;
	res = malloc(size*sizeof(int));

	interval = size/NUM_THREADS;


	int count=0;
	for (int t = 0; t<NUM_THREADS; t++){
		argarr[t].start_arr = input_arr;
		argarr[t].out_arr = res;
		argarr[t].start_pos = interval*count;
		count++;
		argarr[t].end_pos = interval*count-1;
		argarr[t].size=size;

		if(t == NUM_THREADS - 1)
			argarr[t].end_pos += remainder;


		pthread_create(&threads[t], &attr, find_index, (void *)&argarr[t]);
	}

	for(int t=0; t<NUM_THREADS; t++)
		pthread_join(threads[t], &status);

	input_arr = qs_parallel_manager(input_arr, 0, size, size, NUM_THREADS);

	return input_arr;		
}

//The function takes the input arg typcasts it to the datastructure and runs this code for each thread
void *find_index(void *arg){
	long size, start_pos, end_pos;

	data_t *input_arg = (data_t *)arg;
	size = input_arg->size;
	start_pos = input_arg->start_pos;
	end_pos = input_arg->end_pos;

	int *out_arr=input_arg->out_arr, *start_arr=input_arg->start_arr;

	quick_sort(start_arr, start_pos, end_pos);

	pthread_exit(NULL);
}




void quicksort_threeway(int arr[], int left, int right){
	if (right <= left) return;

	int *i, *j;
	i = malloc(sizeof(int));
	j = malloc(sizeof(int));

	partition_three_way(arr, left, right, i, j);

	quicksort_threeway(arr, left, *j);
	quicksort_threeway(arr, left, right);

	free(i);
	free(j);
}

//printing array
void print_arr(int arr[], int size){
	int i;
	for (i = 0; i < size; i++)
		printf("%d ", arr[i]);
	putchar('\n');
}

int main(int argc, char *argv[]){

	start_run_time = get_wall_seconds();
	if (argc != 5){
		printf("Input error!\n./quicksort number_of_elements, intervall, distribution (G for guassian, E for exponential, U for uniform), Number of threads\n");
	} else {

		//INIT ARRAY AND KEEP TRACK OF COST
		int n = atoi(argv[1]);
		printf("n %d\n", n);
		num_size = atoi(argv[2]);
		NUM_THREADS = atoi(argv[4]);
		start_rng_time = get_wall_seconds();
		int *rand_arr = gen_rand_arr(n, *argv[3]);
		end_rng_time = get_wall_seconds();
		//print_arr(rand_arr, n);

		//RUN QUICKSORT ALG
		start_alg_time = get_wall_seconds();
		quick_sort(rand_arr, 0 ,n - 1);
		//print_arr(rand_arr, n);
		free(rand_arr);

		//END TIMER
		end_time = get_wall_seconds();
		total_run_time = end_time - start_run_time;
		total_alg_time = end_time - start_alg_time;
		total_rng_time = end_rng_time - start_rng_time;
		//print_arr(rand_arr, n);

		//RUN THREEWAY QS AND KEEP TRACK OF COST
		index_low = malloc(sizeof(int));
		index_high = malloc(sizeof(int));
		rand_arr = gen_rand_arr(n, *argv[3]);
		//print_arr(rand_arr, n);
		threeway_start_time = get_wall_seconds();
		qs_threeway(rand_arr, 0, n - 1);
		threeway_end_time = get_wall_seconds();
		threeway_tot_time = threeway_end_time - threeway_start_time;
		free(rand_arr);

		//RUN QS PARALLEL
		rand_arr = gen_rand_arr(n, *argv[3]);
		//print_arr(rand_arr, n);
		double start_parallel_time = get_wall_seconds();	
		rand_arr = quick_sort_parallel(n, rand_arr);
		//error check, doesent work with O3 flag
		/*
		int i=0;
		while(arr[i]<=arr[i+1]){
			i++;
			if(i==n-1){
				printf("it was correct\n");
				break;
			}
		}*/
		//print_arr(arr, n);
		double total_parallel_time = get_wall_seconds() - start_parallel_time;


		printf("Total RUN time: %f\nTotal ALG time: %f\nTotal RNG time: %f\nTotal THREEWAY: %f\nTotal PARALLEL: %f\n", total_run_time, total_alg_time, total_rng_time, threeway_tot_time, total_parallel_time);


	}	
	return 0;
}
