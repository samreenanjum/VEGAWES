#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Rdefines.h"

#define EPSILON 0.00001
#define ROUND_CNST 1000


/*-----------------------------------------------------*/
/* ------------ priority queue pqops ------------------*/
#ifndef CS161_PROJLIB
#define CS161_PROJLIB

typedef unsigned long long PriType;

extern PriType newPriority(int CPUusage,int nice);
extern int ageCPUusage(int CPUusage,int nice);
extern int compare(PriType a,PriType b);
extern void initialize(void);
extern void finalize(void);

extern PriType strtopri(char *s);
extern char *pritostr(PriType p);
extern void die(char *s);
extern void warn(char *s);
extern void readCommand(void);
extern void printTop(int jobid,int duration,int nice,int heapsize);

extern char cmdletter;
extern int cmdtick,cmdjobid,cmdduration,cmdnice;
extern int tick;

#define CPU_PER_TICK 100
#define TICKS_UNTIL_AGING 100

#endif

#define FREE(x)  free(x) ; x = NULL            /* core if free'ed pointer is used */
#define LEFT(x)  (2*x)                         /* left child of a node */
#define RIGHT(x) ((2*x)+1)                     /* right child of a node */
#define PARENT(x) (x/2)                        /* parent of a node */
#define SWAP(t,x,y) tmp = x ; x = y ; y = tmp  /* swap to variables */

#define MSGSIZE 128       /* max size of a debug message */
#define FAILED   -1       /* define -1 to indicate failure of some operations */
#define EMPTY     1       /* define  1 to indicate that the heap is empty */ 

/* just in case NULL is not defined */
#ifndef NULL
#define NULL    0
#endif


char messages[MSGSIZE];
FILE *file; 
typedef double priority;


/* define the heap object*/
typedef struct node {
  priority p;
  unsigned int id;
  short sentinel;
} node ;

/* define a structure representing an individual node in the heap, and
 * make it a valid type for convenience */
typedef struct breakpoint {
  priority p;
  unsigned int id;
  double sum;
  int size;
  unsigned int prev;
  unsigned int next;
} breakpoint ;

/* create a global node tmp, for swaping purposes */
node tmp;

/* for convience in function declarations, typedef a pointer to a node
 * as its own type, node_ptr */
typedef node * node_ptr;

/* define a structure representing the heap, and make it a valid type
 * for convenience */
typedef struct binary_heap {
  int heap_size;
  int max_elems;
  node_ptr elements;
} binary_heap ; 

/* function prototypes for functions which operate on a binary heap */ 
extern void        heapify(binary_heap *a,int i);
extern node_ptr    heap_max(binary_heap *a);
extern node        heap_extract_max(binary_heap *a);
extern void        heap_insert(binary_heap *a,node key);
extern void        heap_delete(binary_heap *a,int i);
extern void        heap_increase_key(binary_heap *a,int i,priority p);
extern void        heap_initalize(binary_heap *a,int nodes);
extern void        heap_finalize(binary_heap *a);

/* function prototypes for functions which operate on a node */
extern int         node_find(binary_heap a,unsigned int id);
extern node        node_create(unsigned int id,priority p);
extern breakpoint  breakpoint_create(unsigned int id, priority p, double sum, int size, unsigned int prev, unsigned int next);

/* function prototypes for helper functions */
extern int         compare_priority(node i,node j);
extern void        print_error(char *msg);


void heapify(binary_heap *a,int i) {
  register int l,r,largest;
  
  l = LEFT(i);
  r = RIGHT(i);

  /* check the left child */
  largest = ((l <= a->heap_size && 
        compare_priority(a->elements[l],a->elements[i])) ? l : i);

  /* check the right child */
  if (r <= a->heap_size && 
      compare_priority(a->elements[r],a->elements[largest])) largest = r;

  if (largest != i) { 

    /* swap nodes largest and i, then heapify */

    SWAP(node,a->elements[i],a->elements[largest]);
    heapify(a,largest);
  }
}

/* Function to return the max (first) node of a heap */

node_ptr heap_max(binary_heap *a) { 
  return ((a->heap_size <= 0) ? NULL : &(a->elements[1])); 
}

/* Function to remove the max node from the heap and return it.  The
 * running time is O(lg(n)) since it performs only a costant amount of
 * work on top of the O(lg(n)) of heapify(). Adapted from Introduction
 * to Algorithms (Cormen, Leiserson, Rivest 1990) page 150 */

node heap_extract_max(binary_heap *a) {
  node max;

  max.sentinel = 1;

  /* if there are elements in the heap, make the last item in the heap
   * the first one, shorten the heap by one and call heapify(). */

  if (a->heap_size >= 1) {
    max = a->elements[1];
    a->elements[1] = a->elements[(a->heap_size)--];
    heapify(a,1);
  }

  return max;
}

/* Function to insert an element into the heap, worst case running
 * time is O(lg(n)) on an n element heap, since the path traced from
 * the new leaf to the root has at most length lg(n). This occurs when
 * the new leaf should be the root node.  Adapted from Introduction to
 * Algorithms (Cormen, Leiserson, Rivest 1990) page 150 */

void heap_insert(binary_heap *a,node key) {
  register int i;

  /* if the heap already has the max number of elements we do not
   * allow more elements to be added */
  if (a->heap_size >= a->max_elems) {
    print_error("Heap capacity exceeded, new element not added.");
    return;
  }

  /* increase the heap size to accomidate the new node, and set the
   * inital position of this node to be the last node in the heap */
  i = ++(a->heap_size);

  /* traverse the parth from the leaf to the root to find the a proper
   * place for the new element */
  while (i > 1 && compare_priority(key,a->elements[PARENT(i)])) {
    a->elements[i] = a->elements[PARENT(i)];
    i = PARENT(i);
  }

  /* insert the element at the position that was determined */
  a->elements[i] = key;
}

/* Function to delete a node from the heap. Adapted from Introduction
 * to Algorithms (Cormen, Leiserson, Rivest 1990) page 151 Exercise
 * 7.5-5 */

void heap_delete(binary_heap *a,int i) {
  node deleted;

  /* return with an error if the input is invalid, ie trying to delete
   * elements that are outside of the heap bounds, 1 to heap_size */
  if (i > a->heap_size || i < 1) {    
    sprintf(messages,"heap_delete(): %d, no such element.",i);
    print_error(messages);
    return;
  }

  /* switch the item to be deleted with the last item, and then
   * shorten the heap by one */
  deleted = a->elements[i];
  a->elements[i] = a->elements[(a->heap_size)--];
  
  heapify(a,i);
}


/* Function to increase the key value of a node from in the
 * heap. Adapted from Introduction to Algorithms (Cormen, Leiserson,
 * Rivest 1990) page 151 Exercise 7.5-4 */
void heap_increase_key(binary_heap *a,int i,priority p) {

  /* return with an error if the input is invalid, ie trying to
   * increase elements that are outside of the heap bounds, 1 to
   * heap_size */
  if (i > a->heap_size || i < 1) {
    sprintf(messages,"heap_increase_key(): %d, no such element.",i);
    print_error(messages);
    return;
  }

  /* change and propagate */
  a->elements[i].p = p;
  heapify(a,i);
}

/* function to initalize a given binary heap */
void heap_initalize(binary_heap *a,int nodes) { 

  /* We initalize heap_size to zero, since a newly created heap
   * contains no elements. */
  a->heap_size = 0; 

  /* we set the max elems to the requested number of nodes, and the
   * allocate enough space for this + 1 number of nodes, since the
   * heap is always numbered from 1, but array/pointer accesses are
   * always from 0. */
  a->max_elems = nodes;
  a->elements = (node_ptr)malloc(sizeof(node)*((a->max_elems)+1));

  /* mark the zero'th element of the heap a to be empty, just in case
   * it is every accessed */
  a->elements[0].sentinel = 1;
}


/* function to clean up after we are done with the heap */
void heap_finalize(binary_heap *a) { FREE(a->elements); }

/* function to create a node */
node node_create(unsigned int id,
		 priority p) {
  node n;
  n.id = id;
  n.p = p;
  n.sentinel = 0;
  return n;
}

/* function to create a breakpoint */
breakpoint  breakpoint_create(unsigned int id, priority p, double sum, int size, unsigned int prev, unsigned int next){
	breakpoint b;
	b.id = id;
	b.p = p;
	b.sum = sum;
	b.size = size;
	b.prev = prev;
	b.next = next;
	return b;
}

/* function to compare the priority of two given nodes, this is a
 * wrapper for the given compare routine, since in all heap
 * comparisions, we are only interested in greater than or less than
 * operations */
int compare_priority(node i,node j) {
  if( i.p > j.p ){
	return 1;
  }else{
	if( i.p < j.p ){
		return 0;
	}
	else{
		if (i.id < j.id) return 1;
	}
  }
  return 0;
}

/* function to find if a node is in the heap, O(n) worst case, since
 * we will have to consider every element in a failed search */
int node_find(binary_heap a,unsigned int id) {
  register int i;

  for (i = 1; i<=a.heap_size; i++)
    if (id == a.elements[i].id) return i;
  return FAILED;
}

/* function to print an error message */
void print_error(char *msg) { Rprintf("# ERROR: %s\n",msg); 
}
/* ------------ end priority queue pqops ------------------*/
/*---------------------------------------------------------*/


/* number of observed probes*/
int num_of_probes;

/* minumum region size allowed*/
int min_region_size;

/* constant used in the stop condition*/
double beta;

/* constant used in the third term of the merging condition*/
double alpha;

/* function prototypes for segmentation */ 

/* main function */
void vegawes(double *data, int *markers_start, int *markers_end, double *positions, int *start, int *end, int *size, double *mean, int *label, int *n, double *al, double *be, int *m, double *std, int *n_reg);

/* this function computes the initial trivial segmentation*/
void init_trivial_segmentation(breakpoint *brks, double *data, binary_heap *b);

/* this function updates the priority for the i-th breakpoint*/
double update_priority(breakpoint *brks, int i, int *markers_start, int *markers_end, double *positions);

/* this function deletes the regions having a size lower than the specified min_region_size */
void delete_outlier(breakpoint *brks, int num_of_brks, int *markers_start, int *markers_end, double *positions);

/* Added: Computes the mean of the distances between regions*/
double compute_mean_distances(breakpoint *brks, int *markers_start, int *markers_end);

/* Added: Computes local avg of distance between exons in a region i*/
double compute_local_mean_distance(breakpoint *brks, double *positions, int i);

void vegawes(double *data, int *markers_start, int *markers_end, double *positions, int *start, int *end, int *size, double *mean, int *label, int *n, double *al, double *be, int *m, double *std, int *n_reg){


	int i, j;
	int index;
	double lambda, lambda_gradient, stop_lambda_gradient, tmp_lambda, tmp_lambda2;
	num_of_probes = *n;
  alpha = *al;
	beta=*be;
	min_region_size = -(*m);
	breakpoint  brk_del;
	binary_heap b;
	int node_index;
	node max, *maxptr, tmp_node;
	breakpoint brks[num_of_probes+1];
  
	/* compute the stop condition */
	stop_lambda_gradient = (*std)*beta;

	/* Iiitialize the binary heap*/
	heap_initalize(&b,num_of_probes-1);
	
	
	/*Compute the trivial segmenatation*/
	init_trivial_segmentation(brks, data, &b);

	
	lambda =  (heap_max(&b)->p)-EPSILON;
	lambda_gradient = abs(lambda);
	int dist_ct = 0;
	while( (b.heap_size>0) && (lambda_gradient<=stop_lambda_gradient)){
		
		while( (b.heap_size>0) && (heap_max(&b)->p > lambda) ){
			
			max = heap_extract_max(&b);
			
			if( (max.p > brks[max.id].p) && (brks[max.id].p < lambda) ){			
				max.p = brks[max.id].p;
				heap_insert(&b,max);
				
			}else{
				brk_del = brks[max.id];

				brk_del = brks[max.id];
				brks[brk_del.prev].next = brk_del.next;
				brks[brk_del.next].prev = brk_del.prev;
				brks[brk_del.prev].sum = brks[brk_del.prev].sum + brk_del.sum;
				brks[brk_del.prev].size = brks[brk_del.prev].size + brk_del.size;
				
				tmp_lambda = update_priority(brks, brk_del.prev, markers_start, markers_end, positions);
		
        
				if(brk_del.prev>0){
					node_index = node_find(b, brk_del.prev);					
					tmp_node = b.elements[node_index];
					heap_delete(&b, node_index);
					tmp_node.p = tmp_lambda;
					heap_insert(&b,tmp_node);
				}
				brks[brk_del.prev].p = tmp_lambda;
				
								

				if(brk_del.next<num_of_probes){
				tmp_lambda = update_priority(brks, brk_del.next, markers_start, markers_end, positions);
					node_index = node_find(b, brk_del.next);
					tmp_node = b.elements[node_index];
					heap_delete(&b, node_index);
					tmp_node.p = tmp_lambda;
				
					heap_insert(&b,tmp_node);
				}
				brks[brk_del.next].p = tmp_lambda;
				
			}
			
        
		}
		
    
		/*Update the current lambda value*/
		if( (b.heap_size>0) ){
			lambda_gradient = (float)  lambda -(heap_max(&b)->p);
			lambda = (heap_max(&b)->p)-EPSILON;
		}
	}
	
	delete_outlier(brks, b.heap_size, markers_start, markers_end, positions);
	
	
	for(i=0, j=0; i< num_of_probes;j++){
	
		brk_del = brks[i];
		start[j] =  markers_start[brk_del.id];
		end[j] = markers_end[brk_del.next-1];
		size[j] = brk_del.size;
		mean[j] = (( (double) brk_del.sum) / ( (double) brk_del.size));
		
		label[j] = 0;
 		if(mean[j] < (double) -0.2){
 			label[j] = -1;
 		}else{
 			if(mean[j] > (double) 0.2){
 				label[j] = 1;
 			}
 		}
		i=brk_del.next;
		
	}

	*n_reg = j;
	heap_finalize(&b);
  
  
}



void delete_outlier(breakpoint *brks, int num_of_brks, int *markers_start, int *markers_end, double *positions){

	int i;
	binary_heap s;
	node n;
	breakpoint curr, prev, next;
	double tmp_lambda;


	heap_initalize(&s,num_of_probes);

	for(i=0;i<num_of_probes;){
		heap_insert(&s, node_create(i, (float) -(brks[i].size)));
		i=brks[i].next;
	}
	
	while( (s.heap_size>1) && (heap_max(&s)->p > min_region_size) ){
		n = heap_extract_max(&s);
		curr = brks[n.id];
		/* Check if the current breakpoint is valid*/
		if(curr.prev!=-2){
			
			/* case 1: current breakpoint is b0*/
			if(curr.id==0){
				next = brks[curr.next];

				/*Update the prev/next pointers*/
				brks[curr.id].next = next.next;
				brks[next.id].prev = curr.id;
				/*Update the data*/
				brks[curr.id].sum = curr.sum + next.sum;
				brks[curr.id].size = curr.size + next.size;
				tmp_lambda = update_priority(brks, curr.id, markers_start, markers_end, positions);
				curr.p = tmp_lambda;
				if(next.next!=num_of_probes){
					tmp_lambda = update_priority(brks, next.next, markers_start, markers_end, positions);
					brks[next.next].p = tmp_lambda;
				}
				brks[next.id].prev = -2;
				if(-(curr.size) > min_region_size){
					heap_insert(&s, node_create(curr.id, (float) -(curr.size)));
				}
				
			}else{
				/* If there are only two regions we link b0 and bN*/
				if( s.heap_size==2 ){
					brks[0].next = num_of_probes;
					brks[num_of_probes].prev = 0;
					brks[0].sum = brks[0].sum + curr.sum;
					brks[0].size = brks[0].size + curr.size;
					curr.prev = -2;
				}else{
					/* In this case we must delete the last but one breakpoint*/
					if(curr.next==num_of_probes){
						prev = brks[curr.prev];	
						brks[prev.id].next = num_of_probes;
						brks[num_of_probes-1].prev = prev.id;
						brks[prev.id].sum = prev.sum+curr.sum;
						brks[prev.id].size = prev.size+curr.size;
						tmp_lambda = update_priority(brks, prev.id, markers_start, markers_end, positions);
						brks[prev.id].p = tmp_lambda;
					}else{
						/* Here we consider all remaning cases*/
						next = brks[curr.next];
					
						if( curr.p >= next.p ){
							prev = brks[curr.prev];
							brks[prev.id].next = curr.next;
							brks[next.id].prev = prev.id;
							brks[prev.id].sum = prev.sum + curr.sum;
							brks[prev.id].size = prev.size + curr.size;
							tmp_lambda = update_priority(brks, prev.id, markers_start, markers_end, positions);
							brks[prev.id].p = tmp_lambda;
							tmp_lambda = update_priority(brks, next.id, markers_start, markers_end, positions);
							brks[next.id].p = tmp_lambda;
					
						}else{
							brks[curr.id].next = next.next;
							brks[next.next].prev = curr.id;
							brks[curr.id].sum = curr.sum + next.sum;
							brks[curr.id].size = curr.size + next.size;
							brks[next.id].prev = -2;
							tmp_lambda = update_priority(brks, curr.id, markers_start, markers_end, positions);
							brks[curr.id].p = tmp_lambda;
							if(-(curr.size) > min_region_size){
								heap_insert(&s, node_create(curr.id, (float) -(brks[curr.id].size)));
							}
							if(next.next!=num_of_probes){
								tmp_lambda = update_priority(brks, next.next, markers_start, markers_end, positions);
								brks[next.next].p = tmp_lambda;
							}
						}
					}
				}
			}

		}
	}
}


double update_priority(breakpoint *brks, int i, int *markers_start, int *markers_end, double *positions){
	double lambda, first_term, second_term, prev_mean, curr_mean, next_mean, third_term;
	int prev_size, curr_size, prev_id, dist;
	double tmp, mean_dist;

	prev_id = brks[i].prev;	
	prev_size = brks[prev_id].size;
	prev_mean = (double) (brks[prev_id].sum)/(prev_size);
	
	curr_size = brks[i].size;
	curr_mean = (double) (brks[i].sum)/(curr_size);

 
 /* Compute local average distance between exons */
 
 //Compute the local averages
  double avg_dist_prev = compute_local_mean_distance(brks, positions, prev_id);
  double avg_dist_curr = compute_local_mean_distance(brks, positions, i);
  double diff= (double)abs(avg_dist_prev-avg_dist_curr);


	first_term = (double) (prev_size * curr_size)/(prev_size + curr_size);
	tmp = prev_mean-curr_mean;
	if(tmp<0) tmp=-1*tmp;
	
	second_term = pow(tmp,2);
  
  if(diff==0){third_term=0;}
  else third_term = (log10(diff))*alpha;
  
	lambda = first_term*second_term +third_term;
	return(-lambda);
}

/* Added: function to compute the average distance between exons in a region, given a breakpoint */
double compute_local_mean_distance(breakpoint *brks, double *positions, int i){
  double local_mean_distance = 0, dist_sum = 0, diff = 0;
  double size = brks[i].size; //number of exons
  int j = 0;
  
  //compute the differences between the positions
  for (j = 0 ; j<(size-1)  && j<num_of_probes;j++){   
    dist_sum+=(positions[i+j+1]-positions[i+j]);
  }
  
  //divided by number of probes, or length of the region (position of start and end)
  local_mean_distance = dist_sum/size;
  return local_mean_distance;
}

void init_trivial_segmentation(breakpoint *brks, double *data, binary_heap *b){

	int i;
	double tmp_lambda, data_prev, data_curr ;
	double a, c;
	/*The breakpoint 0 has no prev*/
	brks[0] = breakpoint_create(0, -1, floor(*(data)*ROUND_CNST)/ROUND_CNST, 1, -1, 1);

	/*The breakpoint N has no next/prev and no data it is just a closure*/
	brks[num_of_probes] = breakpoint_create(num_of_probes, -1, -1, -1, -1, -1);
	
	/*Initialize all other nodes in the trivial segmentation*/
	for (i=1; i<num_of_probes; i++){
		data_prev = floor(*(data+(i-1))*ROUND_CNST)/ROUND_CNST;
		data_curr = floor(*(data+(i))*ROUND_CNST)/ROUND_CNST;
		a = data_prev-data_curr;
		if(a<0){
			a= -1*a;
			}
		c = pow(a,2);
		tmp_lambda = 0.5*c;
	
		brks[i] = breakpoint_create(i, -tmp_lambda, floor(*(data+(i))*ROUND_CNST)/ROUND_CNST, 1, i-1, i+1);
		
		heap_insert(b,node_create(i, -tmp_lambda));
	}	
	
}


