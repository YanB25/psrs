#include <stdio.h> 
#include <stdlib.h> 

// Node 
typedef struct node { 
	int data; 
	// Lower values indicate higher priority 
	int priority; 
	struct node* next; 
} Node; 