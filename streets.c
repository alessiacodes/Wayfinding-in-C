#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "streets.h"

/* Struct Definitions */
struct node {
    int id;
    double lat;
    double lon;
    int num_ways;
    int *way_ids; 
};

struct way {
    int id;
    char *name;
    float maxspeed;
    bool oneway;
    int num_nodes;
    int *node_ids; 
};                  

struct ssmap {
    int num_nodes;
    int num_ways;
    struct node **nodes;
    struct way **ways;
};

/* Min Heap Data Structure */
struct minheap_pair{
    int node_id;
    double priority;
};

struct minheap_prio_queue{
    struct minheap_pair **heap_array; 
    int heap_size;
    int heap_capacity;
};

        struct minheap_prio_queue *create_prio_queue(int capacity)
        {   
            struct minheap_prio_queue*minheap_ptr = malloc(sizeof(struct minheap_prio_queue));
                if (minheap_ptr == NULL){
                    free(minheap_ptr);
                    return NULL;
                }
            minheap_ptr->heap_array = malloc(sizeof(struct minheap_pair *) * capacity);
                if (minheap_ptr->heap_array == NULL){
                    free(minheap_ptr->heap_array);
                    return NULL;
                }
            minheap_ptr->heap_capacity = capacity;
        minheap_ptr->heap_size = 0;
        return minheap_ptr;
    }

    void minheap_destroy(struct minheap_prio_queue *minheap)
    {
        for (int i = 0; i < minheap->heap_size; i++){ 
            free((minheap->heap_array)[i]);
        }
        free(minheap->heap_array);
        free(minheap);
    }

    void minheapify(struct minheap_prio_queue *minheap, int parent_index) 
    {   
        int smallest = parent_index;
        int left = 2 * parent_index + 1;  // Left child index
        int right = 2 * parent_index + 2; // Right child index

        if (left < minheap->heap_size - 1 && ((minheap->heap_array)[parent_index])->priority > ((minheap->heap_array)[left])->priority){
            smallest = left;
        }

        if (right < minheap->heap_size - 1 && ((minheap->heap_array)[smallest])->priority > ((minheap->heap_array)[right])->priority){
            smallest = right;
        }

        if (smallest != parent_index){
            struct minheap_pair *hold = (minheap->heap_array)[parent_index];
            (minheap->heap_array)[parent_index] = (minheap->heap_array)[smallest];
            (minheap->heap_array)[smallest] = hold;
            minheapify(minheap, smallest);
        }
    }

    struct minheap_pair* extract_min(struct minheap_prio_queue *minheap)
    {
        // If the heap is empty, return an empty node
        if (minheap->heap_size == 0) {
            return NULL;
        }
        
        else if (minheap->heap_size == 1){
            minheap->heap_size -= 1;
            return (minheap->heap_array)[0];
        }

        else{
            struct minheap_pair* hold = (minheap->heap_array)[0];
            
            // Replace root with last element;
            (minheap->heap_array)[0] = (minheap->heap_array)[minheap->heap_size - 1];
            (minheap->heap_array)[minheap->heap_size - 1] = NULL;
            minheap->heap_size -= 1;
            minheapify(minheap, 0);
            return hold;
        }
    }

    void bubble_up(struct minheap_prio_queue*minheap, int i)
    {
        while (i > 0){
            // check if heap property satisfied, then early return
            if (((minheap->heap_array)[i])->priority >= ((minheap->heap_array)[(i - 1) / 2])->priority){
                return;
            }
            // othw. swap parent with current node
            else{
                struct minheap_pair* hold = (minheap->heap_array)[i];
                // Replace root with last element;
                (minheap->heap_array)[i] = (minheap->heap_array)[(i - 1) / 2];
                (minheap->heap_array)[(i - 1) / 2] = hold;
                i = (i - 1) / 2;
            }
        }
    }

    void add_with_priority(struct minheap_prio_queue *minheap, int node_id, double priority)
    {
        // Check if minheap at capacity
        if (minheap->heap_size == minheap->heap_capacity) {
            printf("Oops! Capacity reached :(\n");
            return;
        }

        struct minheap_pair* new_pair = malloc(sizeof(struct minheap_pair));
            if (new_pair == NULL){
                free(new_pair);
                return;
            }
        
        new_pair->node_id = node_id;
        new_pair->priority = priority;

        // Ensure minheap property is satisfied
        int i = minheap->heap_size;
        minheap->heap_size++;
        (minheap->heap_array)[i] = new_pair;
        bubble_up(minheap, i);
    }

    void decrease_priority(struct minheap_prio_queue*minheap, int node_id, double new_prio)
    {
        // Find which index node_id is located at 
        int index = -1;
        for (int i = 0; i < minheap->heap_size; i++){
            if (((minheap->heap_array)[i])->node_id == node_id){
                index = i;
                break;
            }
        }

        // Check if that node_id actually exists
        if (index == -1){
            printf("Oops! that node doesn't exist :(\n");
            return;
        }

        // Change prio and ensure minheap property satisfied
        ((minheap->heap_array)[index])->priority = new_prio;
        bubble_up(minheap, index);
    }

/* HELPERS */

int 
check_way_association(const struct ssmap * m, int node_id, const char * keyword, int exclusion, int start_index){
    for (int i = start_index; i < (((m->nodes)[node_id])->num_ways); i++){
        int way_id = (((m->nodes)[node_id])->way_ids)[i];

        // See if occurrence of keyword occurs within name
        char *occurence = strstr(((m->ways)[way_id])->name, keyword);

        if(occurence != NULL && i != exclusion){
            return i;
        }
    }
    return -1;
}

void
print_id(bool *first_id, int id)
{
    if(!(*first_id)){
        printf("%d \n", id);
        *first_id = true;
    }
    else{
        printf("%d ", id);
    }
}

int 
find_way_with_nodes(const struct ssmap *m, int node_id1, int node_id2) 
{
    // Dereference node object
    struct node node1 = *((m->nodes)[node_id1]);
    struct node node2 = *((m->nodes)[node_id2]);

    // Search for a way_id associated with both node1 and node2
    for (int i = 0; i < node1.num_ways; i++){
        for (int j = 0; j < node2.num_ways; j++){
            if (node1.way_ids[i] == node2.way_ids[j]){
                return node1.way_ids[i];
            }
        }
    }

    // If search fails, return -1
    return -1;
}

/* Returns whether or not a pair of nodes are validly adjacent (i.e. adhere to oneway) given node1
   and node2. Assume node1 and node2 are in the order of the path. Will print required error
   message as well. */
bool check_adjacency(const struct ssmap *m, int node_id1, int node_id2, int way_id){
    struct way way = *((m->ways)[way_id]);
    bool oneway = way.oneway;

    // Check for possible valid adjacency (error 3)
    bool possible_valid_adj = false;
    int store_first_node = -1; // Store the first node to determine order of pathing for oneway
    for (int i = 0; i < way.num_nodes - 1; i++){
        if ((way.node_ids[i] == node_id1) && (way.node_ids[i + 1] == node_id2)){
            possible_valid_adj = true;
            store_first_node = node_id1;
            break;
        }

        else if ((way.node_ids[i] == node_id2) && (way.node_ids[i + 1] == node_id1)){
            possible_valid_adj = true;
            store_first_node = node_id2;
            break;
        }
    }

    if (possible_valid_adj == false){
        printf("error: cannot go directly from node %d to node %d.\n", node_id1, node_id2);
        return false;
    }
    else{
        if (oneway && store_first_node != node_id1){
            printf("error: cannot go in reverse from node %d to node %d.\n", node_id1, node_id2);
            return false;
        }
    }

    return true;
}

int
is_duplicates(int size, int node_ids[size])
{
    for (int i = 0; i < size; i++){
        for (int j = i + 1; j < size; j++){
            if (node_ids[i] == node_ids[j]){
                return node_ids[i];
            }
        }
    }
    return -1;
}

bool
check_valid_node(const struct ssmap *m, int node_id)
{
    if (node_id >= 0 && node_id < m->num_nodes){
        return true;
    }
    return false;

}

void 
print_sequence(int size, const int *sequence)
{
    // Print everything but last element for formatting purposes.
    for (int i = 0; i < size - 1; i++){
        printf("%d ", sequence[i]);
    }
    printf("%d\n", sequence[size - 1]);
}

/* Official Assignment Functions */

struct ssmap * 
ssmap_create(int nr_nodes, int nr_ways)
{
    struct ssmap *m_ptr = malloc(sizeof(struct ssmap));
        // Handle malloc error
        if (m_ptr == NULL){
            free(m_ptr);
            return NULL;
        }

    // Load Values into map

        // Logic Error Checks
        if (nr_nodes == 0 || nr_ways == 0){
            return NULL;
        }

    m_ptr->num_nodes = nr_nodes;
    m_ptr->num_ways = nr_ways;

    m_ptr->nodes = malloc(sizeof(struct node *) * nr_nodes);
        if (m_ptr->nodes == NULL){
            free(m_ptr->nodes);
            return NULL;
        }
        // Wipe memory
        for (int i = 0; i < nr_nodes; i++){
            (m_ptr->nodes)[i] = NULL;
        }
    m_ptr->ways = malloc(sizeof(struct ways *) * nr_ways);
        if (m_ptr->ways == NULL){
            free(m_ptr->ways);
            return NULL;
        }
        // Wipe memory (Note to TA: I know this is repetitive, but I wasn't able
        // to make a generic "wipe_memory " function without the compiler complaining
        // that I was dereferencing a generic type).
        for (int i = 0; i < nr_ways; i++){
            (m_ptr->ways)[i] = NULL;
        }
    return m_ptr;
}

bool
ssmap_initialize(struct ssmap * m)
{
    return true;
}

void
ssmap_destroy(struct ssmap * m)
{
    // Free Way-Related Memory
    for (int i = 0; i < m->num_ways; i++){
        free((m->ways)[i]->name);
        free((m->ways)[i]->node_ids);
        free((m->ways)[i]);
    }
    free((m->ways));

    // Free Node-related memory
    for (int i = 0; i < m->num_nodes; i++){
        free((m->nodes)[i]->way_ids);
        free((m->nodes)[i]);
    }
    free((m->nodes));

    // Top-level clearing
    free(m);
}

struct way * 
ssmap_add_way(struct ssmap * m, int id, const char * name, float maxspeed, bool oneway, 
              int num_nodes, const int node_ids[num_nodes])
{
    // Load in values into new_way
    struct way *new_way_ptr = malloc(sizeof(struct way)); // ⚠️ this a source of memory leaks unless we readjust so struct has a pointer
        if (new_way_ptr == NULL){
            free(new_way_ptr);
            return NULL;
        }
    struct way new_way;
    new_way.id = id;
    new_way.name = malloc(sizeof(char) * (strlen(name) + 1));
        if (new_way.name == NULL){
            free(new_way.name);
            return NULL;
        }
    strcpy(new_way.name, name);
    new_way.maxspeed = maxspeed;
    new_way.oneway = oneway;
    new_way.num_nodes = num_nodes;
    new_way.node_ids = malloc(sizeof(int) * num_nodes);
        if (new_way.node_ids == NULL){
            free(new_way.node_ids);
            return NULL;
        }

    // Copy values from node_ids to the allocated array
    for (int i = 0; i < num_nodes; i++){
        new_way.node_ids[i] = node_ids[i];
    }
    
    // Add new_way to ssmap
    (m->ways)[id] = new_way_ptr; 
    *new_way_ptr = new_way;
    return new_way_ptr;
}

struct node * 
ssmap_add_node(struct ssmap * m, int id, double lat, double lon, 
               int num_ways, const int way_ids[num_ways])
{
    // Load in values into new_node
    struct node *new_node_ptr = malloc(sizeof(struct node));
        if (new_node_ptr == NULL){
            free(new_node_ptr);
            return NULL;
        }
    struct node new_node;
    new_node.id = id;
    new_node.lat = lat;
    new_node.lon = lon;
    new_node.num_ways = num_ways;
    new_node.way_ids = malloc(sizeof(int) * num_ways);
        if (new_node.way_ids == NULL){
            free(new_node.way_ids);
            return NULL;
        }

    // Copy values from node_ids to the allocated array
    for (int i = 0; i < num_ways; i++){
        new_node.way_ids[i] = way_ids[i];
    }
    
    // Add new_way to ssmap
    (m->nodes)[id] = new_node_ptr; // ⚠️ POSSIBLE ISSUE: Can you set it equal like this?
    *new_node_ptr = new_node;
    return new_node_ptr;
}

void
ssmap_print_way(const struct ssmap * m, int id)
{
    // Check if ID exists
    if (id < 0 || id >= (m->num_ways)){
        printf("error: way %d does not exist.\n", id);
        return;
    }
    else if ((m->ways)[id] == NULL){
        printf("error: way %d does not exist.\n", id);
        return;
    }
    else{
        printf("Way %d: %s\n", id, ((m->ways)[id])->name);
        return;
    }
}

void
ssmap_print_node(const struct ssmap * m, int id)
{
    // Check if ID exists
    if (id < 0 || id >= (m->num_nodes)){
        printf("error: node %d does not exist.\n", id);
        return;
    }
    else if ((m->nodes)[id] == NULL){
        printf("error: node %d does not exist.\n", id);
        return;
    }
    else{
        printf("Node %d: (%.7lf, %.7lf)\n", id, ((m->nodes)[id])->lat, ((m->nodes)[id])->lon);
        return;
    }
}

void 
ssmap_find_way_by_name(const struct ssmap * m, const char * name)
{
    for (int i = 0; i < (m->num_ways); i++){
        // See if occurrence of keyword occurs within name
        char *occurence = strstr(((m->ways)[i])->name, name);

        bool *first_id = malloc(sizeof(bool));
            if (first_id == NULL){
                free(first_id);
                return;
            }
        
        *first_id = true;

        // If keyword does exist, print
        if (occurence != NULL){
            print_id(first_id, i);
        }

        free(first_id);
    }
    printf("\n");
}

void 
ssmap_find_node_by_names(const struct ssmap * m, const char * name1, const char * name2)
{
    bool *first_id = malloc(sizeof(bool));
         if (first_id == NULL){
                free(first_id);
                return;
            }
    *first_id = true;

    if (name2 == NULL){
        for (int i = 0; i < (m->num_nodes); i++){
            if (check_way_association(m, i, name1, -1, 0) != -1){
                print_id(first_id, i);
            }
        }
    }
    else{
        for (int i = 0; i < (m->num_nodes); i++){
            // Store possible way contenders, if these stay -1 we know that the conditions don't hold
            int possible_way_id1 = check_way_association(m, i, name1, -1, 0);
            int possible_way_id2;
            if (possible_way_id1 != -1){
                // If possible way1 does not contain name2 as a keyword, we can proceed searching for name2
                if (strstr(((m->ways)[possible_way_id1])->name, name2) == NULL){
                    possible_way_id2 = check_way_association(m, i, name2, possible_way_id1, 0);
                }

                // Otherwise, we should search for a unique occurence of name1.
                // We know that we know that name2 occurs in the node we just found,
                // So we can proceed searching from where we left off searching for
                // name 1.
                else{
                    possible_way_id2 = check_way_association(m, i, name1, -1, possible_way_id1 + 1);
                }
            }

            if (possible_way_id1 != -1 && possible_way_id2 != -1){
                print_id(first_id, i);
            }

        }
    }

    printf("\n");
    free(first_id);
}

/**
 * Converts from degree to radian
 *
 * @param deg The angle in degrees.
 * @return the equivalent value in radian
 */
#define d2r(deg) ((deg) * M_PI/180.)

/**
 * Calculates the distance between two nodes using the Haversine formula.
 *
 * @param x The first node.
 * @param y the second node.
 * @return the distance between two nodes, in kilometre.
 */
static double
distance_between_nodes(const struct node * x, const struct node * y) {
    double R = 6371.;       
    double lat1 = x->lat;
    double lon1 = x->lon;
    double lat2 = y->lat;
    double lon2 = y->lon;
    double dlat = d2r(lat2-lat1); 
    double dlon = d2r(lon2-lon1); 
    double a = pow(sin(dlat/2), 2) + cos(d2r(lat1)) * cos(d2r(lat2)) * pow(sin(dlon/2), 2);
    double c = 2 * atan2(sqrt(a), sqrt(1-a)); 
    return R * c; 
}

double 
ssmap_path_travel_time(const struct ssmap * m, int size, int node_ids[size])
{
    double sum = 0.0;

    for (int i = 0; i < size - 1; i++){
        // Error 1 check
        if (check_valid_node(m, node_ids[i]) == false){
            printf("error: node %d does not exist.\n", node_ids[i]);
            return -1.0;
        }

        // Also check if the one ahead is invalid to prevent future seg faults
        else if (check_valid_node(m, node_ids[i + 1]) == false){
            printf("error: node %d does not exist.\n", node_ids[i + 1]);
            return -1.0;
        }

        // Error 2 check
        int way_with_nodes = find_way_with_nodes(m, node_ids[i], node_ids[i + 1]); 
        if (way_with_nodes == -1){
            printf("error: there are no roads between node %d and node %d.\n", node_ids[i], node_ids[i + 1]);
            return -1.0;
        }

        // Error 3,4 check (note that the check_adjacency function will print the error message)
        else if (!check_adjacency(m, node_ids[i], node_ids[i + 1], way_with_nodes)){
            return -1.0;
        }

        float maxspeed = ((m->ways)[way_with_nodes])->maxspeed;
        sum += distance_between_nodes((m->nodes)[node_ids[i]], (m->nodes)[node_ids[i + 1]]) / maxspeed;
    }

    // Error 5 Check
    int possible_duplicate = is_duplicates(size, node_ids);
    if (possible_duplicate != -1){
        printf("error: node %d appeared more than once.\n", possible_duplicate);
        return -1.0;
    }
    return (sum * 60);
}


/* Path Create Specific Helpers */

/* Takes in the source node and mutates the two lists to give desired info. Note that size of both arrays
   are implicitly 2 * m->num_ways.*/
void
find_valid_neighbors(const struct ssmap *m, int node_id, int *associated_way_id, int *valid_neighbor_ids)
{
    struct node src = *((m->nodes)[node_id]);

    int skip = 0;

    // Go through each way
    for (int i = 0; i < src.num_ways; i++){
        struct way way = *(m->ways)[src.way_ids[i]];

        int node_is_located = -1;
        // Find where our node is in the way
        for (int j = 0; j < way.num_nodes; j++){
            if (way.node_ids[j] == node_id){
                node_is_located = j;
                break;
            }
        }   

        // node is located at end and not one-way  -> can only check left
        if (node_is_located == way.num_nodes - 1 && !way.oneway && way.num_nodes > 1){
            valid_neighbor_ids[i + skip] = way.node_ids[node_is_located - 1];
            associated_way_id[i + skip] = way.id;
        }

        // One-way or node is located at beginning -> can only check right
        else if ((way.oneway || node_is_located == 0) && way.num_nodes > 1 && node_is_located != (way.num_nodes - 1)){
            valid_neighbor_ids[i + skip] = way.node_ids[node_is_located + 1];
            associated_way_id[i + skip] = way.id;
        }

        // Not end-point and no one-way restriction
        else if (!way.oneway && way.num_nodes > 2){
            // get left and right neighbor
            valid_neighbor_ids[i + skip] = way.node_ids[node_is_located + 1];
            associated_way_id[i + skip] = way.id;
            valid_neighbor_ids[i + skip + 1] = way.node_ids[node_is_located - 1];
            associated_way_id[i + skip + 1] = way.id;
            skip++;
        }

        else{
            valid_neighbor_ids[i + skip] = -2;
        }
    }

    valid_neighbor_ids[src.num_ways + skip] = -1;
}

int *
run_dijkstra(const struct ssmap * m, int start_id, int end_id)
{
    // Run Dijkstra's
    double *dist = malloc(sizeof(double) * m->num_nodes); // Note that index will correspond to node id
         if (dist == NULL){
                free(dist);
                return NULL;
            }
    int *prev = malloc(sizeof(double) * m->num_nodes);
         if (prev == NULL){
                free(prev);
                return NULL;
            }
    dist[start_id] = 0;

    struct minheap_prio_queue * prio_queue_ptr = create_prio_queue(m->num_nodes);
    for (int i = 0; i < m->num_nodes; i++){
        if (i != start_id){
            dist[i] = INFINITY;
            prev[i] = -1;
        }
        add_with_priority(prio_queue_ptr, i, dist[i]);
    }

    while (prio_queue_ptr->heap_size > 0){
        struct minheap_pair* u_pair = extract_min(prio_queue_ptr);
        int u = u_pair->node_id;
        free(u_pair);
        if (u == end_id){
            free(dist);
            minheap_destroy(prio_queue_ptr);
            return prev;
        }
        int *associated_way_id = malloc(sizeof(int) * 2 * m->num_ways);
            if (associated_way_id == NULL){
                free(associated_way_id);
                return NULL;
            }
        int *valid_neighbor_ids = malloc(sizeof(int) * 2 * m->num_ways);
            if (valid_neighbor_ids == NULL){
                free(valid_neighbor_ids);
                return NULL;
            }
        find_valid_neighbors(m, u, associated_way_id, valid_neighbor_ids); 

        int i = 0;
        while (valid_neighbor_ids[i] != -1){
            if (valid_neighbor_ids[i] != -2){ // Only run if neighbor is actually valid
                double alt = dist[u] + distance_between_nodes((m->nodes)[u], (m->nodes)[valid_neighbor_ids[i]]) / (m->ways)[associated_way_id[i]]->maxspeed; 
                if (alt < dist[valid_neighbor_ids[i]]){
                    dist[valid_neighbor_ids[i]] = alt;
                    prev[valid_neighbor_ids[i]] = u;
                    decrease_priority(prio_queue_ptr, valid_neighbor_ids[i], alt);
                }
            }
            i++;
        }
        free(associated_way_id);
        free(valid_neighbor_ids);
    }

    free(dist);
    minheap_destroy(prio_queue_ptr);
    return prev;
}

int 
get_path_to_target(int * prev, int source, int target, int **s)
{
    int size = 0;
    int u = target;
    if (prev[u] != -1 || u == source){
        while (u != -1){
            // insert u at beginning of sequence
            size++;
            int * new_s = malloc(sizeof(int) * (size));
                if (new_s == NULL){
                    free(new_s);
                    return -1;
                }
                
            new_s[0] =  u;
            for (int i = 1; i < size; i++){
                new_s[i] = (*s)[i-1];
            }
            free(*s);
            *s = new_s;
            u = prev[u];
            if ((*s)[0] == source){
                return size;
            }
        }
    }
    return size;
}

void 
ssmap_path_create(const struct ssmap * m, int start_id, int end_id)
{
    int * prev = run_dijkstra(m, start_id, end_id);
    int *s = NULL;
    
    int s_size = get_path_to_target(prev, start_id, end_id, &s);
    free(prev);
    
    // Error Checking
    int error = ssmap_path_travel_time(m, s_size, s);
    if (error != -1){
        print_sequence(s_size, s);
    }
    else{
        printf("error: could not find a path from node %d to node %d.\n", start_id, end_id);
    }

    free(s);
    return;
}

