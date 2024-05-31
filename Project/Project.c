#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>
#include <time.h>
#include <limits.h>

#define MAX_LINE_LENGTH 100 // Maximum length of a line in the input file
#define HASH_SIZE 100000 // Size of the hash table
#define NUM_THREADS 4 //number of threads
#define NUM_PAIRS 10 // number of pairs
#define K 5

// Node structure for adjacency list
typedef struct Node {
    int vertex;
    struct Node* next;
} Node;

// Distance matrix 
typedef struct DistanceEntry {
    int fromNode;
    int toNode;
    int distance;
    struct DistanceEntry* next;
} DistanceEntry;

// Graph structure
typedef struct {
    int numNodes;
    Node* adjacencyList;
    DistanceEntry** distanceMatrix;
} Graph;

// Path structure
typedef struct {
    int length;
    int* nodes;
} Path;

// Function to initialize a graph
Graph* initializeGraph(int numNodes) {
    Graph* graph = (Graph*)malloc(sizeof(Graph));
    if (graph == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }
    graph->numNodes = numNodes;
    graph->adjacencyList = (Node*)malloc(numNodes * sizeof(Node));
    if (graph->adjacencyList == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }
    for (int i = 0; i < numNodes; i++) {
        graph->adjacencyList[i].next = NULL;
    }
    graph->distanceMatrix = (DistanceEntry*)malloc(HASH_SIZE * sizeof(DistanceEntry));
    if (graph->distanceMatrix == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }
    for (int i = 0; i < HASH_SIZE; i++) {
        graph->distanceMatrix[i] = NULL;
    }
    return graph;
}

// Function to add an edge to the graph
void addEdge(Graph* graph, int fromNode, int toNode) {
    // Update the adjacency list
    Node* newNode = (Node*)malloc(sizeof(Node));
    if (newNode == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }
    newNode->vertex = toNode;
    newNode->next = graph->adjacencyList[fromNode].next;
    graph->adjacencyList[fromNode].next = newNode;

    // Updating the distance matrix
    int index = (fromNode * 31 + toNode) % HASH_SIZE;
    DistanceEntry* entry = (DistanceEntry*)malloc(sizeof(DistanceEntry));
    if (entry == NULL) {
        printf("Memory allocation failed.\n");
        exit(1);
    }
    entry->fromNode = fromNode;
    entry->toNode = toNode;

    // Assigns 0 for diagonal vertices and 1 for other vertices
    entry->distance = (fromNode == toNode) ? 0 : 1;
    entry->next = graph->distanceMatrix[index];
    graph->distanceMatrix[index] = entry;
}


// Free memory allocated by the graph
void freeGraph(Graph* graph) {
    free(graph->adjacencyList);
    for (int i = 0; i < HASH_SIZE; i++) {
        DistanceEntry* currentEntry = graph->distanceMatrix[i];
        while (currentEntry != NULL) {
            DistanceEntry* temp = currentEntry;
            currentEntry = currentEntry->next;
            free(temp);
        }
    }
    free(graph->distanceMatrix);
    free(graph);
}

//Initialize the graph with edges after reading it from a file
void initializeGraphFromFile(Graph* graph, const char* filename) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        printf("Error opening file.\n");
        exit(1);
    }

    char line[MAX_LINE_LENGTH];
    // Skiping the first four lines
    for (int i = 0; i < 4; i++) {
        fgets(line, MAX_LINE_LENGTH, file);
    }

    // Reading the file line by line
    int fromNode, toNode;
    while (fscanf(file, "%d %d", &fromNode, &toNode) == 2) {
        addEdge(graph, fromNode, toNode);
    }

    fclose(file);
}

// Function to find the shortest path using Dijkstra's algorithm
Path dijkstra(DistanceEntry* graph[HASH_SIZE], int numNodes, int source, int destination) {
    int distance[numNodes];
    int predecessor[numNodes];
    int visited[numNodes];
    for (int i = 0; i < numNodes; i++) {
        distance[i] = INT_MAX;
        predecessor[i] = -1;
        visited[i] = 0;
    }
    distance[source] = 0;

    for (int count = 0; count < numNodes - 1; count++) {
        int min_distance = INT_MAX;
        int min_index = -1;

        for (int v = 0; v < numNodes; v++) {
            if (!visited[v] && distance[v] < min_distance) {
                min_distance = distance[v];
                min_index = v;
            }
        }

        if (min_index == -1)
            break;

        visited[min_index] = 1;

        DistanceEntry* entry = graph[min_index];
        while (entry != NULL) {
            int v = entry->toNode;
            if (!visited[v] && distance[min_index] != INT_MAX &&
                distance[min_index] + entry->distance < distance[v]) {
                distance[v] = distance[min_index] + entry->distance;
                predecessor[v] = min_index;
            }
            entry = entry->next;
        }
    }

    Path shortest_path;
    shortest_path.length = 0;
    shortest_path.nodes = (int*)malloc(numNodes * sizeof(int));
    int current = destination;
    while (current != -1 && current != source) {
        shortest_path.nodes[shortest_path.length++] = current;
        current = predecessor[current];
    }
    shortest_path.nodes[shortest_path.length++] = source;

    for (int i = 0; i < shortest_path.length / 2; i++) {
        int temp = shortest_path.nodes[i];
        shortest_path.nodes[i] = shortest_path.nodes[shortest_path.length - i - 1];
        shortest_path.nodes[shortest_path.length - i - 1] = temp;
    }

    return shortest_path;
}

//Printing the path
void print_path(Path path) {
    printf("Path length: %d, Nodes: ", path.length);
    for (int i = 0; i < path.length; i++) {
        printf("%d ", path.nodes[i]);
    }
    printf("\n");
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("Usage: %s <filename>\n", argv[0]);
        return 1;
    }

    const char* filename = argv[1];
    int numNodes;

    MPI_Init(&argc, &argv);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (world_rank == 0) {
        printf("Choose the graph you want to use:\n");
        printf("1. CSV Graph (doctorwho.csv)\n");
        printf("2. Text Graph (Email-EuAll.txt)\n");
        printf("3. Text Graph (Email-Enron.txt)\n");


        int choice;
        scanf("%d", &choice);

        switch (choice) {
            case 1:
                numNodes = 100; // Assuming the number of nodes for the CSV graph
                break;
            case 2:
                numNodes = 265214; // Assuming the number of nodes for Text Graph 
                break;
            case 3:
                numNodes = 36692; // Assuming the number of nodes for Text Graph
                break;
            default:
                printf("Invalid choice.\n");
                MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    MPI_Bcast(&numNodes, 1, MPI_INT, 0, MPI_COMM_WORLD);

    Graph* graph = initializeGraph(numNodes);
    initializeGraphFromFile(graph, filename);

    int source, destination;
    int pair_source[NUM_PAIRS], pair_destination[NUM_PAIRS];

    srand(time(NULL) + world_rank); 

    // Generate random pairs of source and destination nodes
    for (int i = 0; i < NUM_PAIRS; i++) {
        pair_source[i] = rand() % numNodes;
        pair_destination[i] = rand() % numNodes;
    }

    // Calculating the number of pairs to process per node
    int pairs_per_node = NUM_PAIRS / world_size;
    int start_index = world_rank * pairs_per_node;
    int end_index = (world_rank + 1) * pairs_per_node;


    // Process pairs assigned to this node
    for (int i = start_index; i < end_index; i++) {
        source = pair_source[i];
        destination = pair_destination[i];
        printf("Node %d calculating shortest paths from node %d to node %d\n", world_rank, source, destination);
        printf("\n");
        #pragma omp parallel for num_threads(NUM_THREADS)
        for (int j = 0; j < K; j++) {
            Path shortest_path = dijkstra(graph->distanceMatrix, graph->numNodes, source, destination);
            printf("Shortest path %d from node %d to node %d: \n", j + 1, source, destination);
            printf("\n");
            print_path(shortest_path);
            
        }
    }

    // Freeing memory
    freeGraph(graph);
    MPI_Finalize();
    return 0;
}