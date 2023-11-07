# FaSTest-Subgraph-Cardinality-Estimation
Supplementary materials for 
**Cardinality Estimation of Subgraph Matching: A Filtering-Sampling Approach** (Submitted to VLDB 2024)

### Dependencies
- Boost Library
- g++ Compiler with C++17 support

### Usage 
Build with the command
```sh
mkdir build
cd build
cmake .. && make  
```

The configuration for the filtering and the upper bound of sample size for the graph sampling can be modified by command line arguments.

|     Argument      |                          Description                           |
|:-----------------:|:--------------------------------------------------------------:|
|        -D         |                        Dataset to solve                        |
|        -A         |               Solve for all queries in [dataset]               |
| -d name, -q name  |            Used for solving one query specifically             |
|       -s K        |        Upper bound for graph sampling is set to Vq * K         |
| --NEIGHBORHOOD NS |        Neighborhood filter condition : Neighbor Safety         |
| --NEIGHBORHOOD NB |   Neighborhood filter condition : ESIC    |
| --NEIGHBORHOOD EB |   Neighborhood filter condition : Edge Bipartite Safety   |
|   --STRUCTURE X   |                   No substructure condition                    |
|   --STRUCTURE 3   |        Substructure filter condition : Triangle Safety         |
|   --STRUCTURE 4   | Substructure filter condition : Triangle and Four-cycle Safety |

Execution example with `yeast` dataset:
```sh
./main/FaSTest -D yeast -A
```
Which will first read queries from `yeast_ans.txt`, and evaluate all queries whose true cardinalities are known. 

The command 
```sh
./main/FaSTest -D yeast -A --NEIGHBORHOOD NB --STRUCTURE 3
```
disables Four-cycle safety, and uses only ESIC (without Edge Bipartite Safety).

To run single query to single data, use:
```sh
./main/FaSTest -D yeast -d ../dataset/yeast/data_graph/yeast.graph -q ../dataset/yeast/query_graph/query_sparse_32_200.graph
```
If the true cardinality is not known, True will show as -1. (ex: query_sparse_32_97.graph)

### Datasets
The datasets and query graphs from [RapidMatch](https://github.com/RapidsAtHKUST/RapidMatch/) is used for evaluation.

#### Input Format 
```
t [#Vertex] [#Edge]
v [ID] [Label] [Degree]
v [ID] [Label] [Degree]
...
e 1235 2586 [EdgeLabel]
```
Edge labels are optional (missing labels are considered as zero).

Example (0-0-1-0 labeled path)
```
t 4 3
v 0 0 1
v 1 0 2
v 2 1 2
v 3 0 1
e 0 1
e 1 2
e 2 3
```

#### True Cardinalities
To assess accuracy, true cardinalities should be provided. We provide the ground truth values for 1,707 out of 1,800 queries for yeast dataset in `dataset/yeast/yeast_ans.txt`, computed with [DAF](https://github.com/SNUCSE-CTA/DAF). 
