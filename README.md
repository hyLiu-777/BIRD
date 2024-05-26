# BIRD
This repo contains codes for paper **"BIRD: Efficient Approximation of Bidirectional Hidden Personalized PageRank"**, PVLDB 2024.

## Abstract
In bipartite graph analysis, similarity measures play a pivotal role in various applications. Among existing metrics, the Bidirectional Hidden Personalized PageRank (BHPP) stands out for its superior query quality. However, the computational expense of BHPP remains a bottleneck. Existing approximation methods either demand significant matrix storage or incur prohibitive time costs. For example, current state-of-the-art methods require over 3 hours to process a single-source BHPP query on the real-world bipartite graph *Orkut*, which contains approximately $3\times10^8$ edges. We introduce **BIRD**, a novel algorithm designed for answering single-source BHPP queries on weighted bipartite graphs. Through meticulous theoretical analysis, we demonstrate that BIRD significantly improves time complexity to $\widetilde{\mathcal{O}}(n)$, as compared to the previous best one, $\widetilde{\mathcal{O}}(m)$, under typical relative error setting and constant failure probability. ($n,m$ denote the number of nodes and edges respectively.) Extensive experiments confirm that BIRD outperforms existing baselines by orders of magnitude in large-scale bipartite graphs. Notably, our proposed method accomplishes a single-source BHPP query on *Orkut* using merely 7 minutes.

## Generating queries.
```python
$ cd data/
$ python genseed.py avito 27736   # python genseed.py data_name |U|
```

## Running 
```shell
$ sh build.sh
$ ./bppr -f ../data/ -g avito -a BDPush -e 0.5 --querynum 100 # -f data path -g graph name -a algorithm 
```
