# BDPush
Code for submission 

"Efficient Approximation of Bidirectional Hidden Personalized
PageRank: Breaking the $\tilde{\mathcal{O}}(m)$ Time Barrier"

## Generating queries.
```python
$ cd data/
$ python genseed.py avito 27736   # python data_name |U|
```

## Running 
```shell
$ sh build.sh
$ ./bppr -f ../data/ -g avito -a BDPush -e 0.5 --querynum 100 # -f data path -g graph name -a algorithm 
```
