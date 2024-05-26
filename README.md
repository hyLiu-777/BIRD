# BIRD
This repo contains codes for paper "BIRD: Efficient Approximation of Bidirectional Hidden Personalized PageRank".

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
