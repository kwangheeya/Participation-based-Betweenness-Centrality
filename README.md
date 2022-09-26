# Participation-based Betweenness Centrality (PBC)


##### Contributors: [Kwang Hee Lee](https://lekwanghee.github.io/) in Database Lab. at KAIST.


Subtitle: Betweenness centrality in directed hypergraphs

Code for the papers
* "[Computing Betweenness Centrality in B-hypergraphs](http://dx.doi.org/10.1145/3132847.3133093)" published in the proceedings of the 2017 ACM on CIKM.
* "[Linearization of Dependency and Sampling for Participation-based Betweenness Centrality in Very Large B-hypergraphs](https://doi.org/10.1145/3375399)" published in ACM TKDD 2020.

## Environment
* Windows 10, Python 3.6.2 64bit

## Requirement
* [numpy](http://www.numpy.org/) package 
* [statistics] package
* [scipy] package
* [collections] package
* [contextlib] package

## Data (B-hypergraph) format
> The data used in this code are stored in folder 'src/data'.

* In data file, the first line must indicate [the number of nodes].

* From the second line, each line indicates each hyperedge.
  * Each hyperedge <{a,b,c},{d}> is represented by a,b,c;d
  * We use semicolon(;) as the separator in each hyperedge
  * Comma(,) is used to distinguish source nodes

> Example data file

```
9
0;1
0;2
0;3
1;4
2,3;5
5,6;7
2,7;4
4,7;8
7;8
```

## Data 'HLMN' reference

We make a B-hypergraph model of "[HLMN](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2964116/bin/msb201056-s5.xml)" data from the paper <Jerby, Livnat, Tomer Shlomi, and Eytan Ruppin. "Computational reconstruction of tissue‐specific metabolic models: application to human liver metabolism." Molecular systems biology 6.1 (2010): 401.>.


## Usage
	python compute_bc.py [filepath] [bc option] [epsilons] [etas]    
    [bc option] = obc | pbc | lpbc | slpbc
    
> Examples: 
```
        python compute_bc.py data/hlmn/hlmn.txt obc
        python compute_bc.py data/hlmn/hlmn.txt pbc
        python compute_bc.py data/hlmn/hlmn.txt lpbc
        python compute_bc.py data/hlmn/hlmn.txt slpbc 0.01
        python compute_bc.py data/hlmn/hlmn.txt slpbc 0.01,0.005
        python compute_bc.py data/hlmn/hlmn.txt slpbc 0.01 0.1
```
