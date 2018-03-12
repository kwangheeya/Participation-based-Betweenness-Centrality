# BBC

Code for the papers
* "[Computing Betweenness Centrality in B-hypergraphs](http://dx.doi.org/10.1145/3132847.3133093)" published in the proceedings of the 2017 ACM on CIKM.
* To be

## Environment
* Windows 10, Python 3.6.2 64bit

## Requirement
* [Numpy](http://www.numpy.org/) package 

## Data format
> The data used in this code are stored in folder 'src/data'.

* In data file, the first line must indicate [the number of nodes].

* From the second line, each line indicates each hyperedge.
  * Each hyperedge <{a,b,c},{d}> is represented by a,b,c;d
  * We use semicolon(;) as the separator in each hyperedge

> Example data file
---
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
---
