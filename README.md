# ShibuChem

## Description

Repository for "CompRet: a comprehensive recommendation framework for chemical synthesis planning with algorithmic enumeration"

## Requirement

### Java

- ChemAxon (liscenced, academic free)

### Python

- rdkit 
- networkx

## Set up

1. Clone this repository

2. Download ChemAxon to your home directory
3. Put liscence certification file (liscence.cxl) under ~/.chemaxon
4. Install anaconda and create environment for extractor

```
$ cd ~
$ wget https://repo.anaconda.com/archive/Anaconda3-2019.07-Linux-x86_64.sh
$ sha256sum Anaconda3-2019.07-Linux-x86_64.sh
$ bash Anaconda3-2019.07-Linux-x86_64.sh
```

5. Create python environment

```
$ conda create -c rdkit -n shibuchem-env rdkit
$ conda activate shibuchem-env
$ conda install tqdm
$ conda install IPython
$ conda install pydot
$ pip install networkx
$ pip install bs4
```

6. Clone scscore repository

```
$ cd CompRet
$ git clone https://github.com/connorcoley/scscore.git
```

7. Prepare your own chemical reaction dataset under ```./reactions``` and building blocks as ```blocks.smi```. ChemAxon provides example set of basic chemical reaction in [their website](https://chemaxon.com/products/reactor/download).

reactions
    |-rxn1.smarts
    |-rxn2.smarts
    |-.
      .
      .


8. Edit run_all.sh and run!