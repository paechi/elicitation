# Monte Carlo Elicitation
This repository contains the code for Monte Carlo Elicitation experiments.
We have borrowed the Polytope.py and Utils.py from https://github.com/federicotoffano/SMMR and minimax 
implementations with linear programming. 
## Install the dependencies

- Create a new virtual/conda environment (I use Python 3.9.6) and install necessary libraries using pip.
```
python3 -m venv venv
source venv/bin/activate
pip3 install -r requirements.txt
```

## Run Monte Carlo elicitation.

```
python3 elicitation.py
```

