# Data Valuation from Data-Driven Optimization

Code and data to replicate the results of the following paper:

[Data Valuation from Data-Driven Optimization](https://arxiv.org/pdf/2305.01775.pdf)

---

**Abstract** 
With the ongoing investment in data collection and communication technology in power systems, data-driven optimization has been established as a powerful tool for system operators to handle stochastic system states caused by weatherand behavior-dependent resources. However, most methods are ignorant to data quality, which may differ based on measurement and underlying privacy-protection mechanisms. This paper addresses this shortcoming by (i) proposing a practical data quality metric based on Wasserstein distance, (ii) leveraging a novel modification of distributionally robust optimization using information from multiple data sets with heterogeneous quality to valuate data, (iii) applying the proposed optimization framework to an optimal power flow problem, and (iv) showing a direct method to valuate data from the optimal solution. We conduct numerical experiments to analyze and illustrate the proposed model and publish the implementation open-source.

---

## Usage

Everything is implemented in [Julia](https://julialang.org/) (Version 1.6.3) and all results can be reproduced by running the code blocks in `experiments_paper.ipynb` using the interactive [Jupyter environment with and IJulia backend](https://github.com/JuliaLang/IJulia.jl). Required packages and their versions used for this implementation are stored in the `Project.toml` file. The project environment is also activated in the notebook. Don't forget to also run `Pkg.instantiate()` during the first run to install and precompile all required packages.

This repository hosts a copy of the [PSDataIO](https://github.com/mieth-robert/PSDataIO) package, which provides custom code to read and handle power system data. This copy, i.e., the folder `PSDataIO`, must remain in the root folder of the project to work properly. 

---

## Citation
```
@article{mieth2023data,
  title={Data Valuation from Data-Driven Optimization},
  author={Mieth, Robert and Morales, Juan M and Poor, H Vincent},
  journal={arXiv preprint arXiv:2305.01775},
  year={2023}
}
```