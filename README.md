# sparseCROM
Sparsity-enabled CROM computes cluster-based reduced-order models (CROM) from compressed data and allows one to find few optimized sensor locations tailored to the specific model. Estimating a CROM from those compressed or few point measurements preserves the model structure and topology as compared to model estimated from the full data. The publication is available on [arXiv](https://arxiv.org/abs/1701.00038).

## Installation

1. Clone this repository to your desktop.
2. Add path to `sparseCROM/src` folder to Matlab search path using `addpath('<path to mds>/sparseCROM/src')`.

## Dependencies
For determining the optimized sensor locations tailored to a specific CROM, the following packages need to be installed.

1. [Sparse Sensor Placement Optimization (SSPOC)](https://github.com/bwbrunton/SSPOC_pub), which sets up the optimization problem. It is sufficient to add the file `SSPOC.m` to the source folder `sparseCROM/src`.

2. The optimization problem is solved using the [cvx](http://cvxr.com) toolbox, which needs to be installed.

## Getting Started

See `examples/example.m` for demonstrating the approach on the period double gyre flow, a simplified model of the gulf stream ocean front. Just execute this file in MatLab and it will generate the plot files in `examples/output`.

## License ([CiteMe OSS](https://github.com/cite-me/oss))

See the [LICENSE file](LICENSE) for details.
