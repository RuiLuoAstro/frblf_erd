# frblf_erd
A Bayesian framework to measure the event rate density of FRB luminosity function

## References

If you would like to use this code to study FRB luminosity function, please cite the paper [Luo et al. 2020, MNRAS, 494, 665](https://ui.adsabs.harvard.edu/abs/2020MNRAS.tmp..667L/abstract) as well as [Luo et al. 2018, MNRAS, 481, 2320](https://ui.adsabs.harvard.edu/abs/2018MNRAS.481.2320L/abstract).

## Dependencies

Python (2.7.x), Numpy (1.14 at least), Scipy (1.0.0 at least), PyMultiNest (see https://github.com/JohannesBuchner/PyMultiNest for more details), Matplotlib

## Verify the mock data using PyMultiNest

``` ./run_simu.sh ``` &emsp;&emsp;&emsp;&emsp;
**Notes: better implement it in cluster where MPI was installed well. The posterior outputs are saved on ./nest_out/simu/**

``` ./draw_sim.sh ``` &emsp;&emsp;&emsp;&emsp;
**Plot the posterior distribution contours of the mock data, which are made on ./plots/simu/**

## Measure the quantified FRB LF with sample

``` ./run_samp.sh ``` &emsp;&emsp;&emsp;&emsp;
**Notes: better implement it in the cluster where MPI was installed well. The posterior outputs are saved on ./nest_out/samp/**

``` ./draw_samp.sh ``` &emsp;&emsp;&emsp;&emsp;
**Plot the posterior distribution contours of the real FRB sample, which are saved on ./plots/samp/**
