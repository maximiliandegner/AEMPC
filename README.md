# Adaptive Economic Model Predictive Control: Performance Guarantees for Nonlinear Systems

## Description

This repository contains the MATLAB code that accompanies the paper 
> Maximilian Degner, Raffaele Soloperto, Melanie N. Zeilinger, John Lygeros, Johannes Köhler, "Adaptive Economic Model Predictive Control: Performance Guarantees for Nonlinear Systems", 2024, arXiv: 2412.13046.


## Prerequisites

- MATLAB (tested with version R2023b) and MPT3 Toolbox
- CasADi (used for online computations, tested with version 3.6.6 and the corresponding IPOPT version)
- YALMIP (used for offline computations)
- MOSEK (used for offline computations, tested with version 10.1.21)

## Usage

Run `main.m` to execute the simulation for comparing our proposed adaptive economic MPC scheme to an economic MPC scheme without adaptation and to an economic MPC scheme with perfect parameter knowledge.

The offline computations include the computation of the tube and the auxiliary feedback and the computation of the terminal cost according to Proposition 2 in our paper (cf. Appendix for more implemtation details).
For running solely the offline computations, you can run `offline_pipeline.m`.

The online simulation itself runs in the file `AE_MPC_journal.m`.

Plotting is performed by the script `figures_journal`.

## License

This project is licensed under the MIT License.

## Citation

If you use this code in your research, please cite our work:

```text
@article{degner2024AEMPC,
  title={Adaptive Economic Model Predictive Control: Performance Guarantees for Nonlinear Systems},
  author={Degner, Maximilian and Soloperto, Raffaele and Zeilinger, Melanie N., Lygeros, John and Köhler, Johannes},
  year={2024}
}
```
  
## Support and Contact

For any questions or issues related to this code, please contact the author:

- Maximilian Degner mdegner(at)ieee(dot)org
