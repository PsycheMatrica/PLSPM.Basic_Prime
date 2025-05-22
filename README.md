# PLSPM.Basic_Prime

### Author:
Gyeongcheol Cho and Heungsun Hwang

## Description:
- The **PLSPM.Basic_Prime** package enables users to estimate and evaluate basic PLSPM models.

## Features:
- Estimate PLSPM model parameters and calculate their standard errors (SE) along with 95% confidence intervals (CI).
- Enable parallel computing for bootstrap sampling.
- Allow users to specify a sign-fixing indicator for each component.
- Provide an option for Dijkstra's correction.
- Handle missing values in the data.

## Installation:
To use this package in MATLAB:
1. Clone or download the repository:
   ```bash
   git clone https://github.com/PsycheMatrica/PLSPM.Basic_Prime.git
   ```
2. Add the package to your MATLAB path:
   ```matlab
    addpath(genpath('PLSPM.Basic_Prime'))
   ```

## Usage:
- For examples on how to use the package, refer to the `Run_Example_BasicPLSPM.m` file. This file demonstrates the implementation of `BasicPLSPM()` using the ACSI dataset.

## Compatibility:
- Tested on MATLAB R2023b.
- Likely compatible with earlier MATLAB versions.

### Citation (APA):
- If you use **PLSPM.Basic_Prime** in your research or publications, please cite it in APA format as follows:

```plaintext
Cho, G & Hwang, H (2024). PLSPM.Basic_Prime: A package for basic partial least squares path modeling [Computer software]. GitHub. https://github.com/PsycheMatrica/PLSPM.Basic_Prime
```
