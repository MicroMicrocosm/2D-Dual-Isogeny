# Two-Dimensional Representation Evaluation via Dual Isogenies

Paper title: **Efficient Two-Dimensional Representation Evaluation via Dual Isogenies with Applications to Isogeny-Based Cryptography**

## Repository Overview

This repository contains the implementations, experiments, and supporting code for the paper. The project is organized into three main subprojects, each focusing on a different isogeny-based cryptographic construction:

- [`poke/`](poke/) contains the POKÉ, INKE, and PIKE implementations, including the variants based on dual HD isogenies.
- [`sqisign-v3/`](sqisign-v3/) contains the SQIsign-v3 implementation with the dual HD-isogeny signing path.
- [`sqisign2dpush/`](sqisign2dpush/) contains the SQIsign2DPush implementation with cached dual-isogeny computations.

## Documentation

Please go to the README in each subfolder for the corresponding implementation's complete documentation, including its background, source-code organization, build instructions, tests, benchmarks, parameter sets, and implementation details:

1. [POKÉ with Dual HD isogenies](poke/README.md)
2. [SQIsign with Dual HD isogenies](sqisign-v3/README.md)
3. [SQIsign2DPush with cached dual isogenies](sqisign2dpush/README.md)

The detailed instructions are maintained in these subproject READMEs. Before building or running experiments, read the README for the subfolder you intend to use and follow its corresponding requirements and commands.

## General Requirements

The subprojects are primarily C/CMake projects. Depending on the selected implementation and target, you may need:

- CMake and a C99/C11-compatible compiler;
- GNU Make or another build tool supported by CMake;
- GMP development files for the POKÉ implementation and some SQIsign configurations;
- Clang for the POKÉ implementation;
- SageMath and the additional tools listed in the relevant subproject README when regenerating precomputed data.

Build directories and configuration options are local to each subproject. Do not assume that a build directory or CMake cache can be shared between different subfolders.

## License and Attribution

Each subproject contains its own license, notices, upstream attribution, and detailed project history. Refer to the corresponding subfolder and its README before redistributing or reusing any component.