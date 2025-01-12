# PGEcore

Welcome to PGEcore, a central repository for scripts that integrate and wrap commonly used bioinformatics tools and bespoke code for common functions. This repository is being developed collaboratively as part of a hackathon (WASABI25) organized by the [PlasmoGenEpi (PGE) group](https://www.plasmogenepi.org/). The scripts here will serve as foundational components for building robust and reusable bioinformatics workflows.

Contents: 
* [Purpose of the Repository](#purpose-of-the-repository)
* [Repository Structure](#repository-structure)
* [How to Contribute](#how-to-contribute)
* [Script Requirements](#script-requirements)

## Purpose of the Repository

This repository aims to:
* Provide a centralized collection of R scripts that wrap around external bioinformatics tools.
* Include bespoke utility functions to simplify and standardize workflow development.
* Serve as a shared resource for bioinformatics pipelines and modular workflows.

## Repository Structure 
```
PGEcore/
├── README.md       # Overview and guidelines
├── scripts/        # Scripts wrapping external tools and bespoke code
├── utils/          # Utility functions used across scripts
├── data/           # Example datasets
├── docs/           # Additional documentation or references
└── .gitignore      # Ignore unnecessary files
```

## How to Contribute 

1. Clone the Repository and checkout the `develop` branch:
```
  git clone https://github.com/PGEgroup/PGEtools.git
  cd PGEtools
  git checkout develop
```

2. Create a branch to develop on:
* Follow the branch naming convention (e.g. feature/your_tool_name`).
```
  git checkout -b <branch name>
```

3. Add your script and documentation:
* See [here](#script-requirements) for more information on this. 

4. Commit Your Changes:
* Write clear and concise commit messages.
```
  git add .
  git commit -m "Added wrapper for ToolX"
  git push origin <branch-name>
```

5. Submit a Pull Request (PR):

* Ensure you have included all of the features outlined in the [script guidlines](#script-guidlines).
* Once you have made your changes, create a PR into the `develop` branch. Someone will review your PR and provide feedback or approve it for merging. Never make any changes to the `main` branch, and please always PR into `develop`.

## How to Add a Script

For an example please see `scripts/FreqEstimationModel_wrapper`.

1. Create a directory under the `scripts/` folder. Make sure to give it an appropriate name. E.g. moire_wrapper.
    ```
    mkdir scripts/<my_dir>
    ```
2. Add the script for running the tool under the directory you just created. 
   Please follow the guidelines below when writting your code.
3. Copy the template from `docs/README.md` into this directory. Fill in the sections of the README for your tool. 
    ```
    cp docs/template_README.md scripts/<my_dir>/README.md
    ```

## Code Guidelines

Each software tool should be wrapped by one module (i.e., R script).

In case-by-case circumstances, multiple modules may be created for one software 
tool, but only if the functionality used and inputs have little overlap.

Within the module, code should be separated into functions according to whether 
it serves to format input, run the software tool, or format output. Multiple 
functions may be written for each purpose at the discretion of the module 
author. For example, if a module takes a common set of input but generates three 
different outputs, a logical set of functions might be read_input, run_tool, 
write_output1, write_output2, and write_output3.

Within the input function(s), contents of tabular inputs should be validated 
with the 
[validate](https://github.com/data-cleaning/validate?tab=readme-ov-file) (see 
[here](https://data-cleaning.github.io/validate/) for examples) package. 
Presence of expected columns, assumptions about column types, and assumptions 
about missingness should all be validated. Other input validation and general 
assertions should be performed with the 
[checkmate](https://mllg.github.io/checkmate/) package.

The main body of the script (i.e., outside of function definitions) should 
contain very little code aside from parsing the shell arguments and calling the 
appropriate sequence of functions depending on the outputs requested.

### Arguments

The script should expose the following as shell arguments:
- Paths for all required and optional input files
- Paths for all required and optional output files
- In the case of software tools that are being used to generate multiple 
  outputs, the module should have flags that control what outputs will be 
  generated
- The number of threads, in the case of tools that support multi-threading
- The random number seed, in the case of tools that use random number 
  generation
- Any other tool parameters the module author deems important to adjust when 
  running the tool on different datasets. If a parameter seems useful to tune 
  during benchmarking and stress testing of the tool, it should probably be 
  exposed as an optional argument.

Arguments should be parsed with the 
[optparse](https://rdrr.io/cran/optparse/f/README.md) package. Specify types and 
include help messages.

### Dependencies

It is encouraged to use the tidyverse packages for all tabular data 
manipulation. Other than that and the packages described above for argument 
parsing and input validation, endeavor to minimize dependencies.

### Documentation

We will use [roxygen2](https://roxygen2.r-lib.org/) to create online 
documentation for each module. Please write docstrings using the roxygen2 format 
for each function in your script.

Otherwise, attempt to make your code as readable as possible and use in-line 
comments to clarify complicated segments.
