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
  git clone https://github.com/PlasmoGenEpi/PGEcore.git
  cd PGEcore
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

## Script Requirements

For an example please see ...

1. Create a directory under the `scripts/` folder. Make sure to give it an appropriate name. E.g. moire_wrapper.
    ```
    mkdir scripts/<my_dir>
    ```
2. Add the script for running the tool under the directory you just created. Make sure to add options for available outputs and any parameters that could be relevant. 
3. Copy the template from `docs/README.md` into this directory. Fill in the sections of the README for your tool. 
    ```
    cp docs/template_README.md scripts/<my_dir>/README.md
    ```
