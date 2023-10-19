# Welcome to Neglected Diagnostics contributing guide

Thank you for investing your time in contributing to our project!

Read our [Code of Conduct](https://github.com/uw-ssec/code-of-conduct/blob/main/CODE_OF_CONDUCT.md) to keep our community approachable and respectable.

## Developing Neglected Diagnostics

Neglected Diagnostics is in its initial phase of development. Follow the below steps for setting up the development environment and contributing to the project.

### Environment Setup

Fork the repository, and create your local version, then follow the installation steps:

1. Create conda environment

```bash
conda create -y -n neglected-diagnostics python=3.10 pip
```

2. Activate conda environment

```bash
conda activate neglected-diagnostics
```

3. Install dependencies

```bash
pip install -e ".[dev,docs]"
```

### Post setup

You should prepare pre-commit, which will help you run the linting and formatting checks configured in the project before every commit.

```bash
pre-commit install # Will install a pre-commit hook into the git repo
```

NOTE: You can also/alternatively run `pre-commit run` (changes only) or `pre-commit run --all-files` to check even without installing the hook.

## Run Neglected Diagnostics

You are now ready to run the Neglected Diagnostics Streamlit application:

```bash
streamlit run src/app/Home.py
```


## Testing

When making changes, please run the following tests:

Unit tests:

```bash
pytest ./tests
```

## Pull Requests

Please follow the below guidelines when you want to raise a Pull Request:

- It may be helpful to review [this tutorial](https://www.dataschool.io/how-to-contribute-on-github/) on how to contribute to open source projects. A typical task workflow is:

    - [Fork](https://docs.github.com/en/get-started/quickstart/fork-a-repo) the code repository specified in the task and [clone](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) it locally.
    - Review the repo's README.md and CONTRIBUTING.md files to understand what is required to run and modify this code. 
    - Create a branch in your local repo to implement the task.
    - Commit your changes to the branch and push it to the remote repo.
    - Create a pull request, adding the task owner as the reviewer.

- Please follow the [Conventional Commits](https://github.com/uw-ssec/rse-guidelines/blob/main/conventional-commits.md) naming for pull request titles.
