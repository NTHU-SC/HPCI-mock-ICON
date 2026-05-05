# ICON Fortran-support library Contribution Guidelines

## Introduction

The ICON Fortran-support library is a collection of supported functionalities commonly shared in code development. Such as exception, initialization, and namelist handling. The library also contains basic functionalities implemented in C with Fortran interfaces. See `README.md` for a list of modules in the library.

## Communication

Use [issues](https://docs.gitlab.com/user/project/issues/) of the repository (fork) as a tool to communicate, track, and obtain information related to work items. Create an issue to report a bug you found, request a feature that needs to be implemented, or to discuss and coordinate with others on a particular topic. Issues can also be employed to address other work items, such as tracking the porting of a module to a different programming language or collecting information regarding an unexpected simulation result.

> **Note:** Prefer GitLab communication over private emails or messages to ensure information is searchable and accessible to a larger audience.

> **Note:** Avoid using [tasks](https://docs.gitlab.com/user/tasks/), as they add an extra level of hierarchy. To keep things simple, use issues only.

### Creating Issues

1. Present issues clearly and include all relevant information so others can easily understand them.

2. Each issue should focus on **one actionable task** (e.g., one bug or feature request). If it becomes too lengthy or complex, break it into smaller issues and link them to the original one.

3. Always assign an issue to someone and communicate this clearly. Mention developers as needed, but keep it to a minimum.

4. Apply as many relevant [labels](https://docs.gitlab.com/user/project/labels) as possible. They help in classifying the issue (e.g. `bug`, `feature request`, `discussion`). Feel free to create new labels if the right ones don't exist.

## Code Contributions

### Coding style

We use [`clang-format`](https://clang.llvm.org/docs/ClangFormat.html) and [`fprettify`](https://github.com/fortran-lang/fprettify) for C and Fortran code styling correspondingly. Install both packages if you wish to contribute to the code.

We provide predefined CMake target `make format` for code formatting using `clang-format` and `fprettify`. Be sure to run `make format` before you commit to the repository.

We also additionally use [`cmake-format`](https://github.com/cheshirekow/cmake_format) for CMake styling. Please also format the CMake files if you made changes in the CMake script.

### General Coding Rules

1. Avoid adding comments about future actions. For example,
    ```fortran
    ! Delete the subroutines below once the module is validated
    ```
    If necessary, include a reference to an issue that provides the progress status (e.g., an issue on the validation of the aforementioned module).
2. Do not add commented-out code to the codebase, as it produces maintenance and development overhead.

### Merge Requests

1. Always **Choose a template** (`bugfix`/`feature`) before writing a merge request description.

2. Follow that instructions in the merge request template and make sure you complete the `Mandatory steps before review` list before requesting a review.

3. Make the merge request title concise (titles become the first line of the commit message when the merge requests are accepted.

4. Please, adhere to the following recommendations for the merge requests descriptions, which will become part of the commit message when the merge request is accepted:

    - use simple English in the active form (e.g. this implements A, updates B);

    - avoid special Markdown symbols and prefer plain ASCII, the message should read well in the terminal;

    - keep it short (excluding details, descriptions are appended to the merge request commit message).

5. The lists of co-authors in merge requests are generated automatically based on the authorship of the commits in the source branches. Please ensure that the commits in the source branch have the correct authorship with the correct email addresses (they can be [automatically-generated private commit emails](https://docs.gitlab.com/user/profile/#use-an-automatically-generated-private-commit-email)). If some commits have the wrong authorship, you can provide the list of co-authors using the following format:
    ```
    Co-authored-by: First-Name Second-Name <email.address@example.de>
    Co-authored-by: Another Name <another.address@example.com>
    ```

6. Check whether you are listed in the `AUTHORS.txt` list. If not, add yourself in the alphabetic order. Your name will be included when your merge request gets merged.

### Versions (Commits) Update in ICON

1. Always use commits from the `master` branch for updating new versions in the ICON repositories. However, you can use commits from other branches for your local development.

2. The `libfortran-support` versions will be irregularly updated directly to the `icon/icon` repository.

3. Each commit that merges back into the ICON repositories must contain a release tag. Please contact the maintainers if you cannot wait for the update from `icon/icon`. The maintainers will issue a release tag.
