# Coding style
We use [`pre-commit`](https://pre-commit.com) hooks to maintain a set of
formatting and linting rules. Although there is a CI job that runs for each
merge request and checks whether the contribution does not break the rules, we
recommend registering the hooks in your local repository clone. This way, each
commit undergoes the formatting and linking checks automatically.

We recommend installing `pre-commit` to a separate Python virtual environment
using `pip`. For example, the following commands install the tool to the user's
home directory:
```bash
python3 -m venv ~/pre-commit
~/pre-commit/bin/python3 -m pip install --upgrade pip
~/pre-commit/bin/python3 -m pip install pre-commit
```

You can now switch to the root of the repository and run the following command
to register the hooks specified in
[`.pre-commit-config.yaml`](/.pre-commit-config.yaml):
```bash
cd libmtime
~/pre-commit/bin/pre-commit install
```

From now on, each commit you make will be checked by a set of formatters and
linters. Normally, the formatting tools are configured to modify the files in
place. This means that if they fail, all you need to do is to accept the
suggested changes and commit them:
```bash
git add .
git commit
```

Note that you will need to register the hooks for each fresh clone of the
repository. Alternatively, you can follow
[these instructions](https://pre-commit.com/#automatically-enabling-pre-commit-on-repositories)
to configure `git` to register hooks automatically for each new clone of a
repository that declares them.
