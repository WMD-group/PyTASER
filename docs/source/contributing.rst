Contributing to PyTASER
=======================

We love your input! We want to make contributing to as easy and
transparent as possible, whether it's:

- Reporting a bug
- Discussing the current state of the code
- Submitting a fix
- Proposing or implementing new features
- Becoming a maintainer

Reporting bugs, getting help, and discussion
--------------------------------------------

Please submit questions and bug reports to the `GitHub issues page <https://github.com/WMD-group/PyTASER/issues>`__.

If you are making a bug report, incorporate as many elements of the
following as possible to ensure a timely response and avoid the
need for followups:

- A quick summary and/or background.
- Steps to reproduce - be specific! **Provide sample code.**
- What you expected would happen, compared to what actually happens.
- The full stack trace of any errors you encounter.
- Notes (possibly including why you think this might be happening,
  or steps you tried that didn't work).

We love thorough bug reports as this means the development team can
make quick and meaningful fixes. When we confirm your bug report,
we'll move it to the GitHub issues where its progress can be
further tracked.

Contributing code modifications or additions through Github
-----------------------------------------------------------

We use github to host code, to track issues and feature requests,
as well as accept pull requests. We maintain a list of all
contributors `here
<https://pytaser.readthedocs.io/en/latest/contributors.html>`__.

Pull requests are the best way to propose changes to the codebase.
Follow the `Github flow
<https://www.atlassian.com/git/tutorials/comparing-workflows/forking-workflow>`__
for more information on this procedure.

The basic procedure for making a PR is:

- Fork the repo and create your branch from main.
- Commit your improvements to your branch and push to your Github fork (repo).
- When you're finished, go to your fork and make a Pull Request. It will
  automatically update if you need to make further changes.

How to Make a **Great** Pull Request
------------------------------------

We have a few tips for writing good PRs that are accepted into the main repo:

- Use the Black Code style for all of your code. Find an example `here
  <https://black.readthedocs.io/en/stable/the_black_code_style/current_style.html>`__.
- Your code should have (4) spaces instead of tabs.
- If needed, update the documentation.
- **Write tests** for new features! Good tests are 100%, absolutely necessary
  for good code. We use the python ``pytest`` framework -- see some of the
  other tests in this repo for examples, or review the `Hitchhiker's guide
  to python <https://docs.python-guide.org/writing/tests/>`__ for some good
  resources on writing good tests.
- Understand your contributions will fall under the same license as this repo.

When you submit your PR, our CI service will automatically run tests.
We welcome good discussion on the best ways to write your code, and the comments
on your PR are an excellent area for discussion.
