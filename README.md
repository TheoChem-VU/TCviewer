# template
This is a template repository for TheoCheM programs. Please use the quick-start guide below to set up your project quickly.

This template features automated documentation building using sphinx-apidoc and publishing to github-pages. The publishing will only be enabled when the repo is made public. Don't forget to set the github-pages source branch to `gh-pages` branch. The template also has automatic testing after every push in any branch. Tests must be written in the test folder at root level. The template also does automatic publishing to PyPI everytime a release is made. For this, create a trusted publisher on the PyPI website pointing to this repo. Don't forget to make a 

# Quick-Start Guide
1. Edit `pyproject.toml` to reflect the new repository, including name of the project, authors and description.
2. Update the `requirements.txt` file to include any packages you need to run your new code. You will need this when building and packaging the project before uploading it to PyPI.
3. The new code should be added in the `src` directory (e.g. `src/pyorb/module1/...`). Anything in that directory will be published on PyPI, so don't include documentation or testing scripts in there.
4. To publish the documentation to github-pages you will have to make the repo public. Don't forget to enable github pages in the repo settings and setting it to publish from the `gh-pages` branch.
5. To publish new releases to PyPI you will have to make the repo public and create a new environment in the settings. Then create a trusted publisher on the PyPI website pointing to the correct repo and environment. If you then make a new release it should automatically be uploaded to the PyPI database and the project can be installed using `pip install ...`
