Edit the setup.py file, as a minimum, update the version tag.
Also check the __init__.py version information and requirements are consistent with development.

The proper source distribution build command is: "python setup.py sdist --formats tar"
This will produce a tar file in the dist folder.

To Deploy
upack the tar file under /dist  (Note: Winzip is getting "pissy" about the gz format.)

cd to dist/alfano-0.xa

Use command "setup install"

Package "alfano" should be deployed to python environment Lib/site-packages/alfano.