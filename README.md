# Blackbox
Blackbox provides a way to generate data which can then be used for building/testing a surrogate model or for any other prupose.
Currently, only airfoil data generation with CST parameterization is ready to use, rest all modules are under developement.

## Documentation
Since the repository is not public, you will have to build the documentation locally only.

- Download the repository. If you are planning to install the package, then clone (not download) the repo in a folder 
  where you usually store packages (``packages`` folder is recommended) since you will be installing 
  in editable mode.
- Open the terminal and ``cd`` into the docs folder in the downloaded (or cloned) repository and run ``make clean html``
  This will create a build folder in the docs folder.
- Go to html folder in the build folder and open the index.html in the browser to view the documentation.
  More details about the project are provided there.

## License
Blackbox is licensed under Apache License, Version 2.0. See `LICENSE` file for more details. 

## Copyright
Copyright 2023 CODE Lab.
