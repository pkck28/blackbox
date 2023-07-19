# Blackbox
Blackbox provides a way to generate aerodynamic/aerostructural data which can then be used for building/testing a surrogate model or for any other purpose. Currently, airfoil data generation using CST and FFD parameterization, and wing data generation using FFD parameterization are ready to use, aerostructural module is under developement.

## Documentation
**Since the repository is not public, you will have to build the documentation locally**.

- Clone (not download) the repo in a folder where you usually store packages (``packages`` folder is recommended)
  since you will be installing in editable mode. Easiest way to clone is to use terminal.
  - Open the terminal in folder where you want to clone the repository and run:
  ```
    git clone https://github.com/ComputationalDesignLab/blackbox
  ```
  It will ask for your username and password of github account.
  - If you get an erorr saying that git is not installed, then run:
  ```
    sudo apt install git
  ```
  and try the clone command once again.
- Open the terminal and ``cd`` into the docs folder in the cloned repository and run ``make clean html``.
  This will create a build folder in the docs folder.
- Open html folder in file explorer and open the index.html in the browser to view the documentation. More details about the project are provided there.

## License
Blackbox is licensed under Apache License, Version 2.0. See `LICENSE` file for more details. 

## Copyright
Copyright 2023 CODE Lab.
