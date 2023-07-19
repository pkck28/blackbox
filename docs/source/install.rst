Installation
============

Follow the below steps for installation:

- You probably have the repository, so can skip this step. If not, then clone the 
  `repository <https://github.com/ComputationalDesignLab/blackbox>`_.
- Open the terminal and ``cd`` into the downloaded repository and run::

    pip install -e .

  This will install the package in editable mode i.e. if you change anything in the 
  source code of the package, then you don't need to re-install the package.

- Since the packge is under constant developement, installing new updates will be very easy.
  Open the terminal and ``cd`` into the cloned repository and run::

    git pull

  It will ask for your github username and password. This will get the latest updates from the  main 
  `repository <https://github.com/ComputationalDesignLab/blackbox>`_, if there are any. If you get an 
  error saying that git is not installed, then run this command in terminal::

    sudo apt install git

  and re-run the ``git pull`` command in the root folder of the cloned repository.
