===============
WHAT TO INSTALL
===============

This is the list of the different packages/software you will have to install in order to 
use the tool.
If you run into any problems, please contact Leandro de Almeida: monolipo.physics@gmail.com.
Any use or modification of this tool can be done, as long as the author is warned
These instructions were adapted from Morisset, C 2017.

git
===

You will need to have installed git software. You can verify if you
already have it by doing in a terminal: ::
   which git

if no link is given, you have to download. Depending on your operating system:

I. Mac OSX: install the package from https://git-scm.com/download/mac
II. Linux: depending on your distribution:

    A. ``sudo yum install git``
    B. ``sudo apt-get install git``

Once done, from a new terminal create a directory dedicated for
this tool with the name LIGHTCURVETOOL, and from it, clone this github files.
::   
   git clone https://github.com/monolipo/Light-Curve-Python-Analisys-Tool.git

This will download the code and instructions.

python 3
======

The most simple way to have all the needed packages is to install a full python distribution using for example Anaconda. If you already have Anaconda installed, you still may need to install the 3.6 version and some packages.

To verify if you have anaconda installed, just look at the answer of ``which python`` command in a terminal. If it points to a directory which name contains "anaconda", it means that you already have anaconda installed .

If you don't have anaconda installed: The anaconda package is available following this link: `https://www.continuum.io/downloads <https://www.continuum.io/downloads>`_. We will use the python 3.6 package.
Once done, from a new terminal, do the following: ::

  conda install pymysql

2. If you already have an anaconda distribution installed for python2:
you only need to install a 3.6 environment, by writing the following in a terminal: ::
   
   conda create --name py3k6 python=3.6 matplotlib scipy numpy ipython h5py astropy pymysql pandas pytest ipykernel

   source activate py3k6

   ipython kernel install --user

In both cases, once you have python 3.6 installed, you still need some
additional libraries. From a terminal, do the following: ::
   pip install atpy
   pip install pillow
   
The last librarie is the kplr tool that you can instal doing the following: ::
   git clone https://github.com/dfm/kplr.git
   cd kplr
   python setup.py install

Test your installation
======================

Once all the above is done, you can open a terminal and go to the directory where we have downloaded the code using git (at the beginning of this page). 

Go to codigo subdirectory and enter: ::
  
   jupyter notebook

This should open a new tab in your web browser. 

Click on ``test_install.ipynb``

A new page appears. You can execute each of the instructions from this page by clicking on the "PLAY" button in the upper part of the page. Or press SHIFT-ENTER.

If everything is OK, it will appear "have fun".
