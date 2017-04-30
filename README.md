PLANET: a Piece-wise LineAr feed-forward NEural network verification Tool
===========================================================================

PLANET is a tool to verify feed-forward neural networks against convex specifications. The node types supported in the Network are ReLU nodes, MaxPool nodes, linear nodes (also called InnerProduct in Caffe), and variants of it that are equivalent from the verification point of view (such as convolutional layer nodes).

This documentation is work in progress, but the tool is already usable. If you have questions, please leave a node under https://github.com/progirep/planet/issues. It shall be noted that verifying neural networks is computationally difficult, hence don't expect verification to scale as well as learning. However, Planet has been successfully used for some case studies.

Part of this repository are scripts to translate models that have been learned with the Caffe framework into input files for the verifier. Note that these scripts currently support only a limited number of layer types and only certain architecture styles -- please look at the example network .prototxt files for supported features.


License
-------

Planet itself is available under GPLv3 license, as it uses the GNU Linear Programming toolkit, which is also under GPLv2.

This repository contains some modified network definition files from the [Caffe](http://caffe.berkeleyvision.org) neural network learning framework, which is licensed under the BSD 2-Clause license. The modifications are available under the same license.

The Planet tool uses a modified version of Minisat v.2.2.0, which itself is available under the MIT license. If you are interested in Minisat, visit http://minisat.se/MiniSat.html and get the solver from there. This repository only contains the Minisat parts that were actually used in Planet.


Preparation:
------------

- Not part of this repository are two scripts to translate a neural network description in the "prototxt" format (as produced by Caffe) into a json file, and a (modified) script to generate a database in "HDF5" format from comma-separated-value files. To obtain those, run the following command from the terminal in the PLANET directory:

    cd tools; wget https://gist.github.com/progirep/fd7d2dc120862faa984a70f503611013/raw/260e1e76cebd0ea58bf1a03b64c3f1e0002fc677/csv_to_hdf5_supervised_classification.py; wget https://raw.githubusercontent.com/vadimkantorov/caffemodel2json/3a8fd443bf1596dad5f517aecdef08a81bf73bfe/caffemodel2json.py
    
- Then, some packages have to be installed in order to build "PLANET" (and to use the python scripts that come with it). For Ubuntu 17.04, the following packages are known to be requirements:
    - libglpk-dev
    - qt5-qmake
    - valgrind

Some scripts in this repository need the root path of the Caffe neural network learning framework, or the root path of the Yices SMT solver. The following environment variables should point to these:

- The $CAFFE environment variable should be set to the Caffe root path.
- The $YICES environment variable should be set to the Yices root path.


Running PLANET on a single specification file
---------------------------------------------

It suffices to run "path/to/planet inputFile.rlv" to start verifying an input file. The result will either be "UNSAT" or "SAT", depending on whether a neural network input has been found for which the network behavior satisfies all given constraints. In the latter case, all node values are printed out as well, and a cumulative error value for the differences between the network values obtained by standard propagation and the linear programming solver is given as well.


Running the Collision Avoidance Benchmarks:
-------------------------------------------
TODO.


Running the MNIST Examples:
-------------------------------------------
- First, run "./testAccuracy.sh" in the casestudies/MNIST directory. The process can be aborted once the first digit results start coming in.
- Then, run "../../tools/rlvCompressor.py testNetwork.rlv > testNetworkB.rlv" to get a slightly optimized version of the benchmark.
- Then, run "./prodNetwork.py testNetworkB.rlv GIVESTRONG 2 <distance>" to obtain an image that is strongly classified as a 2.

