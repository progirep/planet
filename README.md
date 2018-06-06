PLANET: a Piece-wise LineAr feed-forward NEural network verification Tool
===========================================================================

PLANET is a tool to verify feed-forward neural networks against convex specifications. The tool implements the approach described in [this paper](https://arxiv.org/abs/1705.01320).

This documentation is work in progress, but the tool is already usable. If you have questions, please leave a note under https://github.com/progirep/planet/issues. It shall be noted that verifying neural networks is computationally difficult, hence don't expect verification to scale as well as learning. However, Planet has been successfully used for some case studies. The node types supported in the Network are ReLU nodes, MaxPool nodes, linear nodes (also called InnerProduct node in Caffe), and variants of it that are equivalent from the verification point of view (such as convolutional layer nodes).

Part of this repository are scripts to translate models that have been learned with the [Caffe](http://caffe.berkeleyvision.org) framework into input files for the verifier. Note that these scripts currently support only a limited set of layer types and only certain architecture styles -- please look at the example network .prototxt files for supported features.


License
-------

Planet itself is available under GPLv3 license, as it uses the GNU Linear Programming toolkit, which is also under GPLv2.

This repository contains some modified network definition files from the [Caffe](http://caffe.berkeleyvision.org) neural network learning framework, which is licensed under the BSD 2-Clause license. The modifications are available under the same license.

The Planet tool uses a modified version of [Minisat v.2.2.0](http://minisat.se/MiniSat.html), which itself is available under the MIT license. If you are interested in Minisat, visit http://minisat.se/MiniSat.html and get the solver from there. This repository only contains the Minisat parts that were actually used in Planet.


Preparation:
------------

The verifier is written in C++ (using some C++14 features), and the scripts that accompany the verification tool are written in Python 2. We assume a Linux system for the instructions in the documentation. You may be able to get the verification tool and the scripts running on Windows with some modification, but this has not been tested.

Before we can use the tool, some preparation has to be done:

- Not part of this repository are two scripts to translate a neural network description in the "prototxt" format (as produced by Caffe) into a json file, and a (modified) script to generate a database in "HDF5" format from comma-separated-value files. To obtain those, run the following command from the terminal in the PLANET directory:

    cd tools; wget https://gist.github.com/progirep/fd7d2dc120862faa984a70f503611013/raw/260e1e76cebd0ea58bf1a03b64c3f1e0002fc677/csv_to_hdf5_supervised_classification.py; wget https://raw.githubusercontent.com/vadimkantorov/caffemodel2json/3a8fd443bf1596dad5f517aecdef08a81bf73bfe/caffemodel2json.py
    
- Then, some packages have to be installed in order to build "PLANET" (and to use the python scripts that come with it). For Ubuntu 17.04, the following packages are known to be requirements:
    - libglpk-dev
    - qt5-qmake
    - valgrind
    - libltdl-dev
    - protobuf-compiler
    
The tool can then be compiled by running:

    cd src; qmake Tool.pro; makes

Some scripts in this repository need the root path of the Caffe neural network learning framework, or the root path of the Yices SMT solver. The following environment variables should point to these:

- The $CAFFE environment variable should be set to the Caffe root path.
- The $YICES environment variable should be set to the Yices root path.


PLANET specification files
--------------------------

Planet uses a special neural network verification scenario description format as input file. It consists of a neural network description and a specification. The file format is text-based and independent of any neural network framework. While we provide scripts for Caffe to translate networks to input files, it should be possible to write similar scripts for other frameworks, such as TensorFlow.

The input files are also pretty low-level -- nodes in the network are described individually, with all their weights. Every node description and specification part is given on an individual line in the input file. The following lines constitute a simple example specification file that features all node types supported by PLANET.

    Input inX0
    Input inX1
    Input inX2
    ReLU relu1 1.0 -0.5 inX0 -0.5 inX1
    MaxPool maxi inX2 relu1
    Linear lina 0.2 0.3 relu1 0.4 maxi
    Assert <= 0.0 1.0 inX0
    Assert >= 1.0 1.0 inX0
    Assert <= 0.0 1.0 inX1
    Assert >= 1.0 1.0 inX1
    Assert <= 0.0 1.0 inX2
    Assert >= 1.0 1.0 inX2
    Assert >= 0.6 1.0 lina 0.1 inX0
    Assert <= 0.4 1.0 lina 0.1 inX0

The first three lines define Input nodes for the network. The line afterwards defines a ReLU node, with bias 1.0 and two incoming edges from inX0 and inX, each weighted -0.5. The next node is a MaxPool node, called maxi, with two incoming edges from inX2 and relu1. The node afterwards is a linear node called "lina", with a bias of 0.2 and two edges from relu1 and maxi, weighted 0.3 and 0.4.

The lines afterwards are the specification constraints. Every line is of the form "Assert (operation) (constant) (linear parts)", and requires the verification tool to only consider network node valuations for which the verification condition "(constant) (operation) (linear parts)" holds. So the first "Assert" statement requires 0 to be less than or equal to 1.0 times inX0, and the second one requires 1.0 to be greater than or equal to 1.0 times inX0. Together, the constrain inX0 to be in the interval from 0 to 1.

The later constraint lines follow a similar form and restrict the values of the other input nodes. Note that Planet requires all input node values to be bounded - the tool will otherwise abort with an error message. The final two constraints require 1.0 times lina + 0.1 times inX0 to be between 0.4 and 0.6.

It suffices to run "path/to/planet inputFile.rlv" to start verifying an input file. The result will either be "UNSAT" or "SAT", depending on whether a neural network input has been found for which the network behavior satisfies all given constraints. Running Planet on this input file (by just passing its file name as parameter to the tool) yields a lot of debug output and concludes with the following result:

    SAT
    
    Valuation:
    - inX0: 1 / 1
    - inX1: 0.714286 / 0.714286
    - inX2: 0.142857 / 0.142857
    - relu1: 0.142857 / 0.142857
    - maxi: 0.142857 / 0.142857
    - lina: 0.3 / 0.3
    Literals: -4 2 3
    Total error: 1.11022e-16
    
The first line states that there is a node value assignment function for the network that satisfies all "Assert" constraints. The lines below state how it looks like, including the input and output nodes of the network. The "Literals" line tells us the assignment to the SAT solving variables - the SAT solver is used internally in the search process. The final line "Total Error" gives the sum of deviations between the computed node values and the ones returned by the linear programming library employed in Planet. These values should be very low, otherwise the network weights were so high (or low) that the numerical accurracy of the solution process was to low to obtain good results.

In case there is no node value assignment function for the network that satisfies all constraints, the final output of the tool will be "UNSAT".


Translating a Caffe model to an input file for the verifier: 
------------------------------------------------------------
The RLV file format used by Planet is pretty low-level, as it lists the artificial neurons one-by-one. The Caffe deep learning framework on the other hand operates on a layer-by-layer basis and also stores learned network weights in this way. Thus, Caffe network models need to be preprocessed in order to be verified. The planet distribution come with a few Python scripts to make this work.

For this to work, a few requirements need to be met:

1. Currently, the only Caffe layer types supported are "Split", "InnerProduct", "ReLU", "Convolution", "Pooling" (using Max-Pooling), "HDF5Data", "Data", "Reshape", and "Pooling". There may be other layer types, but they must not be backwards-reachable from the "Accuracy" layer (of which there must be exactly one in the training phase).
2. Caffe allows InnerProduct and ReLU layers to be connected to each other (so that when plotting the network architecture, we get a small loop). This is not allowed - rather, the output of an InnerProduct layer can be sent through a ReLU layer.
3. The translation script needs to know the dimensions of all layers, which it cannot always infer from the description of the output network. Whenever the script yields an error due to this problem, a Reshape layer should be added immediately after the data layer to rehape it to the same shape as before. This causes shape information to be written to the file containing the learned network, which is then interpreted by the translator scripts.
4. InnerProduct and ReLU layers currently only operate on one-dimensional layers. If they should be used on n-dimensional layers for n>1 (such as in image recognition applications), they need to be reshaped before the InnerProduct or ReLU layers. They can be reshaped back in the network definition.
5. MaxPool layers assume a three-dimensional input (x-dimension of images, y-dimension of images, and color/feature dimension).

There are two network description files as examples in the repository, which let Caffe output network models that the scripts that come with Planet can process:

- casestudies/vehicleCollision/caffe\_net\_with\_pooling.prototxt 
- casestudies/MNIST/lenet\_train\_test.prototxt 

Training the network can be done as usual with Caffe - by passing a separate .prototxt file with the parameters for learning (learning rate, which network description file is being used, etc.) to the caffe command line utility.


Running the Collision Avoidance Benchmarks:
-------------------------------------------
All files for the collision avoidance benchmark can be found in the folder "casestudies/vehicleCollision". 

The script "generator.py" can be run with "./generator.py" and is generating collision yes/no cases. They are stored in the file "collisions.csv". The GIT repository already contains the collision cases used as benchmarks, but by letting the "generator.py" script run, more lines are added to the CSV file.

Before the CSV file can be used in Caffe, it has to be brought into a compatible formar. The script "build_database.sh" builds the HDF5 version of "collisions.csv"

Next is the training step. Script "./train.sh" starts Caffe on the task. The network structure is defined in "caffe_net_with_pooling.prototxt", while the file "caffe_solver_with_pooling.prototxt" stores Caffe's parameters. The result is put into the "snapshots" directory.

The "verify.sh" script then starts the robustness verification process of the learned network (note that this only really makes sense if the network turned out to have an accurracy of 100 percent). It first translates the network to a JSON file. The JSON file is then used in the script "json_net_against_csv_robustness_verifier.py", which performs the robustness verification. The latter script translates the JSON network description to a .rlv network description as its first step of operation. 

Running the MNIST Examples:
-------------------------------------------

After learning the network (by running the train.sh script), the following steps can be performed:

- First, run "./testAccuracy.sh" in the casestudies/MNIST directory. The process can be aborted once the first digit results start coming in.
- Then, run "../../tools/rlvCompressor.py testNetwork.rlv > testNetworkB.rlv" to get a slightly optimized version of the benchmark.
- Then, run "./prodNetwork.py testNetworkB.rlv GIVESTRONG 2 <distance>" to obtain an image that is strongly classified as a 2.

