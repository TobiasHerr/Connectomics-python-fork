# TE-Causality

This package consists of the full source code used for the Transfer Entropy (TE) based measure of effective connectivity (called Generalized Transfer Entropy, or GTE) published in [Model-free Reconstruction of Neuronal Network Connectivity from Calcium Imaging Signals](http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1002653).

Its main purpose is to calculate pairwise causal influences using one of the following methods:

- *xc*, by finding the peak in the cross-correlogram between the two time series
- *mi*, by finding the lag with the largest Mutual Information
- *gc*, by computing the Granger Causality, based on: C.W.J. Granger, [Investigating Causal Relations by Econometric Models and Cross-Spectral Methods](http://www.jstor.org/stable/1912791) , Econometrica, 1969
- *te-extended*, by computing GTE as defined above.

It also contains four methods of estimating TE and GTE without binning:

- *te-binless-Leonenko*, based on: L.F. Kozachenko and N.N. Leonenko, [Sample Estimate of the Entropy of a Random Vector](http://www.mathnet.ru/php/archive.phtml?wshow=paper&jrnid=ppi&paperid=797&option_lang=eng), Problems of Information Transmission, 1987
- *te-binless-Kraskov*, based on: A. Kraskov et al., [Estimating Mutual Information](http://pre.aps.org/abstract/PRE/v69/i6/e066138), Physical Review E, 2004
- *te-binless-Frenzel*, based on: S. Frenzel and B. Pompe, [Partial Mutual Information for Coupling Analysis of Multivariate Time Series](http://prl.aps.org/abstract/PRL/v99/i20/e204101), Physical Review Letters, 2007
- *te-symbolic* (experimental) based on: M. Staniek and K. Lehnertz, [Symbolic Transfer Entropy](http://link.aps.org/doi/10.1103/PhysRevLett.100.158101), Physical Review Letters, 2008.


## Step by step installation guide

### Ubuntu (precise server)

If you are on a fresh install of ubuntu server you might need some very
basic stuff, like make, rake, gcc, g++ and git, so install them:

    sudo apt-get install make, rake, gcc, g++, git

Now let's start with the real dependencies, first boost:

    sudo apt-get install libboost-dev

Continue with the GSL:

    sudo apt-get install libgsl0-dev

And finally SimKernel. You can either download it or just clone the
repository:

    git clone https://github.com/ChristophKirst/SimKernel.git
    cd SimKernel
    make lib
    sudo make install

Now you should have SimKernel installed as a library at /usr/local. It
is time to install TE-Causality (let's use the challenge branch):

    cd ~
    git clone -b connectomics_challenge https://github.com/dherkova/te-causality.git
    cd te-causality

Since we have skipped installing yaml-cpp (because we don't need that
functionality) we have to undefine something in the file
`te-datainit.h`, change line 37:

    #define ENABLE_YAML_IMPORT_AT_COMPILE_TIME

for:

    #undef ENABLE_YAML_IMPORT_AT_COMPILE_TIME

To directly have outputs compatible with the Challenge we also need to
edit something in the causality files. Here we will only do it for
te-extended. So edit `transferentropy-sim/te-extended.cpp` and make sure
line 42 reads:

    #define FORMAT_TEXT_OUTPUT_FOR_ML_CHALLENGE

We are almost done. We still need to change something in the rakefile so
edit `transferentropy-sim/Rakefile` and change line 32:

    ld_flags_basic = "-lgsl -lgslcblas -lm -lyaml-cpp -L. -lsim -lrt"

for:

    ld_flags_basic = "-lgsl -lgslcblas -lm -L. -lsim -lrt"

We have removed the yaml-cpp dependency. We are  all set! time to compile:

    cd transferentropy-sim
    rake te-extended

Now let's try it out. First download some datasets, let's start with the
small networks, download the data set [small](https://www.kaggle.com/c/connectomics/download/small.tgz)
and store it somewhere (let's call it the challenge folder). Now extract
its contents:

    tar -xvf small.tgz

Also copy the `te-extended` executable to the same folder and create
the following control file `control.txt`:

```c++
// -- SimKernel control file --

size = 100;

FormatOutputForMathematica = False;

// word length
p = 2;
SourceMarkovOrder = p;
TargetMarkovOrder = p;
StartSampleIndex = p;
globalbins = 2;

bins = 2;
binEdges = {-10.0, 0.12, 10.0};

// Data paths
// basedir = "/Users/dherkova/Dropbox/Projects/GTE-Challenge/tests/"; // for UNIX-QSub
basedir = ""; // for local
outputpath = basedir+"";

// Input data
baseFile = "iNet1_Size100_CC03inh";
inputfile = basedir + "fluorescence_"+baseFile+".txt";

// Tag used for the kaggle format
NetworkTag = "test";

samples = 150000;

HighPassFilterQ = True;
InstantFeedbackTermQ = True;
IncludeGlobalSignalQ = True;

// RelativeGlobalConditioningLevelQ = False;
//condList = Table[i*0.025, {i, 1, 25, 1}];
//iCond = Iterator[i, {i, 0, Length[condList]-1}];

conditioningList = {0.05, 0.10};
iC = Iterator[j,{j,0,Length[conditioningList]-1,1}];
GlobalConditioningLevel = conditioningList[[iC]];
CLLabel = ToString[iC];

// Output files
fileindex = ToString[Iteration[]];
outputfile = outputpath+"scores_"+baseFile+"_"+CLLabel+".csv";
outputparsfile = outputpath+baseFile+"_"+CLLabel+"_pars.mx";
```

Pay special attention to the sections in the control file regardin the
folders, input and output files and modify them accordingly. Now if you
run:

    ./te-extended control.txt

Everything should work, the program should do two iterations over 2
different conditioning levels (as defined by the conditioningList
iterator in the control file). If you only want to iterate over the
second conditioning level you can run:

    ./te-extended control.txt 2

If you want more info about how to set the control file read the
[SimKernel documentation](https://github.com/ChristophKirst/SimKernel)
and further down here.

If everything went right you sould get a scores file called
`scores_iNet1_Size100_CC03inh_1.csv` with the scores in Kaggle format.

Now if you use MATLAB you could load the scores with something like:

```matlab
scoresKaggle = dlmread('scores_iNet1_Size100_CC03inh_1.csv',',',1,1);
% Scores should be a complete square matrix, so let's hack it back to
% matrix form
scores = zeros(sqrt(length(scoresKaggle)));
cidx = 1;
for j = 1:length(scores)
    for i = 1:length(scores)
        scores(i,j) = scoresKaggle(cidx);
        cidx = cidx+1;
    end
end
```

To load the network you can just do:

```matlab
networkData = load('network_iNet1_Size100_CC03inh.txt');
N = max(max(networkData(:,1:2)));
network.RS = sparse(networkData(:,1), networkData(:,2), networkData(:,3), N, N);
network.RS(network.RS < 0) = 0;
% No need to be sparse for AUC and ROC computations
network.RS = full(network.RS);
```

And now you can compute the AUC with the provided `calculateROC`
function:

```matlab
[AUC, FPR, TPR, TPRatMark, raw] = calculateROC(network, scores, 'plot', true);
```

Voila! You are done! This should give you an AUC around 0.75.


## Dependencies

### GTE

To compile the GTE binaries based on binned estimates, all you need is a standard C++ compiler, and the following libraries:

- [GNU Scientific Library (GSL)](http://www.gnu.org/s/gsl/)
- [Boost](http://www.boost.org/)
- and Christoph Kirst's [SimKernel](https://github.com/ChristophKirst/SimKernel) package

Please make sure that GSL and Boost are available to your C++ linker, and that the path to the SimKernel files is correctly set in the Rakefile.

### Light scattering

To simulate light scattering, we need to load the spatial positions of each node from a YAML file. You therefore need to have the [yaml-cpp](http://code.google.com/p/yaml-cpp) libraries installed (version 0.5.0 or greater).

### Binless methods

To use the binless estimators, you also need to install Marius Muja's excellent [FLANN](https://github.com/mariusmuja/flann) (Fast Library for Approximate Nearest Neighbors) package and make it available to your linker.



## Input file formats

See below for details on the file formats of the input data to the causality programs. You will need a SimKernel control file, and and some input time series, either directly loaded from disk, or in the form of spike trains. In the latter case, the `te-datainit` library will simulate a fluorescence signal.

Note that minimal sample files are included in the `transferentropy-sim/tests` directory.

### SimKernel control file

All of the reconstruction and signal parameters are set via SimKernel. See the [SimKernel repository](https://github.com/ChristophKirst/SimKernel) for general documentation.

In principle, all parameters like `size` for the number of nodes or `SourceMarkovOrder` for the order of the assumed Markov process are set in this file and can be changed without re-compiling the GTE programs.

Example control file:

```c++
// word length
p = 2; 
SourceMarkovOrder = p;
TargetMarkovOrder = p;
StartSampleIndex = p;

bins = 3;

// input data
size = 100;
spikeindexfile = "data/indices.txt";
spiketimesfile = "data/times.txt";
tauF = 20;

samples = 1000;

FluorescenceModel = "Leogang";
DeltaCalciumOnAP = 50.;
tauCa = 1000.;
noise = 0.03;
saturation = 300;
tauF = 20;

globalbins = 2;
OverrideRescalingQ = False;
HighPassFilterQ = True;
InstantFeedbackTermQ = True;
IncludeGlobalSignalQ = True;
GlobalConditioningLevel = 0.5;

// output files
outputfile = "adjA.mx";
outputparsfile = "pars.mx";
```

This could be a standard setting to simulate a fluorescence signal based on the given spike times, and the calculate the GTE matrix conditioned on the average fluorescence via the parameter `GlobalConditioningLevel`.

One of the main features of SimKernel are iterators. For instance, to calculate GTE for different numbers of samples, one could use code such as the following:

```c++
samplesList = {10,50,100,500,1000,5000,10000,15000};
iS = Iterator[j,{j,0,Length[samplesList]-1,1}];
samples = samplesList[[iS]];
```

Note the double brackets (Mathematica syntax), but that the array index is starting from zero.

### Spike trains

Spike times are loaded simply as line break separated entries in two files of equal length. One containing spike times (SimKernel parameter `spiketimesfile`) and one listing the corresponding node indices (SimKernel parameter `spikeindexfile`).

### Fluorescence time series (or any other continuous signal)

Alternatively, you can simply load an already existing time series (in ASCII format) directly via the SimKernel parameter `inputfile`.

### Network topology

For the purpose of simulating light scattering, we need to know the spatial position of each node. We therefore load a YAML file of the following format:

```yaml
---
size: 3 # number of nodes
cons: 4 # number of connections
notes: "example miniature topology"
createdAt: "Tue 20 Sep 2011 14:55:21"
nodes:
  - id: 1
    pos: [0.516609, 0.187339]
    connectedTo: [2,3]
  - id: 2
    pos: [0.933885, 0.0615908]
    connectedTo: [1]
  - id: 3
    pos: [0.519384, 0.110653]
    connectedTo: [2]
```

### Result

The result of the computation are two files: A parameter file and a file stating the computed effective connectivity matrix. Both files are written in [Mathematica List](http://reference.wolfram.com/mathematica/ref/List.html) format. An example of the connectivity file using `te-extended` might be the following:

```c++
{{{0.000000000000000},{0.142936610901646}},{{0.053183796775799},{0.000000000000000}}}
```

In Mathematica syntax, this corresponds to a 2x2 matrix with zeros on the diagonal. The value of 0.14 is then the computed GTE in bits which is the effective weight of the causal link from the first time series to the second time series.


## Copyright

All of the files (with the exception of `transferentropy-sim/tests/catch.hpp`, taken from [here](https://github.com/philsquared/Catch), which is licensed under the Boost licence) can be copied, modified, used for commercial or non-commercial purpose, as long as you keep the following copyright message in the source files:

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see [here](http://www.gnu.org/licenses/).
