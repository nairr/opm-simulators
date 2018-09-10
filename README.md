# Open Porous Media Simulators and Automatic Differentiation Library [![Build Status](https://travis-ci.org/OPM/opm-simulators.svg?branch=master)](https://travis-ci.org/OPM/opm-simulators)


CONTENT
-------

opm-simulators contains simulator programs for porous media flow. It
also contains a small library for automatic differentiation
built on the Eigen linear algebra package which is used by many of the
simulators to handle the building of Jacobians. The most important parts are:

* flow.cpp (a fully implicit black-oil simulator)
* AutoDiffBlock.hpp (class for AD on vectorized data with sparse jacobians)

LICENSE
-------

The library is distributed under the GNU General Public License,
version 3 or later (GPLv3+).


PLATFORMS
---------

The opm-simulators module is designed to run on Linux platforms. It is
also regularly run on Mac OS X. No efforts have been made to ensure
that the code will compile and run on windows platforms.


REQUIREMENTS
------------

opm-simulators requires opm-output, opm-core, and all their
requirements (see opm-core/README). In addition, opm-simulators
requires the Dune modue dune-istl and Eigen, version 3.1 (has not been
tested with later versions).


DOWNLOADING
-----------

For a read-only download:
git clone git://github.com/OPM/opm-simulators.git

If you want to contribute, fork OPM/opm-simulators on github.


BUILDING
--------

See build instructions at http://opm-project.org/?page_id=36


DOCUMENTATION
-------------

Efforts have been made to document the code with Doxygen.
In order to build the documentation, enter the command

 make doc

in the topmost directory. The class AutoDiffBlock is the most
important and most well-documented.


REPORTING ISSUES
----------------

Issues can be reported in the Git issue tracker online at:

    https://github.com/OPM/opm-simulators/issues

To help diagnose build errors, please provide a link to a build log together
with the issue description.

You can capture such a log from the build using the `script' utility, e.g.:

    LOGFILE=$(date +%Y%m%d-%H%M-)build.log ;
    cmake -E cmake_echo_color --cyan --bold "Log file: $LOGFILE" ;
    script -q $LOGFILE -c 'cmake ../opm-core -DCMAKE_BUILD_TYPE=Debug' &&
    script -q $LOGFILE -a -c 'ionice nice make -j 4 -l 3' ||
    cat CMakeCache.txt CMakeFiles/CMake*.log >> $LOGFILE

The resulting file can be uploaded to for instance gist.github.com.

"First Changes Made by Wessel. This is a check to undeerstand the system"
