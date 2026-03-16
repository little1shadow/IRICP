# Sparse Iterative Closest Point Algorithm #

This repository contains an implementation of the [IRICP]
## Dependencies ##

The dependencies are header-only and are all included in the ext directory. As a consequence, there is nothing to do.
For the record, here is the list of dependencies :

* [nanoflann](https://github.com/jlblancoc/nanoflann)
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main\_Page)

## Usage ##

When the program has been built thanks to CMake and Make, just run `./IRICP -h` in order to see the help.

Other useful lines : 

* `./IRICP -d` runs a demo with the files stored in the media directory
* `./IRICP -i1 /path/to/first/objectFile.obj -i2 ../path/to/second/objectFile.obj -o /path/to/output/directory/ -pl` would be a minimal line in order to use the point-to-plane variant of the sparse ICP.

Notice that you should always specify if you want to use the point-to-plane variant (`-pl`) or the point-to-point variant (`-po`).


##notice##
>* &nbsp; **this code is the ICP algorithm based on errro interval, so the experienment is not satisfied**.<br>

