#!/usr/bin/sh

cd matrix
cmake .
make
make doc
sudo make install

cd ..
cd interp
cmake .
make
make doc
sudo make install

cd ..
cd geom
cmake .
make
make doc
sudo make install
