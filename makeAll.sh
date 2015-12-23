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
cd curve
cmake .
make
make doc
sudo make install

cd ..
cd optim
cmake .
make 
make doc
sudo make install

cd ..
cd root
cmake .
make
make doc
sudo make install

cd ..
cd surface
cmake .
make
make doc
sudo make install

cd ..
cd ode
cmake .
make
make doc
sudo make install
