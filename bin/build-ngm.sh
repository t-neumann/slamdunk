git clone https://github.com/Cibiv/NextGenMap.git

cd NextGenMap
mkdir build
cd build
cmake ..
make -j4

cd ../..
ln -s NextGenMap/bin/ngm-*/ngm
