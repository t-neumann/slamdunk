cd nextgenmap
mkdir -p build
cd build
cmake ..
make -j4

cd ../..
ln -fs nextgenmap/bin/ngm-0.5.1//ngm ngm
