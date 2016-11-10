cd nextgenmap
mkdir -p build
cd build
cmake ..
make -j4

cd ../..
# Get version from version.py
version=`grep "__ngm_version__" ../version.py | cut -d "\"" -f 2`
ln -fs nextgenmap/bin/ngm-${version}/ngm ngm
