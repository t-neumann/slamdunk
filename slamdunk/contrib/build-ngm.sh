git clone https://github.com/Cibiv/NextGenMap.git

cd NextGenMap
#cd nextgenmap
mkdir build
cd build
cmake ..
make -j4

cd ../..
ln -fs NextGenMap/bin/ngm-*/ngm ngm
#ln -fs nextgenmap/bin/ngm-*/ngm ngm

# Update
# git fetch nextgenmap master
# git subtree pull --prefix nextgenmap nextgenmap master --squash