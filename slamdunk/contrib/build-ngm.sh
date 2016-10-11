#git clone https://github.com/Cibiv/NextGenMap.git

cd ../../nextgenmap
mkdir build
cd build
cmake ..
make -j4

cd ../../bin
ln -fs ../nextgenmap/bin/ngm-*/ngm ngm

# Update
# git fetch nextgenmap master
# git subtree pull --prefix nextgenmap nextgenmap master --squash