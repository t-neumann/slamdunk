wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
tar xjvf samtools-1.3.1.tar.bz2
cd samtools-1.3.1
./configure
make
cd ..
rm samtools-1.3.1.tar.bz2
ln -s samtools-1.3.1/samtools