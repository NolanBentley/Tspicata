#From https://github.com/primer3-org/primer3
cd ~/
git clone https://github.com/primer3-org/primer3.git primer3
cd primer3/src
make
make test
#Test
~/primer3/primer3_core ../example
