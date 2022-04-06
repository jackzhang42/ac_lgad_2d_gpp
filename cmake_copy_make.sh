rm -r build
mkdir build
cd build
echo "********************start cmake*********************"
cmake ..
echo "********************start copy**********************"
cp ../dynamic* .
cp ../normal* .
echo "********************start make**********************"
make
