sudo apt-get install libboost-all-dev libusb-1.0.0-dev python-mako doxygen python-docutils cmake build-essential
UHD_REPO_PATH = $PWD
git clone https://github.com/EttusResearch/uhd.git
cd $UHD_REPO_PATH/host
mkdir build
cd build
cmake -DENABLE_PYTHON_API=ON ../
make
make test
sudo make install
sudo /usr/local/lib/uhd/utils/uhd_images_downloader.py
export LD_LIBRARY_PATH=$UHD_REPO_PATH
export UHD_IMAGES_DIR=/usr/local/share/uhd/images
export PYTHONPATH=$UHD_REPO_PATH/host/build/python