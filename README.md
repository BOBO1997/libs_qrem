# build

`python3 setup.py build_ext --inplace`

# install

`python3 setup.py install --record install_record.txt`
or
`cat install_record.txt | xargs rm -rf`
or
`pip install git+https://github.com/BOBO1997/libs_qrem`

# uninstall

`pip uninstall libs_qrem`

# clean

`rm -r dist/ build/ libs_qrem.egg-info/`
`rm libs_qrem/*.cpp`