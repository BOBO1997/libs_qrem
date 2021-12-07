# build

`python3 setup.py build_ext --inplace`

# install

`python3 setup.py install --record install_record.txt`
or
`cat install_record.txt | xargs rm -rf`

# uninstall

`pip uninstall libs_qrem`

# clean

`rm -r dist/ build/ libs_qrem.egg-info/`
`rm libs_qrem/*.cpp`