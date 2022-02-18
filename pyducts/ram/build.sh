python setup.py clean
rm *.so
rm src/epade.c
rm src/ram_main_v4.c
rm src/*.o
rm src/*.mod
cython src/epade.pyx -3
cython src/ram_main_v4.pyx -3
python setup.py build_ext --inplace
