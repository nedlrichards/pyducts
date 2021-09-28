python setup.py clean
rm *.so
rm src/epade.c
rm src/ram_main_v4.c
rm src/*.o
rm src/*.mod
cython src/epade.pyx -3
cython src/ram_main_v4.pyx -3
#cython src/tmp.pyx -3

#gfortran -c -o src/cmplx_roots_sg.f90.o src/cmplx_roots_sg.f90
#gfortran -c -o src/constants.f90.o src/constants.f90
#gfortran -c -o src/pade_coeffs.f90.o src/pade_coeffs.f90
#gfortran -c -o src/pade_coeffs_wrapper.f90.o src/pade_coeffs_wrapper.f90
#gfortran -c -o src/ram_io.f90.o src/ram_io.f90
#gfortran -c -o src/ram_v4.f90.o src/ram_v4.f90
#gfortran -c -o src/ram_v4_wrapper.f90.o src/ram_v4_wrapper.f90
#gfortran -c -o src/tmp_wrapper.f90.o src/tmp_wrapper.f90

#gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I/home/nedrichards/miniconda3/envs/py3/include/python3.9 -I/home/nedrichards/miniconda3/envs/py3/lib/python3.9/site-packages/numpy/core/include -o tmp.so src/tmp.c src/pade_coeffs_wrapper.f90.o src/cmplx_roots_sg.f90.o src/constants.f90.o src/pade_coeffs.f90.o src/ram_io.f90.o src/ram_v4.f90.o src/tmp_wrapper.f90.o -lgfortran

#gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I/home/nedrichards/miniconda3/envs/py3/include/python3.9 -I/home/nedrichards/miniconda3/envs/py3/lib/python3.9/site-packages/numpy/core/include -o epade.so src/epade.c src/pade_coeffs_wrapper.f90.o src/ram_main_v4.c src/cmplx_roots_sg.f90.o src/constants.f90.o src/pade_coeffs.f90.o src/ram_io.f90.o src/ram_v4.f90.o src/ram_v4_wrapper.f90.o -lgfortran

#gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I/home/nedrichards/miniconda3/envs/py3/include/python3.9 -I/home/nedrichards/miniconda3/envs/py3/lib/python3.9/site-packages/numpy/core/include -o ram_main_v4.so src/ram_main_v4.c src/pade_coeffs_wrapper.f90.o src/cmplx_roots_sg.f90.o src/constants.f90.o src/pade_coeffs.f90.o src/ram_io.f90.o src/ram_v4.f90.o src/ram_v4_wrapper.f90.o -lgfortran

python setup.py build_ext --inplace
