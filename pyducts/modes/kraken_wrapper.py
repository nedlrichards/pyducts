import numpy as np
from os import path
import subprocess

at_path = path.abspath(path.expanduser('~/at/bin'))

def run_kraken(file_name):
    """run kraken on an env file"""
    env_file = path.splitext(file_name)[0]

    run_cmd = [path.join(at_path, 'kraken.exe'), env_file]
    # TODO: see if status gives anything
    status = subprocess.run(run_cmd, check=True)
    return

def run_field(flp_file):
    """run kraken on an env file"""
    flp_file = path.splitext(flp_file)[0]

    run_cmd = [path.join(at_path, 'field.exe'), flp_file]
    # TODO: see if status gives anything
    status = subprocess.run(run_cmd, check=True)
    return

def write_env(file_name, fc, z, ssp, bottom_HS =[1600., 1.8]):
    """basic capability env writting"""

    fn = path.abspath(file_name)
    fn_ending = fn.split('.')

    if len(fn_ending) == 1:
        env_file = fn + '.env'
    elif fn_ending[-1] != 'env':
        ValueError('unrecognized file extension')
    else:
        env_file = fn


    # front matter of env file
    title = "'auto_gen'"
    nmedia = 1
    options = "'NVF'"

    # ssp
    z = np.abs(z)
    sigma = 0
    z_end = z[-1]

    # back matter of env file
    bot_options = 'A'
    bottom_c, bottom_rho = bottom_HS

    ps_limits = [0, 3000]  # let kraken bound phase speeds
    num_points = int(10 * np.ceil(z_end * fc / np.min(ssp)))


    with open(env_file, 'w') as writer:
        writer.write("\n".join([title, f'{float(fc):.2f}', str(nmedia), options]))
        writer.write("\n")
        writer.write(f'{num_points} \t {sigma:.2f} \t {z_end:.2f}')
        writer.write("\n")

        np.savetxt(writer, np.array([z, ssp]).T, fmt='%.2f', newline='\t /\n', delimiter='\t')
        writer.write("'" + bot_options + "'" + "\t 0.0" )
        writer.write(f"\n{z_end:.2f} \t {bottom_c:.2f} \t 0.0 \t {bottom_rho:.2f} \t 0.0 /")
        writer.write(f"\n{ps_limits[0]} \t {ps_limits[1]}")
        writer.write(f"\n1000")
        writer.write(f"\n{num_points}")
        writer.write(f"\n{z[0]:.2f} \t {z[-2]:.2f} \t /")
        writer.write(f"\n1 \n {z[-1]:.2f} /")
        writer.write(f"\n")

def read_mod(envpath):
    """read binary mode file output by kraken"""
    envpath = path.splitext(envpath)[0] + '.mod'
    with open(envpath, 'rb') as f:
        lrec = int(4 * np.fromfile(f, dtype='i4', count=1))
        f.seek(4)
        title = f.read(80)
        freq = np.fromfile(f, dtype='f4', count=1)
        nums = np.fromfile(f, dtype='i4', count=3)

        if nums[1] < 0:
            raise Exception('No modes in file')

        f.seek(lrec)
        num=[]
        medium_type = []
        for _ in range(nums[0]):
            num.append(np.fromfile(f, dtype='i4', count=1))
            medium_type.append(f.read(8))

        #sample depths
        f.seek(4 * lrec)
        z = np.fromfile(f, dtype='f4', count = nums[1])

        #number of modes
        f.seek(5 * lrec)
        M = int(np.fromfile(f, dtype='i4', count=1))

        #modes
        phi = np.zeros((M, nums[2]), dtype='complex_')
        for i in range(M):
            rec = f.seek((i+7) * lrec)
            phi_c = np.fromfile(f, dtype='f4', count=2 * nums[2])
            phi[i, :] = phi_c[::2] + 1j * phi_c[1::2]

        #wave numbers
        f.seek((7 + M) * lrec)
        k_c = np.fromfile(f, 'f4', count=2 * M)
        k = k_c[::2] + 1j * k_c[1::2]

    return phi, k, z

def read_shd(shdpath):
    shd_file = path.splitext(shdpath)[0]
    shd_file += '.shd'

    with open(shd_file, 'rb') as f:
        recl = np.fromfile(f, dtype='int32', count=1)[0] # record length in bytes will be 4*recl

        f.seek(2 * 4 * recl) # reposition to end of second record
        Nfreq  = np.fromfile(f, dtype='int32', count=1)[0]
        Ntheta = np.fromfile(f, dtype='int32', count=1)[0]

        Nsx = np.fromfile(f, dtype='int32', count=1)[0]
        Nsy = np.fromfile(f, dtype='int32', count=1)[0]
        Nsz = np.fromfile(f, dtype='int32', count=1)[0]

        Nrz = np.fromfile(f, dtype='int32', count=1)[0]
        Nrr = np.fromfile(f, dtype='int32', count=1)[0]
        freq0 = np.fromfile(f, dtype='float32', count=1)[0]
        atten = np.fromfile(f, dtype='float32', count=1)[0]

        # read receiver postions
        f.seek(8 * 4 * recl)
        z = np.fromfile(f, dtype='float32', count=Nrz)

        f.seek(9 * 4 * recl)
        r = np.fromfile(f, dtype='float32', count=Nrr)
        all_pressure = []

        for irz in np.arange(Nrz - 1):
            recnum = 10 + irz
            status = f.seek(recnum * 4 * recl)

            temp = np.fromfile(f, dtype='float32', count=2 * Nrr)
            pressure = temp[::2] + 1j * temp[1::2]

            all_pressure.append(pressure)
        all_pressure = np.array(all_pressure)

    return z, r, all_pressure
