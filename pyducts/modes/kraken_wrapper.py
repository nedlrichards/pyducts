import numpy as np
from os import path
import subprocess

at_path = path.abspath(path.expanduser('~/at/bin'))

class Kraken:
    """tracks a kraken run"""
    def __init__(self, file_name, r_rcr, z_rcr, z_src=None, c_bounds=[0., 3000.]):
        """all runs will work with the same file name"""
        self.file_name = path.splitext(file_name)[0]
        self.c_bounds = c_bounds
        self.options = "'NVW'"
        self.z_src = np.abs(z_src) if z_src is not None else None
        self.r_rcr = np.array(r_rcr, ndmin=1)
        self.z_rcr = np.array(np.abs(z_rcr), ndmin=1)

    def run_kraken(self):
        """run kraken on an env file"""
        run_cmd = [path.join(at_path, 'kraken.exe'), self.file_name]
        #run_cmd = [path.join(at_path, 'krakenc.exe'), self.file_name]
        # TODO: see if status gives anything
        status = subprocess.run(run_cmd, check=True)

    def write_env(self, fc, z, ssp, bottom_HS=[1600., 1800., 0.5], append=False):
        """basic capability env writting"""
        env_file = self.file_name + '.env'

        # front matter of env file
        title = "'auto_gen'"
        z = np.abs(z)
        sigma = 0

        # check if last recevier is in the bottom HS
        zmax = max(z[-1], self.z_rcr[-1])
        nmedia = 2 if zmax > z[-1] + 0.5 else 1

        # back matter of env file
        bot_options = 'A'
        bottom_c, bottom_rho = bottom_HS[:2]
        if len(bottom_HS) == 3:
            bottom_attn = bottom_HS[2]
        else:
            bottom_attn = 0.0

        # change of dentisty to cgs for kraken
        bottom_rho /= 1000.

        decimation = 15
        num_points = int(decimation * np.ceil(z[-1] * fc / np.min(ssp)))

        write_string = 'a' if append else 'w'

        with open(env_file, write_string) as writer:
            writer.write("\n".join([title, f'{float(fc):.2f}', str(nmedia), self.options]))
            writer.write("\n")
            writer.write(f'{num_points} \t {sigma:.2f} \t {z[-1]:.2f}\n')

            np.savetxt(writer, np.array([z, ssp]).T, fmt='%.2f', newline='\t /\n', delimiter='\t')
            if nmedia == 2:
                writer.write(f'0 \t 0.0 \t {zmax:.2f}')
                writer.write(f"\n{z[-1]:.2f} \t {bottom_c:.2f} \t 0.0 \t {bottom_rho:.2f} \t {bottom_attn:.2f} /")
                writer.write(f"\n{zmax:.2f} \t {bottom_c:.2f} \t 0.0 \t {bottom_rho:.2f} \t {bottom_attn:.2f} /")
            writer.write("\n'" + bot_options + "'" + "\t 0.0" )
            writer.write(f"\n{zmax:.2f} \t {bottom_c:.2f} \t 0.0 \t {bottom_rho:.2f} \t {bottom_attn:.2f} /")
            writer.write(f"\n{self.c_bounds[0]} \t {self.c_bounds[1]}")
            writer.write(f"\n1000")
            writer.write(f"\n{self.z_rcr.size}")
            writer.write(f"\n{self.z_rcr[0]:.2f} \t {self.z_rcr[-1]:.2f} \t /")
            if self.z_src is None:
                writer.write(f"\n1\n {self.z_rcr[-1]:.2f} /")
            else:
                writer.write(f"\n1\n {self.z_src:.2f} /")
            writer.write(f"\n")

    def run_field(self, option='RA', raxis=None, num_modes=9999):
        """run kraken on an env file"""
        if self.z_src is None:
            raise(ValueError("No source depth specified"))
        flp_file = self.file_name + '.flp'

        if raxis is None:
            nprof = 1
        else:
            raxis = np.array(raxis)
            nprof = raxis.size
            rmax = raxis[-1]
            is_uniform = np.diff(raxis)
            is_uniform = np.all(np.abs(is_uniform - is_uniform[0]) < 1e-3)

        with open(flp_file, 'w') as writer:
            writer.write(f"/,\n")
            writer.write("'"+option+"'\n")
            writer.write(f"{num_modes}\n")
            writer.write(f"{nprof}\n")
            if raxis is None:
                writer.write("0.0\n")
            elif is_uniform:
                writer.write(f"0.0\t{rmax/1e3:.2f} /\n")
            else:
                raxis = np.array(raxis, ndmin=2) / 1e3
                np.savetxt(writer, raxis, fmt='%.2f', newline='\t /\n', delimiter=' ')
            if self.r_rcr.size == 1:
                writer.write("1\n")
                writer.write(f"{self.r_rcr[0]/1e3:.2f} /\n")
            else:
                writer.write(f"{self.r_rcr.size}\n")
                writer.write(f"{self.r_rcr[0]/1e3:.2f}\t{self.r_rcr[-1]/1e3:.2f} /\n")
            # assume a single source depth
            writer.write("1\n")
            writer.write(f"{self.z_src:.2f} /\n")
            if self.z_rcr.size == 1:
                writer.write("1\n")
                writer.write(f"{self.z_rcr[0]:.2f} /\n")
            else:
                writer.write(f"{self.z_rcr.size}\n")
                writer.write(f"{self.z_rcr[0]:.2f}\t{self.z_rcr[-1]:.2f} /\n")
            # assume zero array tilt
            if self.z_rcr.size == 1:
                writer.write("1\n")
                writer.write(f"0.0 /\n")
            else:
                writer.write(f"{self.z_rcr.size}\n")
                writer.write(f"0.0\t0.0 /\n")

        run_cmd = [path.join(at_path, 'field.exe'), self.file_name]
        # TODO: see if status gives anything
        status = subprocess.run(run_cmd, check=True)


def read_mod(envpath):
    """read binary mode file output by kraken"""
    envpath = path.splitext(envpath)[0] + '.mod'
    with open(envpath, 'rb') as f:
        lrec = int(4 * np.fromfile(f, dtype='i4', count=1))
        i0 = 0
        f.seek(4)
        title = f.read(80)

        f.seek(i0 + 4)
        t_test = f.read(80)

        freq = np.fromfile(f, dtype='f4', count=1)
        nums = np.fromfile(f, dtype='i4', count=3)

        if nums[1] < 0:
            raise Exception('No modes in file')

        #sample depths
        f.seek(i0 + 4 * lrec)
        z = np.fromfile(f, dtype='f4', count = nums[1])

        #number of modes
        f.seek(i0 + 5 * lrec)
        M = int(np.fromfile(f, dtype='i4', count=1))

        #modes
        phi = np.zeros((M, nums[2]), dtype=np.complex128)
        for i in range(M):
            rec = f.seek(i0 + (i + 7) * lrec)
            phi_c = np.fromfile(f, dtype='f4', count=2 * nums[2])
            phi[i, :] = phi_c[::2] + 1j * phi_c[1::2]

        #wave numbers
        f.seek(i0 + (7 + M) * lrec)
        k_c = np.fromfile(f, 'f4', count=2 * M)
        # Don't use kraken's phase convention
        k = k_c[::2] - 1j * k_c[1::2]

    # match units of denisty
    phi *= np.sqrt(1000.)

    return phi.astype(np.complex128), k.astype(np.complex128), z.astype(np.float64)

def read_shd(shdpath):
    """read pressure from kraken shd file"""
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

        all_pressure = None

        for irz in np.arange(Nrz):
            recnum = 10 + irz
            status = f.seek(recnum * 4 * recl)

            temp = np.fromfile(f, dtype='float32', count=2 * Nrr)
            pressure = temp[::2] + 1j * temp[1::2]

            if all_pressure is None:
                all_pressure = [pressure]
            else:
                all_pressure.append(pressure)
        all_pressure = np.squeeze(np.array(all_pressure, dtype=np.complex128))

    return z.astype(np.float64), r.astype(np.float64), all_pressure

def synthesize_modes(modes_src, modes_rcr, kr, r_axis):
    """synthesize modes into a pressure field"""
    pressure = 1j * np.exp(-1j * pi / 4) * modes_src * modes*rcr \
             * np.exp(1j * kr * r_axis) / np.sqrt(8 * pi * kr * r_axis)
    return pressure
