import subprocess
import shutil
from os import path

TARPMIR_LOCATION = './TarPmiR_Linux'
MIRNA_FILE = 'mirnatest.fasta'
MRNA_FILE = 'cDNA.fasta'
PROBABILITY_CUTOFF = '0.5'
TARPMIR_CONDA_ENV = 'TarPmiR'

def runTarPmiR(conda_env, miRNA_file, mRNA_file, probability_cutoff):
    args = [
        './TarPmiR.py',
        '-a', path.join('.', miRNA_file),
        '-b', path.join('.', mRNA_file),
        '-m', './models/Human_sklearn_0.22.pkl',
        '-p', probability_cutoff
    ]

    conda_exe = shutil.which('conda')
    if conda_exe:
        command = [conda_exe, 'run', '-n', conda_env, 'python'] + args
    else:
        raise RuntimeError('Could not find `conda` on this system')

    try:
        result = subprocess.run(
            command,
            text=True,
            cwd=path.join('/home/nikos/ptixiaki', TARPMIR_LOCATION),
            capture_output=True,
            check=True
        )
        print('TarPmiR stdout:', result.stdout)
    except subprocess.CalledProcessError as e:
        print(f'PID: TarPmiR failed with return code {e.returncode}:\n{e.stderr}')
    except Exception as e:
        print(f'PID: Exception occurred: {e}')

runTarPmiR(TARPMIR_CONDA_ENV, MIRNA_FILE, MRNA_FILE, PROBABILITY_CUTOFF)