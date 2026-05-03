import subprocess
import shutil
from os import path
import os
from multiprocessing import Pool
from functools import partial
from tempfile import NamedTemporaryFile
import re

TARPMIR_LOCATION = './TarPmiR_Linux'
MIRNA_FILE = 'mirnatest.fasta'
MRNA_FILE = 'cDNA.fasta'
PROBABILITY_CUTOFF = '0'
TARPMIR_CONDA_ENV = 'TarPmiR'


def runTarPmiR(conda_env, miRNA_file, mRNA_file, probability_cutoff):
    args = [
        './TarPmiR.py',
        '-a', path.join('.', miRNA_file),
        '-b', path.join('.', mRNA_file),
        '-m', './models/Human_sklearn_0.22.pkl',
        '-p', probability_cutoff
    ]

    condaExe = shutil.which('conda')
    if condaExe:
        command = [condaExe, 'run', '-n', conda_env, 'python'] + args
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


def runTarPmiRWithIds(conda_env, miRNA_file, mRNA_ids, probability_cutoff):
    temp = NamedTemporaryFile(mode='w', suffix='.fasta', delete=False, encoding='utf-8')
    temp.write(mRNA_ids)
    
    args = [
        './TarPmiR.py',
        '-a', path.join('.', miRNA_file),
        '-b', temp.name,
        '-m', './models/Human_sklearn_0.22.pkl',
        '-p', probability_cutoff
    ]

    condaExe = shutil.which('conda')
    if condaExe:
        command = [condaExe, 'run', '-n', conda_env, 'python'] + args
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
        
    return temp.name.split('/')[-1]


def batchTarPmiR(conda_env, miRNA_file, mRNA_file, probability_cutoff, processes):
    def makeChunks(mRNA_file, process_num):
        with open(mRNA_file, mode='r') as file:
            seqs = file.read().split('>')[1:]
            mRNAList = [f'>{seq}' for seq in seqs]
        
        chunkSize = len(mRNAList) // process_num
        chunks = [''.join(mRNAList[i:i+chunkSize]) for i in range(0, len(mRNAList), chunkSize)]
        if len(mRNAList) % process_num != 0:
            chunks[process_num-1] += chunks.pop()
        return chunks
    
    with Pool(processes) as pool:
        tempFiles = pool.map(partial(runTarPmiRWithIds, conda_env, miRNA_file, probability_cutoff=probability_cutoff), makeChunks(mRNA_file, processes))
    
    with open(path.join(TARPMIR_LOCATION, f'{MIRNA_FILE}_{MRNA_FILE}.bp'), 'w') as output:
        for tempFile in tempFiles:
            fileName = path.join(TARPMIR_LOCATION, f'{MIRNA_FILE}_{tempFile}.bp')
            with open(fileName, 'r') as file:
                part = file.read()
            output.write(re.sub(r'(?<=.)hsa', r'\nhsa', part))
            if tempFile != tempFiles[-1]:
                output.write('\n')
            os.remove(fileName)
            


batchTarPmiR(TARPMIR_CONDA_ENV, MIRNA_FILE, path.join(TARPMIR_LOCATION, MRNA_FILE), PROBABILITY_CUTOFF, 3)

# runTarPmiR(TARPMIR_CONDA_ENV, MIRNA_FILE, MRNA_FILE, PROBABILITY_CUTOFF)

# with open(path.join(TARPMIR_LOCATION, MRNA_FILE), mode='r') as file:
#     ids = file.read()
#     runTarPmiRWithIds(TARPMIR_CONDA_ENV, MIRNA_FILE, ids, PROBABILITY_CUTOFF)