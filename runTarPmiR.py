import subprocess
from os import path
import os
from multiprocessing import Pool
from functools import partial
from tempfile import NamedTemporaryFile
import re

PROBABILITY_CUTOFF = '0.5'
MODEL = 'Human_sklearn_0.22.pkl'

def runTarPmiR(miRNA_file, mRNA_file, probability_cutoff=PROBABILITY_CUTOFF, model=MODEL):
    try:
        result = subprocess.run(
            [
                'python',
                './TarPmiR.py',
                '-a', miRNA_file,
                '-b', mRNA_file,
                '-m', 'models/' + model,
                '-p', probability_cutoff
            ],
            text=True,
            capture_output=True,
            check=True
        )
    except subprocess.CalledProcessError as e: # NOTE: add log
        print(f'PID: TarPmiR failed with return code {e.returncode}:\n{e.stderr}')
    except Exception as e: # NOTE: add log
        print(f'PID: Exception occurred: {e}')


def runTarPmiRByText(miRNA_file, mRNA_text, probability_cutoff=PROBABILITY_CUTOFF, model=MODEL):
    temp = NamedTemporaryFile(mode='w', suffix='.fasta', delete=False, encoding='utf-8')
    temp.write(mRNA_text)
    
    runTarPmiR(miRNA_file, temp.name, probability_cutoff, model) # type: ignore
        
    return temp.name.split('/')[-1]


def batchTarPmiR(processes, miRNA_file, mRNA_file, probability_cutoff=PROBABILITY_CUTOFF, model=MODEL):
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
        tempFiles = pool.map(partial(runTarPmiRByText, miRNA_file, probability_cutoff=probability_cutoff, model=model), makeChunks(mRNA_file, processes))
    
    
    with open(path.join('output/', f'{miRNA_file}_{mRNA_file}.bp'), 'w') as output:
        for tempFile in tempFiles:
            fileName = path.join('output/', f'{miRNA_file}_{tempFile}.bp')
            with open(fileName, 'r') as file:
                part = file.read()
            output.write(re.sub(r'(?<=.)hsa', r'\nhsa', part))
            if tempFile != tempFiles[-1]:
                output.write('\n')
            os.remove(fileName)
            


batchTarPmiR(10, 'mirnatest.fasta', 'benchMRNA.fasta')

# batchTarPmiR(10, 'mirnatest.fasta', 'benchMRNA.fasta', '0.5', 'Human_sklearn_0.22.pkl')
# runTarPmiR('mirnatest.fasta', 'benchMRNA.fasta', '0', 'Human_sklearn_0.22.pkl')