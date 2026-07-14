import os
import argparse
import logging
import re
import subprocess
import time
from functools import partial
from multiprocessing import Pool
from os import path
from tempfile import NamedTemporaryFile

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('run.log'),
        logging.StreamHandler()
    ]
)
log = logging.getLogger(__name__)


##### FUNCTIONS #####
def runTarPmiR(miRNA_file:str, mRNA_file:str, probability_cutoff:str, model:str) -> None:
    '''
    Run instance of TarPmiR with args:\n
    miRNA_file: miRNA FASTA file name\n
    mRNA_file: mRNA FASTA file name\n
    probability_cutoff\n
    model: Human_sklearn_0.18.pkl, Human_sklearn_0.19.pkl, Human_sklearn_0.22.pkl, Human.pkl
    '''
    
    try:
        log.info(f'Starting TarPmiR with PID{os.getpid()} with mRNA {mRNA_file}')
        
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
        
        log.info(f'TarPmiR with PID{os.getpid()} completed successfully')
    except subprocess.CalledProcessError as e:
        log.error(f'TarPmiR with PID{os.getpid()} failed with return code {e.returncode}:\n{e.stderr}')
        
    except Exception as e:
        log.exception(f'PID{os.getpid()} exception occurred: {e}')


def runTarPmiRByText(miRNA_file:str, mRNA_text:str, probability_cutoff:str, model:str):
    '''
    Run TarPmiR with mRNA FASTA string instead of file name; creates a temporary file with the FASTA string and passes it to runTarPmiR().\n
    miRNA_file: miRNA FASTA file name\n
    mRNA_text: mRNA FASTA content\n
    probability_cutoff\n
    model: Human_sklearn_0.18.pkl, Human_sklearn_0.19.pkl, Human_sklearn_0.22.pkl, Human.pkl\n
    Returns the name of the temporary file for use in batchTarPmiR.
    '''
    
    temp = NamedTemporaryFile(mode='w', suffix='.fasta', delete=False, encoding='utf-8')
    temp.write(mRNA_text)
    
    runTarPmiR(miRNA_file, temp.name, probability_cutoff, model)
        
    return temp.name.split('/')[-1]


def batchTarPmiR(processes:int, miRNA_file:str, mRNA_file:str, probability_cutoff:str, model:str) -> None:
    '''
    Run multiple processes of runTarPmiRByText() on the same mRNA FASTA for faster processing with the following arguments:\n
    processes: Number of processes to create\n
    miRNA_file: miRNA FASTA file name\n
    mRNA_file: mRNA FASTA file name\n
    probability_cutoff\n
    model: Human_sklearn_0.18.pkl, Human_sklearn_0.19.pkl, Human_sklearn_0.22.pkl, Human.pkl
    '''
    
    log.info(f'Attempting to batch run TarPmiR using miRNA file {miRNA_file} and mRNA file {mRNA_file}')
    def makeChunks(mRNA_file, process_num):
        log.info(f'Splitting FASTA file into {processes} chunks')
        try:
            with open(mRNA_file, mode='r') as file:
                seqs = file.read().split('>')[1:]
                mRNAList = [f'>{seq}' for seq in seqs]
            
            chunkSize = len(mRNAList) // process_num
            chunks = [''.join(mRNAList[i:i+chunkSize]) for i in range(0, len(mRNAList), chunkSize)]
            if len(mRNAList) % process_num != 0:
                chunks[process_num-1] += chunks.pop()
            return chunks
        except Exception as e:
            log.exception(f'Exception occurred during chunk creation: {e}')
            log.error('Aborting')
            os.abort()
    
    
    with Pool(processes) as pool:
        log.info(f'Running TarPmiR with {processes} processes')
        start_time = time.time()
        tempFiles = pool.map(partial(runTarPmiRByText, miRNA_file, probability_cutoff=probability_cutoff, model=model), makeChunks(mRNA_file, processes))
    log.info(f'TarPmiR processes completed, time taken: {(time.time() - start_time)/60:.2f}min')
    
    log.info('Merging files into a single output file')
    with open(path.join('output/', f'{miRNA_file}_{mRNA_file}.bp'), 'w') as output:
        for tempFile in tempFiles:
            fileName = path.join('output/', f'{miRNA_file}_{tempFile}.bp')
            with open(fileName, 'r') as file:
                part = file.read()
            output.write(re.sub(r'(?<=.)hsa', r'\nhsa', part))
            if tempFile != tempFiles[-1]:
                output.write('\n')
            os.remove(fileName)
    log.info('Batch TarPmiR completed')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run TarPmiR in batch mode')
    parser.add_argument('-mrna', type=str, required=True, help='mRNA file')
    parser.add_argument('-mirna', type=str, default='mature_hsa.fa', help='miRNA file')
    parser.add_argument('-pcut', type=str, default='0.5', help='Probability cutoff')
    parser.add_argument('-model', type=str, default='Human_sklearn_0.22.pkl', help='Model filename')
    parser.add_argument('-num', type=int, default=20, help='Number of processes')

    args = parser.parse_args()
    batchTarPmiR(args.num, args.mirna, args.mrna, args.pcut, args.model)