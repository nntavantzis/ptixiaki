import requests
import logging
import argparse
import os
import pandas as pd
from io import StringIO
from runTarPmiR import batchTarPmiR


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('run.log'),
        logging.StreamHandler()
    ]
)
log = logging.getLogger(__name__)


###################### FUNCTIONS ######################

def getTranscriptsFromIDs(idList: 'list[str]') -> str:
    log.info(f'Trying to GET sequences of {len(idList)} IDs')
    query_xml = f"""<?xml version="1.0" encoding="UTF-8"?>
                <!DOCTYPE Query>
                <Query  virtualSchemaName = "default" formatter = "FASTA" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
                    <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
                        <Filter name = "ensembl_transcript_id" value = "{','.join(idList)}"/>
                        <Attribute name = "ensembl_transcript_id" />
                        <Attribute name = "cdna" />
                    </Dataset>
                </Query>"""

    res = requests.get("http://www.ensembl.org/biomart/martservice", params={"query": query_xml})
    log.info('GET succeeded')
    return res.text


###################### MAIN ######################

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True, help='Input VEP TSV file')
    parser.add_argument('-o', '--output', type=str, default='output_disrupted.tsv', help='Output TSV file')
    parser.add_argument('-debug', '--debug', type=int, default=None, help='FOR TESTING ONLY: Number of lines from input file to use, the rest are omitted')
    args = parser.parse_args()

    ##### VARIABLES #####    
    INPUT_FILE = args.input
    OUTPUT_FILE = args.output
    DEBUG = args.debug
    ##### FIXED VARIABLES #####
    MIRNA_FILE = 'mature_hsa.fa'
    PROBABILITY_CUTOFF = '0.5'
    MODEL = 'Human_sklearn_0.22.pkl'
    PROCESS_NUM = 20

    
    log.info('##################### STARTING PROGRAM #####################')
    log.info(f'Loading {INPUT_FILE} to Dataframe')
    df = pd.read_csv(
        INPUT_FILE,
        sep='\t',
        comment='#',
        header=None
    )
    
    log.info('File loaded, formatting Dataframe')
    # with open(INPUT_FILE) as file:
    #     for line in file.readlines():
    #         if line.startswith('#') and not line.startswith('##'):
    #             header = line.lstrip('#').strip().split('\t')
    #             df.columns = header
    #             break
    df.columns = ['Uploaded_variation', 'Location', 'Allele', 'Gene', 'Feature', 'Feature_type', 'Consequence', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'Extra']
    df = df[df['Consequence'].str.contains('3_prime_UTR_variant')]

    if DEBUG:
        df = df.head(DEBUG)
    

    transcriptsToDownload = df['Feature'].drop_duplicates().to_list()
    res = getTranscriptsFromIDs(transcriptsToDownload)
    with open('cDNA.fasta', 'w') as file:
        file.write(res)
        
    log.info('Calling a batch TarPmiR job')
    try:
        batchTarPmiR(PROCESS_NUM, MIRNA_FILE, 'cDNA.fasta', PROBABILITY_CUTOFF, MODEL)
    except:
        log.error('Aborting')
        os.abort()
        
    log.info('Editing Dataframe: Marking IDs where a sequence was found')
    with open('cDNA.fasta', 'r', encoding='utf-8') as file:
        fasta_ids = [line.strip()[1:] for line in file if line.startswith('>')]
    df.loc[:, 'Sequence_exists'] = df['Feature'].isin(fasta_ids)


    log.info('Loading TarPmiR Prediction Output to Dataframe')
    with open(f'output/{MIRNA_FILE}_cDNA.fasta.bp', 'r', encoding='utf-8') as file:
        content = file.read()
        df_pred = pd.read_csv(StringIO(content), sep='\t')
    log.info('File loaded, formatting Dataframe')
    df_pred.columns = ['miRNA','mRNA','binding_site','binding_probability','energy','seed','accessibility','AU_content','PhyloP_Stem','PyloP_Flanking','m/e','number_of_pairings','binding_region_length','longest_consecutive_pairings','position_of_longest_consecutive_pairings','pairings_in_3prime_end','difference_of_pairings_between_seed_and_3prime_end']
    df_pred.loc[:, 'miRNA'] = df_pred['miRNA'].str.split().str[0]
    df_pred[['bs_start', 'bs_end']] = (df_pred['binding_site'].str.split(',', expand=True).astype(int))
    
    
    def findDisrupted(row):
        pos = row['cDNA_position']
        local_df_pred = df_pred[df_pred['mRNA'] == row['Feature']]
        if local_df_pred.empty:
            return None

        if '-' in pos:
            pos1, pos2 = map(int, pos.split('-', 1))
        else:
            pos1 = pos2 = int(pos)

        hits = pd.concat([
            local_df_pred[
                (local_df_pred['bs_start'] <= pos_iter) &
                (local_df_pred['bs_end'] >= pos_iter)
            ]
            for pos_iter in range(pos1, pos2 + 1)
        ], ignore_index=True).drop_duplicates()
        
        if hits.empty:
            return None

        return ','.join(sorted(hits['miRNA'].unique()))
    
    
    log.info('Editing Dataframe: Adding miRNAs that are disrupted to respective IDs')
    df.loc[:, 'Disrupted_miRNA'] = df.apply(findDisrupted, axis=1)
    
    log.info('Exporting output file')
    df.to_csv(f'output/{OUTPUT_FILE}', sep='\t', index=False)
    
    log.info(f'Finished, output written to file {OUTPUT_FILE}')

