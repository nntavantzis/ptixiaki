import pandas as pd
import requests
import re
import logging
from os import path
from io import StringIO

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('run.log'),
        logging.StreamHandler()
    ]
)
log = logging.getLogger(__name__)

VEP_OUTPUT = 'output.tsv'
TARPMIR_LOCATION = './TarPmiR_Linux'
TRANSCRIPT_NUM = 10 # LIMIT TRANSCRIPTS FOR TESTING


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

def runTarPmiR():
    pass # FIXME: PREPEI NA FTIAKSO TO FUNCTION NA DOULEVEI KALA


###################### MAIN ######################

def main():
    log.info('Starting program...')
    log.info(f'Loading file to DataFrame: {VEP_OUTPUT}')
    df = pd.read_csv(
        VEP_OUTPUT,
        sep='\t',
        comment='#',
        header=None
    )
    log.info(f'File \'{VEP_OUTPUT}\' loaded')
    
    # Find the header and add to DataFrame
    with open(VEP_OUTPUT) as file:
        for line in file.readlines():
            if line.startswith('#') and not line.startswith('##'):
                header = line.lstrip('#').strip().split('\t')
                df.columns = header
                break
    df = df[df['Consequence'].str.contains('3_prime_UTR_variant')]

    # Temporary DataFrame subset to use with less transcripts
    df2 = df.head(TRANSCRIPT_NUM).copy()
    
    # DOWNLOAD cDNA SEQUENCE FASTA FOR USE WITH TARPMIR
    transcriptsToDownload = df2['Feature'].drop_duplicates().to_list()
    

    # NOTE: ta parakato einai an kano mono ena tarpmir kai oxi batch
    res = getTranscriptsFromIDs(transcriptsToDownload)
    with open(path.join(TARPMIR_LOCATION, 'cDNA.fasta'), 'w') as file:
        file.write(res)
        
    # FIXME: EDO KANONIKA GINETAI TARPMIR
        
    log.info('Adding to DataFrame: IDs with found FASTA sequences')
    with open('./TarPmiR_Linux/cDNA.fasta', 'r', encoding='utf-8') as file: # FIXME: otan valo to tarpmir kai ta variable file names na allakso to file name edo
        fasta_ids = [line.strip()[1:] for line in file if line.startswith('>')]
    df2.loc[:, 'Sequence_exists'] = df2['Feature'].isin(fasta_ids)



    # FIXME: tha prepei na to ftiakso gia na doulevei me to batch -> na kano merge ola ta results se ena arxeio
    log.info('Loading file to DataFrame: TarPmiR Prediction Output')
    with open('./TarPmiR_Linux/mirnaverified.fasta_cDNA.fasta.bp', 'r', encoding='utf-8') as file: # FIXME: otan valo to tarpmir kai ta variable file names na allakso to file name edo
        content = file.read()
        content = re.sub(r'(?<=.)hsa', r'\nhsa', content)
    df_pred = pd.read_csv(StringIO(content), sep='\t')
    log.info('File loaded')
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
    
    log.info('Adding to DataFrame: Disrupted miRNAs')
    df2.loc[:, 'Disrupted_miRNA'] = df2.apply(findDisrupted, axis=1)
    
    log.info('Exporting modified output file')
    df2.to_csv('./output_disrupted.tsv', sep='\t', index=False) # FIXME: na valo to filename me vasi to variable filename
    
    log.info('Finished')
    
try:
    main()
except Exception as e:
    log.error(e)
