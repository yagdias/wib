#!/usr/bin/env python3

import argparse
import pandas as pd
import os
from ete3 import NCBITaxa


# Definindo os argumentos da linha de comando
parser = argparse.ArgumentParser(description="send blast files path to top5 scripts")
parser.add_argument("--blast_file", type=str, help="blastx path")
parser.add_argument("--sample_name", type=str, help="sample name")
parser.add_argument("--pident", type=str, help="minimal pident to filter sequences")
parser.add_argument("--cut_evalue", type=str, help="minimal evalue to filter sequences")
parser.add_argument("--taxonomy_db_path", type=str, help="path to the NCBI taxonomy database file")
parser.add_argument("--update_db", type=str, help="whether to update the NCBI taxonomy database")

args = parser.parse_args()

# Atribuindo os valores dos argumentos a variáveis
blast_file = args.blast_file
sample_name = args.sample_name
pident = float(args.pident)
cut_evalue = float(args.cut_evalue)
taxonomy_db_path = args.taxonomy_db_path
update_db = args.update_db.lower() == 'true'



# 1 - creating, filtering and sorting dataframe from blast. File name is based on tools. if blast, file name starts with 'b'. If diamond, file name starts with 'd'

def extract_access(x):
    if "|" in x:
        return x.split("|")[1]
    return x

def create_df_from_blast(blast_file, pident, cut_evalue):
    mode = blast_file.split(".")[-2]
    columns = ['qseqid', 'qlen', 'sseqid', 'slen', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'staxids', 'sscinames', 'sskingdoms', 'stitle', 'salltitles', 'nident', 'qcovhsp']
    df = pd.read_csv(blast_file, sep='\t', names=columns)
    df = df[(df['pident'] >= pident) & (df['evalue'] < cut_evalue)].sort_values(['bitscore', 'evalue'], ascending=False)
    df["sseqid"] = df["sseqid"].apply(extract_access).astype(str)
    df["stitle"] = df["sseqid"] + " " + df["stitle"]
    if not mode.startswith("b"):
        df = df.rename(columns={"qcovhsp": "qcovs"})
        df['qcovhsp'] = (df['length'] / df['qlen']) * 100
        df['qcovhsp'] = df['qcovhsp'].astype('int')
    df['taxid'] = df['staxids']
    return df

# 2

def create_top(df):
    top5_list = []
    contig_names = df['qseqid'].unique()
    for contig_name in contig_names:
        filtered_group = df[df['qseqid'] == contig_name].nlargest(1, 'bitscore')
        top5_list.append(filtered_group)
    top5_df = pd.concat(top5_list, ignore_index=True)
    return top5_df

# 3

def get_taxonomy(top5_df, taxonomy_db_path, update_db=True):
    print(f'valor do update_db: {update_db}')
    if not os.path.exists(taxonomy_db_path):
        os.makedirs(taxonomy_db_path, exist_ok=True)
    db_file_path = os.path.join(taxonomy_db_path, 'taxonomy.sqlite')
    ncbi = NCBITaxa(dbfile=db_file_path)
    if update_db:
        ncbi.update_taxonomy_database()

# load top5df

    def get_complete_lineage(taxid):
        try:
            # Obter linhagem e ranks
            lineage = ncbi.get_lineage(taxid)
            lineage_names = ncbi.get_taxid_translator(lineage)
            ranks = ncbi.get_rank(lineage)

            # IDs de superkingdoms esperados
            superkingdom_ids = {
                2: "Bacteria",
                2157: "Archaea",
                2759: "Eukaryota",
                10239: "Viruses"
            }

            # Inicializar dicionário com taxid
            complete_lineage = {'taxid': taxid}

            # Identificar o superkingdom na linhagem
            detected_superkingdom = None
            for tax in lineage:
                if tax in superkingdom_ids:
                    detected_superkingdom = superkingdom_ids[tax]
                    break

            # Atribuir o superkingdom detectado ou Unknown se não encontrado
            complete_lineage['superkingdom'] = detected_superkingdom if detected_superkingdom else "Unknown"

            # Mapear outros ranks relevantes
            rank_mapping = {
                'clade': 'clade',
                'kingdom': 'kingdom',
                'phylum': 'phylum',
                'class': 'class',
                'order': 'order',
                'family': 'family',
                'genus': 'genus',
                'species': 'species'
            }

            for tax in lineage:
                rank = ranks.get(tax, "no_rank")
                if rank in rank_mapping:
                    mapped_rank = rank_mapping[rank]
                    complete_lineage[mapped_rank] = lineage_names.get(tax, "Unknown")

            if 'clade' not in complete_lineage:
                complete_lineage['clade'] = "Unknown"

            print(complete_lineage)
            return complete_lineage
        except Exception as e:
            print(f"Taxid error! {e}")
            return {'taxid': taxid, 'superkingdom': "Unknown"}

    taxon_data = [get_complete_lineage(taxid) for taxid in top5_df['taxid'].unique()]
    taxonomy_df = pd.DataFrame(taxon_data)
    taxonomy_df.fillna("Unknown", inplace=True)
    top5_taxonomy_df = pd.merge(top5_df, taxonomy_df, on='taxid')
    print(top5_taxonomy_df)
    print(top5_taxonomy_df.columns)
    return top5_taxonomy_df

# 5 creating best hits and others best hits

def create_report(top5_taxid_df):

    queries = top5_taxid_df['qseqid'].unique()
    
    # DataFrames para armazenar resultados
    df_best = []
    df_others = []
    accession_desc_list = []

    # Processa cada query
    for query in queries:
        query_df = top5_taxid_df[top5_taxid_df['qseqid'] == query]
        
        # Melhor hit (apenas 1)
        best_hit = query_df.nlargest(1, 'bitscore')
        df_best.append(best_hit)
        
        # Outros hits (até 4 para completar o top 5)
        others = query_df.nlargest(5, 'bitscore').iloc[1:]
        df_others.append(others)
        
        if not others.empty:
            others['qcovs'] = others['qcovs'].astype(str)
            others['qcovhsp'] = others['qcovhsp'].astype(str)
            
            # Formata as informações dos hits adicionais
            formatted_info = "\n".join(
                ["$".join(map(str, row)) for row in others[['sseqid', 'sscinames', 'qcovs', 'qcovhsp']].values]
            )
            accession_desc_list.append([query, formatted_info])
    
    # Concatena resultados em DataFrames
    df_best = pd.concat(df_best).sort_values('bitscore', ascending=False)
    df_others = pd.concat(df_others).fillna('undefined')
    df_others_selected = pd.DataFrame(accession_desc_list, columns=['qseqid', 'other_best_hits_to_consider']).set_index('qseqid')
    
    # Renomeia e ajusta colunas
    df_best = df_best.rename(columns={'sscinames': 'best_hit'}).set_index('qseqid')
    df_top5_taxid_best_hit = top5_taxid_df.set_index('qseqid').join(df_best[['best_hit']], how='inner').reset_index()
    df_top5_taxid_best_hit_others = df_top5_taxid_best_hit.set_index('qseqid').join(df_others_selected, how='outer').reset_index()
    
    # Ajusta e formata valores
    for col in ['pident', 'qcovs', 'qcovhsp']:
        if col in df_top5_taxid_best_hit_others:
            df_top5_taxid_best_hit_others[col] = df_top5_taxid_best_hit_others[col].astype(float).map("{:.2f} %".format)
    
    return df_top5_taxid_best_hit_others.fillna('-')

##

def export_report_to_csv(final_df, file_name):
    final_df.to_csv(f"{file_name}.csv",  index=False)

def export_krona_tsv(final_df, file_name):
    final_df = final_df[final_df['superkingdom'] == "Viruses"]
    krona_df = final_df[['qseqid', 'taxid']].to_csv(f"{file_name}_krona.tsv", sep='\t', index=False, header=False)






dataframe = create_df_from_blast(blast_file, pident, cut_evalue)
top1 = create_top(dataframe)
fulltax = get_taxonomy(top1, taxonomy_db_path, update_db)
fulltax.to_csv(f"{sample_name}.csv", index=False)