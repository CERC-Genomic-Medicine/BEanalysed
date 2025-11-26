#!/usr/bin/env python3


import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from sklearn.utils import resample
import random
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import openpyxl
import matplotlib.lines as mlines
from matplotlib.patches import Patch
import matplotlib.patches as patches
import matplotlib.colors as mcolors
import argparse
from scipy.stats import rankdata
from scipy.stats import binomtest
import sys
import os
import re
import warnings
from itertools import pairwise  # Python 3.10+
from itertools import combinations, product
import itertools
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import percentileofscore
from openpyxl.drawing.image import Image
'''
AUTHOR: Vincent Chapdelaine <vincent.chapdelaine@mcgill.ca>
VERSION: 1.1
YEAR: 2024

Goal :

Create Base Summary

'''

parser = argparse.ArgumentParser(description='This script create a base concatenated summary for analysis performed with MaGeCK ')
parser.add_argument('-k',"--key", metavar='FILE', required=True, dest='keyFile', help="tab delimited key file, with columns ID_guide and Sequence (sequence is the protospacer)")
parser.add_argument('-I',"--input", metavar='FILE', required=True, dest='mageck_files', nargs='+', help="MageCK files per genes")
parser.add_argument('-a', '--annotation_excel', metavar='FILE', dest='vep_excel_file',  nargs='+', required=True, type=str, help='Excel file with predicted consequences')
parser.add_argument('-s', '--sheet', metavar='str', dest='vep_sheet_file', required=True, type=str, help='name of the sheet to be used in the previously stated excel file')
parser.add_argument('-l', '--lib_sheet', metavar='str', dest='vep_lib_sheet', required=False, type=str, default='Library', help='name of the sheet to be used in the previously stated excel file')
parser.add_argument('--Isoform', metavar='str', dest='Isoforms', required=False, type=str, nargs='+', help='Isoforms to select (optional)')
parser.add_argument('--Pick', dest='PickH', action='store_true', help="Select annotation based on VEP\'s PickH (column PICK) (optional)")
parser.add_argument('--e', dest='Control_empty', action='store_true', help="Mark sgRNA with no predicted mutation as Negative Controls")
parser.add_argument('-n', '--negative_Control', metavar='str', dest='Control_N', required=False, nargs='+', type=str, help='list of negative control protein/region (if no negative controls (empty widows/regions) Biological significance ignored)')
parser.add_argument('-p', '--positive_Control', metavar='str', dest='Control_P', required=False, nargs='+', type=str, help='list of positive control protein/region ')
parser.add_argument('-X',"--xvar", dest='variable_names', required=True, help="Comma-separated experimental conditions (e.g., UNT/TREAT,KO/WT) should be reflected in the name of the mageck file * No overlap * eg. treat,treatment")
parser.add_argument('--out', dest='output', default='summary', type=str, help='Name Output File')




variant_consequences_mapping = {
'missense_variant': 'missense',
'intron_variant': 'non_coding',
'downstream_gene_variant': 'non_coding',
'NMD_transcript_variant': 'non-sense',
'upstream_gene_variant': 'non_coding',
'3_prime_UTR_variant': 'non_coding',
'synonymous_variant': 'synonymous',
'non_coding_transcript_exon_variant': 'non_coding',
'splice_region_variant': 'non_coding',
'splice_polypyrimidine_tract_variant': 'non_coding',
'stop_gained': 'non-sense',
'coding_sequence_variant': 'coding', # Should never be alone
'5_prime_UTR_variant': 'non_coding',
'regulatory_region_variant': 'regulatory',
'splice_donor_variant': 'splice',
'splice_acceptor_variant': 'splice',
'non_coding_transcript_variant': 'non_coding',
'splice_donor_region_variant': 'non_coding',
'splice_donor_5th_base_variant': 'non_coding',
'TF_binding_site_variant': 'regulatory',
'start_lost': 'non-sense',
'incomplete_terminal_codon_variant': 'non-sense',
'NA/NotFound' :  'none',
'None' : 'none',
'Empty' : 'empty',
}



def convert_to_int(s):
    try:
        return int(s)
    except ValueError:
        return None


def parse_treatment_levels(xvar):
    return [sorted(group.split('/'), key=len, reverse=True) for group in xvar.split(',')]

def extract_labels_from_filename(filename, levels):
    fname = os.path.basename(filename)
    labels = []
    for group in product(*levels):
        hits = [lvl for lvl in group if re.search(rf'{re.escape(lvl)}', fname)]
        sized = len(group)
        if len(hits) !=sized:
            continue
        else :
            labels.extend(hits)
    return labels

def Transform_MaGeCK(filename, value=None):
    Mage = pd.read_csv(filename,sep='\t',header=0)
    if value :
        returned = pd.DataFrame({
        'id': Mage.id,
        value: [None] * len(Mage.id)
        })
        returned[value] = Mage[value]
    else :
        returned = pd.DataFrame({
        'id': Mage.id,
        'lfc': [None] * len(Mage.id),
        'p_value': [None] * len(Mage.id),
        'fdr': [None] * len(Mage.id)
        })
        returned['pos|rank']=Mage['pos|rank']
        returned['neg|rank']=Mage['neg|rank']
        returned['pos|p-value']=Mage['pos|p-value']
        returned['neg|p-value']=Mage['neg|p-value']
        returned['pos|fdr']=Mage['pos|fdr']
        returned['neg|fdr']=Mage['neg|fdr']
        returned['lfc'] = Mage['pos|lfc']
        returned['p_value'] = [row['neg|p-value'] if row['pos|lfc']<0 else row['pos|p-value'] for _, row in Mage.iterrows()]
        returned['fdr'] = [row['neg|fdr'] if row['pos|lfc']<0 else row['pos|fdr'] for _, row in Mage.iterrows()]
    return returned


def parse_VEP_excel(excel_files, sheet_name, variant_consequences_mapping, dict_IDs, isoforms_selection=None):
    df = pd.concat([pd.read_excel(excel_file, sheet_name=sheet_name, dtype=str).fillna('') for excel_file in excel_files], axis=0)
    replicate_col=sheet_name.split(' ')[2]
    #only Relevant lines
    df=df.loc[[i in set(dict_IDs) for i in df['#Uploaded_variation']],] # Keeping relevant


    ## Verfication
    df['#Uploaded_variation'] = [dict_IDs[ids] for ids in df['#Uploaded_variation']]
    nb_guides_annotated = len(set(df['#Uploaded_variation']))
    already_filtered = nb_guides_annotated == len(df['#Uploaded_variation'])

    ### Isoform Selection
    if (already_filtered and isoforms_selection):
        print(' Warning : The excel annotations seems to already be filtered')
    if isoforms_selection == 'PickH':
        df = df.loc[[i==1 for i in df['PICK']],]
    if isinstance(isoforms_selection, list):
        if set(isoforms_selection) - set(df['Feature']) :
            raise ValueError(f"isoforms specified are not found in the annotations: {set(isoforms_selection) - set(df['Feature'])}")
        else :
            df = df.loc[[j in set(isoforms_selection) for j in df['Feature']],]

    if nb_guides_annotated < len(df['#Uploaded_variation']) : 
        if isoforms_selection :
            raise ValueError(f"Annotation remains ambiguous after isoform selection (This should not happen !) Most likely multiple isoforms were specified for the same genes")
        else :
            raise ValueError(f"Annotations are ambiguous, consider filtering")
    if nb_guides_annotated > len(df['#Uploaded_variation']) :     
        if isoforms_selection :
            raise ValueError(f"Some guides annotation were lost, most likely not some isoforms are missing from the specified isoforms")
        else :
            raise ValueError(f"Annotations are ambiguous, consider filtering")

    for ID, loc, cons, prot_pos, impact, aa, replicate in zip(
        df['ID'],
        df['Location'],
        df['Consequence'],
        df['Protein_position'],
        df['IMPACT'],
        df['Amino_acids'],
        df[replicate_col]
    ):
        mutation=[]
        mutations_AA=""
        ### AminoAccid position
        if '?' in prot_pos or '/' not in aa:
            part = prot_pos.split('-')[0]
            part = part if part != '?' else prot_pos.split('-')[-1]
            POS_AA = convert_to_int(part)
        elif  '/' in aa:
            AA_mutation = aa.split('/')
            if '-' in prot_pos :
                position_initial = int(prot_pos.split('-')[0])
                position_final = int(prot_pos.split('-')[1])
            else :
                position_initial = int(prot_pos)
                position_final = int(prot_pos)
            arr_a = np.frombuffer(AA_mutation[0].encode(), dtype='S1')
            arr_b = np.frombuffer(AA_mutation[1].encode(), dtype='S1')
            diff_idx = np.where(arr_a != arr_b)[0].tolist()
            if diff_idx :
                POS = [iteration + position_initial for iteration in diff_idx]
                mutation = ["_".join([AA_mutation[0][iteration],str(iteration + position_initial)]) for iteration in diff_idx]
                mutations_AA = ",".join(mutation)
                start = min(diff_idx)
                end = max(diff_idx)
                POS_AA = np.round(position_initial + ((start + end) / 2)).astype(int)
            else :
                POS_AA = np.round((position_initial+position_final)/2).astype(int)
        # EStablish consequence
        mapped = [variant_consequences_mapping.get(c, None) for c in cons.split(',')]
        if 'splice' in mapped:
            consequence = 'splice'
            consequence_detail = 'splice'
        elif 'non-sense' in mapped:
            consequence = 'non-sense'
            consequence_detail = 'non-sense'
        elif 'missense' in mapped:
            consequence = 'missense'
            consequence_detail = 'missense'
        elif 'synonymous' in mapped:
            consequence = 'synonymous'
            consequence_detail = 'synonymous'
        elif 'regulatory' in mapped:
            consequence_detail = 'regulatory'
            consequence = None
        elif 'non_coding' in mapped:
            consequence_detail = 'non-coding'
            consequence = None
        elif 'empty' in mapped:
            consequence_detail = 'No predicted Mutation'
            consequence = None
        else :
            consequence_detail = 'N/A'
            consequence = None
        yield ID, consequence, consequence_detail, impact, POS_AA, mutations_AA, replicate


def main():
    args = parser.parse_args()

    #### Read keys of the Assay
    key_file = pd.read_csv(args.keyFile, sep='\t', usecols=['ID_guide', 'Sequence'])
    keys_assay = dict(zip(key_file['Sequence'],key_file['ID_guide']))
    ### Read Annoation
    df = pd.concat([pd.read_excel(excel_file, sheet_name=args.vep_lib_sheet, dtype=str).fillna('') for excel_file in args.vep_excel_file], axis=0)
    ### Verifications
    if len(set(key_file['Sequence']) - set(df['protospacer']))>0:
        raise ValueError('According to keyfile and excel files not all sgRNA guides were annotated')
    df=df.loc[[i in set(key_file['Sequence']) for i in df['protospacer']],] # Keeping relevant

    ### Convert to Key file IDs
    df['old'] = df['ID']
    df['ID'] = [keys_assay[i] for i in df.protospacer]
    dict_ID_ID = dict(zip(df.old, df.ID))
    dict_ID_sgRNA = dict(zip(df.ID, df.protospacer))


    ## Main Annotation
    if args.PickH and args.Isoforms :
        raise ValueError('Annotation selection can not be both pickH and Isoforms')
    if args.PickH :
        isoforms_selection = 'PickH'
    elif args.Isoforms :
        isoforms_selection = args.Isoforms
    else :
        isoforms_selection = None
    VEP= pd.DataFrame(parse_VEP_excel(args.vep_excel_file,args.vep_sheet_file, variant_consequences_mapping, dict_ID_ID , isoforms_selection))
    VEP.columns=['ID', 'consequence', 'consequence_detail', 'impact', 'POS_AA', 'mutations_AA', 'replicate']

    Variant_effect=dict(zip(VEP['ID'],VEP['consequence']))
    Variant_effect_full=dict(zip(VEP['ID'],VEP['consequence_detail']))
    Variant_impact=dict(zip(VEP['ID'],VEP['impact']))
    Variant_Amino_position=dict(zip(VEP['ID'],VEP['POS_AA']))
    Variant_Amino_mutation_list=dict(zip(VEP['ID'],VEP['mutations_AA']))
    replicate_dic=dict(zip(VEP['ID'],VEP['replicate']))

    ### ID empties ###
    empties = df.loc[~df['ID'].isin(set(VEP['ID'])), 'ID']
    Variant_effect.update(dict(zip(empties,itertools.repeat('No predicted Mutation'))))
    Variant_effect_full.update(dict(zip(empties,itertools.repeat('No predicted Mutation'))))
    Variant_impact.update(dict(zip(empties,itertools.repeat(''))))
    Variant_Amino_position.update(zip(empties,itertools.repeat('')))
    Variant_Amino_mutation_list.update(dict(zip(empties,itertools.repeat(''))))
    replicate_dic.update(dict(zip(empties,itertools.repeat(''))))
    print('running')
    with pd.ExcelWriter(f'{args.output}.xlsx') as writer:
        for file in args.mageck_files:
            print(file)
            try :
                data_full = Transform_MaGeCK(file)
            except: 
                raise ValueError('Main input data (-i) does not have a supported format')
            if set(data_full['id']) != set(dict_ID_sgRNA.keys()):
                raise ValueError('Key file and mageck file do not completely overlap')
            data_full.insert(1, 'sgRNA_seq', data_full['id'].map(dict_ID_sgRNA))
            label = "_".join(extract_labels_from_filename(file, parse_treatment_levels(args.variable_names)))
            data_full['proteins']=[i.split('_')[0] for i in data_full['id']]
            data_full['Consequence'] = [Variant_effect[i] for i in data_full['id']]
            data_full['Consequence_Detail'] = [Variant_effect_full[i] for i in data_full['id']]
            data_full['Position']=[Variant_Amino_position[i] for i in data_full['id']]  
            data_full['Mutations']=[Variant_Amino_mutation_list[i] for i in data_full['id']]  
            data_full['replicates_window']=[replicate_dic[i] for i in data_full['id']]  
            data_full['Controls'] = None
            if args.Control_P:
                positive_controls_true = [i in args.Control_P for i in data_full['proteins']]
                data_full.loc[np.array(positive_controls_true),'Controls'] = 'positive_control' 
                data_full.loc[np.array(positive_controls_true),'Consequence'] = f'Positive Control Gene{"s" if len(args.Control_P)!=1 else ""}'
                data_full.loc[np.array(positive_controls_true),'Consequence_Detail'] = f'Positive Control Gene{"s" if len(args.Control_P)!=1 else ""}'
            if args.Control_N:
                negative_controls_true = [i in args.Control_N for i in data_full['proteins']]
                data_full.loc[np.array(negative_controls_true),'Controls'] = 'negative_control' 
                data_full.loc[np.array(negative_controls_true),'Consequence'] = f'Negative Control Gene{"s" if len(args.Control_N)!=1 else ""}'
                data_full.loc[np.array(negative_controls_true),'Consequence_Detail'] = f'Negative Control Gene{"s" if len(args.Control_N)!=1 else ""}'
            if args.Control_empty :
                data_full.loc[[j in empties for j in data_full['id']],'Controls'] = 'negative_control' 
            data_full.to_excel(writer, sheet_name=label, index=False)

if __name__ == "__main__":
    main()
