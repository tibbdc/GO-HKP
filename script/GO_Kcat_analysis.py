import json
import os
import re
import time
import numpy as np
import pandas as pd
import requests

from cobra import Reaction
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
from cobra.util.solver import set_objective


def load_json(file_path):
    with open(file_path, 'r') as f:
        return json.load(f)

def load_csv(file_path, index_col):
    try:
        return pd.read_csv(file_path, index_col=index_col)
    except FileNotFoundError:
        print(f"Error: File not found at {file_path}")
        return None
    
def calculate_range_mean(value):
    if '-' in value:
        start, end = map(float, value.split('-'))
        return (start + end) / 2
    return float(value)

def calculate_mean(lst):
    return sum(lst) / len(lst) if lst else float('nan')


def convert_to_irreversible(model):
    reactions_to_add = []
    coefficients = {}
    for reaction in model.reactions:
        if reaction.lower_bound < 0 and reaction.upper_bound == 0:
            for metabolite in reaction.metabolites:
                coeff = reaction.get_coefficient(metabolite)
                reaction.add_metabolites({metabolite: -2 * coeff})
            reaction.id += "_reverse"
            reaction.upper_bound = -reaction.lower_bound
            reaction.lower_bound = 0
        if reaction.lower_bound < 0 and reaction.upper_bound > 0:
            reverse_reaction = Reaction(reaction.id + "_reverse")
            reverse_reaction.lower_bound = max(0, -reaction.upper_bound)
            reverse_reaction.upper_bound = -reaction.lower_bound
            coefficients[reverse_reaction] = -reaction.objective_coefficient
            reaction.lower_bound = max(0, reaction.lower_bound)
            reaction.upper_bound = max(0, reaction.upper_bound)
            reverse_reaction.notes["reflection"] = reaction.id
            reaction.notes["reflection"] = reverse_reaction.id
            reverse_reaction.add_metabolites({k: -v for k, v in reaction.metabolites.items()})
            reverse_reaction.gene_reaction_rule = reaction.gene_reaction_rule
            reverse_reaction.subsystem = reaction.subsystem
            reverse_reaction.annotation = reaction.annotation
            reactions_to_add.append(reverse_reaction)
    model.add_reactions(reactions_to_add)
    set_objective(model, coefficients, additive=True)

def isoenzyme_split(model):
    for r in model.reactions:
        if " or " in r.gene_reaction_rule:
            base = r.copy()
            parts = r.gene_reaction_rule.split(" or ")
            for i, part in enumerate(parts):
                if i == 0:
                    r.id += f"_num1"
                    r.gene_reaction_rule = part
                else:
                    r_add = base.copy()
                    r_add.id = base.id + f"_num{i+1}"
                    r_add.gene_reaction_rule = part
                    model.add_reactions([r_add])
    for r in model.reactions:
        r.gene_reaction_rule = r.gene_reaction_rule.strip("( )")
    return model


def get_go_term_mean_kcat_by_org_v3_1(norm_model, ns2assc_org, org_type, org_name, org_id, GO_kcat_file, go):
    go_term_mean_kcat_dict = {}

    def get_mean_kcat_list(gene, go_type):
        mean_kcat_list = []
        try:
            ncbigene = int(norm_model.genes.get_by_id(gene).annotation['ncbigene'])
            go_set = ns2assc_org[go_type].get(ncbigene, set())
            print(go_set)
            for eachgo in go_set:
                go_term_mean_kcat = get_kcat_from_GO_1(org_id, org_name, org_type, GO_kcat_file, go, eachgo)
                if go_term_mean_kcat != '-':
                    mean_kcat_list.append(go_term_mean_kcat)
        except Exception as e:
            print(f"Error occurred while processing gene: {gene}. Error: {str(e)}")
        return mean_kcat_list

    def process_reaction(rea):
        gene_list = re.findall(r'\w+', rea.gene_reaction_rule)
        go_term_mean_kcat_dict[rea.id] = {
            'BP': {gene: get_mean_kcat_list(gene, 'BP') for gene in gene_list},
            'MF': {gene: get_mean_kcat_list(gene, 'MF') for gene in gene_list},
            'CC': {gene: get_mean_kcat_list(gene, 'CC') for gene in gene_list}
        }
        go_term_mean_kcat_dict[rea.id]['Total'] = {
            gene: go_term_mean_kcat_dict[rea.id]['BP'][gene] +
                  go_term_mean_kcat_dict[rea.id]['MF'][gene] +
                  go_term_mean_kcat_dict[rea.id]['CC'][gene]
            for gene in go_term_mean_kcat_dict[rea.id]['BP']
        }

    with ThreadPoolExecutor() as executor:
        executor.map(process_reaction, norm_model.reactions)

    return go_term_mean_kcat_dict

def get_kcat_from_GO_1(org_id, org_name, use_type, GO_kcat_file, go, go_term_id):
    GO_term_kcat_df = pd.read_csv(GO_kcat_file)

    if use_type == 'org_name':
        df = GO_term_kcat_df[GO_term_kcat_df['Organism'].str.contains(org_name)]
    elif use_type == 'org_id':
        GO_term_kcat_df['Organism_ID'] = GO_term_kcat_df['Organism'].str.extract(r'\[(\d+)\]')
        df = GO_term_kcat_df[GO_term_kcat_df['Organism_ID'] == org_id]
    else:
        df = GO_term_kcat_df

    df = df[df['Kcat'].notnull()]

    if len(df[df['GO_term'] == go_term_id]) > 0:
        kcat_values = df[df['GO_term'] == go_term_id]['Kcat'].values
        kcat_list = [calculate_range_mean(value) for eachlist in kcat_values for value in eachlist.split(';')]
        return calculate_mean(kcat_list)
    elif go_term_id in go:
        child_terms = go[go_term_id].get_all_children()
        kcats = []
        for child in child_terms:
            if len(df[df['GO_term'] == child]) > 0:
                child_vals = df[df['GO_term'] == child]['Kcat'].values
                kcats += [calculate_range_mean(val) for each in child_vals for val in each.split(';')]
        return calculate_mean(kcats) if kcats else '-'
    else:
        print(f"GO Term with ID {go_term_id} is obsolete or not found.")
        return '-'


_go_kcat_cache = {'BP': {}, 'MF': {}, 'CC': {}}

def get_kcat_from_GO_2(org_name, use_type, GO_kcat_file, go, go_term_id, go_type):
    if go_term_id in _go_kcat_cache[go_type]:
        return _go_kcat_cache[go_type][go_term_id]

    df = pd.read_csv(GO_kcat_file)
    if use_type == 'org_name':
        df = df[df['Organism'].str.contains(org_name)]
    df = df[df['Kcat'].notnull()]

    matched = df[df['GO_term'] == go_term_id]
    if not matched.empty:
        kcats = [calculate_range_mean(k) for v in matched['Kcat'] for k in v.split(';')]
        result = calculate_mean(kcats)
        _go_kcat_cache[go_type][go_term_id] = result
        return result

    if go_term_id in go:
        children = go[go_term_id].get_all_children()
        kcats = []
        for ct in children:
            matched = df[df['GO_term'] == ct]
            for v in matched['Kcat']:
                kcats += [calculate_range_mean(k) for k in v.split(';')]
        if kcats:
            result = calculate_mean(kcats)
            _go_kcat_cache[go_type][go_term_id] = result
            return result

    _go_kcat_cache[go_type][go_term_id] = '-'
    return '-'

def get_go_term_mean_kcat_by_org_v3_2(norm_model, result, org_type, org_name, GO_kcat_file, go):
    def get_mean_kcat_list(gene, go_type):
        try:
            uniprot_ids = norm_model.genes.get_by_id(gene).annotation['uniprot']
            go_terms = result.get(go_type, {}).get(uniprot_ids, [])
            return [k for g in go_terms if (k := get_kcat_from_GO_2(org_name, org_type, GO_kcat_file, go, g, go_type)) != '-']
        except Exception:
            return []

    def process_reaction(rea):
        genes = re.findall(r'\w+', rea.gene_reaction_rule)
        return rea.id, {
            go_type: {gene: get_mean_kcat_list(gene, go_type) for gene in genes}
            for go_type in ['BP', 'MF', 'CC']
        }

    result_dict = {}
    with ThreadPoolExecutor(max_workers=16) as ex:
        futures = {ex.submit(process_reaction, r): r.id for r in norm_model.reactions}
        for fut in tqdm(as_completed(futures), total=len(futures)):
            rid, kcat_data = fut.result()
            kcat_data['Total'] = {
                g: kcat_data['BP'].get(g, []) + kcat_data['MF'].get(g, []) + kcat_data['CC'].get(g, [])
                for g in set(kcat_data['BP']) | set(kcat_data['MF']) | set(kcat_data['CC'])
            }
            result_dict[rid] = kcat_data

    return result_dict


def process_data(go_term_mean_kcat_dict, go_type, clc_type):
    result = {}
    for rid, entry in go_term_mean_kcat_dict.items():
        vals = [v for vs in entry[go_type].values() for v in vs if v != '-']
        if vals:
            result[rid] = {
                'max': np.max(vals),
                'mean': np.mean(vals),
                'median': np.median(vals)
            }[clc_type]
    return result

def extract_kcat_GO_terms_from_protein_info(protein_info, GO_kcat_outfile):
    GO_pattern = r'GO:\d+'

    if not os.path.isfile(GO_kcat_outfile) or os.path.getsize(GO_kcat_outfile) == 0:
        with open(GO_kcat_outfile, 'w') as file:
            file.write("GO_term,Organism,Kcat\n")

    with open(GO_kcat_outfile, 'a') as file:
        for pkb_id, protein in protein_info.items():
            kcat_values = []
            for substrate, substrate_kcat_info in protein['Properties']['Kcat'].items():
                for value in substrate_kcat_info:
                    kcat_values.append(value['value'])
            if kcat_values:
                for entry, entry_info in protein['Consist'].items():
                    if 'Gene Ontology (GO)' in entry_info:
                        gene_terms = re.findall(GO_pattern, entry_info['Gene Ontology (GO)'])
                        for term in gene_terms:
                            org = f"[{protein['Organism (ID)']}]{protein['Organism']}"
                            file.write(f"{term},{org},{';'.join(kcat_values)}\n")

def get_uniprot_sequence_by_gene(gene_name):
    url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene_name}&format=fasta"
    try:
        response = requests.get(url, headers={"Accept": "text/plain"}, timeout=30)
        if response.ok:
            return response.text.strip()
    except Exception as e:
        print(f"Error fetching sequence for {gene_name}: {str(e)}")
    return None
