"""
"""


from seq_script import extract_seq, convert_pam
from .pages_utils import ANNOTATIONS_DIR, JOBID_ITERATIONS_MAX, JOBID_MAXLEN, RESULTS_DIR, VALID_CHARS
from app import (
    URL, 
    app, 
    operators,
    current_working_directory, 
    app_main_directory, 
    DISPLAY_OFFLINE,
    ONLINE,
    exeggutor
)

from dash.exceptions import PreventUpdate
from dash.dependencies import Input, Output, State
from typing import List, Tuple
from datetime import datetime

import dash_bootstrap_components as dbc
import dash_html_components as html
import dash_core_components as dcc

import collections
import subprocess
import filecmp
import random
import string
import glob
import os


# Allowed mismatches and bulges
# ONLINE = False  # NOTE change to True for online version, False for offline
if ONLINE:
    MAX_MMS = 7  # max allowed mismatches
    MAX_BULGES = 3  # max allowed bulges
else:
    # NOTE modify value for increasing/decreasing max mms or bulges available on 
    # Dropdown selection
    MAX_MMS = 7  # max allowed mismatches
    MAX_BULGES = 3  # max allowed bulges
# mismatches, bulges and guides values
AV_MISMATCHES = [{"label": i, "value": i} for i in range(MAX_MMS)]
AV_BULGES = [{"label": i, "value": i} for i in range(MAX_BULGES)]
AV_GUIDE_SEQUENCE = [{"label": i, "value": i} for i in range(15, 26)]


def split_filter_part(filter_part: str) -> Tuple:
    """Recover filtering operator.

    ...

    Paramters
    --------
    filter_part : str
        Filter
    
    Returns
    -------
    Tuple
    """

    if not isinstance(filter_part, str):
        raise TypeError(f"Expected {str.__name__}, got {type(filter_part).__name__}")
    for operator_type in operators:
        for operator in operator_type:
            if operator in filter_part:
                name_part, value_part = filter_part.split(operator, 1)
                name = name_part[(name_part.find("{") + 1):name_part.rfind("}")]
                value_part = value_part.strip()
                v0 = value_part[0]
                if (v0 == value_part[-1] and v0 in ("'", '"', '`')):
                    value = value_part[1:-1].replace(f"\\{v0}", v0)
                else:
                    try:
                        value = float(value_part)
                    except:
                        value = value_part
                # word operators need spaces after them in the filter string,
                # but we don't want these later
                return name, operator_type[0].strip(), value
    return [None] * 3


# load example data
@app.callback(
    [Output('text-guides', 'value'),
     Output('available-cas', 'value'),
     Output('available-pam', 'value'),
     Output('available-genome', 'value'),
     Output('checklist-variants', 'value'),
     Output('mms', 'value'),
     Output('dna', 'value'),
     Output('rna', 'value')],
    [Input('load-example-button', 'n_clicks')]
)
def load_example_data(load_button_click: int) -> List[str]:
    """Load data for CRISPRme example run.

    ...

    Parameters
    ----------
    load_button_click : int
        Click on "Load Example" button

    Returns
    -------
    List
        Example parameters
    """

    return [
        "CTAACAGTTGCTTTTATCAC", 
        "SpCas9", 
        "20bp-NGG-SpCas9", 
        "hg38", 
        ["1000G"], 
        "6", 
        "2", 
        "2",
    ]


# Job submission and results URL definition
@app.callback(
    [Output('url', 'pathname'),
     Output('url', 'search')],
    [Input('submit-job', 'n_clicks')],
    [State('url', 'href'),
     State('available-cas', 'value'),
     State('available-genome', 'value'),
     State('checklist-variants', 'value'),
     State('checklist-annotations', 'value'),
     State('vcf-dropdown', 'value'),
     State('annotation-dropdown', 'value'),
     State('available-pam', 'value'),
     State('radio-guide', 'value'),
     State('text-guides', 'value'),
     State('mms', 'value'),
     State('dna', 'value'),
     State('rna', 'value'),
     State('checklist-mail', 'value'),
     State('example-email', 'value'),
     State('job-name', 'value')
     ]
)
def change_url(
    n: int, 
    href: str, 
    nuclease: str, 
    genome_selected: str, 
    ref_var: List, 
    annotation_var: List, 
    vcf_input: str, 
    annotation_input: str, 
    pam: str, 
    guide_type: str, 
    text_guides: List[str], 
    mms: int, 
    dna: int, 
    rna: int, 
    adv_opts: List, 
    dest_email: str, 
    job_name: str
) -> Tuple[str, str]:
    """Launch the targets search and generates the input files for 
    post-processing operations, and results visualization.

    It manages the input data given by the user in the main webpage of CRISPRme 
    and run the search, notify the user by sending an email when the job is
    completed, and produce the link to the webpage used to visualize the results.

    ** Further details **
    
    Perform checks on input parameters' consistency.

    To each received job is assigned a different identifier (or job name). This 
    allows to easily recognize different job submissions. The job IDs consist
    in alpha-numeric strings of 10 characters (A-z 0-9). The IDs are randomly 
    generated. If the generated ID is already assigned to some other job, 
    compute another ID. Every 7 iterations, add +1 to ID length (up to 20 chars 
    as max length). Once generated the ID, create the job directory. Within the
    job directory, create the `queue.txt` file (for job queueing).

    If the input parameters match those of an already processed search, the 
    current job ID is modified to match that of the available analysis (even if
    completed/currently submitted/in queue). Update the email address associated 
    to the job and reset the 3 days availability of the results.

    The current policy of CRISPRme allows up to 2 jobs to run concurrently. The 
    others are put in queue.

    ...

    Parameters
    ----------
    n : int
        Clicks
    href : str
        URL 
    nuclease : str
        Selected nuclease
    genome_selected : str
        Selected genome
    ref_var : List
        Reference variants
    annotation_var : str
        Annotation variants
    vcf_input : str
        Input VCF
    annotation_input : str
        Annotation file
    pam : str
        Selected PAM
    guide_type : str
        RNA guide type
    text_guides : str
        Input guides
    mms : int
        Number of mismatches 
    dna : int
        Number of DNA bulges
    rna : int
        Number of RNA bulges
    adv_opts : List
        Selected advanced options
    dest_email : str
        User mail address
    job_name : str
        Submitted job ID

    Returns
    -------
    Tuple[str, str]
        URL to retrieve CRISPRme analysis results
    """

    if n is not None:
        if not isinstance(n, int):
            raise TypeError(f"Expected {int.__name__}, got {type(n).__name__}")
    if not isinstance(href, str):
        raise TypeError(f"Expected {str.__name__}, got {type(href).__name__}")
    if not isinstance(nuclease, str):
        raise TypeError(
            f"Expected {str.__name__}, got {type(nuclease).__name__}"
        )
    if genome_selected is not None:
        if not isinstance(genome_selected, str):
            raise TypeError(
                f"Expected {str.__name__}, got {type(genome_selected).__name__}"
            )
    if not isinstance(ref_var, list):
        raise TypeError(
            f"Expected {list.__name__}, got {type(ref_var).__name__}"
        )
    if pam is not None:
        if not isinstance(pam, str):
            raise TypeError(
                f"Exepcted {str.__name__}, got {type(pam).__name__}"
            )
    if text_guides is not None:
        if not isinstance(text_guides, str):
            raise TypeError(
                f"Expected {str.__name__}, got {type(text_guides).__name__}"
            )
    if job_name is not None:
        if not isinstance(job_name, str):
            raise TypeError(
                f"Expected {str.__name__}, got {type(job_name).__name__}"
            )

    if n is None:
        raise PreventUpdate  # do not update the page
    # job start
    print("Launching JOB")
    # ---- Check input. If fails, give simple input
    if (genome_selected is None) or (not genome_selected):
        genome_selected = "hg38_ref"  # use hg38 by default
    if( pam is None) or (not pam):
        pam = "20bp-NGG-SpCas9"  # use Cas9 PAM
        guide_seqlen = 20
    else:
        for c in pam.split("-"):
            if "bp" in c:  # use length specified in PAM
                guide_seqlen = int(c.replace("bp", ""))  
    if (text_guides is None) or (not text_guides):
        text_guides = "A" * guide_seqlen
    elif guide_type != "GS":
        text_guides = text_guides.strip()
        if (
            not all(
                [
                    len(c) == len(text_guides.split("\n")[0]) 
                    for c in text_guides.split("\n")
                ]
            )
        ):
            text_guides = select_same_len_guides(text_guides)
    # remove Ns from guides
    guides_tmp = "\n".join(
        [guide.replace("N", "") for guide in text_guides.split("\n")]
    )
    text_guides = guides_tmp.strip()
    # ---- Generate random job ids
    id_len = 10
    for i in range(JOBID_ITERATIONS_MAX):
        # get already assigned job ids
        assigned_ids = [
            d for d in os.listdir(
                os.path.join(current_working_directory, RESULTS_DIR)
            )
            if (
                os.path.isdir(
                    os.path.join(current_working_directory, RESULTS_DIR, d)
                ) and not d.startswith(".")  # avoid hidden files/directories
            )
        ]
        job_id = "".join(
            random.choices(string.ascii_uppercase + string.digits, k=id_len)
        )
        if job_id not in assigned_ids:  # suitable job id
            break
        if i > 7:
            i = 0  # restart
            id_len += 1  # increase ID length
            if id_len > JOBID_MAXLEN:  # reached maximum length
                break
    if job_name:
        assert isinstance(job_name, str)
        job_id = f"{job_name}_{job_id}"
    result_dir = os.path.join(current_working_directory, RESULTS_DIR, job_id)
    # create results directory
    cmd = f"mkdir {result_dir}"
    code = subprocess.call(cmd, shell=True)
    if code != 0:
        raise ValueError(f"An error occurred while running {cmd}")
    # NOTE test command for queue
    cmd = f"touch {os.path.join(current_working_directory, RESULTS_DIR, job_id, 'queue.txt')}"
    subprocess.call(cmd, shell=True)
    # ---- Set search parameters
    # ANNOTATION CHECK
    gencode_name = "gencode.protein_coding.bed"
    annotation_name = ".dummy.bed"  # to proceed without annotation file
    if "EN" in annotation_var:
        annotation_name = "encode+gencode.hg38.bed"
        if "MA" in annotation_var:
            annotation_name = "".join(
                [
                    f"{annotation_name}+",
                    "".join(annotation_input.split(".")[:-1]),
                    ".bed"
                ]
            )
            annotation_dir = os.path.join(
                current_working_directory, ANNOTATIONS_DIR
            )
            annotation_tmp = os.path.join(
                annotation_dir, f"ann_tmp_{job_id}.bed"
            )
            cmd = f"cp {os.path.join(annotation_dir, annotation_name)} {annotation_tmp}"
            code = subprocess.call(cmd, shell=True)
            if code != 0:
                raise ValueError(f"An error occurred while running {cmd}")
            annotation_input_tmp = f"{os.path.join(annotation_dir, annotation_input)}.tmp"
            cmd = f"awk '$4 = $4\"_personal\"' {os.path.join(annotation_dir, annotation_input)} > {annotation_input_tmp}"
            code = subprocess.call(cmd, shell=True)
            if code != 0:
                raise ValueError(f"An error occurred while running {cmd}")
            cmd = f"mv {annotation_input_tmp} {os.path.join(annotation_dir, annotation_input)}"
            code = subprocess.call(cmd, shell=True)
            if code != 0:
                raise ValueError(f"An error occurred while running {cmd}")
            cmd = f"tail -n +1 {os.path.join(annotation_dir, annotation_input)} >> {annotation_tmp}"
            code = subprocess.call(cmd, shell=True)
            if code != 0:
                raise ValueError(f"An error occurred while running {cmd}")
            cmd = f"mv {annotation_tmp} {os.path.join(annotation_dir, annotation_name)}"
            code = subprocess.call(cmd, shell=True)
            if code != 0:
                raise ValueError(f"An error occurred while running {cmd}")
    elif "MA" in annotation_var:
        annotation_input_tmp = f"{os.path.join(annotation_dir, annotation_input)}.tmp"
        cmd = f"awk '$4 = $4\"_personal\"' {os.path.join(annotation_dir, annotation_input)} > {annotation_input_tmp}"
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise ValueError(f"an error occurred while running {cmd}")
        annotation_name = f"{annotation_input}.tmp"
    if "EN" not in annotation_var:
        cmd = f"rm -rf {os.path.join(annotation_dir, '.dummy.bed')}"
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise ValueError(f"An error occurred while running {cmd}")
        cmd = f"touch {os.path.join(annotation_dir, '.dummy.bed')}"
        code = subprocess.call(cmd, shell=True)
        if code != 0:
            raise ValueError(f"An error occurred while running {cmd}")
        gencode_name = '.dummy.bed'

    # GENOME TYPE CHECK
    ref_comparison = False
    genome_type = 'ref'  # Indicates if search is 'ref' or 'both'
    if len(ref_var) > 0:
        ref_comparison = True
        genome_type = 'both'
    # if '+' in genome_selected:
    #     genome_type = 'var'
    # if ref_comparison:
    #     genome_type = 'both'
    # generate_index_ref = False
    # generate_index_enr = False
    search_index = True
    # search = True
    # annotation = True
    # report = True
    genome_selected = genome_selected.replace(' ', '_')
    # genome_ref = genome_selected.split('+')[0]              # + char to separate ref and enr, eg Human_Genome_ref+Human_Genome_1000_genome_project
    # if genome_ref == genome_selected:
    #    ref_comparison = False
    genome_ref = genome_selected
    # NOTE Indexed genomes names are PAM + _ + bMax + _ + genome_selected

    # VCF CHECK
    if genome_type == 'ref':
        sample_list = None
        # dictionary_directory = None
    # else:
    #    sample_list = current_working_directory + \
    #        'samplesID/samples_' + genome_selected + '.txt'
    #    dictionary_directory = current_working_directory + \
    #        'dictionaries/dictionary_' + genome_selected

    # print('valore del vcf', ref_var)
    sample_list = []
    with open(result_dir+"/.list_vcfs.txt", 'w') as vcf_file:
        if len(ref_var) == 0:
            vcf_folder = "_"
            vcf_file.write(vcf_folder+"\n")
        if '1000G' in ref_var:
            # vcf_folder = current_working_directory+"/VCFs/hg38_1000G/"
            vcf_folder = "hg38_1000G/"
            vcf_file.write(vcf_folder+"\n")
            # sample_list.append(current_working_directory +                               'samplesIDs/hg38_1000G.samplesID.txt')
            sample_list.append('hg38_1000G.samplesID.txt')
            # dictionary_directory = current_working_directory + \
            #     'dictionaries/dictionary_VCFs_1000_genome_project'
        if 'HGDP' in ref_var:
            # vcf_folder = current_working_directory+"/VCFs/hg38_HGDP/"
            vcf_folder = "hg38_HGDP/"
            vcf_file.write(vcf_folder+"\n")
            # sample_list.append(current_working_directory +'samplesIDs/hg38_HGDP.samplesID.txt')
            sample_list.append('hg38_HGDP.samplesID.txt')
            # dictionary_directory = current_working_directory + \
            #     'dictionaries/dictionary_VCFs_1000_genome_project'
        if "PV" in ref_var:
            # genome_selected = data_personal_genome[sel_cel_genome['row']
            #                                       ]["Personal Genomes"]
            # vcf_folder = current_working_directory+"/VCFs/" + vcf_input
            vcf_folder = vcf_input
            vcf_file.write(vcf_folder+"\n")
            # genome_ref = genome_selected.split("+")[0].replace(" ", "_")
            # sample_list = current_working_directory + 'samplesID/' + vcf_input
            # sample_list.append(current_working_directory + "/samplesID/" + vcf_input + "/" + [f for f in listdir(
            #     current_working_directory + 'samplesID/' + vcf_input) if isfile(join(current_working_directory + 'samplesID/' + vcf_input, f))][0])
            # sample_list.append(current_working_directory + "/samplesIDs/" + vcf_input + ".samplesID.txt")
            sample_list.append(vcf_input + ".samplesID.txt")
            # dictionary_directory = current_working_directory + \
            #     'dictionaries/'+genome_ref+'+'+vcf_input

    with open(result_dir+"/.samplesID.txt", 'w') as sampleID_file:
        for ele in sample_list:
            sampleID_file.write(ele+"\n")
    # for ele in sample_list:
    #     os.system(f"tail -n +2 {ele} >> {result_dir}/sampleID.txt")
    # os.system(f"sed -i 1i\"#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tGENDER\" {result_dir}/sampleID.txt")

    send_email = False
    if adv_opts is None:
        adv_opts = []

    if 'email' in adv_opts and dest_email is not None and len(dest_email.split('@')) > 1 and dest_email.split('@')[-1] != '':
        send_email = True
        with open(result_dir + '/email.txt', 'w') as e:
            e.write(dest_email + '\n')
            # e.write(URL + '/load?job=' + job_id + '\n')
            e.write(''.join(href.split('/')[:-1]
                            ) + '/load?job=' + job_id + '\n')
            e.write(datetime.utcnow().strftime("%m/%d/%Y, %H:%M:%S") + '\n')
            # e.write('Job done. Parameters: etc etc')
            e.close()
    else:
        dest_email = "_"

    pam_len = 0
    with open(current_working_directory + 'PAMs/' + pam + '.txt') as pam_file:
        pam_char = pam_file.readline()
        index_pam_value = pam_char.split(' ')[-1]
        if int(pam_char.split(' ')[-1]) < 0:
            end_idx = int(pam_char.split(' ')[-1]) * (-1)
            pam_char = pam_char.split(' ')[0][0: end_idx]
            pam_len = end_idx
            pam_begin = True
        else:
            end_idx = int(pam_char.split(' ')[-1])
            pam_char = pam_char.split(' ')[0][end_idx * (-1):]
            pam_len = end_idx
            pam_begin = False

    if guide_type == 'GS':
        text_sequence = text_guides
        # print(text_sequence)
        # exit()
        # Extract sequence and create the guides
        guides = []
        # extracted_seqs = list()
        # for lines in text_sequence:
        #     print('linea', lines)
        for name_and_seq in text_sequence.split('>'):
            if '' == name_and_seq:
                continue
            name = name_and_seq[:name_and_seq.find('\n')]
            seq = name_and_seq[name_and_seq.find('\n'):]
            # seq = seq.strip().split()
            # seq = ''.join(seq)
            seq = seq.strip()
            # name, seq = name_and_seq.strip().split('\n')
            if 'chr' in seq:
                # extracted_seq = extract_seq.extractSequence(
                #         name, seq, genome_ref.replace(' ', '_'))
                for single_row in seq.split('\n'):
                    if '' == single_row:
                        continue
                    pieces_of_row = single_row.strip().split()
                    seq_to_extract = pieces_of_row[0]+":" + \
                        pieces_of_row[1]+"-"+pieces_of_row[2]
                    extracted_seq = extract_seq.extractSequence(
                        name, seq_to_extract, genome_ref.replace(' ', '_'))
                    guides.extend(convert_pam.getGuides(
                        extracted_seq, pam_char, guide_seqlen, pam_begin))
            else:
                seq = seq.split()
                seq = ''.join(seq)
                extracted_seq = seq.strip()
                guides.extend(convert_pam.getGuides(
                    extracted_seq, pam_char, guide_seqlen, pam_begin))
            # print('extracted seq', extracted_seq)
            # guides.extend(convert_pam.getGuides(
            #     extracted_seq, pam_char, len_guide_sequence, pam_begin))
        guides = list(set(guides))
        if not guides:
            guides = 'A'*guide_seqlen
        text_guides = '\n'.join(guides).strip()
    # print(text_guides, 'and', guides, 'and', pam_char)
    # exit()
    text_guides = text_guides.upper()
    new_test_guides = list()
    for guide in text_guides.split('\n'):
        if len(guide) == guide_seqlen:
            new_test_guides.append(guide)
    if not new_test_guides:
        new_test_guides.append('A'*guide_seqlen)
    text_guides = '\n'.join(new_test_guides)
    for g in text_guides.split('\n'):
        for c in g:
            if c not in VALID_CHARS:
                text_guides = text_guides.replace(c, '')
    if len(text_guides.split('\n')) > 100:  # set limit to 100 guides per run in the website
        text_guides = '\n'.join(text_guides.split('\n')[:100]).strip()
    # len_guides = len(text_guides.split('\n')[0])
    len_guides = guide_seqlen

    # if 'A'*len_guide_sequence in text_guides:
    #     error_in_seq_extract = True
    #     return '/index', ''
    #     error = dbc.Alert(
    #         "The input spacer(s) or sequence is not a valid input, please check your input and retry.", color="danger")
    #     visibility = {'display': 'true'}
    #     return '/index', '', error, visibility

    # Adjust guides by adding Ns to make compatible with Crispritz
    if (pam_begin):
        pam_to_file = pam_char + ('N' * len_guides) + ' ' + index_pam_value
        pam_to_indexing = pam_char + ('N' * 25) + ' ' + index_pam_value
    else:
        pam_to_file = ('N' * len_guides) + pam_char + ' ' + index_pam_value
        pam_to_indexing = ('N' * 25) + pam_char + ' ' + index_pam_value

    save_pam_file = open(result_dir + '/.pam.txt', 'w')
    save_pam_file.write(pam_to_file)
    save_pam_file.close()
    pam = result_dir + '/.pam.txt'

    guides_file = result_dir + '/.guides.txt'
    if text_guides is not None and text_guides != '':
        save_guides_file = open(result_dir + '/.guides.txt', 'w')
        if (pam_begin):
            text_guides = 'N' * pam_len + \
                text_guides.replace('\n', '\n' + 'N' * pam_len)
        else:
            text_guides = text_guides.replace(
                '\n', 'N' * pam_len + '\n') + 'N' * pam_len
        save_guides_file.write(text_guides)
        save_guides_file.close()

    # if (int(dna) == 0 and int(rna) == 0):
    #     search_index = False
    max_bulges = rna
    if (int(dna) > int(rna)):
        max_bulges = dna

    if (search_index):
        search = False

    # Check if index exists, otherwise set generate_index to true
    # genome_indices_created = [f for f in listdir(
        # current_working_directory + 'genome_library') if isdir(join(current_working_directory + 'genome_library', f))]
    genome_idx_list = []
    if genome_type == 'ref':
        genome_idx = pam_char + '_' + str(max_bulges) + '_' + genome_selected
        genome_idx_list.append(genome_idx)
    else:
        if "1000G" in ref_var:
            genome_idx = pam_char + '_' + \
                str(max_bulges) + '_' + genome_selected + \
                '+hg38_1000G'
            genome_idx_list.append(genome_idx)
        if 'HGDP' in ref_var:
            genome_idx = pam_char + '_' + \
                str(max_bulges) + '_' + genome_selected + \
                '+hg38_HGDP'
            genome_idx_list.append(genome_idx)
        if "PV" in ref_var:
            genome_idx = pam_char + '_' + \
                str(max_bulges) + '_' + genome_selected + '+' + vcf_input
            genome_idx_list.append(genome_idx)
    genome_idx = ','.join(genome_idx_list)
    # genome_idx = ''
    # genome_idx_ref = ''
    # for gidx in genome_indices_created:
    # if pam_char in gidx and genome_selected in gidx:
    # if int(gidx.split('_')[1]) >= int(max_bulges):
    # if '+' in gidx:
    # genome_idx = gidx
    # genome_idx_ref = gidx.split('+')[0]

    # if genome_idx == '':
    # if int(max_bulges) > 0:
    # generate_index_enr = True
    # with open(result_dir + '/pam_indexing.txt', 'w+') as pam_id_file:
    # pam_id_file.write(pam_to_indexing)
    # genome_idx = pam_char + '_' + str(max_bulges) + '_' + genome_selected
    # if genome_idx_ref == '':
    # if int(max_bulges) > 0:
    # generate_index_ref = True
    # with open(result_dir + '/pam_indexing.txt', 'w+') as pam_id_file:
    # pam_id_file.write(pam_to_indexing)
    # genome_idx_ref = pam_char + '_' + \
    # str(max_bulges) + '_' + genome_selected.split('+')[0]

    # if ref_var != '':
    #     generate_index_enr = False
    # genome_idx = pam_char + '_' + '2' + '_' + genome_selected
    # genome_idx_ref = genome_idx.split('+')[0]

    # Create .Params.txt file
    with open(result_dir + '/.Params.txt', 'w') as p:
        p.write('Genome_selected\t' + genome_selected + '\n')
        p.write('Genome_ref\t' + genome_ref + '\n')
        if search_index:
            p.write('Genome_idx\t' + genome_idx + '\n')
        else:
            p.write('Genome_idx\t' + 'None\n')
        p.write('Pam\t' + pam_char + '\n')
        p.write('Max_bulges\t' + str(max_bulges) + '\n')
        p.write('Mismatches\t' + str(mms) + '\n')
        p.write('DNA\t' + str(dna) + '\n')
        p.write('RNA\t' + str(rna) + '\n')
        p.write('Annotation\t' + str(annotation_name) + '\n')
        p.write('Nuclease\t' + str(nuclease) + '\n')
        # p.write('Gecko\t' + str(gecko_comp) + '\n')
        p.write('Ref_comp\t' + str(ref_comparison) + '\n')
        p.close()

    # 4) Check if input parameters (mms, bulges, pam, guides, genome) are the same as a previous search
    all_result_dirs = [f for f in listdir(
        current_working_directory + 'Results') if isdir(join(current_working_directory + 'Results', f))]
    all_result_dirs.remove(job_id)
    # all_result_dirs.remove('test')
    for check_param_dir in all_result_dirs:
        if os.path.exists(current_working_directory + 'Results/' + check_param_dir + '/.Params.txt'):
            # if os.path.exists(current_working_directory + 'Results/' + check_param_dir + '/log.txt'):
            # with open(current_working_directory + 'Results/' + check_param_dir + '/log.txt') as log:
            # if ('Job\tDone' in log.read()):
            if (filecmp.cmp(current_working_directory + 'Results/' + check_param_dir + '/.Params.txt', result_dir + '/.Params.txt')):
                guides1 = open(current_working_directory + 'Results/' +
                               check_param_dir + '/.guides.txt').read().split('\n')
                guides2 = open(current_working_directory + 'Results/' +
                               job_id + '/.guides.txt').read().split('\n')
                if (collections.Counter(guides1) == collections.Counter(guides2)):
                    if os.path.exists(current_working_directory + 'Results/' + check_param_dir + '/log.txt'):
                        adj_date = False
                        with open(current_working_directory + 'Results/' + check_param_dir + '/log.txt') as log:
                            log_content = log.read().strip()
                            if ('Job\tDone' in log_content):
                                adj_date = True
                                log_content = log_content.split('\n')
                                new_date = subprocess.Popen(
                                    ['echo $(date)'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                                out, err = new_date.communicate()
                                rewrite = '\n'.join(
                                    log_content[:-1]) + '\nJob\tDone\t' + out.decode('UTF-8').strip()
                        if adj_date:
                            with open(current_working_directory + 'Results/' + check_param_dir + '/log.txt', 'w+') as log:
                                log.write(rewrite)
                                # Send mail
                            if send_email:
                                with open(current_working_directory + 'Results/' + job_id + '/email.txt', 'w+') as e:
                                    e.write(dest_email + '\n')
                                    e.write(
                                        ''.join(href.split('/')[:-1]) + '/load?job=' + job_id + '\n')
                                    e.write(datetime.utcnow().strftime(
                                        "%m/%d/%Y, %H:%M:%S") + '\n')
                                # Send mail with file in job_id dir with link to job already done, note that job_id directory will be deleted
                                # subprocess.call(['python ' + app_main_directory + 'send_mail.py ' +
                                #                  current_working_directory + 'Results/' + job_id], shell=True)

                        elif send_email:
                            # Job is not finished, add this user email to email.txt and when job is done send to both
                            if os.path.exists(current_working_directory + 'Results/' + check_param_dir + '/email.txt'):
                                with open(current_working_directory + 'Results/' + check_param_dir + '/email.txt', 'a+') as e:
                                    e.write('--OTHEREMAIL--')
                                    e.write(dest_email + '\n')
                                    e.write(
                                        ''.join(href.split('/')[:-1]) + '/load?job=' + job_id + '\n')
                                    e.write(datetime.utcnow().strftime(
                                        "%m/%d/%Y, %H:%M:%S") + '\n')
                            else:
                                with open(current_working_directory + 'Results/' + check_param_dir + '/email.txt', 'w+') as e:
                                    e.write(dest_email + '\n')
                                    e.write(
                                        ''.join(href.split('/')[:-1]) + '/load?job=' + job_id + '\n')
                                    e.write(datetime.utcnow().strftime(
                                        "%m/%d/%Y, %H:%M:%S") + '\n')

                        subprocess.call(
                            ['rm -r ' + current_working_directory + 'Results/' + job_id], shell=True)
                        return '/load', '?job=' + check_param_dir
                    else:
                        # We may have entered a jobdir that was in queue
                        if os.path.exists(current_working_directory + 'Results/' + check_param_dir + '/queue.txt'):
                            if send_email:
                                if os.path.exists(current_working_directory + 'Results/' + check_param_dir + '/email.txt'):
                                    with open(current_working_directory + 'Results/' + check_param_dir + '/email.txt', 'a+') as e:
                                        e.write('--OTHEREMAIL--')
                                        e.write(dest_email + '\n')
                                        e.write(
                                            ''.join(href.split('/')[:-1]) + '/load?job=' + job_id + '\n')
                                        e.write(datetime.utcnow().strftime(
                                            "%m/%d/%Y, %H:%M:%S") + '\n')
                                else:
                                    with open(current_working_directory + 'Results/' + check_param_dir + '/email.txt', 'w+') as e:
                                        e.write(dest_email + '\n')
                                        e.write(
                                            ''.join(href.split('/')[:-1]) + '/load?job=' + job_id + '\n')
                                        e.write(datetime.utcnow().strftime(
                                            "%m/%d/%Y, %H:%M:%S") + '\n')
                            return '/load', '?job=' + check_param_dir

    '''
    command = app_main_directory + 'PostProcess/./submit_job.final.sh ' + current_working_directory + 'Results/' + job_id + ' ' + current_working_directory +  'Genomes/' + genome_selected + ' ' + current_working_directory + 'Genomes/' + genome_ref + ' ' + current_working_directory + 'genome_library/' + genome_idx + (
        ' ' + pam + ' ' + guides_file + ' ' + str(mms) + ' ' + str(dna) + ' ' + str(rna) + ' ' + str(search_index) + ' ' + str(search) + ' ' + str(annotation) + (
            ' ' + str(report) + ' ' + str(gecko_comp) + ' ' + str(ref_comparison) + ' ' + current_working_directory +  'genome_library/' + genome_idx_ref + ' ' + str(send_email) + ' '  + current_working_directory +  'Annotations/' + annotation_file +
            ' ' + genome_type + ' ' + app_main_directory + ' ' + str(dictionary_directory) + ' ' + str(sample_list) + ' ' + str(generate_index_ref) + ' ' + str(generate_index_enr) + ' ' + current_working_directory))
    '''
    # merge default is 3 nt wide
    merge_default = 3
    print(
        f"Submitted JOB {job_id}. The stdout is redirected in log_verbose.txt and stderr is redirected in log_error.txt")
    command = f"{app_main_directory}/PostProcess/./submit_job_automated_new_multiple_vcfs.sh {current_working_directory}/Genomes/{genome_ref} {result_dir}/.list_vcfs.txt {guides_file} {pam} {current_working_directory}/Annotations/{annotation_name} {result_dir}/.samplesID.txt {max([int(dna), int(rna)])} {mms} {dna} {rna} {merge_default} {result_dir} {app_main_directory}/PostProcess {8} {current_working_directory} {current_working_directory}/Annotations/{gencode_name} {dest_email} 1> {result_dir}/log_verbose.txt 2>{result_dir}/log_error.txt"
    # with open(f"{result_dir}/log_verbose.txt", 'w') as log_verbose:
    # log_verbose = open(f"{result_dir}/log_verbose.txt", 'w')
    # , stdout=log_verbose)
    exeggutor.submit(subprocess.run, command, shell=True)
    # subprocess.run(command, shell=True, stdout=log_verbose)
    return '/load', '?job=' + job_id


# Check input presence
@ app.callback(
    [Output('submit-job', 'n_clicks'),
     Output('modal', 'is_open'),
     Output('available-genome', 'className'),
     Output('available-pam', 'className'),
     Output('text-guides', 'style'),
     Output('mms', 'className'),
     Output('dna', 'className'),
     Output('rna', 'className'),
     Output('len-guide-sequence-ver', 'className'),
     Output('warning-list', 'children')],
    [Input('check-job', 'n_clicks'),
     Input('close', 'n_clicks')],
    [State('available-genome', 'value'),
     State('available-pam', 'value'),
     State('radio-guide', 'value'),
     State('text-guides', 'value'),
     State('mms', 'value'),
     State('dna', 'value'),
     State('rna', 'value'),
     State("modal", "is_open")]
)
# len_guide_seq, active_tab ,
def checkInput(n, n_close, genome_selected, pam, guide_type, text_guides, mms, dna, rna, is_open):
    '''
    Checks the presence and correctness of the input fields, changing their border to red if it's missing and displaying a Modal element with
    a list of the missing inputs. The callback is triggered when the user clicks on the 'Submit' button or when the modal is closed
    (by the 'Close' button or by clicking on-screen when the Modal element is open)

    ***Args***

    + [**n**] **check-job** (*n_clicks*): button that starts the job after checking the correctness of the input
    + [**n_close**] **close** (*n_clicks*): button that closes the opened modal
    + [**genome_selected**] **available-genome** (*value*): string of the selected genome from the Dropdown. `None` if not selected, '' if selected
    and then deleted
    + [**pam**] **available-pam** (*value*): string of the selected PAM from the Dropdown. `None` if not selected, '' if selected
    and then deleted
    + [**text_guides**] **text-guides** (*value*): string of the input guides from the Textarea. `None` if not selected, '' if selected
    and then deleted
    + [**mms**] **mms** (*value*): int of the selected mismatch value. `None` if not selected, '' if selected
    and then deleted
    + [**dna**] **dna** (*value*): int of the selected DNA bulge value. `None` if not selected, '' if selected
    and then deleted
    + [**rna**] **rna** (*value*): int of the selected RNA bulge value. `None` if not selected, '' if selected
    and then deleted
    + [**len_guide_seq**] **len-guide-sequence-ver** (*value*): int value of the length of the guides (available when 'Sequence' tab is active). `None` if not selected, '' if selected
    and then deleted
    + [**active_tab**] **tabs** (*active_tab*): string of the ID of the active tab ('Guide' or 'Sequence')
    + [**is_open**] **modal** (*is_open*): True if the modal is displayed, false otherwise

    ***Returns***

    + **submit-job** (*n_clicks*): reset the click counter (`None`) if some inputs are missing, else put to `1` in order to trigger the `changeUrl()`
    function and proceed with the job
    + **modal** (*is_open*): `True` if some inputs are missing, `False` otherwise
    + **available-genome** (*className*): string containing the name of the css class for the red-border ('missing-input'), indicating a missing input. `None` if
    input is ok
    + **available-pam** (*className*): string containing the name of the css class for the red-border ('missing-input'), indicating a missing input. `None` if
    input is ok
    + **text-guides** (*style*): dictionary for the style element (NOTE not updated if input is missing)
    + **mms** (*className*): string containing the name of the css class for the red-border ('missing-input'), indicating a missing input. `None` if
    input is ok
    + **dna** (*className*): string containing the name of the css class for the red-border ('missing-input'), indicating a missing input. `None` if
    input is ok
    + **rna** (*className*): string containing the name of the css class for the red-border ('missing-input'), indicating a missing input. `None` if
    input is ok
    + **len-guide-sequence-ver** (*className*): string containing the name of the css class for the red-border ('missing-input'), indicating a missing input. `None` if
    input is ok
    + **warning-list** (*children*): html.Div for the Modal component, listing the missing inputs
    '''
    print("Check input for JOB")
    if n is None:
        raise PreventUpdate
    if is_open is None:
        is_open = False

    classname_red = 'missing-input'
    genome_update = None
    pam_update = None
    text_update = {'width': '300px', 'height': '30px'}
    mms_update = None
    dna_update = None
    rna_update = None
    len_guide_update = None
    update_style = False
    miss_input_list = []

    if genome_selected is None or genome_selected == '':
        genome_update = classname_red
        update_style = True
        miss_input_list.append('Genome')
    if genome_selected is None or genome_selected == '':
        genome_selected = 'hg38_ref'
    genome_ref = genome_selected
    if pam is None or pam == '':
        pam_update = classname_red
        update_style = True
        miss_input_list.append('PAM')
    # if text_guides is None or text_guides == '':
        # text_update = {'width':'450px', 'height':'160px','border': '1px solid red'}
        # update_style = True
        # miss_input_list.append('crRNA sequence(s)')
    if mms is None or str(mms) == '':
        mms_update = classname_red
        update_style = True
        miss_input_list.append('Allowed Mismatches')
    if dna is None or str(dna) == '':
        dna_update = classname_red
        update_style = True
        miss_input_list.append('Bulge DNA size')
    if rna is None or str(rna) == '':
        rna_update = classname_red
        update_style = True
        miss_input_list.append('Bulge RNA size')
    # if (len_guide_seq is None or str(len_guide_seq) == ''): # and ('sequence-tab' in active_tab)
    #    len_guide_update = classname_red
    #    update_style = True
    #    miss_input_list.append('gRNA length')
    if pam is None or pam == '':
        pam = '20bp-NGG-SpCas9'
        len_guide_sequence = 20
    else:
        for elem in pam.split('-'):
            if 'bp' in elem:
                len_guide_sequence = int(elem.replace('bp', ''))

    no_guides = False
    if text_guides is None or text_guides == '':
        text_guides = 'A'*len_guide_sequence
        no_guides = True
        # text_guides = 'GAGTCCGAGCAGAAGAAGAA\nCCATCGGTGGCCGTTTGCCC'
    elif guide_type != 'GS':
        text_guides = text_guides.strip()
        if (not all(len(elem) == len(text_guides.split('\n')[0]) for elem in text_guides.split('\n'))):
            text_guides = selectSameLenGuides(text_guides)

    pam_len = 0
    with open(current_working_directory + 'PAMs/' + pam + '.txt') as pam_file:
        pam_char = pam_file.readline()
        index_pam_value = pam_char.split(' ')[-1]
        if int(pam_char.split(' ')[-1]) < 0:
            end_idx = int(pam_char.split(' ')[-1]) * (-1)
            pam_char = pam_char.split(' ')[0][0: end_idx]
            pam_len = end_idx
            pam_begin = True
        else:
            end_idx = int(pam_char.split(' ')[-1])
            pam_char = pam_char.split(' ')[0][end_idx * (-1):]
            pam_len = end_idx
            pam_begin = False

    if guide_type == 'GS':
        text_sequence = text_guides
        # print(text_sequence)
        # exit()
        # Extract sequence and create the guides
        guides = []
        # extracted_seqs = list()
        # for lines in text_sequence:
        #     print('linea', lines)
        for name_and_seq in text_sequence.split('>'):
            if '' == name_and_seq:
                continue
            name = name_and_seq[:name_and_seq.find('\n')]
            seq = name_and_seq[name_and_seq.find('\n'):]
            # seq = seq.strip().split()
            # seq = ''.join(seq)
            seq = seq.strip()
            # name, seq = name_and_seq.strip().split('\n')
            if 'chr' in seq:
                # extracted_seq = extract_seq.extractSequence(
                #         name, seq, genome_ref.replace(' ', '_'))
                for single_row in seq.split('\n'):
                    if '' == single_row:
                        continue
                    pieces_of_row = single_row.strip().split()
                    seq_to_extract = pieces_of_row[0]+":" + \
                        pieces_of_row[1]+"-"+pieces_of_row[2]
                    extracted_seq = extract_seq.extractSequence(
                        name, seq_to_extract, genome_ref.replace(' ', '_'))
                    guides.extend(convert_pam.getGuides(
                        extracted_seq, pam_char, len_guide_sequence, pam_begin))
            else:
                seq = seq.split()
                seq = ''.join(seq)
                extracted_seq = seq.strip()
                guides.extend(convert_pam.getGuides(
                    extracted_seq, pam_char, len_guide_sequence, pam_begin))
            # print('extracted seq', extracted_seq)
            # guides.extend(convert_pam.getGuides(
            #     extracted_seq, pam_char, len_guide_sequence, pam_begin))
        guides = list(set(guides))
        if not guides:
            guides = 'A'*len_guide_sequence
            no_guides = True
        text_guides = '\n'.join(guides).strip()
    # print(text_guides, 'and', guides, 'and', pam_char)
    # exit()
    text_guides = text_guides.upper()
    new_test_guides = list()
    for guide in text_guides.split('\n'):
        guide = guide.replace('N', '')
        if len(guide) == len_guide_sequence:
            new_test_guides.append(guide)
    if not new_test_guides:
        new_test_guides.append('A'*len_guide_sequence)
        no_guides = True
    text_guides = '\n'.join(new_test_guides)
    for g in text_guides.split('\n'):
        for c in g:
            if c not in VALID_CHARS:
                text_guides = text_guides.replace(c, '')
    # set limit to 100 guides per run in the website
    if len(text_guides.split('\n')) > 1000000000:
        text_guides = '\n'.join(text_guides.split('\n')[:1000000000]).strip()
    # len_guides = len(text_guides.split('\n')[0])
    len_guides = len_guide_sequence

    if no_guides:
        text_update = {'width': '300px',
                       'height': '30px', 'border': '1px solid red'}
        update_style = True
        miss_input_list.append(
            'Input at least one correct guide, correct guides must have the length requested for the selected PAM sequence (e.g., 20bp, 21bp, etc)')

    miss_input = html.Div(
        [
            html.P('The following inputs are missing:'),
            html.Ul([html.Li(x) for x in miss_input_list]),
            html.P('Please fill in the values before submitting the job')
        ]
    )

    if not update_style:
        print('All good')
        return 1, False, genome_update, pam_update, text_update, mms_update, dna_update, rna_update, len_guide_update, miss_input
    return None, not is_open, genome_update, pam_update, text_update, mms_update, dna_update, rna_update, len_guide_update, miss_input


@ app.callback(
    Output('fade-len-guide', 'is_in'),
    [Input('tabs', 'active_tab')],
    [State('fade-len-guide', 'is_in')]
)
def resetTab(current_tab, is_in):
    '''
    Manages the fading of the Dropdown for the guide length when the tab 'Sequence' is active.

    ***Args***

    + [**current_tab**] **tabs** (*active_tab*): string of the ID of the current active tab
    + [**is_in**] **fade-len-guide** (*is_in*): True if the Dropdown guide length element is displayed, False otherwise

    ***Returns***

    + **fade-len-guide** (*is_in*): True in order to show the Dropdown guide length element, False to hide it
    '''
    if current_tab is None:
        raise PreventUpdate

    if current_tab == 'guide-tab':
        return False
    return True

# Email validity


@ app.callback(
    Output('example-email', 'style'),
    [Input('example-email', 'value')]
)
def checkEmailValidity(val):
    '''
    Checks email validity (i.e. an '@' is present) and changes the border to green or red accordingly.

    ***Args***

    + [**val**] **example-email** (*value*): string of the email

    ***Returns***

    + **example-email** (*style*): dictionary of the style for the Input email: change border to red (if no '@' in **val**) or to green
    '''
    if val is None:
        raise PreventUpdate

    if '@' in val:
        return {'border': '1px solid #94f033', 'outline': '0'}
    return {'border': '1px solid red'}

#################################################
# Fade in/out email


@ app.callback(
    Output('example-email', 'disabled'),
    [Input('checklist-mail', 'value')]
)
def disabled_mail(checklist_value):
    # print('value', checklist_ value)
    if 'email' not in checklist_value:
        return True
    elif 'email' in checklist_value:
        return False


@ app.callback(
    Output('job-name', 'disabled'),
    [Input('checklist-job-name', 'value')]
)
def disable_job_name(checklist_value):
    # print('value', checklist_value)
    if 'job_name' not in checklist_value:
        return True
    elif 'job_name' in checklist_value:
        return False

    # @app.callback(
    #     [Output("fade", "is_in"),
    #      Output('example-email', 'disabled')],
    #     [Input("checklist-mail", "value")],
    #     [State("fade", "is_in")],
    # )
    # def fade_toggle(selected_options, is_in):
    #     '''
    #     Manages the fading of the InputBox for the email when the option 'Notify me by email' is selected/deselected.

    #     ***Args***

    #     + [**selected_options**] **checklist-mail** (*value*): list of IDs of the options checked
    #     + [**is_in**] **fade** (*is_in*): True if the Input email element is displayed, False otherwise

    #     ***Returns***

    #     + **fade** (*is_in*): True in order to show the Input email element, False to hide it
    #     '''
    #     if selected_options is None:
    #         # return False
    #         return True, True
    #     if 'email' in selected_options:
    #         return True, False
    #     # return False
    #     return True, True

    # @app.callback(
    #     [Output('div-browse-PV', 'style')],
    #     [Input('radio-genome', 'value')]
    # )
    # def changePlaceholderGuideTextArea(value):
    #     # print(value)
    #     if value == 'PV':
    #         return [{'visibility': 'visible'}]
    #     else:
    #         return [{'visibility': 'visible'}]


@ app.callback(
    [Output('vcf-dropdown', 'disabled'),
     Output('vcf-dropdown', 'value')],
    [Input('checklist-variants', 'value')]
)
def change_disabled_vcf_dropdown(checklist_value):
    if 'PV' in checklist_value:
        return False, ""
    else:
        return True, ""


@ app.callback(
    [Output('annotation-dropdown', 'disabled'),
     Output('annotation-dropdown', 'value')],
    [Input('checklist-annotations', 'value')]
)
def change_disabled_vcf_dropdown(checklist_value):
    if 'MA' in checklist_value:
        return False, ""
    else:
        return True, ""


# @app.callback(
#     [Output('div-browse-annotation', 'style')],
#     [Input('checklist-annotations', 'value')]
# )
# def changePlaceholderGuideTextArea(value):
#     # print(value)
#     if value == 'MA':
#         return [{'visibility': 'visible'}]
#     else:
#         return [{'visibility': 'visible'}]


@ app.callback(
    [Output('available-pam', 'options')],
    [Input('available-cas', 'value')]
)
def changePlaceholderGuideTextArea(value):
    all_options = availablePAM()
    correct_options = []
    for option in all_options:
        if value == option['label'].split('.')[0].split('-')[2]:
            correct_options.append(
                {'label': option['label'], 'value': option['value']})
    return [correct_options]


@ app.callback(
    [Output('text-guides', 'placeholder')],
    [Input('radio-guide', 'value')]
)
def changePlaceholderGuideTextArea(value):
    if value == 'IP':
        return ['GAGTCCGAGCAGAAGAAGAA\nCCATCGGTGGCCGTTTGCCC']
    elif value == 'GS':
        return ['>sequence1\nAAGTCCCAGGACTTCAGAAGagctgtgagaccttggc\n>sequence_bed\nchr1 11130540 11130751\nchr1 1023000 1024000']


def select_same_len_guides(guides: str) -> str:
    """If the user provides guides of different lengths, compute the length of 
    the first given guide and keep only those with the same length.

    ...

    Parameters
    ----------
    guides : str
        Guides

    Returns
    -------
    str
        Selected guides
    """

    if not isinstance(guides, str):
        raise TypeError(f"Expected {str.__name__}, got {type(guides).__name__}")
    length = len(guides.split("\n")[0])
    same_len_guides = [
        guide for guide in guides.split("\n") if len(guide) == length
    ]
    same_len_guides = "\n".join(same_len_guides).strip()
    return same_len_guides


def availablePAM():
    '''
    Returns a list of dictionaries of the available PAMs in the 'PAM' directory.

    ***Returns***

    + **pam_file** (*list* of {'label': pam, 'value': pam}): list containing a series of dictionaries, one for each PAM file found in
        the 'PAM' directory. Used as input parameter for the 'options' element of a Dash Droplist
    '''
    onlyfile = [
        f for f in os.listdir(os.path.join(current_working_directory, "PAMs"))
        if (
            not f.startswith(".") and 
            os.path.isfile(
                os.path.join(current_working_directory, "PAMs", f)
            )
        )
    ]
    # removed .txt for better visualization
    onlyfile = [x.replace('.txt', '') for x in onlyfile]
    pam_file = []
    for pam_name in onlyfile:
        if 'tempPAM' in pam_name:  # Skip the temp pam used for updating dictionaries
            pass
        else:
            pam_file.append({'label': pam_name, 'value': pam_name})
    return pam_file


def availableCAS():
    onlyfile = [
        f for f in os.listdir(os.path.join(current_working_directory, "PAMs"))
        if (
            not f.startswith(".") and 
            os.path.isfile(
                os.path.join(current_working_directory, "PAMs", f)
            )
        )
    ]
    # removed .txt for better visualization
    onlyfile = [x.replace('.txt', '') for x in onlyfile]
    cas_file = []
    cas_set = set()
    for cas_name in onlyfile:
        if 'tempPAM' in cas_name:  # Skip the temp pam used for updating dictionaries
            pass
        else:
            # if 'Cas' in cas_name:
            cas_prot = cas_name.split('.')[0].split('-')[2]
            cas_set.add(cas_prot)
    for ele in sorted(cas_set):
        cas_file.append({'label': ele, 'value': ele})
    return cas_file


# @app.callback(
#     [Output('checklist-variants', 'value')],
#     [Input('available-genome', 'value')]
# )
# def reset_radio_genome(value):
#     return ['']


@ app.callback(
    [Output('checklist-variants', 'options'),
     Output('vcf-dropdown', 'options')],
    [Input('available-genome', 'value')]
)
def changeVariantsChecklistState(genome_value):
    if genome_value is not None:
        checklist_variants_options = []
        checklist_variants_options.append({'label': ' plus 1000 Genome Project variants',
                                           'value': '1000G', 'disabled': False})
        checklist_variants_options.append({'label': ' plus HGDP variants',
                                           'value': 'HGDP', 'disabled': False})
        checklist_variants_options.append({'label': ' plus personal variants*',
                                           'value': 'PV', 'disabled': True})
    personal_vcf = get_more_VCF(genome_value)
    return [checklist_variants_options, personal_vcf]


def availableGenomes():
    '''
    Returns a list of dictionaries of the available genomes in the 'Genomes' directory.

    ***Returns***

    + **gen_dir** (*list* of {'label': genome, 'value': genome}): list containing a series of dictionaries, one for each directory (genome) found in
    the 'Genomes' directory. Used as input parameter for the 'options' element of a Dash Droplist
    '''
    onlydir = [f for f in os.listdir(current_working_directory + 'Genomes')
               if os.path.isdir(os.path.join(current_working_directory + 'Genomes', f))]
    onlydir = [x.replace('_', ' ') for x in onlydir]
    gen_dir = []
    for dir in onlydir:
        if "+" not in dir and 'None' not in dir:
            gen_dir.append({'label': dir, 'value': dir})
    return gen_dir


def get_more_annotations():
    annotation_dir = glob.glob(current_working_directory + 'Annotations/*.bed')
    annotation_list = []

    for elem in annotation_dir:
        if 'encode' not in elem and 'None' not in elem and 'dummy' not in elem and 'tmp' not in elem:
            annotation_list.append({'label': elem.strip().split(
                '/')[-1], 'value': elem.strip().split('/')[-1]})

    return annotation_list


def get_more_VCF(genome_value):
    onlydir = [f for f in os.listdir(current_working_directory + 'VCFs')
               if os.path.isdir(os.path.join(current_working_directory + 'VCFs', f))]
    vcf_dir = []
    genome_value = genome_value.replace(" ", "_")
    for dir in onlydir:
        if 'hg38_HGDP' not in dir and 'hg38_1000G' not in dir and 'None' not in dir and genome_value in dir:
            vcf_dir.append({'label': dir, 'value': dir})
    return vcf_dir


def indexPage():
    '''
    Creates the layout of the main page ('/'). The creation of the main page is put under a function in order to reload the genome and pam dropdown
    if a new genome is added, simply by reloading the page.

    ***Returns***

    + **index_page** (*list*): list of html, dcc and dbc components for the layout.
    '''

    final_list = list()

    introduction_content = html.Div(
        [
            # html.Div('CRISPRme is a web application, also available offline or command-line for comprehensive off-target assessment. It integrates human genetic variant datasets with orthogonal genomic annotations to predict and prioritize CRISPR-Cas off-target sites at scale. The method considers both single-nucleotide variants (SNVs) and indels, accounts for bona fide haplotypes, accepts spacer:spacer mismatches and bulges, and is suitable for population and personal genome analyses.'),
            html.Div('CRISPRme is a web application, also available offline or command line, for comprehensive off-target assessment. It integrates human genetic variant datasets with orthogonal genomic annotations to predict and prioritize CRISPR-Cas off-target sites at scale. The method considers both single-nucleotide variants (SNVs) and indels, accounts for bona fide haplotypes, accepts spacer:protospacer mismatches and bulges, and is suitable for population and personal genome analyses.'),
            html.Div(['Check out our preprint on bioRxiv ', html.A(
                'here!', target='_blank', href='https://www.biorxiv.org/content/10.1101/2021.05.20.445054v1')]),
            html.Div(['CRISPRme offline version can be downloaded from ', html.A(
                'Github', target='_blank', href='https://github.com/pinellolab/CRISPRme')]),
            html.Br()
        ]
    )

    white_space_line = html.Div(style={'border-right': 'solid 1px white'})

    modal = html.Div(
        [
            dbc.Modal(
                [
                    dbc.ModalHeader("WARNING! Missing inputs"),
                    dbc.ModalBody(
                        'The following inputs are missing, please select values before submitting the job', id='warning-list'),
                    dbc.ModalFooter(
                        dbc.Button("Close", id="close",
                                   className="modal-button")
                    ),
                ],
                id="modal",
                centered=True
            ),
        ]
    )

    tab_guides_content = html.Div(
        [
            html.H4('Select gRNA'),
            dcc.RadioItems(id="radio-guide",
                           options=[
                               {'label': ' Input individual spacer(s)',
                                'value': 'IP'},
                               {'label': ' Input genomic sequence(s)',
                                'value': 'GS'},
                           ],
                           value='IP',
                           #    style={'margin-bottom': '10px'}
                           ),
            dcc.Textarea(id='text-guides', placeholder='GAGTCCGAGCAGAAGAAGAA\nCCATCGGTGGCCGTTTGCCC', style={
                         'width': '300px', 'height': '30px'}),
            dbc.FormText(
                'Spacer must be provided as a DNA sequence without a PAM. A maximum of 100 spacer sequences can be provided . If using the sequence extraction feature, only the first 100 spacer sequences (starting from the top strand) will be extracted.*', color='secondary')
        ],
        style={'width': '300px'}  # NOTE same as text-area
    )

    cas_protein_content = html.Div(
        [
            html.H4('Select Cas protein'),
            html.Div(
                dcc.Dropdown(options=availableCAS(
                ), clearable=False, id='available-cas', style={'width': '300px'})
            )
        ]
    )

    pam_content = html.Div(
        [
            html.H4('Select PAM'),
            html.Div(
                dcc.Dropdown(
                    options=[], clearable=False, id='available-pam', style={'width': '300px'})
            )
        ],
        # style={'flex': '0 0 100%',
        #        'margin-top': '10%'}
    )

    personal_data_management_content = html.Div(
        [
            html.Br(),
            html.A(html.Button("Personal Data Management", id='add-genome', style={'display': DISPLAY_OFFLINE}),
                   href=URL + '/genome-dictionary-management', target='', style={'text-decoration': 'none', 'color': '#555'})
        ]
    )

    genome_content = html.Div(
        [
            html.H4('Select genome'),
            html.Div(
                # style = {'width':'75%'})
                dcc.Dropdown(options=availableGenomes(),
                             clearable=False, id="available-genome"),
                style={'width': '300px'}
            ),
            html.Div(
                dcc.Checklist(options=[
                    {'label': ' plus 1000 Genome Project variants',
                     'value': '1000G', 'disabled': True},
                    {'label': ' plus HGDP variants',
                     'value': 'HGDP', 'disabled': True},
                    {'label': ' plus personal variants*',
                     'value': 'PV', 'disabled': True}
                ],
                    id='checklist-variants', value=[])
            ),
            html.Div(
                dcc.Dropdown(
                    options=[], id='vcf-dropdown', style={'width': '300px'}, disabled=True),
                id='div-browse-PV',
                # style={'visibility': 'hidden'},
            ),
        ]
    )

    thresholds_content = html.Div(
        [
            html.H4('Select thresholds'),
            html.Div(
                [
                    html.P('Mismatches'),
                    dcc.Dropdown(options=AV_MISMATCHES, clearable=False, id='mms', style={
                        'width': '60px'})
                ],
                style={'display': 'inline-block',
                       "margin-right": "20px"}
            ),
            html.Div(
                [
                    html.P(
                        [
                            'DNA', html.Br(), 'Bulges'
                        ]
                    ),
                    dcc.Dropdown(options=AV_BULGES, clearable=False, id='dna', style={
                        'width': '60px'})
                ],
                style={'display': 'inline-block',
                       "margin-right": "20px"}
            ),
            html.Div
            (
                [
                    html.P(
                        [
                            'RNA', html.Br(), 'Bulges'
                        ]
                    ),
                    dcc.Dropdown(options=AV_BULGES, clearable=False, id='rna', style={
                        'width': '60px'}),
                ],
                style={'display': 'inline-block'}
            ),
        ],
        style={'margin-top': '10%'}
    )

    annotation_content = html.Div(
        [
            html.H4('Select annotation'),
            html.Div(
                dcc.Checklist(options=[
                    {'label': ' ENCODE cCREs + GENCODE gene',
                     'value': 'EN'},
                    {'label': ' Personal annotations*',
                     'value': 'MA', 'disabled': True},
                ],
                    id='checklist-annotations', value=['EN'])
            ),
            html.Div(
                dcc.Dropdown(
                    options=[i for i in get_more_annotations()], id='annotation-dropdown', style={'width': '300px'}, disabled=True),
                id='div-browse-annotation',
                # style={'visibility': 'hidden'},
            )
        ]
    )

    mail_content = html.Div(
        [
            dcc.Checklist(options=[{'label': ' Notify me by email', 'value': 'email', 'disabled': False}],
                          id='checklist-mail', value=[]),
            # dbc.Fade(
            dbc.FormGroup(dbc.Input(type="email", id="example-email", placeholder="name@mail.com",
                                    className='exampleEmail', disabled=True, style={'width': '300px'}))
            #     # id='fade', is_in=False, appear=False
            #     id='fade', is_in=False, appear=True
            # )
        ]
    )

    job_name_content = html.Div(
        [
            dcc.Checklist(options=[{'label': ' Job name', 'value': 'job_name',
                                    'disabled': False}], id='checklist-job-name', value=[]),
            dbc.FormGroup(dbc.Input(type="text", id="job-name", placeholder="my_job",
                                    className='jobName', disabled=True, style={'width': '300px'}))
        ]
    )

    submit_content = html.Div(
        [
            html.Button('Submit', id='check-job',
                        style={'background-color': '#E6E6E6','width': '260px'}),
            html.Button('', id='submit-job', style={'display': 'none'}),
        ]
    )

    example_content = html.Div(
        [
            html.Button('Load Example', id='load-example-button',
                        style={'background-color': '#E6E6E6','width': '260px'}),
        ]
    )

    terms_and_conditions_content = html.Div(
        [
            html.Div('By clicking submit you are agreeing to the'),
            html.Div(html.A('Terms and Conditions!', target='_blank',
                            href='https://github.com/pinellolab/CRISPRme/blob/main/LICENSE'))
        ]
    )
    # insert introduction in the layout
    final_list.append(introduction_content)
    final_list.append(
        html.Div(
            [
                dbc.Row(
                    [
                        dbc.Col(
                            [
                                modal,
                                dbc.Row(
                                    dbc.Col(tab_guides_content)
                                ),
                                dbc.Row(
                                    dbc.Col(cas_protein_content)
                                ),
                                dbc.Row(
                                    dbc.Col(pam_content)
                                )
                                # html.Div(
                                #     [
                                #         tab_guides_content,
                                #         cas_protein_content,
                                #         pam_content
                                #     ],
                                #     id='column-one-step-1',
                                #     # style={'flex': '0 0 30%', 'tex-align': 'center'}
                                #     # style={'tex-align': 'center'}
                                #     # style={'margin': '1px'}
                                #     # style={'margin-left': '-60px'}
                                # ),
                            ],
                            width="auto"
                        ),
                        dbc.Col(
                            [
                                dbc.Row(
                                    dbc.Col(genome_content)
                                ),
                                dbc.Row(
                                    dbc.Col(thresholds_content)
                                ),
                                html.Br(),
                                dbc.Row(
                                    dbc.Col(example_content)
                                )
                                # html.Div(
                                #     [
                                #         genome_content,
                                #         # html.Br(),
                                #         # html.Br(),
                                #         thresholds_content
                                #     ],
                                #     id='column-two-step-2',
                                #     # style={'flex': '0 0 30%', 'tex-align': 'center'}
                                #     # style={'tex-align': 'center'}
                                #     # style={'margin': '1px'}
                                #     # style={'margin-left': '-100px'}
                                # )
                            ],
                            width="auto"
                        ),
                        dbc.Col(
                            [
                                dbc.Row(
                                    annotation_content
                                ),
                                dbc.Row(
                                    mail_content
                                ),
                                dbc.Row(
                                    job_name_content
                                ),
                                # dbc.Row(
                                #     dbc.Col(example_content)
                                # ),
                                html.Br(),
                                dbc.Row(
                                    dbc.Col(submit_content)
                                ),
                                dbc.Row(
                                    terms_and_conditions_content
                                )
                                # html.Div(
                                #     [
                                #         annotation_content,
                                #         html.Br(),
                                #         html.Br(),
                                #         mail_content,
                                #         job_name_content,
                                #         html.Div(example_content, style={
                                #                  'margin-left': '20%'}),
                                #         html.Br(),
                                #         html.Div(submit_content, style={
                                #                  'margin-left': '30%'}),

                                #         terms_and_conditions_content
                                #         # test update in main
                                #     ],
                                #     id='column-three-step-3',
                                #     # style={'flex': '0 0 30%', 'tex-align': 'center'}
                                #     # style={'tex-align': 'center'}
                                #     # style={'margin-right': '40px'}
                                # )
                            ],
                            width="auto"
                        )
                    ],
                    justify="center",
                    style={'margin-bottom': '1%'}
                )
                # html.Div(
                #     [

                # white_space_line,

                # white_space_line,

                #     ],
                #     id='div-steps',
                #     className='flex-div-steps',
                #     # style={'margin-left': '10%'}
                # ),
            ],
            style={'background-color': 'rgba(157, 195, 230, 0.39)', 'border-radius': '5px',
                   #    'border': '1px solid black', 'margin-left': '5%', 'margin-right': '5%'},
                   'border': '1px solid black'},
            # style={'background-color': 'rgba(157, 195, 230, 0.39)', 'border-radius': '5px',
            #        'border': '1px solid black'},
            id='steps-background'
        )
    )
    final_list.append(html.Br())
    final_list.append(
        html.P('*The offline version of CRISPRme can be downloaded from GitHub and offers additional functionalities, including the option to input personal data (such as genetic variants, annotations, and/or empirical off-target results) as well as custom PAMs and genomes. There is no limit on the number of spacers, mismatches, and/or bulges used in the offline search.'))
    # final_list.append(html.P(
    #     '[1] Cancellieri, Samuele, et al. \"Crispritz: rapid, high-throughput, and variant-aware in silico off-target site identification for crispr genome editing.\" Bioinformatics (2019).'))
    # final_list.append(
    #     html.P(['Download CRISPRitz here: ', html.A('InfOmics/CRISPRitz', href='https://github.com/InfOmics/CRISPRitz',
    #                                                 target="_blank"), ' or ', html.A('Pinellolab/CRISPRitz', href='https://github.com/pinellolab/CRISPRitz', target="_blank")])
    # )
    index_page = html.Div(final_list, style={'margin': '1%'})
    return index_page
