""" """

from typing import List, Dict, Tuple
from pysam.utils import SamtoolsError

import warnings
import random
import pysam
import json
import sys
import os

# dismiss warning messages
warnings.filterwarnings("ignore")
# set random seed
random.seed(a=None, version=2)
# input target file selection criteria
CRITERIA = {"CFD", "CRISTA", "fewest"}

class RadarChartInputArgs:    
    def __init__(self, args: List[str]):
        self._validate_arg_count(args)  # validate number of input arguments
        self._guide_fname = args[0]  # guides file
        self._report_fname = args[1]  # annotated report file
        self._retrieve_annotation_files(args[2])  # retrieve annotation files
        self._retrieve_annotation_colnames(args[3])  # retrieve annotation colnames
        self._output_dir = args[4]  # output results folder
        self._criterion = args[5]  # targets selection criteria
        self._validate_arguments()  # validate input args
    
    def _validate_arg_count(self, args: List[str]) -> None:
        if len(args) != 6:
            raise ValueError(f"Too many/few input arguments ({len(args)})")
        
    def _retrieve_annotation_files(self, annotations: str) -> None:
        self._annotations = [fname for fname in annotations.split(",")]

    def _retrieve_annotation_colnames(self, annotation_colnames: str) -> None:
        self._annotation_colnames = [colname for colname in annotation_colnames.split(",")]
        
    def _validate_arguments(self):        
        for fname in [self._guide_fname, self._report_fname]:  
            if not os.path.isfile(fname):  # check guides and report existance
                raise FileNotFoundError(f"Cannot find {fname}")
        for ann_fname in self._annotations:  # check annotation files existance
            if not os.path.isfile(ann_fname):
                raise FileNotFoundError(f"Cannot find annotation file: {ann_fname}")
        if not os.path.isdir(self._output_dir):  # check output folder
            raise FileNotFoundError(f"Cannot find output directory: {self._output_dir}")
        # Validate selection criteria
        if self._criterion not in CRITERIA:
            raise ValueError(f"Forbidden criterion: {self._criterion}")

    @property
    def guide_fname(self) -> str:
        return self._guide_fname
    
    @property
    def report_fname(self) -> str:
        return self._report_fname
    
    @property
    def annotations(self) -> List[str]:
        return self._annotations
    
    @property
    def annotation_colnames(self) -> List[str]:
        return self._annotation_colnames
    
    @property
    def output_dir(self) -> str:
        return self._output_dir
    
    @property
    def criterion(self) -> str:
        return self._criterion
    

def parse_commandline(args: List[str]) -> RadarChartInputArgs:
    # parse command line input arguments
    return RadarChartInputArgs(args)

def _remove_colname_id(term: str, colnames: List[str]) -> str:
    for cname in colnames:
        term = term.replace(f"__{cname}", "")
    return term

def retrieve_annotation_terms(annotations: List[str], colnames: List[str]) -> List[str]:
    terms = set()  # annotation terms set
    if len(annotations) == 1 and annotations[0] == "empty.txt":  # no annotation
        return list(terms) 
    for annotation in annotations:  # iterate over each input annotation file
        annotation_bed = pysam.TabixFile(annotation)
        try:
            terms = terms.union({_remove_colname_id(line.strip().split()[3], colnames) for line in annotation_bed.fetch()})
        except (SamtoolsError, Exception) as e:
            raise OSError(f"Failed retrieving annotation terms on {annotation}") from e
    return sorted(terms)

def retrieve_guides(guides_fname: str) -> List[str]:
    try:
        with open(guides_fname, mode="r") as infile:
            return [line.strip() for line in infile]
    except (IOError, Exception) as e:
        raise OSError(f"Failed retrieving guides") from e
    

def create_guide_dict(annotation_terms: List[str]) -> Dict[int, Dict[str, int]]:
    return {
        i : {term: 0 for term in ["General"] + annotation_terms}
        for i in range(15)
    }

def create_motif_dict(guidelen: int) -> Dict[int, Dict[str, List[int]]]:
    return {
        i: {e: [0] * guidelen for e in ["A", "C", "G", "T", "RNA", "DNA"]}
        for i in range(15)
    }

def _update_guide_dict(guide_dict: Dict[int, Dict[str, int]], fields: List[str], colnames: List[str], idx: int) -> Dict[int, Dict[str, int]]:
    guide_dict[idx]["General"] += 1
    if fields[14] != "NA":  # annotation available for target
        for annotation in set(fields[14].strip().split(",")):  # retrieve annotations
            annotation = _remove_colname_id(annotation, colnames)  # remove colname id
            guide_dict[idx][annotation] += 1
    return guide_dict


def _update_motif_dict(motif_dict: Dict[int, Dict[str, List[int]]], guide: str, fields: List[str], idx: int) -> Dict[int, Dict[str, List[int]]]:
    algndguide, algndsequence = fields[1], fields[2]  # aligned sequences
    mm, bulge = int(fields[8]), int(fields[9])  # mismatches and bulges
    bulge_type = fields[0]  # bulge type (DNA, RNA, X)
    if bulge_type == "DNA":
        for i, nt in enumerate(algndguide[bulge:]):
            if nt == "-":  # bulge in guide
                motif_dict[idx][bulge_type][i] += 1
        for i, nt in enumerate(algndsequence[bulge:]):
            if nt.islower():  # mismatch
                motif_dict[idx][nt.upper()][i] += 1
            if guide[0] != "N" and guide[i] == "N":
                motif_dict[idx][nt.upper()][i] += 1
        for i, ns in enumerate(guide):
            if ns != "N":
                break
            motif_dict[idx][algndsequence[i].upper()][i] += 1
    elif bulge_type in ["RNA", "X"]:
        for i, nt in enumerate(algndsequence):
            if nt.islower(): 
                motif_dict[idx][nt.upper()][i] += 1
            elif nt == "-":
                motif_dict[idx][bulge_type][i] += 1
            if guide[i] == "N":  # in pam
                motif_dict[idx][nt.upper()][i] += 1      
    return motif_dict  


def fill_dictionary(report_fname: str, guide: str, annotation_colnames: List[str], guide_dict: Dict[int, Dict[str, int]], motif_dict: Dict[int, Dict[str, List[int]]]) -> Tuple[Dict[int, Dict[str, int]], Dict[int, Dict[str, List[int]]]]:
    try:
        with open(report_fname, mode="r") as infile:
            for line in infile:
                if line.startswith("#"):  # skip header
                    continue
                fields = line.strip().split()  # split report fields
                if guide not in fields[15]:  # skip if not current guide
                    continue
                for i in range(int(fields[8]) + int(fields[9]), 15):
                    guide_dict = _update_guide_dict(guide_dict, fields, annotation_colnames, i)
                    motif_dict = _update_motif_dict(motif_dict, guide, fields, i)
    except (IOError, Exception) as e:
        raise OSError(f"An error occurred while reading off-target file: {report_fname}") from e
    return guide_dict, motif_dict


def main() -> None:
    inargs = parse_commandline(sys.argv[1:])  # read command line arguments
    annterms = retrieve_annotation_terms(inargs.annotations, inargs.annotation_colnames)  # retrieve annotation terms
    guides = retrieve_guides(inargs.guide_fname)  # retrieve guides
    for guide in guides:   # construct guide and motif dict for each input guide
        guide_dict = create_guide_dict(annterms)
        motif_dict = create_motif_dict(len(guide))
        guide_dict, motif_dict = fill_dictionary(inargs.report_fname, guide, inargs.annotation_colnames, guide_dict, motif_dict)
        try:
            with open(os.path.join(inargs.output_dir, f".guide_dict_{guide}_{inargs.criterion}.json"), mode="w") as outfile:
                json.dump(guide_dict, outfile)
        except (IOError, Exception) as e:
            raise OSError(f"Failed writing guide dictionary JSON file") from e
        try:
            with open(os.path.join(inargs.output_dir, f".motif_dict_{guide}_{inargs.criterion}.json"), mode="w") as outfile:
                json.dump(motif_dict, outfile)
        except (IOError, Exception) as e:
            raise OSError(f"Failed writing motif dictionary JSON file") from e


if __name__ == "__main__":
    main()
