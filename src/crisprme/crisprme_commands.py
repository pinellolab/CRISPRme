"""Commands classes definiton.

CRISPRme provides 5 commands:
* complete-search
* gnomAD-converter
* targets-integration
* web-interface
* generate-personal-card

For each available command, we define the corresponding class. Each class stores
the arguments required to run the corresponding CRISPRme functionality.
"""


from crisprme.utils import exception_handler

from typing import List

import os


class CRISPRmeCommand(object):
    """CRISPRme command base class.

    ...

    Attributes
    ----------
    _threads : int
    _verbose : bool
    _debug : bool

    Methods
    -------
    None
    """

    def __init__(self, threads: int, verbose: bool, debug: bool) -> None:
        self._threads = threads
        self._verbose = verbose
        self._debug = debug

    def _get_threads(self):
        return self._threads
    
    @property
    def threads(self):
        return self._get_threads()

    def _get_verbose(self):
        return self._verbose
    
    @property
    def verbose(self):
        return self._get_verbose()

    def _get_debug(self):
        return self._debug

    @property
    def debug(self):
        return self._get_debug()


class CompleteSearch(CRISPRmeCommand):
    """Complete search command class. The class extends `CRISPRmeCommand` class.

    ...

    Attributes
    ----------
    _ref_genome : str
    _search_index : str
    _genome_index : str
    _pam_seq : str
    _bmax : int
    _mm : int
    _bdna : int
    _brna : int
    _annotation_file : str
    _nuclease : str
    _ref_comparison : bool
    _outdir : str
    _guides : List

    Methods
    -------
    write_params_file()
    """

    def __init__(
        self, 
        threads : int,
        verbose : bool,
        debug : bool,
        genome: str,
        ref_genome: str,
        search_index: bool,
        genome_index: str,
        vcf: str, 
        guides: str,
        pam_seq: str,
        pam_full: str,
        pam_start: int,
        bmax: int,
        mm: int,
        bdna: int,
        brna: int,
        annotation_file: str,
        samples_file: str,
        nuclease: str,
        ref_comparison: bool,
        outname: str,
        outdir: str,
        merge_thresh: int,
    ) -> None:
        super().__init__(threads, verbose, debug)  # initialize parent class
        self._genome = genome
        self._ref_genome = ref_genome
        self._search_index = search_index
        self._genome_index = genome_index
        self._vcf = vcf
        self._guides = guides
        self._pam_seq = pam_seq
        self._pam_full = pam_full
        self._pam_start = pam_start
        self._bmax = bmax
        self._mm = mm
        self._bdna = bdna
        self._brna = brna
        self._annotation_file = annotation_file
        self._samples_file = samples_file
        self._nuclease = nuclease
        self._ref_comparison = ref_comparison
        self._outname = outname
        self._outdir = outdir
        self.merge_thresh = merge_thresh

    def set_guides(self, guides: List[str]) -> None:
        """Set _guides attribute.

        ...

        Parameters
        ----------
        guides : List
            Guides

        Returns
        -------
        None 
        """
        if not isinstance(guides, list):
            exception_handler(
                TypeError,
                f"Expected {list.__name__}, got {type(guides).__name__}",
                True
            )
        self._guides = guides

    def set_mail(self, mail_address: str) -> None:
        """Set mail_address (never used on command-line version).

        ...

        Parameters
        ----------
        mail_address : str
            Mail address

        Returns
        -------
        None
        """
        if not isinstance(mail_address, str):
            exception_handler(
                TypeError,
                f"Expected {str.__name__}, got {type(mail_address).__name__}",
                True
            )
        self._mail_address = mail_address

    def _get_outname(self):
        return self._outname
    
    @property 
    def outname(self):
        return self._get_outname()

    def _get_vcf(self):
        return self._vcf

    @property
    def vcf(self):
        return self._get_vcf()

    def _get_genome(self):
        return self._genome

    @property
    def genome(self):
        return self._get_genome()

    def _get_PAM(self):
        return self._pam_seq, self._pam_full, self._pam_start

    @property
    def pam(self):
        return self._get_PAM()

    def write_params_file(self) -> None:
        """Write complete search Paramaters in a TXT file.
        
        ...

        Parameters
        ----------
        self

        Returns
        -------
        None
        """

        try:
            with open(
                os.path.join(self._outdir, ".Params.txt"), mode="w"
            ) as handle:
                handle.write(
                    f"Genome_selected\t{self._ref_genome.replace(' ', '_')}\n"
                )
                handle.write(f"Genome_ref\t{self._ref_genome}\n")
                if self._search_index:
                    handle.write(f"Genome_idx\t{self._genome_index}\n")
                else:
                    handle.write(f"Genome_idx\tNone\n")
                handle.write(f"Pam\t{self._pam_seq}\n")
                handle.write(f"Max_bulges\t{self._bmax}\n")
                handle.write(f"Mismatches\t{self._mm}\n")
                handle.write(f"DNA\t{self._bdna}\n")
                handle.write(f"RNA\t{self._brna}\n")
                handle.write(f"Annotation\t{self._annotation_file}\n")
                handle.write(f"Nuclease\t{self._nuclease}\n")
                handle.write(f"Ref_comp\t{self._ref_comparison}\n")
        except:
            # aleays better to trace this kind of errors ;)
            exception_handler(OSError, "Unable to write '.Params.txt'", True)
        finally:
            handle.close()  # close channel


class GnomADConverter(CRISPRmeCommand):
    """gnomAD VCF converter command class. The class extends `CRISPRmeCommand` 
    class.

    ...

    Attributes
    ----------
    _vcf : str
    _samples : str

    Methods
    -------
    """

    def __init__(
        self, threads: int, verbose: bool, debug: bool, vcf: str, samples: str
    ) -> None:
        super().__init__(threads, verbose, debug)
        self._vcf = vcf
        self._samples = samples

    def _get_vcf(self):
        return self._vcf
    
    @property
    def vcf(self):
        return self._get_vcf()

    def _get_samples(self):
        return self._samples

    @property
    def samples(self):
        return self._get_samples()


class TargetsIntegration(CRISPRmeCommand):
    """Targets integration command class. The class extends `CRISPRmeCommand` 
    class.

    ...

    Attributes
    ----------
    _vcf : str
    _samples : str

    Methods
    -----
    """

    def __init__(
        self, 
        threads: int, 
        verbose: bool, 
        debug: bool, 
        targets_file: str, 
        empirical_data: str, 
        outdir: str
    ) -> None:
        super().__init__(threads, verbose, debug)
        self._targets_file = targets_file
        self._empirical_data = empirical_data
        self._outdir = outdir

    def _get_targets_file(self):
        return self._targets_file

    @property
    def targets_file(self):
        return self._get_targets_file()

    def _get_empirical_data(self):
        return self._empirical_data
    
    @property
    def empirical_data(self):
        return self._get_empirical_data()

    def _get_outdir(self):
        return self._outdir
    
    @property
    def outdir(self):
        return self._get_outdir()
        

class WebInterface(CRISPRmeCommand):
    """Web interface command class. The class extends `CRISPRmeCommand` 
    class.

    ...

    Attributes
    ----------
    _help_web : bool
    
    Methods
    -----
    """
    
    def __init__(
        self, threads: int, verbose: bool, debug: bool, help_web: bool
    ) -> None:
        super().__init__(threads, verbose, debug)
        self._help_web = help_web
    
    def _get_help_web(self):
        return self._help_web

    @property
    def help_web(self):
        return self._get_help_web()


class GeneratePersonalCard(CRISPRmeCommand):
    """Generate personal card command class. The class extends `CRISPRmeCommand` 
    class.

    ...

    Attributes
    ----------
    _inputdir : str
    _guide_seq : str
    _sample_id : str
    
    Methods
    -----
    """

    def __init__(
        self, 
        threads: int, 
        verbose: bool, 
        debug: bool, 
        inputdir: str, 
        guide_seq: str,
        sample_id: str
    ) -> None:
        super().__init__(threads, verbose, debug)
        self._inputdir = inputdir
        self._guide_seq = guide_seq
        self._sample_id = sample_id

    def _get_inputdir(self):
        return self._inputdir

    @property
    def inputdir(self):
        return self._get_inputdir()
    
    def _get_guide_seq(self):
        return self._guide_seq

    @property
    def guide_seq(self):
        return self._get_guide_seq()

    def _get_sample_id(self):
        return self._sample_id

    @property
    def sample_id(self):
        return self._get_sample_id()