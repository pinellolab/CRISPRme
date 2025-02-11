"""Display CRISPRme Manual page. 

The CRISPRme manual can be accessed by the user through the tooglebar placed on
top of the CRISPRme web app.

The manual contains the instructions on CRISPRme web app usage, and the 
explanantion of the different results and plots returned/displayed throughout
the webapp. 
"""

from .pages_utils import ASSETS_DIR, MANUAL_IMGS
from app import app_directory

import dash_html_components as html

import base64  # for decoding upload content
import os


def load_html_image(image_fname: str, width: str, height: str) -> html.Ul:
    """Load an image and encode it for HTML display.

    This function loads an image from the specified file, encodes it in base64 format,
    and creates an HTML Img element with the specified width and height.

    Args:
        image_fname: The path to the image file.
        width: The width of the image.
        height: The height of the image.

    Returns:
        An HTML Ul element containing the encoded image.
    """
    style = {
        "padding-left": "500px",
        "padding-top": "15px",
        "padding-bottom": "15px",
    }  # image styling
    return html.Ul(
        html.Img(
            src="data:image/png;base64,{}".format(
                base64.b64encode(open(image_fname, mode="rb").read()).decode(),
            ),
            style={"height": height, "width": width, "border": "2px solid #555"},
        ),
        style=style,
    )


def about_div() -> html.Div:
    """Create the about section of the help page.

    This function creates an HTML div containing information about CRISPRme,
    including its availability as a web app and command-line package, and a brief
    overview of the required input parameters for an online search.

    Returns:
        An HTML div containing the about section.
    """
    return html.Div(
        [
            html.H3("About"),
            html.P(
                [
                    "CRISPRme is available as both an online web app (",
                    html.A(
                        "http://crisprme.di.univr.it",
                        href="http://crisprme.di.univr.it",
                        target="_blank",
                    ),
                    (
                        ") and a standalone command-line package. To perform an "
                        "online search, users need to provide: gRNA spacer(s), "
                        "Cas protein, PAM sequence, genome build (with or without "
                        "genetic variants from 1000G, HGDP, and/or personal "
                        "datasets), and thresholds for mismatches and RNA/DNA bulges."
                    ),
                ]
            ),
        ]
    )


def homepage_intro() -> html.P:
    """Create the homepage introduction section of the help page.

    This function creates an HTML div containing an introductory paragraph about
    CRISPRme's  interface and customization options, along with an image of the
    main page.

    Returns:
        An HTML P element containing the introduction.
    """
    return html.P(
        [
            html.P(
                (
                    "Performing a search on CRISPRme is quick and intuitive, thanks "
                    "to its user-friendly graphical interface. Users can customize "
                    "off-target searches with various parameters to suit their "
                    "specific needs."
                ),
            ),
            load_html_image(
                os.path.join(app_directory, ASSETS_DIR, MANUAL_IMGS, "mainpage.png"),
                "50%",
                "auto",
            ),
        ]
    )


def homepage_spacer_() -> html.Li:
    """Create the spacer description section for the homepage.

    This function creates an HTML list item containing a description of the gRNA
    spacer sequence, its representation in CRISPRme, and input requirements.

    Returns:
        An HTML Li element containing the spacer description.
    """
    return html.Li(
        [
            html.Strong("Spacer(s): "),
            (
                "The guide RNA (gRNA) spacer sequence is a 20-nucleotide "
                "region that matches the genomic protospacer target, directing "
                "Cas protein binding in the presence of a protospacer adjacent "
                "motif (PAM). In CRISPRme, the spacer sequence is represented "
                "as DNA (rather than RNA) for easy comparison with the aligned "
                "protospacer sequence. CRISPRme accepts up to 100 gRNA spacers, "
                "each listed on a separate line and of the same length. The "
                "input spacer sequence should not include the PAM."
            ),
        ]
    )


def homepage_sequences_() -> html.Li:
    """Create the genomic sequences description section for the homepage.

    This function creates an HTML list item containing a description of the
    supported genomic sequence input formats (BED and FASTA), their
    interpretation in CRISPRme, and input requirements.

    Returns:
        An HTML Li element containing the genomic sequences description.
    """
    return html.Li(
        [
            html.Strong("Genomic sequence(s): "),
            (
                "CRISPRme also supports genomic coordinates in BED format "
                "(chromosome# start end) or DNA sequences in FASTA format "
                "(up to 1,000 characters). BED file coordinates are treated "
                "as 0-based, and in the online version, CRISPRme extracts "
                "the first 100 possible spacer sequences within these regions, "
                "starting from the positive strand. For FASTA input, each entry "
                "must be preceded by a >header."
            ),
        ]
    )


def homepage_pam_() -> html.Li:
    """Create the PAM description section for the homepage.

    This function creates an HTML list item containing a description of the
    protospacer adjacent motif (PAM), its role in CRISPRme, and input
    requirements.

    Returns:
        An HTML Li element containing the PAM description.
    """
    return html.Li(
        [
            html.Strong("PAM sequence: "),
            (
                "The protospacer adjacent motif (PAM) is a short (~2-6 nucleotide) "
                "DNA sequence adjacent to the protospacer, essential for Cas protein "
                "binding to its target. CRISPRme supports a predefined set of PAMs, "
                "requiring users to select one for their search. The software "
                "accommodates both 3' PAMs (e.g., SpCas9) and 5' PAMs (e.g., Cas12a)."
            ),
        ]
    )


def homepage_step1() -> html.Li:
    """Create the step 1 description section for the homepage.

    This function creates an HTML list item containing a description of the first
    step in using CRISPRme, which involves selecting the spacer, Cas protein,
    and PAM sequence.  It combines the spacer, sequence, and PAM descriptions.

    Returns:
        An HTML Li element containing the step 1 description.
    """
    return html.Li(
        [
            html.Strong("STEP 1: Spacer, Cas Protein and PAM selection"),
            html.Ul(
                [
                    homepage_spacer_(),  # spacer description
                    homepage_sequences_(),  # sequences description
                    homepage_pam_(),  # pam description
                ],
                style={"padding-left": "15px"},
            ),
        ]
    )


def homepage_genomes_() -> html.Li:
    """Create the genome builds description section for the homepage.

    This function creates an HTML list item containing a description of the
    supported genome builds in CRISPRme, including default builds, variant
    integration options, and considerations for RNA-targeting strategies.

    Returns:
        An HTML Li element containing the genome builds description.
    """
    return html.Li(
        [
            html.Strong("Genome builds: "),
            (
                "CRISPRme supports genome builds based on FASTA files from UCSC, "
                "allowing users to work with various references, including "
                "transcriptomes, non-human genomes, and cancer genomes. By "
                "default, the hg38 genome build (including mitochondrial DNA) "
                "is available, with the option to incorporate variants from "
                "1000G and/or HGDP. Adding personal variants is only supported "
                "in the local offline and command-line versions. For RNA-targeting "
                "strategies, users can either provide a custom transcriptome or "
                "use a variant-enriched genome, though the latter may miss "
                "off-targets at splice junctions."
            ),
        ]
    )


def homepage_thresholds_() -> html.Li:
    """Create the search thresholds description section for the homepage.

    This function creates an HTML list item containing a description of the
    search thresholds in CRISPRme, including mismatch and bulge tolerances for
    both the web tool and command-line version.

    Returns:
        An HTML Li element containing the search thresholds description.
    """
    return html.Li(
        [
            html.Strong("Search thresholds: "),
            (
                "CRISPRme enables users to define tolerance levels for mismatches, "
                "DNA bulges, and RNA bulges when identifying potential off-targets. "
                "The web tool supports up to 6 mismatches and up to 2 RNA/DNA bulges, "
                "which can be consecutive (NN--NN) or interleaved (NN-N-NN). In the "
                "command-line version, these thresholds are unrestricted and can "
                "be adjusted based on available computational resources."
            ),
        ]
    )


def homepage_baseediting_() -> html.Li:
    """Create the base editing thresholds description section for the homepage.

    This function creates an HTML list item containing a description of the base
    editing thresholds in CRISPRme, including the window for base editing
    susceptibility and the output file generated.

    Returns:
        An HTML Li element containing the base editing description.
    """
    return html.Li(
        [
            html.Strong("Base editing thresholds (optional): "),
            (
                "When a base editor is selected as the Cas protein, CRISPRme "
                "allows users to define a window for base editing susceptibility. "
                "The “Window start” and “Window stop” dropdowns, constrained by "
                "the length of the input guide, specify the region where the "
                "selected “Nucleotide” should be identified within potential "
                "on- and off-targets. The tool generates a comprehensive output "
                "file indicating the base editing susceptibility of candidate "
                "off-target sites."
            ),
        ]
    )


def homepage_step2() -> html.Li:
    """Create the step 2 description section for the homepage.

    This function creates an HTML list item containing a description of the
    second step in using CRISPRme, which involves selecting the genome build
    and configuring search thresholds. It combines the genomes, thresholds, and
    base editing descriptions.

    Returns:
        An HTML Li element containing the step 2 description.
    """
    return html.Li(
        [
            html.Strong("STEP 2: Genome selection and threshold configuration"),
            html.Ul(
                [
                    homepage_genomes_(),  # genome assemblies description
                    homepage_thresholds_(),  # thresholds description
                    homepage_baseediting_(),  # base editing options description
                ],
                style={"padding-left": "15px"},
            ),
        ]
    )


def homepage_annotation_() -> html.Li:
    """Create the functional annotation description section for the homepage.

    This function creates an HTML list item containing a description of the
    functional annotation options in CRISPRme, including the use of ENCODE and
    GENCODE annotations, and the possibility of adding custom annotations in the
    offline version.

    Returns:
        An HTML Li element containing the functional annotation description.
    """
    return html.Li(
        [
            html.Strong("Functional annotation (optional): "),
            (
                "To evaluate the potential impact of off-target activity, CRISPRme "
                "provides functional annotations for both coding and non-coding regions. "
                "These annotations are derived from ENCODE, which includes candidate "
                "cis-regulatory elements, and GENCODE, which covers protein-coding "
                "genes, untranslated regions, and introns. In the offline version, "
                "users can integrate custom genome annotations, such as cell-type-specific "
                "chromatin marks or experimentally identified off-target sites, "
                "by uploading BED files."
            ),
        ]
    )


def homepage_email_() -> html.Li:
    """Create the email notification description section for the homepage.

    This function creates an HTML list item containing a description of the
    optional email notification feature in CRISPRme, which notifies the user
    upon job completion.

    Returns:
        An HTML Li element containing the email notification description.
    """
    return html.Li(
        [
            html.Strong("Email notification (optional): "),
            "If provided, the server notifies the user via email upon job completion.",
        ]
    )


def homepage_jobname_() -> html.Li:
    """Create the job name description section for the homepage.

    This function creates an HTML list item containing a description of the
    optional job name feature in CRISPRme, which allows users to provide a prefix
    for the unique job ID.

    Returns:
        An HTML Li element containing the job name description.
    """
    return html.Li(
        [
            html.Strong("Job name (optional): "),
            (
                "If a job name is provided, it will be used as a prefix for the "
                "unique job ID, making it easier to identify a specific search "
                "(e.g., my_job_G05B8KHU0H)."
            ),
        ]
    )


def homepage_step3() -> html.Li:
    """Create the step 3 description section for the homepage.

    This function creates an HTML list item containing a description of the
    third step in using CRISPRme, which involves optional settings for
    functional annotations, email notification, and job naming. It combines the
    annotation, email, and job name descriptions.

    Returns:
        An HTML Li element containing the step 3 description.
    """
    return html.Li(
        [
            html.Strong("STEP 3 Annotation(s), email notification, and job name"),
            html.Ul(
                [
                    homepage_annotation_(),  # annotation description
                    homepage_email_(),
                    homepage_jobname_(),
                ],
                style={"padding-left": "15px"},
            ),
        ]
    )


def homepage_div() -> html.Div:
    """Create the homepage section of the help page.

    This function creates an HTML div containing a detailed description of the
    CRISPRme homepage, including the input parameters and steps involved in
    performing a search. It combines the descriptions for each step and input
    element.

    Returns:
        An HTML Div element containing the homepage section.
    """
    return html.Div(
        [
            html.H3("Homepage"),
            homepage_intro(),  # homepage introduction
            homepage_step1(),  # step 1 description
            homepage_step2(),  # step 2 description
            homepage_step3(),  # step 3 description
        ]
    )


def loadpage_intro() -> html.P:
    """Create the job status page introduction section of the help page.

    This function creates an HTML div containing an introductory paragraph about
    CRISPRme's job status page, explaining what happens after submitting a search
    and how to access the results. It also includes an image of the job status
    page.

    Returns:
        An HTML P element containing the introduction.
    """
    return html.P(
        [
            html.P(
                (
                    "After selecting the desired inputs, clicking the Submit button "
                    "initiates the search. A progress page will display the search "
                    'status, and once completed, a "View Results" link will appear at '
                    "the bottom of the status report page."
                ),
            ),
            load_html_image(
                os.path.join(app_directory, ASSETS_DIR, MANUAL_IMGS, "jobpage.png"),
                "50%",
                "auto",
            ),
        ]
    )


def loadpage_div() -> html.Div:
    """Create the job status section of the help page.

    This function creates an HTML div containing a description of the CRISPRme
    job status page, including an introduction and an image of the page.

    Returns:
        An HTML Div element containing the job status section.
    """
    return html.Div(
        [
            html.H3("Job Status"),
            loadpage_intro(),  # job status page introduction
        ]
    )


def resultpage_intro() -> html.P:
    """Create the results page introduction section of the help page.

    This function creates an HTML div containing an introductory paragraph about
    CRISPRme's results page, explaining the summary table, CFD specificity score,
    on/off-target counts, and download functionality. It also includes an image
    of the result summary table.

    Returns:
        An HTML P element containing the introduction.
    """
    return html.P(
        [
            html.P(
                [
                    (
                        "CRISPRme summarizes the results in a table, presenting "
                        "each gRNA's CFD specificity score along with the number "
                        "of on-targets and off-targets identified in both the "
                        "reference and variant genomes. These targets are grouped "
                        "based on mismatches and bulges. Notably, the CFD specificity "
                        "score was originally designed for searches with up to 3-4 "
                        "mismatches; as the number of mismatches increases, the "
                        "score declines non-linearly. Therefore, comparisons between "
                        "searches with different mismatch/bulge thresholds or variant "
                        "datasets should be made with caution. At the top of the "
                        "page, a summary table provides an overview of the nuclease, "
                        "CFD specificity score, and target counts across mismatch and "
                        'bulge categories. In the top left corner, a "Download '
                        'General Table" button allows users to export the table '
                        "as a text file."
                    ),
                ]
            ),
            load_html_image(
                os.path.join(
                    app_directory, ASSETS_DIR, MANUAL_IMGS, "resultsummary.png"
                ),
                "50%",
                "auto",
            ),
        ]
    )


def resultspage_table_() -> html.P:
    """Create the results table description section for the results page.

    This function creates an HTML unordered list containing a description of the
    results table columns, including the CFD score, off-targets in the reference
    genome, and off-targets in the variant genome.

    Returns:
        An HTML P element containing the results table description.
    """
    return html.P(
        [
            html.Strong("Results Table"),
            html.Ul(
                [
                    html.Li(
                        [
                            html.Strong("CFD: "),
                            "The Off-Target Cutting Frequency Determination (CFD) "
                            "Score quantifies a guide RNA's affinity for off-target "
                            "sites, indicating the likelihood of unintended cleavage "
                            "events.",
                        ]
                    ),
                    html.Li(
                        [
                            html.Strong(
                                "Off-Targets Reference (0 - n Mismatches + Bulges): "
                            ),
                            "Displays the number of potential off-target sites for the "
                            "guide in the reference genome, grouped by mismatch and bulge "
                            "values.",
                        ]
                    ),
                    html.Li(
                        [
                            html.Strong(
                                "Off-Targets Variant (0 - n Mismatches + Bulges): "
                            ),
                            "Displays the number of potential off-target sites for the "
                            "guide in the variant genome, grouped by mismatch and bulge "
                            "values.",
                        ]
                    ),
                ],
                style={"padding-left": "15px"},
            ),
        ]
    )


def resultspage_custom_ranking_() -> html.P:
    """Create the results page section of the help page.

    This function creates an HTML div containing a detailed description of the
    CRISPRme results page, including an introduction, a description of the results
    table, and descriptions of the interactive reports available for each guide.

    Returns:
        An HTML Div element containing the results page section.
    """
    return html.P(
        [
            html.P(html.Strong("Custom ranking")),
            html.P(
                (
                    "This report allows users to filter and rank potential "
                    "off-targets based on the number of mismatches and bulges, "
                    "CFD score, Risk Score (increase in CFD score due to genetic "
                    "variants), or a combination of these factors."
                ),
            ),
            load_html_image(
                os.path.join(
                    app_directory, ASSETS_DIR, MANUAL_IMGS, "customranking.png"
                ),
                width="50%",
                height="auto",
            ),
        ]
    )


def resultspage_mmbsummary_fields_() -> html.P:
    """Create the fields description for the Summary by Sample section.

    This function creates an HTML paragraph containing a description of the fields
    present in the Summary by Sample table on the results page.  These fields
    include sample demographics, target counts, and PAM creation information.

    Returns:
        An HTML P element containing the field descriptions.
    """
    return html.P(
        [
            html.P("Fields:"),
            html.Ul(
                [
                    html.Li(
                        [
                            html.Strong("Bulge Type: "),
                            "Specifies the bulge type of the targets",
                        ]
                    ),
                    html.Li(
                        [
                            html.Strong("Bulge Size: "),
                            "Indicates the size of the bulge in the targets.",
                        ]
                    ),
                    html.Li(
                        [
                            html.Strong("Mismatches: "),
                            "Indicates the number of mismatches in the targets.",
                        ]
                    ),
                    html.Li(
                        [
                            html.Strong("Targets in Reference: "),
                            (
                                "Displays the number of targets found in the reference "
                                "genome for the specified mismatch/bulge combination."
                            ),
                        ]
                    ),
                    html.Li(
                        [
                            html.Strong("Targets in Sample: "),
                            (
                                "Displays the number of targets found in the "
                                "variant genome for the specified mismatch/bulge "
                                "combination, with each target associated with at "
                                "least one sample."
                            ),
                        ]
                    ),
                    html.Li(
                        [
                            html.Strong("PAM Creation: "),
                            "Indicates the number of potential PAMs created by the addition of variants.",
                        ]
                    ),
                    html.Li(
                        [
                            html.Strong("Show Targets: "),
                            (
                                "Opens a new page displaying all the targets for "
                                "the selected row, as shown in the following image:"
                            ),
                        ]
                    ),
                ],
            ),
            load_html_image(
                os.path.join(
                    app_directory, ASSETS_DIR, MANUAL_IMGS, "summarybyguide_targets.png"
                ),
                width="50%",
                height="auto",
            ),
        ]
    )


def resultspage_mmbsummary_() -> html.P:
    """Create the Summary by Mismatches/Bulges section for the results page.

    This function creates an HTML paragraph describing the Summary by
    Mismatches/Bulges report, which categorizes targets based on type, mismatch
    count, and bulge size. It includes an image of the report and a description
    of its fields.

    Returns:
        An HTML P element containing the Summary by Mismatches/Bulges section.
    """
    return html.P(
        [
            html.P(html.Strong("Summary by Mismatches/Bulges")),
            (
                "This report presents a matrix that categorizes targets based on "
                'target type, mismatch count, and bulge size. "X" targets '
                'contain only mismatches, "DNA" targets include DNA bulges '
                '(with or without mismatches), and "RNA" targets include RNA '
                "bulges (with or without mismatches)."
            ),
            html.P(
                load_html_image(
                    os.path.join(
                        app_directory, ASSETS_DIR, MANUAL_IMGS, "summarybyguide.png"
                    ),
                    width="50%",
                    height="auto",
                )
            ),
            resultspage_mmbsummary_fields_(),  # mm+bulges table fields
        ]
    )


def resultspage_summarysample_fields_() -> html.P:
    """Create the fields description for the Summary by Mismatches/Bulges section.

    This function creates an HTML paragraph containing a description of the fields
    present in the Summary by Mismatches/Bulges table on the results page. These
    fields provide information about bulge type, size, mismatches, target counts
    in reference and variant genomes, PAM creation, and a link to view specific
    targets.

    Returns:
        An HTML P element containing the field descriptions.
    """
    return html.P(
        [
            html.P("Fields:"),
            html.Ul(
                [
                    html.Li([html.Strong("Gender: "), "Sample's gender"]),
                    html.Li(
                        [
                            html.Strong("Population: "),
                            "The population to which the sample belongs.",
                        ]
                    ),
                    html.Li(
                        [
                            html.Strong("Super Population: "),
                            "The superpopulation to which the sample belongs.",
                        ]
                    ),
                    html.Li(
                        [
                            html.Strong("Targets in sample: "),
                            "The number of targets in the variant genome generated by that sample.",
                        ]
                    ),
                    html.Li(
                        [
                            html.Strong("Targets in Population: "),
                            "The number of targets in the variant genome generated by all the samples in the population.",
                        ]
                    ),
                    html.Li(
                        [
                            html.Strong("Targets in Super Population: "),
                            "The number of targets in the variant genome generated by all populations.",
                        ]
                    ),
                    html.Li(
                        [
                            html.Strong("PAM Creation: "),
                            "Number of potential PAMs generated by the addition of variants.",
                        ]
                    ),
                    html.Li(
                        [
                            html.Strong("Show Targets: "),
                            "Opens a new page to display all the targets for the selected "
                            "row, as illustrated in the following image:",
                        ]
                    ),
                ],
            ),
            load_html_image(
                os.path.join(
                    app_directory,
                    ASSETS_DIR,
                    MANUAL_IMGS,
                    "summarybysample_targets.png",
                ),
                width="50%",
                height="auto",
            ),
        ]
    )


def resultspage_summarysample_() -> html.P:
    """Create the Summary by Sample section for the results page.

    This function creates an HTML paragraph describing the Summary by Sample report,
    which displays targets associated with each sample in the provided VCFs. It
    includes an image of the report and a description of its fields.

    Returns:
        An HTML P element containing the Summary by Sample section.
    """
    return html.P(
        [
            html.P(html.Strong("Summary by Sample")),
            "This page displays all the samples from the VCFs and enables users "
            "to extract and visualize targets associated with each sample.",
            html.P(
                load_html_image(
                    os.path.join(
                        app_directory, ASSETS_DIR, MANUAL_IMGS, "summarybysample.png"
                    ),
                    width="50%",
                    height="auto",
                )
            ),
            resultspage_summarysample_fields_(),  # summary by samples table fields
        ]
    )


def resultspage_genomicregions_() -> html.P:
    """Create the Query Genomic Regions section for the results page.

    This function creates an HTML paragraph describing the Query Genomic Regions
    report, which allows users to retrieve targets overlapping a specific genomic
    region. It includes an image of the report.

    Returns:
        An HTML P element containing the Query Genomic Regions section.
    """
    return html.P(
        [
            html.P(html.Strong("Query Genomic Regions")),
            "This page allows users to retrieve targets that overlap a specific "
            "genomic region, making it easy to quickly assess potential off-targets "
            "within a regulatory element or coding region.",
            html.P(
                load_html_image(
                    os.path.join(
                        app_directory, ASSETS_DIR, MANUAL_IMGS, "summarybyposition.png"
                    ),
                    width="50%",
                    height="auto",
                )
            ),
        ]
    )


def resultspage_graphicalreports_list_() -> html.P:
    """Create the list of graphical reports descriptions for the results page.

    This function creates an HTML paragraph containing a description of each
    graphical report available on the results page, including stem plots, bar
    plots, radar charts, and motif logos. It also includes images illustrating
    these reports.

    Returns:
        An HTML P element containing the graphical reports descriptions.
    """
    return html.P(
        [
            html.Li(
                "A stem plot demonstrates how genetic variants impact predicted "
                "off-target potential. The arrow connecting the red (reference "
                "allele off-target) and blue (alternative allele off-target) dots "
                "indicates the increase in predicted cleavage potential caused by "
                "the variant."
            ),
            load_html_image(
                os.path.join(app_directory, ASSETS_DIR, MANUAL_IMGS, "lolliplot.png"),
                width="50%",
                height="auto",
            ),
            html.Li(
                "Bar plots illustrate the distribution of candidate off-targets "
                "across super-populations, categorized by the number of mismatches "
                "and bulges."
            ),
            html.Li(
                "A radar chart displays the relative specificity of the analyzed "
                "guide for each functional region, based on annotations from GENCODE "
                "and ENCODE. A larger area on the chart indicates a gRNA with more "
                "potential off-targets within annotated regions, which may suggest "
                "an undesirable outcome. A summary table accompanies the chart, providing "
                "the count and percentage of off-targets with each annotation."
            ),
            html.Li(
                "A motif logo summarizing the frequency of mismatches and bulges (b) "
                "at each base of the protospacer and PAM."
            ),
            load_html_image(
                os.path.join(
                    app_directory, ASSETS_DIR, MANUAL_IMGS, "barplotradarmotif.png"
                ),
                width="50%",
                height="auto",
            ),
        ]
    )


def resultspage_graphicalreports_() -> html.P:
    """Create the list of graphical reports descriptions for the results page.

    This function creates an HTML paragraph containing a description of each
    graphical report available on the results page, including stem plots, bar
    plots, radar charts, and motif logos. It also includes images illustrating
    these reports.

    Returns:
        An HTML P element containing the graphical reports descriptions.
    """
    return html.P(
        [
            html.P(html.Strong("Graphical Reports")),
            "This page generates several graphical reports for each selected sgRNA:",
            resultspage_graphicalreports_list_(),  # graphical reports
        ]
    )


def resultspage_personalrisk_list_() -> html.P:
    """Create the fields description for the Personal Risk Card section.

    This function creates an HTML paragraph containing a description of the fields
    present in the Personal Risk Card table on the results page. These fields
    include counts of personal, PAM creation, and private off-targets.

    Returns:
        An HTML P element containing the field descriptions.
    """
    return html.P(
        [
            html.Li(
                [
                    html.Strong("Personal: "),
                    "The count of all candidate variant off-targets for the selected "
                    "sample, including both unique and shared variants.",
                ]
            ),
            html.Li(
                [
                    html.Strong("PAM creation: "),
                    "The count of instances where a genetic variant in the sample "
                    "creates a new PAM.",
                ]
            ),
            html.Li(
                [
                    html.Strong("Private: "),
                    "The count of candidate variant off-targets found exclusively "
                    "in the selected sample.",
                ]
            ),
            load_html_image(
                os.path.join(
                    app_directory, ASSETS_DIR, MANUAL_IMGS, "personalcard.png"
                ),
                width="50%",
                height="auto",
            ),
        ]
    )


def resultspage_personalrisk_() -> html.P:
    """Create the Personal Risk Card section for the results page.

    This function creates an HTML paragraph describing the Personal Risk Card report,
    which summarizes potential off-target effects of a specific gRNA in an
    individual due to genetic variants. It highlights the importance of
    variant-aware off-target assessments.  It includes a description of the report
    contents, including plots and tables.

    Returns:
        An HTML P element containing the Personal Risk Card section.
    """
    return html.P(
        [
            html.P(html.Strong("Personal Risk Card")),
            "CRISPRme offers a dedicated page for generating Personal Risk Cards "
            "reports, which summarize potential off-target editing by a specific "
            "gRNA in an individual due to genetic variants. This feature is especially "
            "useful for identifying and investigating private off-targets. The "
            "report includes two dynamically generated plots that depict all "
            "candidate variant off-targets for the sample, including those unique "
            "to the individual and those shared with others. These plots highlight "
            "how genetic variants influence predicted off-target cleavage potential, "
            "emphasizing the importance of variant-aware off-target assessments, "
            "as demonstrated by CRISPRme. The report also includes two tables: a "
            "summary table (Table 1, top) and detailed information about each "
            "extracted candidate off-target (Table 2, bottom), with the following "
            "columns:",
            resultspage_personalrisk_list_(),
        ]
    )


def resultspage_tabs_() -> html.P:
    """Create the interactive reports tabs description section for the results page.

    This function creates an HTML paragraph describing the six interactive reports
    available on the CRISPRme results page: Custom Ranking, Summary by
    Mismatches/Bulges, Summary by Sample, Query Genomic Region, Graphical Reports,
    and Personal Risk Cards. It lists each report and provides a brief
    explanation of their purpose.

    Returns:
        An HTML P element containing the interactive reports tabs description.
    """
    return html.P(
        [
            html.P(
                "Additionally, six interactive reports are generated for each guide "
                "and available for download: Custom Ranking, Summary by Mismatches/Bulges, "
                "Summary by Sample, Query Genomic Region, Graphical Reports, and Personal "
                "Risk Cards."
            ),
            html.Ul(
                [
                    resultspage_custom_ranking_(),  # custom rank page description
                    resultspage_mmbsummary_(),  # summary by mm+bulges description
                    resultspage_summarysample_(),  # summary by sample description
                    resultspage_genomicregions_(),  # genomic regions queries description
                    resultspage_graphicalreports_(),  # graphical reports description
                    resultspage_personalrisk_(),  # personal risk description
                ],
            ),
        ]
    )


def resultspage_div() -> html.Div:
    """Create the results page section of the help page.

    This function creates an HTML div containing a detailed description of the
    CRISPRme results page, including an introduction, a description of the results
    table, and descriptions of the interactive reports available for each guide.

    Returns:
        An HTML Div element containing the results page section.
    """
    return html.Div(
        [
            html.H3("Results Page"),
            resultpage_intro(),  # results page introduction
            resultspage_table_(),  # results table description
            resultspage_tabs_(),
        ]
    )


def helpPage() -> html.Div:
    """Create the main help page for the CRISPRme web application.

    This function assembles the various sections of the help page, including
    the about section, homepage description, job status description, and results
    page description. It returns an HTML Div element containing the complete
    help page content.

    Returns:
        An HTML Div element representing the help page.
    """
    # construct manual page for dipsly on website/gui
    html_divs = [
        about_div(),  # about section
        homepage_div(),  # homepage description section
        loadpage_div(),  # job status description section
        resultspage_div(),  # results page description section
    ]
    return html.Div(html_divs, style={"margin": "1%"})
