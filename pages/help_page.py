"""Display CRISPRme Manual page. 

The CRISPRme manual can be accessed by the user through the tooglebar placed on
top of the CRISPRme web app.

The manual contains the instructions on CRISPRme web app usage, and the 
explanantion of the different results and plots returned/displayed throughout
the webapp. 
"""


from .results_page_utils import ASSETS_DIR
from app import app_main_directory

import dash_html_components as html

import base64  # for decoding upload content
import os


def helpPage() -> html.Div:
    """Construct the Manual web page.

    ...

    Parameters
    ----------
    None

    Returns
    -------
    html.Div
    """

    # begin web page construction
    final_list = []
    final_list.append(
        html.Div(
            [
                html.H3("About"),
                html.P(
                    [
                        "CRISPRme is available as an online web app at ",
                        html.A(
                            "http://crisprme.di.univr.it",
                            href='http://crisprme.di.univr.it', 
                            target="_blank"
                        ),
                        str(
                            " a standalone command line package. The required "
                            "inputs to perform an online search are: gRNA "
                            "spacer(s), Cas protein, PAM sequence, genome build "
                            "with or without the inclusion of genetic variants "
                            "(1000G, HGDP and/or personal variants), and thresholds "
                            "of mismatches and RNA/DNA bulges."
                        )
                    ]
                )
            ]
        )
    )
    # main page description
    final_list.append(html.H3("Main Page"))
    final_list.append(
        html.Div(
            [
                html.P(
                    [
                        str(
                            "A search on CRISPRme can be performed in three "
                            "simple steps thanks to the user-friendly user "
                            "interface. Several options are available to "
                            "personalize a search."
                        ),
                        html.Img(
                            src="data:image/png;base64,{}".format(
                                base64.b64encode(
                                    open(
                                        os.path.join(
                                            app_main_directory, 
                                            ASSETS_DIR, 
                                            "main_page.png"
                                        ),
                                        mode="rb"
                                    ).read()
                                ).decode(),
                            ),
                            width="100%",
                            height="auto"
                        )
                    ]
                ),
                # describe analysis steps
                html.Ul(
                    [
                        # first step of CRISPRme analysis
                        html.Li(
                            [
                                "STEP 1: Spacer, Cas Protein and PAM selection",
                                html.Ul(
                                    [
                                        html.Li(
                                            str(
                                                "Spacer(s): The guide RNA (gRNA) "
                                                "spacer sequence matches the "
                                                "genomic target protospacer "
                                                "sequence (typically 20 nucleotides) " 
                                                "and directs Cas protein binding "
                                                "to the protospacer in the presence "
                                                "of a protospacer adjacent motif "
                                                "(PAM). The spacer sequence is "
                                                "represented as DNA (rather than "
                                                "RNA) in CRISPRme to allow easy "
                                                "comparison to the aligned protospacer "
                                                "sequence. CRISPRme accepts a set "
                                                "of gRNA spacer(s), one per line, "
                                                "each with the same length (max "
                                                "100 sequences). The input spacer "
                                                "sequence should not include PAM."
                                            )
                                        ),
                                        html.Li(
                                            str(
                                                "Genomic sequence(s): CRISPRme can "
                                                "alternatively take as input a set "
                                                "of genomic coordinates in BED format "
                                                "(chromosome# start end) or DNA sequences " 
                                                "in FASTA format (max 1000 characters). "
                                                "The BED file coordinates will be treated "
                                                "as 0-based and CRISPRme (online version) "
                                                "will extract the first 100 possible spacer "
                                                "sequences within these coordinates starting " 
                                                "with the positive strand. To use this type of "
                                                "input, the user must delimit each entry with "
                                                "a >header."
                                            )
                                        ),
                                        html.Li(
                                            str(
                                                "PAM sequence: The PAM is a short "
                                                "(∼2-6 nucleotide) DNA sequence "
                                                "adjacent to the protospacer necessary "
                                                "for the Cas protein to bind to a "
                                                "specific DNA target. CRISPRme "
                                                "supports a set of PAMs and users "
                                                "must select one of them in order "
                                                "to perform the search. The software "
                                                "supports both 3’ (e.g. SpCas9) and 5’ "
                                                "(e.g. Cas12a) PAM sequences."
                                            )
                                        ),
                                    ], 
                                    style={'padding': '15px'}
                                )
                            ]
                        ),
                        # second step of CRISPRme analysis
                        html.Li(
                            [
                                "STEP 2: Genome selection and threshold configuration",
                                html.Ul(
                                    [
                                        html.Li(
                                            str(
                                                "Genome builds: The genome builds "
                                                "are based on FASTA files from UCSC. "
                                                "The hg38 and hg19 genomic builds are "
                                                "available with the option to incorporate "
                                                "variants from 1000G, HGDP, and/or personal "
                                                "variants in the search. The option to add "
                                                "personal variants is enabled only for the "
                                                "local offline and command line versions."
                                            )
                                        ),
                                        html.Li(
                                            str(
                                                "Search thresholds: CRISPRme allows "
                                                "users to specify the number of "
                                                "mismatches, DNA and RNA bulges "
                                                "tolerated in enumerating potential "
                                                "off-targets. The web-tool allows up "
                                                "to 6 mismatches and up to 2 RNA/DNA "
                                                "bulges (which can be consecutive "
                                                "(NN--NN) or interleaved (NN-N-NN)). "
                                                "However, for the command line version, "
                                                "these thresholds can be set freely and "
                                                "depend only on the available computational "
                                                "resources."
                                            )
                                        ),
                                    ], 
                                    style={'padding': '15px'}
                                )
                            ]
                        ),
                        # third step of CRISPRme analysis
                        html.Li(
                            [
                                "STEP 3 Annotation(s), email notification, and job name",
                                html.Ul(
                                    [
                                        html.Li(
                                            str(
                                                "Functional annotations: To assess "
                                                "the potential impact of off-target "
                                                "activity, CRISPRme provides a set "
                                                "of functional annotations for coding "
                                                "and non-coding regions. The annotations "
                                                "are based on files obtained from the "
                                                "Encyclopedia of DNA Elements (ENCODE) "
                                                "containing candidate cis regulatory "
                                                "elements21 and from GENCODE25 containing "
                                                "annotations for protein coding genes, "
                                                "transcribed but untranslated regions, and "
                                                "introns. In the offline versions of CRISPRme, "
                                                "users can add custom genome annotations, such "
                                                "as cell-type specific chromatin marks or "
                                                "off-target sites nominated by in vitro and/or "
                                                "cellular assays as simple BED files."
                                            )
                                        ),
                                        html.Li(
                                            str(
                                                "Email notification: If an email "
                                                "address is provided, the user will "
                                                "receive a notification upon the job "
                                                "completion."
                                            )
                                        ),
                                        html.Li(
                                            str(
                                                "Job name: If a string is provided, "
                                                "a prefix will be added to the unique "
                                                "job id to facilitate the identification "
                                                "of a particular search e.g. "
                                                "my_job_G05B8KHU0H."
                                            )
                                        ),
                                    ], 
                                    style={'padding': '15px'}
                                )
                            ]
                        )
                    ], 
                    style={'padding': '15px'}
                )
            ]
        )
    )
    # report status description
    final_list.append(
        html.P(
            [
                str(
                    "After selecting the desired inputs, clicking the Submit "
                    "button starts the search. After the submission, a new page "
                    "will show the status and progress."
                )
            ]
        )
    )
    final_list.append(
        html.Div(
            [
                html.Img(
                    src="data:image/png;base64,{}".format(
                        base64.b64encode(
                            open(
                                os.path.join(
                                    app_main_directory,
                                    ASSETS_DIR,
                                    "helpPage",
                                    "load_page.png"
                                ),
                                mode="rb",
                            ).read()
                        ).decode()
                    ),
                    width="100%",
                )
            ]
        )
    )
    final_list.append(
        html.P(
            [
                str(
                    "Upon completion of the job, a link \"View Results\" will "
                    "appear to view the results at the bottom of the status "
                    "report page."
                ),
                html.P(
                    html.Img(
                        src="data:image/png;base64,{}".format(
                            base64.b64encode(
                                open(
                                    os.path.join(
                                        app_main_directory,
                                        ASSETS_DIR,
                                        "helpPage",
                                        "jobDone.png"
                                    ),
                                    mode="rb",
                                ).read()
                            ).decode()
                        ),
                        width="100%"
                    )
                )
            ]
        )
    )
    # result page description
    final_list.append(html.H3("Result Page"))
    final_list.append(
        html.P(
            [
                str(
                    "CRISPRme summarizes the results in a table highlighting for "
                    "each gRNA its CFD specificity score and the count of on-targets "
                    "and off-targets found in the reference and variant genomes grouped "
                    "by number of mismatches and bulges. Of note, the CFD specificity "
                    "score was initially proposed for searches of up to 3 or 4 mismatches; "
                    "as the number of mismatches increase, the score decreases non-linearly. "
                    "Importantly, these scores should be compared with caution between "
                    "searches with different numbers of mismatches/bulges and/or different "
                    "genetic variant datasets.\nAt the top of the page, the user can find "
                    "a summary table reporting the nuclease, the CFD specificity score "
                    "and the number of targets in each category of mismatches and bulges. "
                    "In the top left corner there is a \"Download General Table\" button "
                    "allowing the download of the table as a text file."
                ),
                html.P(
                    html.Img(
                        src="data:image/png;base64,{}".format(
                            base64.b64encode(
                                open(
                                    os.path.join(
                                        app_main_directory,
                                        ASSETS_DIR,
                                        "resultPage",
                                        "resultsSummary.png"
                                    ),
                                    mode="rb"
                                ).read()
                            ).decode()
                        ),
                        width="100%"
                    )
                ),
                html.Ul(
                    [
                        html.Li(
                            str(
                                "CFD: Off-Target Cutting Frequency Determination "
                                "Score, calculates how much is the affinity of "
                                "the guides with the off-targets, basically tells "
                                "you the likelihood of the guide to perform cut in "
                                "off-target regions."
                            )
                        ),
                        html.Li(
                            str(
                                "Off-Targets Reference (0 - n Mismatches + Bulges): "
                                "shows how many possible Off-Targets the guide can "
                                "have in the Reference Genome. Targets are also grouped "
                                "by Mismatch + Bulge value."
                            )
                        ),
                        html.Li(
                            str(
                                "Off-Targets Variant (0 - n Mismatches + Bulges): "
                                "shows how many possible Off-Targets the guide can "
                                "have in the Variant Genome. Targets are also grouped "
                                "by Mismatch + Bulge value."
                            )
                        )
                    ], 
                    style={"padding": "15px"}
                ),
                html.P(
                    str(
                        "In addition, for each guide, six different interactive "
                        "reports are generated and are available to be downloaded: "
                        "Custom Ranking, Summary by Mismatches/Bulges, Summary by Sample, "
                        "Query Genomic Region, Graphical Reports and Personal Risk Cards."
                    )
                ),
                html.Ul(
                    [
                        html.Li(
                            [
                                html.Span(
                                    "Custom ranking: ", style={"color": "red"}
                                ), 
                                str(
                                    "In this report, users can filter and rank "
                                    "potential off-targets based on number of "
                                    "mismatches and/or bulges, CFD score, Risk "
                                    "Score (increase in CFD score due to genetic "
                                    "variants), or a combination of them."
                                )
                            ]
                        ),
                        html.Img(
                            src="data:image/png;base64,{}".format(
                                base64.b64encode(
                                    open(
                                        os.path.join(
                                            app_main_directory,
                                            ASSETS_DIR,
                                            "resultPage",
                                            "customRanking.png"
                                        ),
                                        mode="rb",
                                    ).read()
                                ).decode()
                            ),
                            width="100%"
                        ),
                        html.Li(
                            [
                                html.Span(
                                    "Summary by Mismatches/Bulges: ", 
                                    style={"color": "red"}
                                ), 
                                str(
                                    "This report shows a matrix separating targets "
                                    "into subgroups based on the type of target, "
                                    "mismatch count and bulge size. “X” targets "
                                    "contain only mismatches, “DNA” targets contain "
                                    "DNA bulges (and may contain mismatches), and \"RNA\" "
                                    "targets contain RNA bulges (and may contain mismatches)"
                                )
                            ]
                        ),
                        html.Img(
                            src="data:image/png;base64,{}".format(
                                base64.b64encode(
                                    open(
                                        os.path.join(
                                            app_main_directory,
                                            ASSETS_DIR,
                                            "resultPage",
                                            "summaryByGuide.png"
                                        ),
                                        mode="rb"
                                    ).read()
                                ).decode()
                            ),
                            width="100%"
                        ),
                        html.Ul(
                            [
                                html.Li(
                                    "Bulge Type: type of bulge of the targets, can be X (no bulge), DNA or RNA."
                                ),
                                html.Li(
                                    "Bulge Size: size of the bulge present in the targets."
                                ),
                                html.Li(
                                    "Mismatches: number of mismatches present in the targets."
                                ),
                                html.Li(
                                    "Targets in Reference: number of targets found in the Reference Genome for that combination of mismatch/bulge."
                                ),
                                html.Li(
                                    "Targets in Variant: number of targets found in the Variant Genome for that combination of mismatch/bulge. Each of these targets is associated with at least one sample."
                                ),
                                html.Li(
                                    "PAM Creation: number of possible created PAMs due to variants addition."
                                ),
                                html.Li(
                                    "Show Targets: open a new page to display all the targets of the row of interest as in the following image:"),
                                html.Img(
                                    src="data:image/png;base64,{}".format(
                                        base64.b64encode(
                                            open(
                                                os.path.join(
                                                    app_main_directory,
                                                    ASSETS_DIR,
                                                    "resultPage",
                                                    "summaryByGuide_show_targets.png"
                                                ),
                                                mode="rb",
                                            ).read()
                                        ).decode()
                                    ),
                                    width="100%"
                                )
                            ]
                        ),
                        html.Li(
                            [
                                html.Span(
                                    "Summary by Sample: ", style={"color": "red"}
                                ), 
                                str(
                                    "This page shows all the samples present in "
                                    "the VCFs and allows users to extract and "
                                    "visualize targets related to each sample."
                                )
                            ]
                        ),
                        html.Img(
                            src="data:image/png;base64,{}".format(
                                base64.b64encode(
                                    open(
                                        os.path.join(
                                            app_main_directory,
                                            ASSETS_DIR,
                                            "resultPage",
                                            "summaryBySamples.png"
                                        ),
                                        mode="rb"
                                    ).read()
                                ).decode()
                            ),
                            width="100%"
                        ),
                        html.Ul(
                            [
                                html.Li("Gender: the sample gender"),
                                html.Li(
                                    "Population: population which the sample belong to"
                                ),
                                html.Li(
                                    "Super Population: continent which the sample belong to"
                                ),
                                html.Li(
                                    "Targets in Variant: number of targets found in the Variant Genome that are generated by that sample"
                                ),
                                html.Li(
                                    "Targets in Population: number of targets found in the Variant Genome that are generated by all the sample of the population"
                                ),
                                html.Li(
                                    "Targets in Super Population: number of targets found in the Variant Genome that are generated by all the populations"
                                ),
                                html.Li(
                                    "PAM Creation: number of possible created PAMs due to variants addition"
                                ),
                                html.Li(
                                    "Show Targets: open a new page to display all the targets of the row of interest as in the following image:"),
                                html.Img(
                                    src="data:image/png;base64,{}".format(
                                        base64.b64encode(
                                            open(
                                                os.path.join(
                                                    app_main_directory,
                                                    ASSETS_DIR,
                                                    "resultPage",
                                                    "summaryBySamples_show_targets.png"
                                                ),
                                                mode="rb"
                                            ).read()
                                        ).decode()
                                    ),
                                    width="100%"
                                )
                            ]
                        ),
                        html.Li(
                            [
                                html.Span(
                                    "Query Genomic Region: ", 
                                    style={"color": "red"}
                                ), 
                                str(
                                    "This page allows the user to retrieve "
                                    "targets overlapping a specific genomic "
                                    "region, for example to quickly assess "
                                    "potential off-targets in a given regulatory "
                                    "element or coding region."
                                )
                            ]
                        ),
                        html.Img(
                            src="data:image/png;base64,{}".format(
                                base64.b64encode(
                                    open(
                                        os.path.join(
                                            app_main_directory,
                                            ASSETS_DIR,
                                            "resultPage",
                                            "summaryByPosition.png"
                                        ),
                                        mode="rb"
                                    ).read()
                                ).decode()
                            ),
                            width="100%"
                        ),
                        html.Li(
                            [
                                html.Span(
                                    "Graphical Reports: ", style={"color": "red"}
                                ), 
                                str(
                                    "This page creates several graphical reports "
                                    "for each selected sgRNA."
                                )
                            ]
                        ),
                        html.Li(
                            str(
                                "A stem plot shows how genetic variants affect "
                                "predicted off-target potential. The arrow "
                                "connecting the red (reference allele off-target) " 
                                "and blue (alternative allele off-target) dots "
                                "shows the increase in predicted cleavage potential "
                                "due to the variant."
                            )
                        ),
                        html.Img(
                            src="data:image/png;base64,{}".format(
                                base64.b64encode(
                                    open(
                                        os.path.join(
                                            app_main_directory,
                                            ASSETS_DIR,
                                            "resultPage",
                                            "lolli_plot.png"
                                        ),
                                        mode="rb"
                                    ).read()
                                ).decode()
                            ),
                            width="100%"
                        ),
                        html.Li(
                            str(
                                "Bar plots depict  how candidate off-targets are "
                                "distributed across super-populations based on the "
                                "number of mismatches and bulges."
                            )
                        ),
                        html.Li(
                            str(
                                "A radar chart with the relative specificity of "
                                "the analyzed guide for each functional region based "
                                "on annotations from GENCODE and ENCODE. A larger area "
                                "in the chart represents a gRNA with more potential "
                                "off-targets falling in annotated regions, possibly "
                                "representing an undesirable outcome. A summary table "
                                "provides the count and percentage of off-targets with "
                                "a given annotation."
                            )
                        ),
                        html.Li(
                            str(
                                "A motif logo summarizing the frequency of mismatches "
                                "and bulges (b) for each base of the protospacer + PAM."
                            )
                        ),
                        html.Img(
                            src="data:image/png;base64,{}".format(
                                base64.b64encode(
                                    open(
                                        os.path.join(
                                            app_main_directory,
                                            ASSETS_DIR,
                                            "resultPage",
                                            "barplot_radar_motif.png"
                                        ),
                                        mode="rb"
                                    ).read()
                                ).decode()
                            ),
                            width="100%"
                        ),
                        html.Li(
                            [
                                html.Span(
                                    "Personal Risk Cards: ", 
                                    style={"color": "red"}
                                ), 
                                str(
                                    "CRISPRme provides a dedicated page to generate "
                                    "reports called Personal Risk Cards that summarize "
                                    "potential off-target editing by a particular gRNA "
                                    "in a given individual due to genetic variants. This "
                                    "feature is particularly useful for retrieving and "
                                    "investigating private off-targets. The report "
                                    "contains two dynamically generated plots depicting "
                                    "all the candidate variant off-targets for the sample "
                                    "including those non-unique to the individual and those "
                                    "that are unique to the individual. These plots highlight "
                                    "how the introduction of genetic variants can change the "
                                    "predicted off-target cleavage potential, thereby "
                                    "demonstrating the importance of variant-aware off-target "
                                    "assessment as in CRISPRme. The report also contains two "
                                    "tables, consisting of a summary (Table 1, top) and "
                                    "information on each extracted candidate off-target "
                                    "(Table 2, bottom) with the following columns:"
                                )
                            ]
                        ),
                        html.Img(
                            src="data:image/png;base64,{}".format(
                                base64.b64encode(
                                    open(
                                        os.path.join(
                                            app_main_directory,
                                            ASSETS_DIR,
                                            "resultPage",
                                            "personalCard_top.png"
                                        ),
                                        mode="rb"
                                    ).read()
                                ).decode()
                            ),
                            width="100%"
                        ),
                        html.Img(
                            src="data:image/png;base64,{}".format(
                                base64.b64encode(
                                    open(
                                        os.path.join(
                                            app_main_directory,
                                            ASSETS_DIR,
                                            "resultPage",
                                            "personalCard_bottom.png"
                                        ),
                                        mode="rb",
                                    ).read()
                                ).decode()
                            ),
                            width="100%"
                        ),
                        html.Ul(
                            [
                                html.P(
                                    "Table 1:"),
                                html.Li(
                                    str(
                                        "Personal, count of all the candidate "
                                        "variant off-targets for the selected sample "
                                        "(including variants unique and non-unique to "
                                        "the individual)"
                                    )
                                ),
                                html.Li(
                                    str(
                                        "PAM creation, count of all the instances "
                                        "where a genetic variant in the sample "
                                        "introduces a new PAM."
                                    )
                                ),
                                html.Li(
                                    str(
                                        "Private, count of all the candidate "
                                        "variant off-targets uniquely found in "
                                        "the selected sample."
                                    )
                                ),
                            ]
                        )                        
                    ],
                    style={"padding": "15px"},
                )
            ]
        )
    )
    page = html.Div(final_list, style={"margin": "1%"})
    return page
