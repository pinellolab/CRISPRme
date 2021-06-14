
import dash
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
import dash_daq as daq
from dash_html_components import Li
import dash_table
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import base64  # for decoding upload content
import io  # for decoding upload content
from app import app_main_directory


def helpPage():
    final_list = []
    final_list.append(
        html.Div(
            [
                html.H3('About'),
                html.P(
                    ['CRISPRme is available as an online web app at ',
                     html.A('http://crisprme.di.univr.it',
                            href='http://crisprme.di.univr.it', target='_blank'),
                     ' a standalone command line package. The required inputs to perform an online search are: gRNA spacer(s), Cas protein, PAM sequence, genome build with or without the inclusion of genetic variants (1000G, HGDP and/or personal variants), and thresholds of mismatches and RNA/DNA bulges.'
                     ]
                )
            ]
        )
    )

    final_list.append(html.H3('Main Page'))

    final_list.append(
        html.Div([
            html.P(
                [
                    'A search on CRISPRme can be performed in three simple steps thanks to the user-friendly user interface. Several options are available to personalize a search.',
                 html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(open(app_main_directory +
                                                                                      '/assets/main_page.png', 'rb').read()).decode()), width="100%", height="auto")
                 ]
            ),
            html.Ul(
                [
                    html.Li(
                        [
                            'STEP 1: Spacer, Cas Protein and PAM selection',
                            html.Ul(
                                [
                                    # html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(
                                    #     open(app_main_directory+'/assets/helpPage/guides.png', 'rb').read()).decode()), width='30%'),
                                    html.Li(
                                        'Spacer(s): The guide RNA (gRNA) spacer sequence matches the genomic target protospacer sequence (typically 20 nucleotides) and directs Cas protein binding to the protospacer in the presence of a protospacer adjacent motif (PAM). The spacer sequence is represented as DNA (rather than RNA) in CRISPRme to allow easy comparison to the aligned protospacer sequence. CRISPRme accepts a set of gRNA spacer(s), one per line, each with the same length (max 100 sequences). The input spacer sequence should not include PAM.'),
                                    # html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(
                                    #     open(app_main_directory+'/assets/helpPage/sequence.png', 'rb').read()).decode()), width='40%'),
                                    html.Li('Genomic sequence(s): CRISPRme can alternatively take as input a set of genomic coordinates in BED format (chromosome# start end) or DNA sequences in FASTA format (max 1000 characters). The BED file coordinates will be treated as 0-based and CRISPRme (online version) will extract the first 100 possible spacer sequences within these coordinates starting with the positive strand. To use this type of input, the user must delimit each entry with a >header.'),
                                    # html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(
                                    #     open(app_main_directory+'/assets/helpPage/nuclease.png', 'rb').read()).decode()), width='30%'),
                                    html.Li(
                                        'PAM sequence: The PAM is a short (∼2-6 nucleotide) DNA sequence adjacent to the protospacer necessary for the Cas protein to bind to a specific DNA target. CRISPRme supports a set of PAMs and users must select one of them in order to perform the search. The software supports both 3’ (e.g. SpCas9) and 5’ (e.g. Cas12a) PAM sequences.'),
                                    # html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(
                                    #     open(app_main_directory+'/assets/helpPage/pam.png', 'rb').read()).decode()), width='30%'),
                                    # html.Li(
                                    #     'PAM: here you can select a Protospacer Adjacent Motif for the specified Cas protein.'),
                                ], style={'padding': '15px'}
                            )
                        ]
                    ),
                    # html.Li(
                    #     [
                    #         'STEP 1: Genome and PAM selection',
                    #         html.Ul(
                    #             [
                    #                 html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(
                    #                     open(app_main_directory+'/assets/helpPage/genome.png', 'rb').read()).decode()), width='30%'),
                    #                 html.Li(
                    #                     'Genome: here you can select a genome from the ones present, some genomes have also the variant version enriched (indicated with a \'+\' symbol) with genetic variant from the 1000genome project'),
                    #                 html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(
                    #                     open(app_main_directory+'/assets/helpPage/pam.png', 'rb').read()).decode()), width='30%'),
                    #                 html.Li(
                    #                     'PAM: here you can select a Protospacer Adjacent Motif for a specific Cas protein.'),
                    #             ], style={'padding': '15px'}
                    #         )
                    #     ]
                    # ),
                    html.Li(
                        [
                            'STEP 2: Genome selection and threshold configuration',
                            html.Ul(
                                [
                                    # html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(
                                    #     open(app_main_directory+'/assets/helpPage/genome.png', 'rb').read()).decode()), width='40%'),
                                    html.Li(
                                        'Genome builds: The genome builds are based on FASTA files from UCSC. The hg38 and hg19 genomic builds are available with the option to incorporate variants from 1000G, HGDP, and/or personal variants in the search. The option to add personal variants is enabled only for the local offline and command line versions.'),
                                    # html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(
                                    #     open(app_main_directory+'/assets/helpPage/thresholds.png', 'rb').read()).decode()), width='20%'),
                                    html.Li('Search thresholds: CRISPRme allows users to specify the number of mismatches, DNA and RNA bulges tolerated in enumerating potential off-targets. The web-tool allows up to 6 mismatches and up to 2 RNA/DNA bulges (which can be consecutive (NN--NN) or interleaved (NN-N-NN)). However, for the command line version, these thresholds can be set freely and depend only on the available computational resources.'),
                                    # html.Li(
                                    #     'Bulge DNA size: size of bubbles tolerated on the DNA sequence (can be consecutive(AA--AA) or interleaved(AA-A-AA)).'),
                                    # html.Li(
                                    #     'Bulge RNA size: size of bubbles tolerated on the RNA sequence (can be consecutive(AA--AA) or interleaved(AA-A-AA))'),
                                ], style={'padding': '15px'}
                            )
                        ]
                    ),
                    # html.Li(
                    #     [
                    #         'STEP 2: Guide insertion and threshold configuration',
                    #         html.Ul(
                    #             [
                    #                 html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(
                    #                     open(app_main_directory+'/assets/helpPage/guides.png', 'rb').read()).decode()), width='40%'),
                    #                 html.Li(
                    #                     'Guides: a list of crRNAs sequences, consisting in 1 or more sequences (max 1000 sequences) to search on the genome'),
                    #                 html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(
                    #                     open(app_main_directory+'/assets/helpPage/sequence.png', 'rb').read()).decode()), width='40%'),
                    #                 html.Li('Sequence: one or more genetic sequences (max 1000 characters), each sequence MUST BE separated with the header \'>name\'. The sequence can be also submitted with a ' +
                    #                         'chromosome range, also provided with an header. The region will be extracted from the Genome selected in STEP 1'),
                    #                 html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(
                    #                     open(app_main_directory+'/assets/helpPage/crRNA.png', 'rb').read()).decode()), width='20%'),
                    #                 html.Li(
                    #                     'Allowed mismatches: number of tolerated mismatches in a target'),
                    #                 html.Li(
                    #                     'Bulge DNA size: size of bubbles tolerated on the DNA sequence (can be consecutive(AA--AA) or interleaved(AA-A-AA)).'),
                    #                 html.Li(
                    #                     'Bulge RNA size: size of bubbles tolerated on the RNA sequence (can be consecutive(AA--AA) or interleaved(AA-A-AA))'),
                    #                 # html.Img(src = 'data:image/png;base64,{}'.format(base64.b64encode(open(app_main_directory+'/assets/helpPage/crRNA.png', 'rb').read()).decode()), width='20%' ),
                    #                 html.Li(
                    #                     'crRNA length: available only when a genetic sequence is given as input, represents the length of the guides (without PAM) that you want to extract from the sequence.')
                    #             ], style={'padding': '15px'}
                    #         )
                    #     ]
                    # ),
                    html.Li(
                        [
                            'STEP 3 Annotation(s), email notification, and job name',
                            html.Ul(
                                [
                                    # html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(
                                    #     open(app_main_directory+'/assets/helpPage/advOpt.png', 'rb').read()).decode()), width='40%'),
                                    html.Li('Functional annotations: To assess the potential impact of off-target activity, CRISPRme provides a set of functional annotations for coding and non-coding regions. The annotations are based on files obtained from the Encyclopedia of DNA Elements (ENCODE) containing candidate cis regulatory elements21 and from GENCODE25 containing annotations for protein coding genes, transcribed but untranslated regions, and introns. In the offline versions of CRISPRme, users can add custom genome annotations, such as cell-type specific chromatin marks or off-target sites nominated by in vitro and/or cellular assays as simple BED files.'),
                                    html.Li(
                                        'Email notification: If an email address is provided, the user will receive a notification upon the job completion.'),
                                    html.Li(
                                        'Job name: If a string is provided, a prefix will be added to the unique job id to facilitate the identification of a particular search e.g. my_job_G05B8KHU0H.'),
                                ], style={'padding': '15px'}
                            )
                        ]
                    )
                    # html.Li(
                    #     [
                    #         'STEP 3 (Advanced options): Select various comparisons',
                    #         html.Ul(
                    #             [
                    #                 html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(
                    #                     open(app_main_directory+'/assets/helpPage/advOpt.png', 'rb').read()).decode()), width='40%'),
                    #                 html.Li(
                    #                     'Compare your results with the GeCKO v2 library: selected by default, compares the results of your guides with the results obtained in a previous search with guides from the well-known GeCKO library.'),
                    #                 html.Li('Compare your results with the corresponding reference genome: selected by default when an enriched genome is chosen, compares the results with the respective reference genome to evaluate differences when variant are added.'),
                    #                 html.Li(
                    #                     'Notify me by email: if selected, let you insert an email to receive a notification when your job is terminated.'),
                    #             ], style={'padding': '15px'}
                    #         )
                    #     ]
                    # )
                ], style={'padding': '15px'}
            )
        ])
    )

    final_list.append(
        html.P(
            ['After selecting the desired inputs, clicking the Submit button starts the search . After the submission, a new page will show the status and progress.']
        )
    )

    final_list.append(
        html.Div(
            [
                html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(
                    open(app_main_directory+'/assets/helpPage/load_page.png', 'rb').read()).decode()), width='100%')
            ]
        )
    )

    final_list.append(
        html.P(
            [
                # 'After the submission, the status of the search will be displayed in a new page',
                # html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(
                #     open(app_main_directory+'/assets/waitPage/loadPage.png', 'rb').read()).decode()), width='100%'),
                'Upon completion of the job, a link “View Results” will appear to view the results at the bottom of the status report page.',
                html.P(html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(
                    open(app_main_directory+'/assets/helpPage/jobDone.png', 'rb').read()).decode()), width='100%'))
            ]
        )
    )

    final_list.append(html.H3('Result Page'))
    final_list.append(
        html.P(
            [
                'CRISPRme summarizes the results in a table highlighting for each gRNA its CFD specificity score and the count of on-targets and off-targets found in the reference and variant genomes grouped by number of mismatches and bulges. Of note, the CFD specificity score was initially proposed for searches of up to 3 or 4 mismatches; as the number of mismatches increase, the score decreases non-linearly. Importantly, these scores should be compared with caution between searches with different numbers of mismatches/bulges and/or different genetic variant datasets,\n      \
                At the top of the page, the user can find a summary table reporting the nuclease, the CFD specificity score and the number of targets in each category of mismatches and bulges. In the top left corner there is a “Download General Table”’ button allowing the download of the table as a text file.',
                html.P(html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(open(app_main_directory +
                                                                                            '/assets/resultPage/resultsSummary.png', 'rb').read()).decode()), width='100%')),
                html.Ul(
                    [
                        html.Li('CFD: Off-Target Cutting Frequency Determination Score, calculates how much is the affinity of the guides with the off-targets, basically tells you the likelihood of the guide to perform cut in off-target regions.'),
                        # html.Li('Doench 2016: On-Target Efficacy Scoring (Azimuth 2.0), it’s a trained machine learning model that gives you the likelihood of on-target activity for the selected guide.'),
                        # html.Li(
                        #     'On-Targets Reference: shows how many possible On-Targets the guide can have in the Reference Genome.'),
                        # html.Li(['Samples in Class 0 - 0+ - 1 - 1+: shows the number of samples grouped by Sample Class:',
                        #          html.Ul([
                        #              html.Li(
                        #                  'Class 0: Samples that does not have any On-Targets'),
                        #              html.Li(
                        #                  'Class 0+: Samples that have a subset of the Reference Genome On-Targets'),
                        #              html.Li(
                        #                  'Class 1: Samples that have the same On-Targets as the Reference Genome'),
                        #              html.Li(
                        #                  'Class 1+: Samples that creates at least a new On-Target, that is not present in the Reference Genome')
                        #          ])
                        #          ]),
                        html.Li('Off-Targets Reference (0 - n Mismatches + Bulges): shows how many possible Off-Targets the guide can have in the Reference Genome. Targets are also grouped by Mismatch + Bulge value.'),
                        html.Li('Off-Targets Variant (0 - n Mismatches + Bulges): shows how many possible Off-Targets the guide can have in the Variant Genome. Targets are also grouped by Mismatch + Bulge value.')
                    ], style={'padding': '15px'}
                ),
                # html.P(html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(open(app_main_directory+
                #     '/assets/resultPage/populationDistribution.png', 'rb').read()).decode()), width='100%')),
                # 'The Show Target Distribution in Populations button opens a section where informations about the number of targets found in each Superpopulation (EAS, EUR, AFR, AMR, SAS) are provided by means of a barplot for each Mismatch + Bulge value. ',
                html.P('In addition, for each guide, six different interactive reports are generated and are available to be downloaded:  Custom Ranking, Summary by Mismatches/Bulges, Summary by Sample, Query Genomic Region, Graphical Reports and Personal Risk Cards:'),
                html.Ul(
                    [
                        html.Li([html.Span('Custom ranking: ', style={
                                'color': 'red'}), 'In this report, users can filter and rank potential off-targets based on number of mismatches and/or bulges, CFD score, Risk Score (increase in CFD score due to genetic variants), or a combination of them.']),
                        html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(open(app_main_directory +
                                                                                             '/assets/resultPage/customRanking.png', 'rb').read()).decode()), width='100%'),


                        html.Li([html.Span('Summary by Mismatches/Bulges: ', style={
                                'color': 'red'}), 'This report shows a matrix separating targets into subgroups based on the type of target, mismatch count and bulge size. “X” targets contain only mismatches, “DNA” targets contain DNA bulges (and may contain mismatches), and “RNA” targets contain RNA bulges (and may contain mismatches)']),
                        html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(open(app_main_directory +
                                                                                             '/assets/resultPage/summaryByGuide.png', 'rb').read()).decode()), width='100%'),
                        html.Ul(
                            [
                                html.Li(
                                    'Bulge Type: type of bulge of the targets, can be X (no bulge), DNA or RNA.'),
                                html.Li(
                                    'Bulge Size: size of the bulge present in the targets.'),
                                html.Li(
                                    'Mismatches: number of mismatches present in the targets.'),
                                html.Li(
                                    'Targets in Reference: number of targets found in the Reference Genome for that combination of mismatch/bulge.'),
                                html.Li(
                                    'Targets in Variant: number of targets found in the Variant Genome for that combination of mismatch/bulge. Each of these targets is associated with at least one sample.'),
                                html.Li(
                                    'PAM Creation: number of possible created PAMs due to variants addition.'),
                                html.Li(
                                    'Show Targets: open a new page to display all the targets of the row of interest as in the following image:'),
                                html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(open(app_main_directory +
                                                                                                     '/assets/resultPage/summaryByGuide_show_targets.png', 'rb').read()).decode()), width='100%'),
                            ]
                        ),


                        html.Li([html.Span('Summary by Sample: ', style={
                                'color': 'red'}), 'This page shows all the samples present in the VCFs and allows users to extract and visualize targets related to each sample.']),
                        html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(open(app_main_directory +
                                                                                             '/assets/resultPage/summaryBySamples.png', 'rb').read()).decode()), width='100%'),
                        html.Ul(
                            [
                                html.Li(
                                    'Gender: the sample gender'),
                                html.Li(
                                    'Population: population which the sample belong to'),
                                html.Li(
                                    'Super Population: continent which the sample belong to'),
                                html.Li(
                                    'Targets in Variant: number of targets found in the Variant Genome that are generated by that sample'),
                                html.Li(
                                    'Targets in Population: number of targets found in the Variant Genome that are generated by all the sample of the population'),
                                html.Li(
                                    'Targets in Super Population: number of targets found in the Variant Genome that are generated by all the populations'),
                                html.Li(
                                    'PAM Creation: number of possible created PAMs due to variants addition'),
                                # html.Li(
                                #     'Class: Sample Class (0 - 0+ - 1 - 1+) associated with the sample'),
                                html.Li(
                                    'Show Targets: open a new page to display all the targets of the row of interest as in the following image:'),
                                html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(open(app_main_directory +
                                                                                                     '/assets/resultPage/summaryBySamples_show_targets.png', 'rb').read()).decode()), width='100%'),
                            ]
                        ),

                        html.Li([html.Span('Query Genomic Region: ', style={
                                'color': 'red'}), 'This page allows the user to retrieve targets overlapping a specific genomic region, for example to quickly assess potential off-targets in a given regulatory element or coding region.']),
                        html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(open(app_main_directory +
                                                                                             '/assets/resultPage/summaryByPosition.png', 'rb').read()).decode()), width='100%'),
                        # html.Ul(
                        #     [
                        #         html.Li(
                        #             'Position: chromosome relative position of the first letter of the guide'),
                        #         html.Li(
                        #             'Best Target: best target found in that position'),
                        #         html.Li(
                        #             'Min Mismatch: minimum number of mismatches present in the targets in that position'),
                        #         html.Li(
                        #             'Min Bulge: minimum number of bulges present in the targets in that position'),
                        #         html.Li(
                        #             'Bulge: number of bulges present in the targets in that position'),
                        #         html.Li(
                        #             'Targets in Cluster by Mismatch Value: Matrix showing the distribution of the targets grouped by mismatch/bulge count'),
                        #     ]

                        # ),
                        # html.Li([html.Span('Graphical Reports: ', style={
                        #         'color': 'red'}), 'This page shows graphics about a specific guide, including genomic annotation and motif logos. The main feature introduced is the possibility to visualize graphical reports at individual level.']),
                        # html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(open(app_main_directory +
                        #                                                                      '/assets/resultPage/summaryByGraphic_population.png', 'rb').read()).decode()), width='100%'),
                        # html.Li(
                        #     'Select a Mismatch and Bulge Value: generate graphics with the specified mismatch+bulge value'),
                        # html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(open(app_main_directory +
                        #                                                                      '/assets/resultPage/summaryByGraphic_sample.png', 'rb').read()).decode()), width='100%'),
                        # html.Li(
                        #     'Select Individual Data: generate individual data, by selecting Super Population, Population and Sample'),

                        # html.Li([html.Span('Personal Risk Cards: ', style={
                        #         'color': 'red'}), 'This page shows at individual level the most important data for a given sample.']),
                        # html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(open(app_main_directory +
                        #                                                                      '/assets/resultPage/personalCard.png', 'rb').read()).decode()), width='100%'),
                        # html.Ul(
                        #     [
                        #         html.Li(
                        #             'Plot showing the difference in terms of CFD score in a specific genetic region, when accounting for variants in targets in which the selected sample appear.'),

                        #         html.Li(
                        #             'Plot showing the difference in terms of CFD score in a specific genetic region, when accounting for variants in targets unique to the selected sample.'),
                        #         html.Li(
                        #             'Set of tables reporting the personal information for the selected sample. \n \
                        #             The first table reports personal targets, PAM creation target and private targets. \n \
                        #             The second table lists for each target the crRNA and DNA sequences, the position and cluster position, \
                        #             the chromosome, direction, mismatches, bulge size and total. It also reports the CFD score for the reference and variant target, \
                        #              the annotation and eventually the new PAM generation due to substitution or insertion/deletion (not showed in the figure).')
                        #     ]
                        # )
                        html.Li([html.Span('Graphical Reports: ', style={
                                'color': 'red'}), 'This page creates several graphical reports for each selected sgRNA.']),
                        # html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(open(app_main_directory +
                        #                                                                      '/assets/resultPage/summaryByGraphic_population.png', 'rb').read()).decode()), width='100%'),
                        html.Li(
                            'A stem plot shows how genetic variants affect predicted off-target potential. The arrow connecting the red (reference allele off-target) and blue (alternative allele off-target) dots shows the increase in predicted cleavage potential due to the variant.'),
                        html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(open(app_main_directory +
                                                                                             '/assets/resultPage/lolli_plot.png', 'rb').read()).decode()), width='100%'),
                        html.Li(
                            'Bar plots depict  how candidate off-targets are distributed across super-populations based on the number of mismatches and bulges '),
                        html.Li(
                            'A radar chart with the relative specificity of the analyzed guide for each functional region based on annotations from GENCODE and ENCODE. A larger area in the chart represents a gRNA with more potential off-targets falling in annotated regions, possibly representing an undesirable outcome. A summary table provides the count and percentage of off-targets with a given annotation.'),
                        html.Li(
                            'A motif logo summarizing the frequency of mismatches and bulges (b) for each base of the protospacer + PAM.'),
                        html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(open(app_main_directory +
                                                                                             '/assets/resultPage/barplot_radar_motif.png', 'rb').read()).decode()), width='100%'),


                        html.Li([html.Span('Personal Risk Cards: ', style={
                                'color': 'red'}), 'CRISPRme provides a dedicated page to generate reports called Personal Risk Cards that  summarize potential off-target editing by a particular gRNA in a given individual due to genetic variants. This feature is particularly useful for retrieving and investigating private off-targets. \
                                    The report contains two dynamically generated plots depicting all the candidate variant off-targets for the sample including those non-unique to the individual and those that are unique to the individual. These plots highlight how the introduction of genetic variants can change the predicted off-target cleavage potential, \
                                    thereby demonstrating the importance of variant-aware off-target assessment as in CRISPRme. The report also contains two tables, consisting of a summary (Table 1, top) and information on each extracted candidate off-target (Table 2, bottom) with the following columns:']),
                        html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(open(app_main_directory + \
                                                                                             '/assets/resultPage/personalCard_top.png', 'rb').read()).decode()), width='100%'),
                        html.Img(src='data:image/png;base64,{}'.format(base64.b64encode(open(app_main_directory + \
                                                                                             '/assets/resultPage/personalCard_bottom.png', 'rb').read()).decode()), width='100%'),
                        # html.Ul(
                        #     [
                        #         html.Li(
                        #             'Plot showing the difference in terms of CFD score in a specific genetic region, when accounting for variants in targets in which the selected sample appear.'),

                        #         html.Li(
                        #             'Plot showing the difference in terms of CFD score in a specific genetic region, when accounting for variants in targets unique to the selected sample.'),
                        #         html.Li(
                        #             'Set of tables reporting the personal information for the selected sample. \n \
                        #             The first table reports personal targets, PAM creation target and private targets. \n \
                        #             The second table lists for each target the crRNA and DNA sequences, the position and cluster position, \
                        #             the chromosome, direction, mismatches, bulge size and total. It also reports the CFD score for the reference and variant target, \
                        #              the annotation and eventually the new PAM generation due to substitution or insertion/deletion (not showed in the figure).')
                        #     ]
                        # )
                        html.Ul(
                            [
                                html.P(
                                    'Table 1:'),
                                html.Li(
                                    'Personal, count of all the candidate variant off-targets for the selected sample (including variants unique and non-unique to the individual)'),
                                html.Li(
                                    'PAM creation, count of all the instances where a genetic variant in the sample introduces a new PAM.'),
                                html.Li(
                                    'Private, count of all the candidate variant off-targets uniquely found in the selected sample.'),
                            ]
                        )
                    ], style={'padding': '15px'}
                )

            ]
        )
    )

    # final_list.append(
    #     html.Div(
    #         [
    #             html.H3('Browser Compatibility'),
    #             html.Div([
    #                 dash_table.DataTable(
    #                     data=[{'OS': 'Linux', 'V': 'Ubuntu 18.04.2 LTS', 'Ch': '79.0', 'S': 'n/a', 'ME': 'n/a', 'F': '71.0'},
    #                           {'OS': 'MacOS', 'V': 'Mojave', 'Ch': ' 79.0',
    #                            'S': '13.0.4', 'ME': 'n/a', 'F': '71.0'},
    #                           {'OS': 'Windows', 'V': '10', 'Ch': '79.0', 'S': 'n/a', 'ME': '44.18362.449.0', 'F': '71.0'}],

    #                     columns=[{'id': 'OS', 'name': 'Operative System'}, {'id': 'V', 'name': 'Version'}, {'id': 'Ch', 'name': 'Chrome'},
    #                              {'id': 'S', 'name': 'Safari'}, {'id': 'ME', 'name': 'Microsoft Edge'}, {
    #                                  'id': 'F', 'name': 'Firefox'}
    #                              ],

    #                     style_cell={
    #                         'textAlign': 'center',
    #                         'width': '20'
    #                     },

    #                     style_data_conditional=[
    #                         {
    #                             'if': {'row_index': 'odd'},
    #                             'backgroundColor': 'rgb(248, 248, 248)'
    #                         }
    #                     ],
    #                     style_header={
    #                         'backgroundColor': 'rgb(230, 230, 230)',
    #                         'fontWeight': 'bold'
    #                     }
    #                 )
    #             ],   style={'display': 'inline-block', 'width': '48%'})
    #         ]
    #     )
    # )
    page = html.Div(final_list, style={'margin': '1%'})
    return page
