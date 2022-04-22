"""Display the contact page on CRISPRme website.

Contacts list:
    - Samuele Cancellieri, PhD 
    - Yingqi Linda
    - Daniel Bauer, MD, PhD <bauer@bloodgroup.tch.harvard.edu>
    - Nicola Bombieri, PhD
    - Rosalba Giugno, PhD <rosalba.giugno@univr.it>
    - Luca Pinello, PhD <lpinello@mgh.harvard.edu>


Contacts list for bug reports and maintainers:
    - Samuele Cancellieri, PhD <samuele.cancellieri@univr.it>
    - Manuel Tognon <manuel.tognon@univr.it>
    - Luca Pinello, PhD <lpinello@mgh.harvard.edu>

CRISPRme developers:
    - Samuele Cancellieri, PhD <samuele.cancellieri@univr.it>
    - Elia Dirupo
    - Francesco Masillo
    - Manuel Tognon <manuel.tognon@univr.it>
    
Current CRISPRme maintainers:
    - Samuele Cancellieri, PhD <samuele.cancellieri@univr.it>
    - Manuel Tognon <manuel.tognon@univr.it>
    - Luca Pinello, PhD <lpinello@mgh.harvard.edu>
"""


from .pages_utils import GITHUB_LINK

from typing import List

import dash_html_components as html


# the following list could be changed during time
# TODO: maintain the following list up to date as much as possible
INFOMICS_DEVS = [
    "Samuele Cancellieri",
    "Francesco Masillo",
    "Manuel Tognon",
    "Nicola Bombieri",
    "Rosalba Giugno",
]
PINELLOLAB_DEVS = ["Luca Pinello"]
BAUERLAB_DEVS = ["Linda Yingqi", "Daniel Bauer"]


def contact_page() -> List:
    """Create the layout of CRISPRme contacts page.

    ...

    Parameters
    ----------
    None

    Returns
    -------
    List
        Contacts webpage
    """

    infomics_p = str(
        "InfOmics Lab, Department of Computer Science, University of Verona, Italy "
    )
    infomics_p = ". ".join([", ".join([dev for dev in INFOMICS_DEVS]), infomics_p])
    pinellolab_p = str(
        "Molecular Pathology Unit and Center for Cancer Research, "
        "Massachusetts General Hospital, Charlestown, MA, USA; Department of "
        "Pathology, Harvard Medical School, Boston, MA, USA; Broad Institute "
        "of MIT and Harvard, Cambridge, MA, USA "
    )
    pinellolab_p = ". ".join(
        [", ".join([dev for dev in PINELLOLAB_DEVS]), pinellolab_p]
    )
    bauerlab_p = str(
        "Division of Hematology/Oncology, Boston Children's Hospital, "
        "Department of Pediatric Oncology, Dana-Farber Cancer Institute, "
        "Harvard Stem Cell Institute, Broad Institute, Department of "
        "Pediatrics, Harvard Medical School, Boston, MA, USA "
    )
    bauerlab_p = ". ".join([", ".join([dev for dev in BAUERLAB_DEVS]), bauerlab_p])
    f = []
    f.append(
        html.Div(
            [
                html.H3("CONTACTS"),
                html.P("CRISPRme was developed by:"),
                html.Ul(
                    [
                        html.Li(
                            [
                                infomics_p,
                                html.A(
                                    "(https://infomics.github.io/InfOmics/)",
                                    href="https://infomics.github.io/InfOmics/",
                                    target="_blank",
                                ),
                            ]
                        ),
                        html.Li(
                            [
                                pinellolab_p,
                                html.A(
                                    "(http://pinellolab.org/)",
                                    href="http://pinellolab.org/",
                                    target="_blank",
                                ),
                            ]
                        ),
                        html.Li(
                            [
                                bauerlab_p,
                                html.A(
                                    "(http://bauerlab.org/)",
                                    href="http://bauerlab.org/",
                                    target="_blank",
                                ),
                            ]
                        ),
                    ],
                    style={"padding": "15px"},
                ),
                html.P("Please send any comment or bug to:"),
                html.Ul(
                    [
                        html.Li("rosalba DOT giugno AT univr DOT it"),
                        html.Li("lpinello AT mgh DOT harvard DOT edu"),
                        html.Li("bauer AT bloodgroup DOT tch DOT harvard DOT edu"),
                    ],
                    style={"padding": "15px"},
                ),
                html.Div(
                    [
                        html.P("Alternatively, please open an issue on GitHub: "),
                        html.A(
                            f"{GITHUB_LINK}/issues",
                            href=f"{GITHUB_LINK}/issues",
                            target="_blank",
                        ),
                    ]
                ),
            ],
            style={"margin-left": "1%"},
        )
    )
    return f
