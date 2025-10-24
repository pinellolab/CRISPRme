"""Display the contact page on CRISPRme website.

Contacts list:
    - Manuel Tognon, PhD <manuel.tognon@univr.it> <mtognon@mgh.harvard.edu>
    - Daniel Bauer, MD, PhD <bauer@bloodgroup.tch.harvard.edu>
    - Rosalba Giugno, PhD <rosalba.giugno@univr.it>
    - Luca Pinello, PhD <lpinello@mgh.harvard.edu>


Contacts list for bug reports and maintainers:
    - Manuel Tognon, PhD <manuel.tognon@univr.it> <mtognon@mgh.harvard.edu>
    - Luca Pinello, PhD <lpinello@mgh.harvard.edu>

CRISPRme developers:
    - Samuele Cancellieri, PhD <samuele.cancellieri@univr.it>
    - Elia Dirupo
    - Francesco Masillo
    - Manuel Tognon, PhD <manuel.tognon@univr.it> <mtognon@mgh.harvard.edu>
    
Current CRISPRme maintainers:
    - Manuel Tognon, PhD <manuel.tognon@univr.it> <mtognon@mgh.harvard.edu>
    - Luca Pinello, PhD <lpinello@mgh.harvard.edu>
"""

from .pages_utils import GITHUB_LINK

from typing import List

import dash_html_components as html


# NOTE: the following list could be changed during time
INFOMICS_DEVS = ["Manuel Tognon", "Nicola Bombieri", "Rosalba Giugno"]
INFOMICS_AFF = "Department of Computer Science, University of Verona, Italy"
PINELLOLAB_DEVS = ["Luca Pinello"]
PINELLOLAB_AFF = (
    "Molecular Pathology Unit and Center for Cancer Research, Massachusetts " 
    "General Hospital, Charlestown, MA, USA; Department of Pathology, Harvard " 
    "Medical School, Boston, MA, USA; Broad Institute of MIT and Harvard, " 
    "Cambridge, MA, USA"
)
BAUERLAB_DEVS = ["Linda Lin", "Daniel Bauer"]
BAUERLAB_AFF = (
     "Division of Hematology/Oncology, Boston Children's Hospital, Department of " 
     "Pediatric Oncology, Dana-Farber Cancer Institute, Harvard Stem Cell " 
     "Institute, Broad Institute, Department of Pediatrics, Harvard Medical " 
     "School, Boston, MA, USA"
)


def format_devs_list(devs_list: List[str]) -> str:
    """Format a list of developers' names into a comma-separated string.

    This function takes a list of strings representing developers' names and
    joins them into a single string, separated by commas.

    Args:
        devs_list: A list of strings, where each string is a developer's name.

    Returns:
        A string containing the comma-separated list of developer names.
    """
    return ",".join(devs_list)

def lab_devs(labinfo: str, url: str) -> html.Li:
    """Format lab information and URL for display.

    This function creates a list item containing the lab information and a
    clickable link to the provided URL.

    Args:
        labinfo: The lab information string.
        url: The URL of the lab's website.

    Returns:
        An html.Li element containing the formatted lab information and link.
    """
    return html.Li([f"{labinfo} - ", html.A(f"{url}", href=url, target="_blank")])

def devs() -> html.Ul:
    """Create an HTML unordered list of CRISPRme developers.

    This function formats the developers' information and affiliations,
    including links to their respective labs' websites.

    Returns:
        An html.Ul element containing the formatted list of developers.
    """
    # developers format
    infomics = f"{format_devs_list(INFOMICS_DEVS)}. {INFOMICS_AFF}"
    pinellolab = f"{format_devs_list(PINELLOLAB_DEVS)}. {PINELLOLAB_AFF}"
    bauerlab = f"{format_devs_list(BAUERLAB_DEVS)}. {BAUERLAB_AFF}"
    return html.Ul(
        [
            lab_devs(infomics, "https://infomics.github.io/InfOmics/"),
            lab_devs(pinellolab, "http://pinellolab.org/"),
            lab_devs(bauerlab, "http://bauerlab.org/"),
        ],
        style={"padding": "15px"},
    )

def emails() -> html.Ul:
    """Create an HTML unordered list of email contacts.

    This function generates an HTML list of email addresses with clickable
    mailto links for each contact person.

    Returns:
        An html.Ul element containing the formatted list of email addresses.
    """
    rgmail = "rosalba.giugno@univr.it"
    lpmail = "lpinello@mgh.harvard.edu"
    dbmail = "bauer@bloodgroup.tch.harvard.edu"
    mtmail = "manuel.tognon@univr.it"
    return html.Ul(
        [
            html.Li(["Rosalba Giugno - ", html.A(f"{rgmail}", href=f"mailto:{rgmail}")]),
            html.Li(["Luca Pinello - ", html.A(f"{lpmail}", href=f"mailto:{lpmail}")]),
            html.Li(["Daniel Bauer - ", html.A(f"{dbmail}", href=f"mailto:{dbmail}")]),
            html.Li(["Manuel Tognon - ", html.A(f"{mtmail}", href=f"mailto:{mtmail}")]),
        ],
        style={"padding": "15px"},
    )


def contact_page() -> List:
    """Create the layout for the contact page.

    This function generates the HTML structure for the contact page,
    including the list of developers, contact emails, and a link
    to the GitHub issue tracker.

    Returns:
        A list of HTML elements representing the contact page layout.
    """
    return html.Div(
        [
            html.H3("Contacts"),
            html.P("CRISPRme developers:"),
            devs(),
            html.P("Please send any comment or bug to:"),
            emails(),
            html.P("Alternatively, please open an issue on GitHub: "),
            html.A(f"{GITHUB_LINK}/issues", href=f"{GITHUB_LINK}/issues", target="_blank")
        ],
        style={"margin-left": "1%"},
    )
