"""Send email to notify the user that CRISPRme analysis was completed and
the corresponding results are available to be visualized and explored.

Structure of email.txt
    - Contact
    - link to CRISPRme job
    - Job submission date

TODO: avoid shell and call send_mail() in other python scripts
TODO: add run parameters to mail (job date + other params)
"""

import os
import sys
import ssl
import smtplib
from email.message import EmailMessage

PSW = "tlgylyubfalteffs"


# from .pages_utils import MAIL_SUBJECT, MAIL_SENDER, SSL_PORT
# import pages_utils


def send_mail() -> None:
    """Send the email to notify the user that his/her CRISPRme is
    complete.

    ...

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    with open(sys.argv[1] + "/email.txt", "r") as e:
        # message building
        all_content = e.read().strip().split("--OTHEREMAIL--")
        for em in all_content:
            # read mail content from file
            msg = EmailMessage()
            em = em.strip().split("\n")
            msg["To"] = em[0]
            job_link = em[1]
            date_submission = em[2]
            msg["Subject"] = "CRISPRme - Job completed"

            msg["From"] = "crisprme.job@gmail.com"
            content_email = (
                "The requested job is completed, visit the following link "
                + job_link
                + " to view the report."
            )

            # TODO add Parameters section with date and other parameters
            msg.set_content(content_email)

            print("EMAIL SENT")
            # gmail settings
            port = 465  # For SSL (used for gmail account)
            # Create a secure SSL context
            context = ssl.create_default_context()
            with smtplib.SMTP_SSL("smtp.gmail.com", port, context=context) as server:
                server.login("crisprme.job@gmail.com", PSW)
                server.send_message(msg)


# call the above funtion
# TODO: avoid calling through bash


def main():
    """Call the above function and send mails to notify user about
    CRSPRme job completion.
    """

    # try:
    send_mail()
    # sys.stderr.write("Mail service not available")
    # except:
    #     sys.stderr.write("No email request by the user")


# entry point
if __name__ == "__main__":
    main()
