"""Send email to notify the user that CRISPRme analysis was completed and
the corresponding results are available to be visualized and explored.

Structure of email.txt
    - Contact
    - link to CRISPRme job
    - Job submission date

TODO: avoid shell and call send_mail() in other python scripts
TODO: add run parameters to mail (job date + other params)
"""

from .pages_utils import MAIL_SUBJECT, MAIL_SENDER, SSL_PORT

from email.message import EmailMessage

import smtplib
import ssl
import sys
import os


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

    try:
        # in future the folowing 2 lines should be removed
        job_directory = sys.argv[1]
        assert os.path.exists(job_directory)
        with open(os.path.join(job_directory), "email.txt", mode="r") as handle:
            # build mail message
            mail_content = handle.read().strip().split("--OTHEREMAIL--")
            for em in mail_content:
                msg = EmailMessage()  # create the email object
                em = em.strip().split("\n")  # get mail fields
                msg["To"] = em[0]
                job_link = em[1]
                submission_date = em[2]
                msg["Subject"] = MAIL_SUBJECT
                msg["From"] = MAIL_SENDER
                mail_body = "Your CRISPRme job is completed.\n\n"
                mail_body += f"visit the following link to view the report:\n{job_link}"
                msg.set_content(mail_body)
                print("EMAIL SENT")
                # univr settings
                # server = smtplib.SMTP(host="smtp.univr.it", port=25)
                # server.ehlo_or_helo_if_needed()
                # server.send_message(msg, from_addr='crisprme-job@crisprme.di.univr.it')
                # server.quit()
                # gmail settings
                port = SSL_PORT  # For SSL (used for gmail account)
                # Create a secure SSL context
                context = ssl.create_default_context()
                with smtplib.SMTP_SSL(
                    "<INSERT SMTP SERVER>", port, context=context
                ) as server:
                    server.login("<USERNAME>", "<PASSWORD>")
                    server.send_message(msg)
    except OSError as e:
        raise e


# call the above funtion
# TODO: avoid calling through bash
def main():
    """Call the above function and send mails to notify user about
    CRSPRme job completion.
    """

    try:
        send_mail()
    except:
        sys.stderr.write("No email request by the user")


# entry point
if __name__ == "__main__":
    main()
