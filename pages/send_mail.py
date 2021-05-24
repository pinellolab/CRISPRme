# import for send mail
# Structure of email.txt
# destinatario
# job link
# date submitted job
# parameters (not yet implemented)

# argv[1] is job directory, eg Results/72C1MNXDWF
# argv[2] is mail config file

import sys
import smtplib
from email.message import EmailMessage
import ssl


def sendMail():
    with open(sys.argv[1] + '/email.txt', 'r') as e:
        # message building
        all_content = e.read().strip().split('--OTHEREMAIL--')
        for em in all_content:

            msg = EmailMessage()
            em = em.strip().split('\n')
            msg['To'] = em[0]
            job_link = em[1]
            date_submission = em[2]
            msg['Subject'] = 'CRISPRme - Job completed'

            # msg['From'] = 'crisprme-job@crisprme.di.univr.it'
            msg['From'] = 'crisprme.job@gmail.com'
            content_email = 'The requested job is completed, visit the following link ' + \
                job_link + ' to view the report.'

            # TODO add Parameters section with date and other parameters
            msg.set_content(content_email)

            print('EMAIL SENT')

            # gmail settings
            # port = 465  # For SSL
            # # Create a secure SSL context
            # context = ssl.create_default_context()
            # # with smtplib.SMTP_SSL("smtp.gmail.com", port, context=context) as server:
            #     # insert data to login into server

            #     server.send_message(msg)


# disabled call until fixed with user personal mail server
# sendMail()
