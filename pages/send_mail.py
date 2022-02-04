# import for send mail
# Structure of email.txt
# destinatario
# job link
# date submitted job
# parameters (not yet implemented)

# argv[1] is job directory, eg Results/72C1MNXDWF

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

            msg['From'] = '<SENDER OF RESULT MAIL>'
            content_email = 'The requested job is completed, visit the following link ' + \
                job_link + ' to view the report.'

            # TODO add Parameters section with date and other parameters
            msg.set_content(content_email)

            print('EMAIL SENT')

            # univr settings
            # server = smtplib.SMTP(host="smtp.univr.it", port=25)
            # server.ehlo_or_helo_if_needed()
            # server.send_message(msg, from_addr='crisprme-job@crisprme.di.univr.it')
            # server.quit()

            # gmail settings
            port = 465  # For SSL (used for gmail account)
            # Create a secure SSL context
            context = ssl.create_default_context()
            with smtplib.SMTP_SSL("<INSERT SMTP SERVER>", port, context=context) as server:
                server.login("<USERNAME>", '<PASSWORD>')
                server.send_message(msg)


# function call
try:
    sendMail()
except:
    print('NO mail requested by user')
