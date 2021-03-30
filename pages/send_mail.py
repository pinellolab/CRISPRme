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
            msg['to'] = em[0]
            job_link = em[1]
            date_submission = em[2]
            msg['Subject'] = 'CRISPRitz - Job completed'

            msg['From'] = 'admin@crispritz.di.univr.it'
            content_email = 'The requested job is completed, visit the following link ' + \
                job_link + ' to view the report.'

            # TODO add Parameters section with date and other parameters
            msg.set_content(content_email)

            context = ssl.SSLContext(ssl.PROTOCOL_TLS)

            server = smtplib.SMTP("smtp.univr.it")
            server.send_message(msg, from_addr='admin@crispritz.di.univr.it')

            # server = smtplib.SMTP('smtp.univr.it',25)
            # server = smtplib.SMTP('smtp-mail.outlook.com', 587)

            #server = smtplib.SMTP_SSL("smtp.live.com",587)
            # for example:
            #server = smtplib.SMTP_SSL("smtp.libero.it", port=465)
            # #start connection
            # server.ehlo()
            # server.starttls(context=context)
            # server.ehlo()
            # #login and send message
            # # server.login("test.cri@hotmail.com", "univrCrispritz")
            # server.login('admin@crispritz.di.univr.it')

            # server.send_message(msg)
            # #close connection
            # server.quit()
