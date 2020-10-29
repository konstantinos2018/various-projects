#!/usr/bin/env python
# -*- coding: utf-8 -*-
# === IMPORTS ===
import cgi;
import cgitb;cgitb.enable()
import sys, os
import time, datetime
import smtplib as smt
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

print "Content-Type: text/html"
print "" #use this double quote print statement to add a blank line in the script
print '<meta charset="utf-8">'

# The following works for GMAIL SENDER only
# Message info
sender = 'put-sender-mail-here'
password = 'put-sender-password-here'
receiver = 'put-receiver-mail-here'

msg = MIMEMultipart('alternative')
msg['Subject'] = u'subject-here'.encode('utf-8')
msg['From'] = u'sender-name <{0}>'.format(sender).encode('utf-8')
msg['To'] = u'receiver-name <{0}>'.format(receiver).encode('utf-8')


# Define intro word
t = datetime.datetime.now()
if (t.hour >= 23) and (t.hour <= 1):
  in_word = u'Καληνύχτα'
elif (t.hour > 1) and (t.hour <= 5):
  in_word = u'Καλό ξημέρωμα'
elif (t.hour > 5) and (t.hour < 12):
  in_word = u'Καλημέρα'
elif (t.hour >= 12) and (t.hour <= 16):
  in_word = u'Καλησπέρα'
elif (t.hour > 16) and (t.hour <= 19):
  in_word = u'Καλό απόγευμα'
elif (t.hour > 19) and (t.hour < 23 ):
  in_word = u'Καλό βράδυ'

# Compute time period until critical timestamp
t_leave = datetime.datetime(2020, 12, 25, 15, 0, 0, 0)
t_diff = t_leave - t

# convert time to minutes
t_mins = t_diff.days*24*60 + t_diff.seconds/60


# Create message
html = u"""\
<html>
  <head></head>
  <body>
    <p>{0},</p>
    <p>
       Hey, hi! <b>my friend</b>;<br>
       These are example emojis &#128536;&#128536;&#128536; &#128150; &#128522;<br><br>
	   {1} minutes left
    </p>
    <p>Regards,<br>
    Κώστας</p>
  </body>
</html>
""".format(in_word, t_mins)

message = MIMEText(html.encode('utf-8'), 'html', 'utf-8')
msg.attach(message)

# Send message
# Get current datetime
print '<p><b>Second</b>: {0}</p><br>'.format(time.ctime())

# Send mail
server = smt.SMTP('smtp.gmail.com', 587)
server.ehlo()
server.starttls()
server.ehlo()
server.login(sender, password)

server.sendmail(sender, receiver, msg.as_string())
server.quit()
