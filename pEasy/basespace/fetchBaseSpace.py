#!/usr/bin/env python


'''
pip install pycurl
git clone git@github.com:basespace/basespace-python-sdk.git
'''

import sys
import BaseSpacePy

from BaseSpacePy.api.BaseSpaceAPI import BaseSpaceAPI

# REST server information and user access_token
BaseSpaceUrl = 'https://api.basespace.illumina.com/'
version = 'v1pre3'
client_key = '2e5ab9db34184f11996ab37b6de02b71'
client_secret = '2ef67eb9d75c45b28ee66fe09f60b5a7'
AppSessionId = 'test_session'
accessToken = "2ae66c45785b430b9e749737cd6cc61c"

# First, create a client for making calls for this user session
BSapi = BaseSpaceAPI(client_key, client_secret, BaseSpaceUrl, version, AppSessionId, AccessToken=accessToken)


#if __name__=="__main__":
#    print >> sys.stderr, "Yippie!"