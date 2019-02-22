import sys
import ChartReport as ChartReport
sys.path.append('../')

import logging
import traceback as tb
# import suds.metrics as metrics
from tests import *
from suds import *
from suds.client import Client
from datetime import datetime

errors = 0
#
# setup_logging()
#
# logging.getLogger('suds.client').setLevel(logging.DEBUG)

url = 'http://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl'

print 'url=%s' % url

#
# create a service client using the wsdl.
#
client = Client(url)
client.wsdl.services[0].setlocation(
	'https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/')

#
# print the service (introspection)
#
# print client


# authenticate user email
print client.service.authenticate('akankshita.dash@u.nus.edu')

# #add a list
inputIds = '31741_at,31734_at,32696_at,37559_at,41400_at,' \
		   '35985_at,39304_g_at,41438_at,35067_at,32919_at,' \
		   '35429_at,36674_at,967_g_at,36669_at,39242_at,39573_at,' \
		   '39407_at,33346_r_at,40319_at,2043_s_at,1788_s_at'
idType = 'AFFYMETRIX_3PRIME_IVT_ID'
listName = 'make_up'
listType = 0
print client.service.addList(inputIds, idType, listName, listType)
print client.service.getAllPopulationNames()
print client.service.getSpecies()
print client.service.setCategories(['KEGG_PATHWAY,GOTERM_CC_DIRECT,GOTERM_MF_DIRECT,GOTERM_BP_DIRECT'])

# getChartReport
thd = 0.1
count = 2
print client.service.getChartReport(thd,count)
# print client.service.getGeneReportCategories()
# print client.service.getListReport()
#
# # getTermClusterReport
# overlap = 3
# initialSeed = 3
# finalSeed = 3
# linkage = 0.5
# kappa = 20
# print client.service.getTermClusterReport(overlap, initialSeed, finalSeed, linkage, kappa)
