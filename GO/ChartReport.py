from suds.client import Client
import os, pandas as pd
def DAVIDenrich(listF, idType, bgF='', resF='', bgName = 'Background1', email = 'akankshita.dash@u.nus.edu', listName='List1', category = '', thd=0.1, ct=2):
    inputListIds = []
    inputBgIds = []
    if len(listF) > 0 and os.path.exists(listF):
        inputListIds = ','.join(open(listF).read().split('\n'))
        print 'List loaded.'
    else:
        print 'Empty list.'

    flagBg = False
    if len(bgF) > 0 and os.path.exists(bgF):
        inputBgIds = ','.join(open(bgF).read().split('\n'))
        flagBg = True
        print 'Use file background.'
    else:
        print 'Use default background.'

    client = Client('http://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl')
    client.wsdl.services[0].setlocation('https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/')
    print 'User Authenticating',client.service.authenticate(email)

    listType = 0
    print 'Percentage mapped(list):', client.service.addList(inputListIds,idType,listName,listType)
    if flagBg:
        listType = 1
        print 'Percentage mapped(background):', client.service.addList(inputBgIds,idType,bgName,listType)

    print 'Use categories:', client.service.setCategories(category)
    chartReport = client.service.getChartReport(thd,ct)
    chartRow = len(chartReport)
    print 'Total chart records:',chartRow
    
    if len(resF) == 0 or not os.path.exists(resF):
        if flagBg:
            resF = listF[:-4] + '_withBG_chartReport.txt'
        else:
            resF = listF[:-4] + '_chartReport.txt'
    with open(resF, 'w') as fOut:
        fOut.write('Category\tTerm\tCount\t%\tPvalue\tGenes\tList Total\tPop Hits\tPop Total\tFold Enrichment\tBonferroni\tBenjamini\tFDR\n')
        for row in chartReport:
            rowDict = dict(row)
            categoryName = str(rowDict['categoryName'])
            termName = str(rowDict['termName'])
            listHits = str(rowDict['listHits'])
            percent = str(rowDict['percent'])
            ease = str(rowDict['ease'])
            Genes = str(rowDict['geneIds'])
            listTotals = str(rowDict['listTotals'])
            popHits = str(rowDict['popHits'])
            popTotals = str(rowDict['popTotals'])
            foldEnrichment = str(rowDict['foldEnrichment'])
            bonferroni = str(rowDict['bonferroni'])
            benjamini = str(rowDict['benjamini'])
            FDR = str(rowDict['afdr'])
            rowList = [categoryName,termName,listHits,percent,ease,Genes,listTotals,popHits,popTotals,foldEnrichment,bonferroni,benjamini,FDR]
            fOut.write('\t'.join(rowList)+'\n')
        print 'write file:', resF, 'finished!'


email = raw_input("Enter authentication e-mail: ")
gene_path = raw_input("Enter file path: " )
background = raw_input("Enter background path: " )
clusters = [6]
num_clusters='6'
for i in range(1,int(num_clusters)+1):
    list_genes = gene_path + '/AccNums/AccNum'+str(i)+'.txt'
    print(list_genes)
    DAVIDenrich(listF = list_genes,bgF = background, idType = 'REFSEQ_MRNA',  email=email,listName = list_genes, category = 'KEGG_PATHWAY,GOTERM_CC_DIRECT,GOTERM_MF_DIRECT,GOTERM_BP_DIRECT')

for num_clusters in clusters:
    common_go_file = gene_path[:-2]+'GO_'+str(num_clusters)+'.txt'
    for i in range(1,int(num_clusters)+1):
        file = gene_path+'/AccNums/AccNum'+str(i)+'_withBG_chartReport.txt'
        df = pd.read_csv(file, sep='\t')
        # df_write = df_write.append(['Cluster'+str(i)], ignore_index=True)
        empty = pd.Series(['-'] * 13, index=['Category', 'Term', 'Count', '%', 'Pvalue', 'Genes', 'List Total',
                                         'Pop Hits', 'Pop Total', 'Fold Enrichment', 'Bonferroni', 'Benjamini',
                                         'FDR'])
        df_bp = df[df.Category == 'GOTERM_BP_DIRECT'].head(10).append(empty, ignore_index=True)
        df_cc = df[df.Category == 'GOTERM_CC_DIRECT'].head(10).append(empty, ignore_index=True)
        df_mf = df[df.Category == 'GOTERM_MF_DIRECT'].head(10).append(empty, ignore_index=True)
        df_kegg = df[df.Category == 'KEGG_PATHWAY'].head(10).append(empty, ignore_index=True)
        header = ['Cluster' + str(i)] + ['*'] * 12
        print(header)
        df_write = pd.DataFrame(columns=header)
        df_write.to_csv(common_go_file, sep='\t', index=False, mode='a')
        df_write = df_bp.append(df_cc.append(df_mf.append(df_kegg)))
        df_write.to_csv(common_go_file, sep='\t', index=False, mode='a')