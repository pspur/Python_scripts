from timeit import default_timer as timer
import gzip
import json
import elasticsearch
#from elasticsearch import helpers
from pprint import pprint as pp

# only genes which are present in sample db, used to search entire test space quickly
gl = ['GPR153','DRAXIN','TTLL10','MRPL20','TNFRSF4','FBXO44','HNRNPCL3','AGTRAP','ICMT','ACTRT2', 
 'PER3','SLC45A1','MAD2L2','H6PD','MEGF6','ENO1','THAP3','ARHGEF16','TNFRSF8','CTNNBIP1',
 'TNFRSF14','MIIP','PRDM16','TMEM88B','TPRG1L','RERE','TARDBP','PLOD1','GABRD','ACAP3','RER1',
 'CAMTA1','HES4','VWA1','CPSF3L','CPTP','CORT','SCNN1D','ACOT7','ATAD3B','PUSL1','TNFRSF18',
 'CDK11A','DFFA','CCNL2','LZIC','CALML6','KCNAB2','TAS1R3','TNFRSF25','CHD5','SDF4','SLC2A5',
 'PRAMEF11','MFN2','NMNAT1','AURKAIP1','CLSTN1','VPS13D','MMEL1','SLC2A7','ZBTB48','C1orf158',
 'PARK7','FBXO6','PLEKHG5','ATAD3A','KLHL21','NADK','PANK4','SMIM1','ANGPTL7','FAM132A',
 'TMEM240','PHF13','AADACL4','LRRC47','ERRFI1','RNF207','UBE2J2','ATAD3C','GPR157','DHRS3',
 'SLC25A33','MASP2','ESPN','RBP7','PRAMEF12','AJAP1','TNFRSF1B','ANKRD65','ISG15','RNF223',
 'MORN1','PLCH2','NPPA','PGD','C1orf174','MTOR','TP73','TMEM201','MTHFR','TNFRSF9','PRAMEF1',
 'TTC34','PRKCZ','CFAP74','SKI','GNB1','VAMP3','PEX10','KIAA2013','WRAP73','APITD1-CORT',
 'FBXO2','PLEKHN1','SLC35E2','SSU72','UTS2','SAMD11','AADACL3','NPPB','MMP23B','OR4F16',
 'SLC35E2B','PERM1','SRM','KIF1B','MXRA8','B3GALT6','CASZ1','CEP104','UBIAD1','DVL1','TAS1R1',
 'DFFB','PIK3CD','CA6','DNAJC11','C1orf159','OR4F5','HES2','NOC2L','TMEM52','PEX14','C1orf127',
 'EXOSC10','HES5','C1orf167','UBE4B','CCDC27','HES3','RPL22','CDK11B','NOL9','MIB2','SPSB1',
 'FAM213B','AGRN','CLCN6','KLHL17','NPHP4']

not_exonic = ['splicing','ncRNA','UTR5','UTR3','intronic','upstream','downstream','intergenic',
              'upstream;downstream','exonic;splicing','UTR5;UTR3',]

# fake familial relationships for testing
fr = {'f01': ['2157','1406','2162'],
      'f02': ['2208','2214','2221'],}


def build_list_from_refflat(rf):
    gene_list = []
    with gzip.open(rf,'rt') as fin:
        for line in fin:
            cols = line.split('\t')
            if cols[1].startswith('NM_'):
                gene_list.append(cols[0])

    gene_list = list(set(gene_list))
    return(gene_list)


def query_by_gene(gene):
    query_template = """
        {
            "size": 1000,
            "query": {
                "nested" : {
                    "path" : "refGene",
                    "query" : {
                        "bool" : {
                            "must" : [
                                { "match" : {"refGene.refGene_symbol" : "%s"} }
                            ]
                        }
                    }
                }
            }
        }  
    """ %(gene)
    return(query_template)


def analyse_gene(res):
    comp_het = { fid:{'mom':[],'dad':[],'child':[]} for fid in fr.keys() }
    results = {}
    variant_id_map = {}
    
    for doc in res['hits']['hits']:
        if doc['_source']['Func_refGene'] in not_exonic:
            continue

        variant_id_map[doc['_source']['Variant']] = doc['_id']
        samples = doc['_source']['sample']
        sam_list = [sample['sample_ID'] for sample in samples]
        sid_gt = {d['sample_ID']: d['sample_GT'] for d in samples}
        for fid,family in fr.items():
            try:
                momgt = sid_gt[family[0]]
                dadgt = sid_gt[family[1]]
                childgt = sid_gt[family[2]]
            except Exception as e:
                continue
            
            if momgt == '0/0' and dadgt == '0/0' and '1' in childgt:
                if doc['_id'] in results:
                    results[doc['_id']].setdefault('denovo',[]).append(fid)
                else:
                    results[doc['_id']] = {'denovo':[fid]}
                #continue
            elif momgt == '0/1' and dadgt == '0/1' and childgt == '1/1':
                if doc['_id'] in results:
                    results[doc['_id']].setdefault('hom_recess',[]).append(fid)
                else:
                    results[doc['_id']] = {'hom_recess':[fid]}
                #continue
                
            if momgt == '0/1':
                comp_het[fid]['mom'].append(doc['_source']['Variant'])
            if dadgt == '0/1':
                comp_het[fid]['dad'].append(doc['_source']['Variant'])
            if childgt == '0/1':
                comp_het[fid]['child'].append(doc['_source']['Variant'])
            #gene = doc['_source']['refGene'][0]['refGene_symbol']
    
    
    for fid,family in comp_het.items():
        #print(family)
        if len(set(family['mom'] + family['dad'])) <= 2:
            continue
        for var1 in family['mom']:
            if var1 in family['dad']:
                continue
            for var2 in family['dad']:
                if var2 in family['mom']:
                    continue
                if var1 in family['child'] and var2 in family['child']:
                    #print(var1,var2)
                    if variant_id_map[var1] in results:
                        results[variant_id_map[var1]].setdefault('comp_het',[]).append({fid:var2})
                    else:
                        results[variant_id_map[var1]] = {'comp_het':[{fid:var2}]}
                    if variant_id_map[var2] in results:
                        results[variant_id_map[var2]].setdefault('comp_het',[]).append({fid:var1})
                    else:
                        results[variant_id_map[var2]] = {'comp_het':[{fid:var1}]}
                        
    #pp(results)
    #print()
    return(results)    
            

def create_update_body(trio_output):
    comphet_string = ''
    denovo_string = ''
    homrecess_string = ''
    for var_type,var_info in trio_output.items():
        if var_type == 'comp_het':
            comphet_string = '"comp_het" : [' 
            for e in var_info:
                for family_id,assoc_var in e.items():
                    comphet_string += '{{"family" : "{}",\n"assocvar" : "{}"}}'.format(family_id,assoc_var)
                    
            comphet_string = comphet_string.replace('}{','},{')
            comphet_string += ']'

        elif var_type == 'denovo':
            denovo_string = '"denovo" : {}'.format(str(var_info).replace("'",'"'))
        
        elif var_type == 'hom_recess':
            homrecess_string = '"hom_recess" : {}'.format(str(var_info).replace("'",'"'))
    
    update_string = ','.join(filter(None,[comphet_string,denovo_string,homrecess_string]))
    
    body = """
        {
            "doc" : {
                %s
            }
        }"""%(update_string)
    #print(body)
    #print(json.loads(body))
    return(json.loads(body))


def main():
    gene_list = build_list_from_refflat('refFlat.txt.gz')
    start = timer()
    es = elasticsearch.Elasticsearch(host='199.109.192.6', port='9200')
    #for gene in gene_list:
    for gene in gl:
    #for gene in ['MEGF6']:
        gene_query = query_by_gene(gene)
        res = es.search(index='sim_sam',doc_type='sim',body=gene_query)
        if int(res['hits']['total']) > 0:
            results = analyse_gene(res)
            if results:
                #pp('final results')
                #pp(results)
                for esid in results:
                    print(esid)
                    es.update(index='sim_sam', doc_type='sim', id=esid,
                              body=create_update_body(results[esid]))
                    print(create_update_body(results[esid]))
                    print()
    #ok, result = helpers.bulk(es, assign_varstatus(es,f), chunk_size=1000)
    end = timer()
    #print(end-start)

if __name__=='__main__':
    main()
