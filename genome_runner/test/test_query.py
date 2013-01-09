import test_data
from nose.tools import with_setup,eq_
from genome_runner import query


#_Enrichment_Par = namedtuple("Enrichment_Par","a b A B Background n flt genome genome_fn organism obs background")

#_Enrichment["A","B","nA","nB","","expected","p_value","obsprox","expprox","pybed_p_value",
#"pybed_expected","jaccard_observed","jaccard_p_value",
#"jaccard_expected","proximity_p_value","kolmogorov_p_value","hypergeometric_p_value"])

#run_enrichments(id, f, gfeatures,background, niter, name, score, strand,organism,run)
data = None
def setup():
    global data
    data = test_data.get_test_data()
def teardown():
    pass

@with_setup(setup,teardown)
def test_run_enrichments():
    d = data
    results = query.run_enrichments(0,d['foi_duplicate_test'],[d['foi_duplicate_test']],
                                    "",10,None,None,None,'hg19',['pvalue'])
    # creates a new namedtuple with the 'correct' values to be compared.
    # Values that are not passed into _replace will be assumed to be correct. 
    # A: Name of the features of interest
    # B: Nme of the genomic feature
    # nA: the number of Feature of Interest
    # nB: the number of Genomic Features
    answers = results[0]._replace(nA = 3,nB=3,p_value=0.0,pybed_expected='NA')
    chk_results(results[0],answers)

@with_setup(setup,teardown)
def test_overlap_run_enrichments():
    d = data
    results = query.run_enrichments(0,d['foi_overlap_test'],[d['gf_overlap_test']],"",
                                    10,None, None, None, 'hg19',['pvalue'])
    answers = results[0]._replace(nA = 22,nB=5,observed=17,p_value=0.0,pybed_expected='NA')
    chk_results(results[0],answers)

def chk_results(results,answers):
    '''Takes in two query._Enrichment tuples.  results contains the results from the run.
    answers contains what the actual results should be. gives a description of which results
    failed.
    '''
    r = results
    a = answers
    for i in range(len(results)):
        eq_(r[i],a[i],"{} should be: {} instead is: {}".format(r._fields[i],a[i],r[i]))
