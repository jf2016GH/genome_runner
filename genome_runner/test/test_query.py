import test_data
from nose.tools import with_setup,eq_
from nose.tools import assert_almost_equal as almost_eq
from genome_runner import query


#_Enrichment_Par = namedtuple("Enrichment_Par","a b A B Background n flt genome genome_fn organism obs background")

#_Enrichment["A","B","nA","nB","observed","expected","p_value","obsprox","expprox","pybed_p_value",
#"pybed_expected","jaccard_observed","jaccard_p_value",
#"jaccard_expected","proximity_p_value","kolmogorov_p_value","hypergeometric_p_value"])

#run_enrichments(id, f, gfeatures,background, niter, name, score, strand,organism,run)

# available tests are : pvalue,pybedtool,jaccard,kolmogorov,proximity,hypergeometric
d = None
def setup():
    global d
    d = test_data.get_test_data()
def teardown():
    pass

@with_setup(setup,teardown)
def test_runenrichment_foiandgfcount():
    ''' Tests whether only one hit is returned per FOI
    '''
    results = query.run_enrichments(0,d['foi_duplicate_test'],[d['foi_duplicate_test']],
                                    "",10,None,None,None,'hg19',['pvalue'])
    answers = results[0]._replace(nA = 3,nB=3,observed=3)

@with_setup(setup,teardown)
def test_runenrichment_sameregions():
    # Runs the test:
    # 0 is the id of the run, d is a dictionary (pass in name of test file created in test_data.py to run),
    # third argument is a list of test files to run,background, number of iterations, 
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
def test_runenrichment_diffchromosome():
    ''' Tests whether  a GFs on diff chrom than FOI are not counted as hits'''
    results = query.run_enrichments(0,d['same_chrom_2'],[d['gf_chrom2_2']],"",
                                    10,None, None, None, 'hg19',['pvalue','jaccard','proximity'])
    answers = results[0]._replace(nA=6,nB=1, observed=0,jaccard_observed=0,obsprox="NA")
    chk_results(results[0],answers)

@with_setup(setup,teardown)
def test_runenrichment_chromspecific():
    ''' Tests where GF hits are only returned for FOI on the same chromosome
    '''
    results = query.run_enrichments(0,d['diff_chrom_2'],[d['gf_chrom2_2']],"",
                                    10,None, None, None, 'hg19',['pvalue','jaccard','proximity'])
    answers = results[0]._replace(nA=7, nB=1,observed=1,obsprox=300.666666667)
    chk_results(results[0],answers)

@with_setup(setup,teardown)
def test_runenrichment_strandspecific():
    ''' Tests strand specific tests '''

    # negative strand only
    results = query.run_enrichments(0,d['foi_strand_3'],[d['gf_strand_3']],"",
                                    10,None, None, "-", 'hg19',['pvalue'])
    answers = results[0]._replace(nA=5,nB=1,observed=2)
    chk_results(results[0],answers)

    # positive strand only
    results = query.run_enrichments(0,d['foi_strand_3'],[d['gf_strand_3']],"",
                                    10,None, None, "+", 'hg19',['pvalue'])
    answers = results[0]._replace(nA=5,nB=1,observed=3)
    chk_results(results[0],answers)

    # negative strand only for proximity
    results = query.run_enrichments(0,d['foi_strand_prox'],[d['gf_strand_prox']],"",
                                    10,None, None, "-", 'hg19',['proximity'])
    answers = results[0]._replace(nA=1,nB=1,obsprox=201)
    chk_results(results[0],answers)

    # positive strand only for proximity
    results = query.run_enrichments(0,d['foi_strand_prox'],[d['gf_strand_prox']],"",
                                    10,None, None, "+", 'hg19',['proximity'])
    answers = results[0]._replace(nA=1,nB=1,obsprox=101)
    chk_results(results[0],answers)



def chk_results(results,answers):
    '''Takes in two query._Enrichment tuples.  results contains the results from the run.
    answers contains what the actual results should be. gives a description of which results
    failed
.    '''
    r = results
    a = answers
    for i in range(len(results)):
        if is_number(r[i]) and is_number(a[i]):
            almost_eq(float(r[i]),float(a[i]),7,"{} should be: {} instead is: {}".format(r._fields[i],a[i],r[i]))
        else:
            eq_(r[i],a[i],"{} should be: {} instead is: {}".format(r._fields[i],a[i],r[i]))


def is_number(s):
    if s==None: return None 
    try:
        float(s)
        return True
    except ValueError:
        return False
