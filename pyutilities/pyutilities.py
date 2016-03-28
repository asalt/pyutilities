import re
import time
import sys
import os
import re
import shutil
from datetime import datetime
from collections import defaultdict, namedtuple
from pprint import pprint
from configparser import ConfigParser
import bs4
from bs4 import BeautifulSoup
from Bio import Entrez, Medline
sys.modules['BeautifulSoup'] = bs4
import click
from pyutilities.utils import make_iterable

__version__ = '1.2'
__config__ = '.pyutilities.ini'

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
class Config(object):

    def __init__(self):
        self.email = None
        self.outfile = '-'

pass_config = click.make_pass_decorator(Config, ensure=True)

class CaseConfigParser(ConfigParser):
    def optionxform(self, optionstr):
        return optionstr

parser = CaseConfigParser()
def make_configfile(path=None):
    """Make a configfile with necessary sections.
    Default places it in os.path.expanduser home directory"""
    if path is None:
        path = os.path.expanduser('~')
    parser['user'] = {'email': 'None'}
    with open(os.path.join(path, __config__), 'w') as configfile:
        parser.write(configfile)

def get_parser(path=None):
    """Get the parser for the configfile"""
    if path is None:
        path = os.path.expanduser('~')
    configfile = os.path.join(path, __config__)
    if not os.path.isfile(configfile):
        make_configfile(path)
    parser.read(configfile)
    return parser

def get_stored_email(path=None):
    """Get email from the config file"""
    parser = get_parser(path)
    email = parser['user']['email']
    return email

def update_email(email, path=None):
    """update the email in the config file"""
    parser = get_parser(path)
    parser['user']['email'] = email
    if path is None:
        path = os.path.expanduser('~')
    with open(os.path.join(path, __config__), 'w') as configfile:
        parser.write(configfile)

def check_email(email):
    """Check validity of an input email, then once valid
    update the email in the config file and return the email."""
    pat = re.compile(r"(^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+$)")
    while not pat.match(email):
        email = input('Enter your email address : ')
    return email

@pass_config
def set_email(config, path=None):
    """Set email for NCBI eutilities"""
    email = check_email(config.email)
    if email != get_stored_email(path):
        update_email(email)
    Entrez.email = email

def pmid_fetcher(pmid_list):
    """Fetch pubmed info for a given list of
    pubmed IDs"""
    handle = Entrez.efetch(db='pubmed', id=pmid_list, rettype='medline')
    records = list(Medline.parse(handle))
    return records

def gene_to_pmid(genelist):
    """Returns a dictionary that maps
    gene -> pmid"""
    gene_pmids = dict()
    for gene in genelist:
        handle = Entrez.elink(db='pubmed', dbfrom='gene', id=gene)
        result = handle.read()
        soup = BeautifulSoup(result, 'xml')
        linkset = soup.eLinkResult.LinkSet  # based on ncbi xml structure
        links = [x.Id.text for x in linkset.LinkSetDb.find_all('Link')]
        gene_pmids[gene] = links
    return gene_pmids

def pubmed_to_pubmedinfo(genepmids):
    """Takes in a dictionary with a genes that
    point to a list of pmids. Returns a dictionary
    of genes that point to a list of pubmed info
    as a namedtuple
    """
    PubMedRecord = namedtuple('PubMedRecord',
                              'pmid, title, journal, date, authors, keywords, abstract')
    genelit = defaultdict(list)  # this will get big fast..probably
    for gene in genepmids:  # iterate through each gene and pass each pmid list
        pmids = genepmids.get(gene)  # pmid list
        if pmids:  # if we have a nonempty list
            pmidrecords = pmid_fetcher(pmids)  # returns a list of medline records
            for r in pmidrecords:
                if not isinstance(r.get('OT'), list):
                    r['OT'] = list('')
                record=PubMedRecord(pmid=r.get('PMID'),
                                    title=r.get('TI'),
                                    journal=r.get('JT'),
                                    date=r.get('EDAT')[:10], # don't care about
                                    authors=r.get('AU', ''), # H:M:S
                                    keywords=r.get('OT'),
                                    abstract = r.get('AB', ''))
                genelit[gene].append(record)
    return genelit

def pop_nokeywords(genelit, filter_kws):
    pat = re.compile('|'.join(filter_kws), re.I)
    for gene in genelit:
        genelit[gene] = [x for x in genelit[gene] if pat.search(x.abstract) or pat.search(x.title)]
    return genelit

def read_stepwise(genelit, genedict=None):
    term = shutil.get_terminal_size((80,20))
    columns = term.columns
    user = ''
    for gene in genelit:
        if 'q' in user.lower():
            break
        if genedict:
            gene_abbrev = genedict.get(gene)
        else:
            gene_abbrev = gene
        tot = len(genelit[gene])
        for ix, record in enumerate(genelit[gene]):
            print('Gene search : {}'.format(gene_abbrev))
            print('Record {} of {}'.format(ix+1, tot))
            print(record.pmid, record.title,)
            print(record.journal, record.date)
            pprint(' '.join(record.authors))
            pprint(record.abstract, width=columns-5)
            print('='*24,'\n\n')
            if ix != tot:
                user = input('Enter to continue or q to quit:')
            if 'q' in user.lower():
                break

@click.group()
@click.option('-e', '--email', default=lambda: get_stored_email(),
        help='Email to use for eutilities access')
@click.option('-o', 'outfile', type=click.File('w'), default='-')
@pass_config
def cli(config, email, outfile):
    config.email = email
    config.outfile = outfile

@cli.command(context_settings=CONTEXT_SETTINGS)
@click.option('-k', '--keyword', multiple=True,
              help='Keyword to use for filtering publications. Regular expression not allowed (yet). Multiple keywords allowed.')
@click.argument('geneids', nargs=-1)
@pass_config
def gene_pmid(config, geneids, keyword,):
    """Search for pubmed IDs from an input of geneids and an (optional) search term."""
    #click.echo(geneids))
    #print(keyword)
    set_email(config.email)
    gene_pmids = gene_to_pmid(geneids)
    genelit = pubmed_to_pubmedinfo(gene_pmids)

    if keyword:
        genelit = pop_nokeywords(genelit, keyword)
    read_stepwise(genelit)

@make_iterable(0, str)
def gene_pmid_api(geneids, keywords=None, email=None):
    """An API for gene_pmid
    geneids is one or more geneids to query.
    keywords is a list of keywords to filter by"""
    if email is None:
        email = get_stored_email()
    if email:
        Entrez.email = email
    gene_pmids = gene_to_pmid(geneids)
    genelit = pubmed_to_pubmedinfo(gene_pmids)

    if keyword:
        genelit = pop_nokeywords(genelit, keyword)
    return genelit
