import pandas as pd
from gql import gql, Client
from gql.transport.aiohttp import AIOHTTPTransport
import os


def clinvar_columns():
    keyD = {'AF_ESP': 'af_esp',
            'AF_EXAC': 'af_exac',
            'AF_TGP': 'af_tgp',
            'ALLELEID': 'allele_id',
            'CLNDISDB': 'crossrefs',
            'CLNDN': 'trait',
            'CLNHGVS': 'hgvs_id',
            'CLNREVSTAT': 'clinvar_review_status',
            'CLNSIG': 'clinvar_significance',
            'CLNSIGCONF': 'clinvar_confidence',
            'CLNSIGSCV': 'clinvar_significance_source',
            'CLNVC': 'clinvar_variant_type',
            'CLNVCSO': 'clinvar_sequence_ontology_id',
            'CLNVI': 'other_ids',
            'GENEINFO': 'gene_ids',
            'MC': 'molecular_consequence',
            'ORIGIN': 'germline/acquired',
            'RS': 'rsid',
            'ONC': 'oncogenicity',
            'ONCDISDB': 'oncogenicity_crossref',
            'ONCDN': 'onco_gene_name',
            'ONCREVSTAT': 'onco_review_status',
            'ONCSCV': 'onco_significance_source'}
    return (keyD)


def get_clinvar_variants(path, outdir, chrom, start, end):
    statement = f'tabix {path} {chrom}:{start}-{end} \
    > {outdir}/{chrom}_{start}_{end}.tsv'
    os.system(statement)
    tab = pd.read_csv(f"{outdir}/{chrom}_{start}_{end}.tsv",
                      sep="\t", header=None)
    keyD = clinvar_columns()
    DD = dict()
    keys = set()
    for ID, pos, seg in zip(tab[2], tab[1], tab[7]):
        D = dict()
        ss = seg.split(";")
        for s in ss:
            key, val = s.split("=")
            D[keyD[key]] = val
            keys.add(key)
        D['pos'] = pos
        DD[str(ID)] = D
    df = pd.DataFrame(DD).T
    df['clinvar_variation_id'] = df.index.values
    return (df)



async def get_gnomad_variants(chrom, start, end):
    # Setup for the API call

    # This makes a "transport point" which sends the GraphQL request over HTTP
    # to the API endpoint
    transport = AIOHTTPTransport(url="https://gnomad.broadinstitute.org/api",
                                 ssl=False)

    # And the interface with the GraphQL API
    client = Client(transport=transport, fetch_schema_from_transport=True)

    # Column names for results table
    cols = ['variant_id', 'chrom', 'pos', 'ref', 'alt', 'AC', 'AN', 'AF',
            'rsids',
            'flags', 'gene_symbol', 'consequence', 'lof', 'lof_flags',
            'lof_filter', 'hgvsc', 'hcvsp']

    cols2 = ['variant_id', 'chrom', 'pos', 'ref', 'alt',
             'clinical_significance', 'clinvar_variation_id',
             'gold_stars', 'review_status', 'major_consequence']

    # database query string
    string = f"""
    query VariantsInRegion {{
      region(chrom: "{chrom}", start: {start}, stop: {end},
      reference_genome: GRCh38) {{
        variants(dataset: gnomad_r4) {{
          variant_id chrom pos ref alt
          exome {{ ac an af }}
          genome {{ ac an af }}
          rsids
          flags
          transcript_consequence {{
            gene_symbol
            transcript_id
            consequence_terms
            lof
            lof_flags
            lof_filter
            hgvsc
            hgvsp
            sift_prediction
            polyphen_prediction
          }}
        }}
        clinvar_variants {{
          variant_id chrom pos ref alt
          clinvar_variation_id
          clinical_significance
          gold_stars
          review_status
          major_consequence
        }}
      }}
    }}

    """
    # Execute the query on the transport
    query = gql(string)
    result = await client.execute_async(query)

    # Extract the results from the JSON
    rows = []
    for v in result['region']['variants']:
        if v['exome']:
            row = [v['variant_id'],
                   v['chrom'],
                   v['pos'],
                   v['ref'],
                   v['alt'],
                   v['exome']['ac'],
                   v['exome']['an'],
                   v['exome']['af'],
                   ",".join(v['rsids']),
                   ",".join(v['flags']),
                   v['transcript_consequence']['gene_symbol'],
                   ",".join(v['transcript_consequence']['consequence_terms']),
                   v['transcript_consequence']['lof'],
                   v['transcript_consequence']['lof_flags'],
                   v['transcript_consequence']['lof_filter'],
                   v['transcript_consequence']['hgvsc'],
                   v['transcript_consequence']['hgvsp']]
        else:
            row = [v['variant_id'],
                   v['chrom'],
                   v['pos'],
                   v['ref'],
                   v['alt'],
                   v['genome']['ac'],
                   v['genome']['an'],
                   v['genome']['af'],
                   ",".join(v['rsids']),
                   ",".join(v['flags']),
                   v['transcript_consequence']['gene_symbol'],
                   ",".join(v['transcript_consequence']['consequence_terms']),
                   v['transcript_consequence']['lof'],
                   v['transcript_consequence']['lof_flags'],
                   v['transcript_consequence']['lof_filter'],
                   v['transcript_consequence']['hgvsc'],
                   v['transcript_consequence']['hgvsp']]
        rows.append(row)

    # Make a dataframe
    dfv = pd.DataFrame(rows, columns=cols)

    # The clinvar variants are separate for some reason
    rows = []
    for v in result['region']['clinvar_variants']:
        row = [v['variant_id'],
               v['chrom'],
               v['pos'],
               v['ref'],
               v['alt'],
               v['clinical_significance'],
               v['clinvar_variation_id'],
               v['gold_stars'],
               v['review_status'],
               v['major_consequence']]
        rows.append(row)
    dfc = pd.DataFrame(rows, columns=cols2)

    df = dfv.merge(dfc, 'outer')

    # Assign an ID to the region
    df['ID'] = f"{chrom}:{start}-{end}"

    return (df)
