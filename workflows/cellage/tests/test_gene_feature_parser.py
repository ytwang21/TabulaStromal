from scripts.common import parse_feature


def test_longest_verified_celltype_prefix_preserves_underscores():
    c,g=parse_feature('Type_A_pan_GENE_WITH_UNDERSCORE',['Type','Type_A','Type_A_pan'])
    assert c=='Type_A_pan' and g=='GENE_WITH_UNDERSCORE'
    assert parse_feature('Unmatched_GENE',['Type_A'])==(None,None)
