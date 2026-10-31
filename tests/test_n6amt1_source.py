from pathlib import Path
import sys
import unittest
import xml.etree.ElementTree as ET
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from audit_n6amt1_source import source_attributes,require_identity


class SourceAuditGuards(unittest.TestCase):
    def test_missing_age_is_not_imputed_from_tissue_or_label(self):
        node=ET.fromstring('<SAMPLE><SAMPLE_ATTRIBUTES><SAMPLE_ATTRIBUTE><TAG>genotype</TAG><VALUE>CTL</VALUE></SAMPLE_ATTRIBUTE><SAMPLE_ATTRIBUTE><TAG>tissue</TAG><VALUE>testis</VALUE></SAMPLE_ATTRIBUTE></SAMPLE_ATTRIBUTES></SAMPLE>')
        attrs=source_attributes(node);self.assertNotIn('age',attrs);self.assertEqual(attrs['genotype'],'CTL')

    def test_contact_fields_excluded(self):
        node=ET.fromstring('<SAMPLE><SAMPLE_ATTRIBUTES><SAMPLE_ATTRIBUTE><TAG>biomaterial_provider</TAG><VALUE>private</VALUE></SAMPLE_ATTRIBUTE></SAMPLE_ATTRIBUTES></SAMPLE>')
        self.assertEqual(source_attributes(node),{})

    def test_wrong_biosample_rejected(self):
        experiment=ET.fromstring('<EXPERIMENT accession="SRX1"/>');sample=ET.fromstring('<SAMPLE><IDENTIFIERS><EXTERNAL_ID namespace="BioSample">SAMN2</EXTERNAL_ID></IDENTIFIERS></SAMPLE>')
        with self.assertRaises(ValueError):require_identity('GSM1',{'GSM1':dict(BioSample='SAMN1',experiment='SRX1')},experiment,sample)


if __name__=='__main__':unittest.main()
