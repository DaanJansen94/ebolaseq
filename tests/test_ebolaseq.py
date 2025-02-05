import unittest
from ebolaseq.core import get_reference_length, standardize_country, get_location

class TestEbolaSeq(unittest.TestCase):
    def test_reference_length(self):
        """Test reference genome length retrieval"""
        self.assertEqual(get_reference_length('Zaire'), 18959)
        self.assertEqual(get_reference_length('Sudan'), 18875)
        self.assertEqual(get_reference_length('Unknown'), 18959)  # Default value

    def test_country_standardization(self):
        """Test country name standardization"""
        self.assertEqual(standardize_country('DRC'), 'Democratic_Republic_of_the_Congo')
        self.assertEqual(standardize_country('Guinea'), 'Republic_of_Guinea')
        self.assertEqual(standardize_country('Conakry'), 'Guinea')
        self.assertEqual(standardize_country('Unknown'), 'Unknown')

    def test_location_parsing(self):
        """Test location parsing from GenBank records"""
        # Mock GenBank record
        class MockFeature:
            def __init__(self, qualifiers):
                self.type = "source"
                self.qualifiers = qualifiers

        class MockRecord:
            def __init__(self, features):
                self.features = features

        # Test with country qualifier
        record1 = MockRecord([MockFeature({'country': ['Guinea']})])
        self.assertEqual(get_location(record1), 'Republic_of_Guinea')

        # Test with geo_loc_name qualifier
        record2 = MockRecord([MockFeature({'geo_loc_name': ['Conakry']})])
        self.assertEqual(get_location(record2), 'Guinea')

        # Test with no location info
        record3 = MockRecord([MockFeature({})])
        self.assertEqual(get_location(record3), 'Unknown')

if __name__ == '__main__':
    unittest.main() 