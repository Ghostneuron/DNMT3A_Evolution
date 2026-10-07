import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from human_mouse_promoter_synteny import (
    forward_query_coordinate,
    transformed_strand,
)


class HumanMousePromoterSyntenyTests(unittest.TestCase):
    def test_forward_query_coordinates(self):
        self.assertEqual(forward_query_coordinate(1000, "+", 100), 100)
        self.assertEqual(forward_query_coordinate(1000, "-", 100), 899)

    def test_strand_transformation(self):
        self.assertEqual(transformed_strand("-", "+"), "-")
        self.assertEqual(transformed_strand("-", "-"), "+")


if __name__ == "__main__":
    unittest.main()
