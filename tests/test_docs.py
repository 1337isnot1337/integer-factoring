"""Unit tests for the offline documentation-link checker."""
import importlib.util
import unittest
from pathlib import Path

PATH = Path(__file__).resolve().parents[1] / "tools" / "check_docs.py"
SPEC = importlib.util.spec_from_file_location("check_docs", PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class LinkCheckerTests(unittest.TestCase):
    def setUp(self):
        self.files = {"README.md", "research/proof.md", "research/INDEX.md"}

    def test_valid_links_and_fences(self):
        markdown = (
            "[proof](research/proof.md#proof)\n"
            "[directory](research/)\n"
            "[external](https://example.com/missing.md)\n"
            "[section](#methods)\n"
            "[relative](./README.md)\n"
            "[linked file]: <research/proof.md>\n"
            "~~~md\n[not a real link](missing.md)\n~~~\n"
        )
        self.assertEqual(MODULE.link_targets("README.md", markdown, self.files), [])

    def test_reports_missing_file(self):
        self.assertEqual(
            MODULE.link_targets("research/INDEX.md", "[missing](absent.md)", self.files),
            ["research/INDEX.md:1: missing research/absent.md (from absent.md)"],
        )

    def test_rejects_escaping_repository(self):
        self.assertEqual(
            MODULE.link_targets("README.md", "[missing](../elsewhere.md)", self.files),
            ["README.md:1: outside repository: ../elsewhere.md"],
        )

    def test_percent_encoded_file(self):
        files = self.files | {"research/data table.md"}
        self.assertEqual(
            MODULE.link_targets("research/INDEX.md", "[data](data%20table.md)", files), []
        )


if __name__ == "__main__":
    unittest.main()
