import macromaker.utils as utils
import pytest

class TestUtils(object):
    def test_get_file_extension_with_a_regular_filename(self):
        assert utils.get_file_extension("file.tar") == "tar"

    def test_(self):
        pass