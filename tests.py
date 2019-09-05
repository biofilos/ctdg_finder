import unittest
from ctdg_finder import load_config

class TestConfig(unittest.TestCase):
    config = load_config("./config.json")
    def test_load_config(self):
        self.assertIsInstance(self.config, dict)
     

if __name__ == "__main__":
    unittest.main()
    