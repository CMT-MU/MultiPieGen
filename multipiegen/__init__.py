__version__ = "0.1.0"

import os
import importlib.resources as res
import multipiegen

__top_dir__ = str(res.files(multipiegen).parent) + os.sep
__bin_dir__ = os.path.join(__top_dir__, "binary_data") + os.sep
