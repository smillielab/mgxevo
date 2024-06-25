import pkg_resources
from Bio.Phylo import read as read_tree

from . import majority_snp_utils as snp_util

script_path = pkg_resources.resource_filename('mgxevo', 'tl/tree/')