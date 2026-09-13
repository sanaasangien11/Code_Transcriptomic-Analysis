import pandas as pd
from tqdm import tqdm
import os


if not os.path.exists(f"tauschii_data/hisat2_index/tauschi.8.ht2"):
 chroms = ",".join([f"/nam-99/ablage/nam/ss_rnaseq/tauschii_data/CHROM/{x}" for x in os.listdir(f"/nam-99/ablage/nam/ss_rnaseq/tauschii_data/CHROM")])
 os.system(f"hisat2-build -p 7 {chroms} /nam-99/ablage/nam/ss_rnaseq/tauschi_data/hisat2_index/tauschi")

