from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from Bio.Seq import Seq
HTTP = HTTPRemoteProvider()

configfile: "config.json"
