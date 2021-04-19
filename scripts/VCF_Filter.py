#!/usr/bin/env python


import click
import vcfpy
from os import path
import warnings




@click.command()
@click.option("-vi", help="Input VCF File", type=str, default=None, show_default=True)
@click.option("-vo", help="Output VCF File. If not given, will be printed in terminal", type=str, default=None, show_default=True)
@click.option("-mmc", help="Minor allele minimum count for each strand", type=int, default=5, show_default=True)
@click.option("-mmf", help="Minor allele minimum fraction of D4", type=float, default=0.03, show_default=True)
@click.option("-mmfs", help="Minor allele minimum strand fraction of D4", type=float, default=0.02, show_default=True)
@click.option("-mff", help="Minor allele fixation fraction of D4", type=float, default=0.99, show_default=True)
def run(vi, vo, mmc, mmf, mmfs, mff):
    """
    Filter VCF based on given parameter 
    """
    mff = 0.99
    if not vi:
        exit("VCF file path not given. Exiting . . . . . .")
    if not path.isfile(vi):
        exit("Given file path doesn't exist or it is not a file. Exiting . . . . .")
    try:
        reader = vcfpy.Reader.from_path(vi)
        reader.header.add_filter_line(vcfpy.OrderedDict([
            ('ID', 'min_alt_count'), ('Description', f'Minimum {mmc} from each atrand')]))
        reader.header.add_filter_line(vcfpy.OrderedDict([
            ('ID', 'min_alt_frac'), ('Description', f'Minimum {mmf} of total DP4')]))
        reader.header.add_filter_line(vcfpy.OrderedDict([
            ('ID', 'min_alt_fication_frac'), ('Description', f'Minimum {mff} of total DP4')]))
        writer = vcfpy.Writer.from_path('/dev/stdout', reader.header) if not vo else vcfpy.Writer.from_path(vo, reader.header)

        for record in reader:
            base_sum = sum(record.INFO["DP4"])
            if((record.INFO["DP4"][2] >= mmc) & (record.INFO["DP4"][2] >= mmfs*(record.INFO["DP4"][2]+record.INFO["DP4"][3])) & 
            (record.INFO["DP4"][3] >= mmc) & (record.INFO["DP4"][3] >= mmfs*(record.INFO["DP4"][2]+record.INFO["DP4"][3])) & 
            (record.INFO["DP4"][2]+record.INFO["DP4"][3] >= mmf* base_sum) ):
                if record.INFO["DP4"][2]+record.INFO["DP4"][3] >= mff* base_sum:
                    record.INFO["DP4"][0]=0
                    record.INFO["DP4"][1]=0
                    record.INFO["AF"] = 1.0
                writer.write_record(record)
    except:
        exit("Please check given file. Is it in proper format??. Exiting . . . . .")


if __name__=='__main__':
    warnings.filterwarnings("ignore")
    run()