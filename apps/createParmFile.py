import hapi as hp
import numpy as np
import json

from tempfile import TemporaryDirectory
from struct import calcsize, pack, unpack
from pathlib import Path
from contextlib import contextmanager
from os import chdir as os_chdir

@contextmanager
def chDir(dest):
    cwd = Path.cwd()
    try:
        os_chdir(dest)
        yield
    finally:
        os_chdir(cwd)

def _getMolecToIsoNum() -> dict:
    ans = dict()
    molec = list({v[-1] for _,v in hp.ISO_ID.items()})
    for m in molec:
        ans[m] = [v[0] for _,v in hp.ISO.items() if v[4]==m]
    return ans

def hapiDownload(json_input: dict):
    """
    Using input for HalIR to download and create a file with the parameters
    needed for the calulation uses Hapi
    link
    """

    _entries = ['molec_id', 'local_iso_id', 'nu', 'sw', 'a', 'gamma_air',
                'gamma_self', 'elower', 'n_air', 'delta_air',
                'global_upper_quanta', 'global_lower_quanta',
                'local_upper_quanta', 'local_lower_quanta', 'ierr', 'iref',
                'line_mixing_flag', 'gp', 'gpp']

    struct_header = '6s12siffiP'
    struct_hline = 'i'*2+'d'*8+'16s'*4+'6s'+'12s'+'2s'+'d'*4

    isoTxtToNum = {v[1]: v[0] for _, v in hp.ISO.items()}
    molecToIsoNum = _getMolecToIsoNum()

    with TemporaryDirectory() as tmpDir:
        with chDir(tmpDir):
            rootDir = Path(json_input['project']['rootDir'])
            hp.db_begin(json_input['project']['pname'])
            for comp in json_input['composition']:
                filename = rootDir/Path(comp['molec']+'.hpar')
                if filename.exists():
                    print("Warning file exists what to do")

                if comp['isotop'] == 'Natural':
                    molec = comp['molec']
                    molecs = molecToIsoNum[molec]
                    hp.fetch_by_ids(molec,
                            molecs,
                            *json_input['sampleEnv']['ROI'])

                else:
                    molec = comp['molec']
                    molecs = [isoTxtToNum[comp['isotop']]]
                    hp.fetch_by_ids(molec,
                            molecs,
                            *json_input['sampleEnv']['ROI'])

                line_data = np.array(
                        hp.getColumns(molec, _entries)).transpose()

                with open(filename, "wb") as outf:
                   mbin = pack(struct_header, molec.encode(),
                              comp['isotop'].encode(), len(molecs),
                              json_input['sampleEnv']['ROI'][0], json_input['sampleEnv']['ROI'][1], len(line_data), 0)
                   outf.write(mbin)
                   mbin = pack('i'*len(molecs), *molecs)
                   outf.write(mbin)
                   for hl in line_data:
                       mbin = pack(struct_hline, int(hl[0]),
                                  int(hl[1]), float(hl[2]),
                                  float(hl[3]), float(hl[4]),
                                  float(hl[5]), float(hl[6]),
                                  float(hl[7]), float(hl[8]),
                                  float(hl[9]), hl[10].encode(),
                                  hl[11].encode(), hl[12].encode(),
                                  hl[13].encode(), hl[14].encode(),
                                  hl[15].encode(), hl[16].encode(),
                                  float(hl[17]), float(hl[18]),
                                  hp.abundance(int(hl[0]), int(hl[1])),
                                  hp.molecularMass(int(hl[0]), int(hl[1])))
                       outf.write(mbin)

def main(infile: Path) -> None:
    """
    Parse the HalIR input file and download the needed data to the project dir

    input: infile is a json file
    """
    with open(infile) as json_input:
        input_dict = json.load(json_input)

    for key in ['project', 'sampleEnv', 'composition']:
        if key not in input_dict['input']:
            raise KeyError(f"The field: \"{key}\" is missing from input file \"{infile.name}\" in dir \"{infile.parent}\"")

    hapiDownload(input_dict['input'])

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", "-i", type=str,
            help="JSON formated input file, according to HalIR documentation")
    parser.add_argument("--listMolec", "-lM", action='store_true',
            help="List molecules in Hapi database")

    args = parser.parse_args()
    if args.input:
        infile = Path(args.input)
        if infile.exists() and infile.is_file:
            main(infile)
    elif args.listMolec:
        mol = list({v[-1] for _,v in hp.ISO_ID.items()})
        print(len(mol))
        molec = {n:v for n,v in enumerate(list({v[-1] for _,v in hp.ISO_ID.items()}))}
        for np in zip(*[iter(molec.items())]*5):
            for p in np:
                print(f'{p[0]} {p[1]}', end=", ")
            print("")
