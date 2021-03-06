# -*- coding: utf-8 -*-
import hapi as hp
import numpy as np
from tempfile import TemporaryDirectory
from contextlib import contextmanager
from struct import calcsize, pack, unpack
import os


@contextmanager
def chDir(dest):
    try:
        cwd = os.getcwd()
        os.chdir(dest)
        yield
    finally:
        os.chdir(cwd)


def _getMolecToIsoNum():
    ans = {}
    molec = _getMolecs()
    for m in molec:
        ans[m] =\
            [v[0] for k, v in hp.ISO.items() if v[4] == m]
    return ans


def _getMolecs():
    ans = []
    for k, v in hp.ISO_ID.items():
        if v[-1] not in ans:
            ans.append(v[-1])
    return ans


def hitranDownload(projDict, sampleDict):
    """ Facilitates download of parameters to HalIR \
        using Hapi. And saves them as binary file """

    _entries = ['molec_id', 'local_iso_id', 'nu', 'sw', 'a', 'gamma_air',
                'gamma_self', 'elower', 'n_air', 'delta_air',
                'global_upper_quanta', 'global_lower_quanta',
                'local_upper_quanta', 'local_lower_quanta', 'ierr', 'iref',
                'line_mixing_flag', 'gp', 'gpp']

    struct_header = '6s12siffP'
    struct_hline = 'i'*2+'f'*8+'16s'*4+'6s'+'12s'+'2s'+'f'*4

    isoTxtToNum =\
        {v[1]: int(str(k[0])+str(k[1])) for k, v in hp.ISO.items()}
    molecToIsoNum = _getMolecToIsoNum()
    components = sampleDict['comp']
    # Download the files
    with TemporaryDirectory() as tmpDir:
        with chDir(tmpDir):
            print(tmpDir)
            wDir = projDict['pdir']
            hp.db_begin(projDict['pname'])
            for comp in components:
                # print(comp)
                # print(self.molecToIsoNum[comp['molec']])
                wDir = projDict['pdir']
                filename = os.path.join(wDir, comp['molec']+'.hpar')
                if os.path.exists(filename):
                    with open(filename, 'rb') as headfil:
                        head = unpack(struct_header,
                                      headfil.read(calcsize(struct_header)))
                    # if (head[0].decode().strip() == comp['molec'])
                    if(head[0].decode().strip('\x00') == comp['molec'] and
                       head[1].decode().strip('\x00') == comp['isotop'] and
                       head[2] == len(molecToIsoNum[comp['molec']]) and
                       head[3] == sampleDict['ROI'][0] and
                       head[4] == sampleDict['ROI'][1]):
                        continue

                if comp['isotop'] == 'Natural':
                    molec = comp['molec']
                    molecs = molecToIsoNum[molec]
                    hp.fetch_by_ids(molec,
                                    molecs,
                                    *sampleDict['ROI'])

                else:
                    molec = comp['molec']
                    molecs = isoTxtToNum[comp['isotop']]
                    hp.fetch_by_ids(molec,
                                    molecs,
                                    *sampleDict['ROI'])

                line_data = np.array(
                            hp.getColumns(molec, _entries)).transpose()

                with open(filename, 'wb') as outf:
                    bin = pack(struct_header, molec.encode(),
                               comp['isotop'].encode(), len(molecs),
                               sampleDict['ROI'][0], sampleDict['ROI'][1], 0)
                    outf.write(bin)
                    bin = pack('i'*len(molecs), *molecs)
                    outf.write(bin)
                    for hl in line_data:
                        bin = pack(struct_hline, int(hl[0]),
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
                        outf.write(bin)


def parseDictToHalIR(projDict, sampleDict):
    """ Parse the projDict and sampleDict to HalIR input str """
    pass
