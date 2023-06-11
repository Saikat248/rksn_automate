from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import pybel
import glob

def tanimoto_calc(xyzfile1, xyzfile2):
    xyz1 = next(pybel.readfile('xyz', xyzfile1))
    xyz2 = next(pybel.readfile('xyz', xyzfile2))
    smi1 = xyz1.write(format='smi').split()[0].strip()
    smi2 = xyz2.write(format='smi').split()[0].strip()
    mol1 = Chem.MolFromSmiles(smi1)
    mol2 = Chem.MolFromSmiles(smi2)
    fp1 = Chem.RDKFingerprint(mol1)
    fp2 = Chem.RDKFingerprint(mol2)

    # fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 3, nBits=2048)
    # fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 3, nBits=2048)
    # s = round(DataStructs.TanimotoSimilarity(fp1,fp2),3)
    s = round(DataStructs.FingerprintSimilarity(fp1,fp2),3)
    return s

conf_lst = glob.glob('conf*.xyz')

for i in range(0, len(conf_lst)):
    for j in range(i+1, len(conf_lst)):
        print(conf_lst[i], conf_lst[j])
        s = tanimoto_calc(conf_lst[i], conf_lst[j])
        print(s)


# print(tanimoto_calc('conf86.xyz', 'conf96.xyz'))