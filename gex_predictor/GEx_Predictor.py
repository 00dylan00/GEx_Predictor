import torch
from rdkit import Chem
from rdkit.Chem.MolStandardize  import rdMolStandardize
from signaturizer import Signaturizer
import sys
import numpy as np 
from pathlib import Path


project_root = Path(__file__).parent.parent.resolve()  
if str(project_root) not in sys.path:
    sys.path.append(str(project_root))

from model.models import GenomicExpressionNet2

# Ensure the model directory is in the Python path
model_dir = project_root / "model"
if str(model_dir) not in sys.path:
    sys.path.append(str(model_dir))

class GEx_Predictor():
    BEST_FOLD = 4
    ALL_FOLDS = [1, 2, 3, 4, 5]

    def __init__(self, all_folds = False):
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.all_folds = all_folds        
        print(f"[INFO] Using device: {self.device}")
        self.predictor = self._load_all_models()
        
    
    def _load_single_model(self, fold_idx: int=4):
        pt_path = project_root / "model" / f"fold_{fold_idx}.pt"
        print(f"[INFO] Loading fold {fold_idx} from {pt_path} ...")
        model = torch.load(pt_path, map_location=self.device, weights_only=False)
        model.eval()

        return model

    def _load_all_models(self,):
        # We can predict using our best fold or all the folds
        folds = self.ALL_FOLDS if self.all_folds else [self.BEST_FOLD]
        models = [self._load_single_model(i) for i in folds]
        label = f"all {len(models)} folds" if self.all_folds else "best fold only"
        print(f"[INFO] {len(models)} model(s) loaded ({label}).")
        return models
    
    def _standardize_smiles(self, input_smiles):
        standardize = []
        for smile in input_smiles:
            try:
                mol = Chem.MolFromSmiles(smile)
                if mol is None:
                    return None

                mol = rdMolStandardize.Cleanup(mol)

                mol = rdMolStandardize.LargestFragmentChooser().choose(mol)
                mol = rdMolStandardize.Uncharger().uncharge(mol)

                # Add tautomer canonicalization if needed
                mol = rdMolStandardize.TautomerEnumerator().Canonicalize(mol)

                Chem.SanitizeMol(mol)

                standardize.append(Chem.MolToSmiles(mol))

            except Exception as e:
                print(f"Could not process {input_smiles}: {e}")
                return None

        return standardize

    def _get_global_Signature(self, input_smiles):

        # We parse our standarized smiles to GLOBAL Signature
        signature = 'GLOBAL'
        sign = Signaturizer(signature)
        results = sign.predict(input_smiles)

        return results.signature[:]
    
    def _predict_single_model(self, model: torch.nn.Module, x: torch.Tensor):
        with torch.no_grad():
            return model(x).cpu().numpy()
        
    def predict(self, input_smile, input_type = "SMILES"):

        print(f"[INFO] Starting prediction — input_type={input_type}")
        print(f"ensemble={f'{len(self.predictor)} folds' if self.all_folds else 'best fold'}")
               
        if input_type.upper() == "SMILES":
            standardized_smiles = self._standardize_smiles(input_smile)
            if standardized_smiles is None:
                raise ValueError("[ERROR] Failed to standardize input SMILES.")
        else:
            raise NotImplementedError("Only 'SMILES' input_type is implemented.")

        signatures_to_predict = self._get_global_Signature(standardized_smiles)

        print("[STEP] Converting signatures to tensor and performing prediction...")
        x = torch.tensor(signatures_to_predict, dtype=torch.float32, device=self.device)

        fold_preds = [self._predict_single_model(m, x) for m in self.predictor]
        preds = np.mean(fold_preds, axis=0)

        print("[INFO] Prediction completed successfully.")
        return preds