# views.py
from django.shortcuts import render
from rdkit import Chem
from rdkit.Chem import Draw
import base64
from io import BytesIO

def generate_molecule_image(smiles, size=(300, 300)):
    try:
        # Remove any whitespace and newlines from SMILES
        smiles = smiles.strip()
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        # Generate 2D depiction
        img = Draw.MolToImage(mol, size=size)
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        return base64.b64encode(buffered.getvalue()).decode()
    except Exception as e:
        print(f"Error generating molecule image: {e}")
        return None

def lipases_view(request):
    # Simplified SMILES representations for the molecules
    # These are simplified versions that RDKit can handle better
    palm_oil = "CCCCCCCCCCCCCCCCCC(=O)OCC(COC(=O)CCCCCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCCCC"
    stearic_acid = "CCCCCCCCCCCCCCCCCC(=O)O"
    cbs_product = "CCCCCCCCCCCCCCCCCC(=O)OCC(COC(=O)CCCCCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCCCC"
    distearyl_ether = "CCCCCCCCCCCCCCCCCCOCCCCCCCCCCCCCCCCCC"

    # Generate images with error handling
    context = {}
    molecules = {
        'palm_oil_img': palm_oil,
        'stearic_acid_img': stearic_acid,
        'cbs_img': cbs_product,
        'distearyl_img': distearyl_ether
    }

    for key, smiles in molecules.items():
        img = generate_molecule_image(smiles)
        if img:
            context[key] = img
        else:
            context[key] = ''  # Provide empty string if molecule generation fails
            print(f"Failed to generate image for {key}")

    return render(request, 'pages/lipases.html', context)