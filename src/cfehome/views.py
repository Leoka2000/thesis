from django.shortcuts import render
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
import base64
from io import BytesIO

def home_view(request):
    # Create an RDKit molecule from a SMILES string
    smiles = "CCO"  # Ethanol, for example
    molecule = Chem.MolFromSmiles(smiles)
    
    # Generate an image of the molecule
    img = Draw.MolToImage(molecule)
    
    # Convert the image to base64
    buffered = BytesIO()
    img.save(buffered, format="PNG")
    img_base64 = base64.b64encode(buffered.getvalue()).decode("utf-8")
    
    # Pass the image to the template
    context = {
        'molecule_image': img_base64,
        'smiles': smiles
    }
    
    return render(request, "pages/home.html", context)