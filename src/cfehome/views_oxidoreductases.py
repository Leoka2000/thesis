from django.shortcuts import render
from rdkit import Chem
from rdkit.Chem import Draw
import base64
from io import BytesIO

def oxidoreductases_view(request):
    # Define the SMILES representations for reagents and products
    reagents_smiles = ["C1=CC(=CC=C1O)O", "O=O"]  # 4-benzenediol (1,2-benzenediol) and O₂
    products_smiles = ["C1=CC(=CC=C1[O])O", "O"]   # Benzosemiquinone and H₂O
    
    # Set default for no reaction available
    reaction_available = False
    reagents_images = []
    products_images = []

    # Check if 'specific_enzyme' and 'chemical_group' are in request GET parameters
    selected_enzyme = request.GET.get('specific_enzyme', '').lower()
    selected_chemical = request.GET.get('chemical_group', '').lower()
    
    # Check for the specific enzyme and chemical group combination
    if selected_enzyme == "laccase" and selected_chemical == "carbonyl":
        reaction_available = True
        # Generate images for reagents and products
        for smi in reagents_smiles:
            mol = Chem.MolFromSmiles(smi)
            img = Draw.MolToImage(mol)
            buffer = BytesIO()
            img.save(buffer, format="PNG")
            encoded_image = base64.b64encode(buffer.getvalue()).decode("utf-8")
            reagents_images.append(encoded_image)

        for smi in products_smiles:
            mol = Chem.MolFromSmiles(smi)
            img = Draw.MolToImage(mol)
            buffer = BytesIO()
            img.save(buffer, format="PNG")
            encoded_image = base64.b64encode(buffer.getvalue()).decode("utf-8")
            products_images.append(encoded_image)

    context = {
        "reaction_available": reaction_available,
        "reagents_images": reagents_images,
        "products_images": products_images,
    }
    
    return render(request, "pages/oxidoreductases.html", context)
