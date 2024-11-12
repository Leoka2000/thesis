
from django.shortcuts import render
from rdkit import Chem
from rdkit.Chem import Draw
import base64
from io import BytesIO

def oxidoreductases_view(request):
    # Define reagents and products as SMILES
    reagents_smiles = ["C1=CC(=CC=C1O)O", "O=O"]  # 4-benzenediol (1,2-benzenediol) and O2
    products_smiles = ["C1=CC(=CC=C1[O])O", "O"]  # benzosemiquinone and water

    # Generate molecule images and convert them to base64
    def smiles_to_base64(smiles):
        mol = Chem.MolFromSmiles(smiles)
        img = Draw.MolToImage(mol, size=(300, 300))
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        return base64.b64encode(buffered.getvalue()).decode("utf-8")

    # Convert reagents and products to base64 images
    reagents_images = [smiles_to_base64(smiles) for smiles in reagents_smiles]
    products_images = [smiles_to_base64(smiles) for smiles in products_smiles]

    # Pass images to the template
    context = {
        "reagents_images": reagents_images,
        "products_images": products_images,
    }
    return render(request, "pages/oxidoreductases.html", context)