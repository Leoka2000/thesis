from django.shortcuts import render
from rdkit import Chem
from rdkit.Chem import Draw
import base64
from io import BytesIO

def home_view(request):
    # Define the SMILES strings for glucose monomers
    glucose_monomer_smiles = 'C(C1C(C(C(O1)O)O)O)O'
    glucose_monomer = Chem.MolFromSmiles(glucose_monomer_smiles)

    # Generate the image for each monomer and the combined form
    images = []
    for mol in [glucose_monomer, glucose_monomer]:  # To illustrate the two starting glucose monomers
        buffer = BytesIO()
        img = Draw.MolToImage(mol, size=(200, 200))
        img.save(buffer, format="PNG")
        images.append(base64.b64encode(buffer.getvalue()).decode("utf-8"))

    # Illustrate the polymerization (simplified view of two glucose monomers linked)
    # Create a placeholder molecule representing the linked form for visualization purposes
    polymer_smiles = 'C(C1C(C(C(O1)O)O)O)OC(C2C(C(C(O2)O)O)O)O'
    polymer_molecule = Chem.MolFromSmiles(polymer_smiles)
    polymer_image = Draw.MolToImage(polymer_molecule, size=(300, 300))

    # Convert to base64 for the template
    buffer = BytesIO()
    polymer_image.save(buffer, format="PNG")
    polymer_img_str = base64.b64encode(buffer.getvalue()).decode("utf-8")

    return render(request, "pages/home.html", {
        "monomer_images": images,
        "polymer_image": polymer_img_str
    })
