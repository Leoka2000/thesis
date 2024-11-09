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


def lipases_view(request):
    return render(request, "pages/lipases.html", {})


def hydrolases_view(request):
    # Define the SMILES strings
    pet_monomer_smiles = 'O=C(OCC)C1=CC=CC=C1C(=O)OCC'
    pet_polymer_smiles = 'O=C(OCC)C1=CC=CC=C1C(=O)OCCC(=O)C2=CC=CC=C2C(=O)OCCC(=O)C3=CC=CC=C3C(=O)OCC'
    vanillin_smiles = 'COC1=CC=C(C=C1)C=O'
    cyclic_acetal_smiles = 'O1CCOC1'
    propadienol_smiles = 'C=CCO'  # SMILES for 1,3-propadienol

    # Generate the images throughtb the .Chem AND BytesIO
    pet_monomer = Chem.MolFromSmiles(pet_monomer_smiles)
    buffer = BytesIO()
    pet_image = Draw.MolToImage(pet_monomer, size=(200, 200))
    pet_image.save(buffer, format="PNG")
    pet_img_str = base64.b64encode(buffer.getvalue()).decode("utf-8")

    # Generate and encode the image for the PET polymer (4 monomers)
    pet_polymer = Chem.MolFromSmiles(pet_polymer_smiles)
    buffer = BytesIO()
    pet_polymer_image = Draw.MolToImage(pet_polymer, size=(400, 400))
    pet_polymer_image.save(buffer, format="PNG")
    pet_polymer_img_str = base64.b64encode(buffer.getvalue()).decode("utf-8")

    # # Generate the images throughtb the .Chem AND BytesIOfor Vanillin
    vanillin_molecule = Chem.MolFromSmiles(vanillin_smiles)
    buffer = BytesIO()
    vanillin_image = Draw.MolToImage(vanillin_molecule, size=(200, 200))
    vanillin_image.save(buffer, format="PNG")
    vanillin_img_str = base64.b64encode(buffer.getvalue()).decode("utf-8")

    # # Generate the images throughtb the .Chem AND BytesIO Cyclic Acetal
    cyclic_acetal_molecule = Chem.MolFromSmiles(cyclic_acetal_smiles)
    buffer = BytesIO()
    cyclic_acetal_image = Draw.MolToImage(cyclic_acetal_molecule, size=(200, 200))
    cyclic_acetal_image.save(buffer, format="PNG")
    cyclic_acetal_img_str = base64.b64encode(buffer.getvalue()).decode("utf-8")

    #  # Generate the images throughtb the .Chem AND BytesIO1,3-Propadienol
    propadienol_molecule = Chem.MolFromSmiles(propadienol_smiles)
    buffer = BytesIO()
    propadienol_image = Draw.MolToImage(propadienol_molecule, size=(200, 200))
    propadienol_image.save(buffer, format="PNG")
    propadienol_img_str = base64.b64encode(buffer.getvalue()).decode("utf-8")

    return render(request, "pages/hydrolases.html", {
        "pet_image": pet_img_str,
        "pet_polymer_image": pet_polymer_img_str,
        "vanillin_image": vanillin_img_str,
        "cyclic_acetal_image": cyclic_acetal_img_str,
        "propadienol_image": propadienol_img_str
    })


def transferases_view(request):
    return render(request, "pages/transferases.html", {})

def oxidoreductases_view(request):
    return render(request, "pages/oxidoreductases.html", {})